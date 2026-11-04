# Pseudo-nCounter conversion and frozen ENET validation in public cohorts
rm(list = ls()); options(stringsAsFactors = FALSE); set.seed(20261009)
suppressPackageStartupMessages({library(data.table); library(dplyr); library(GSVA); library(pROC); library(readxl)})
path_env <- function(key) {value <- Sys.getenv(key, unset = ""); if (!nzchar(value)) stop("Missing environment variable: ", key); value}
retro_bulk_file <- path_env("RETRO_BULK_FILE")
retro_nc_file <- path_env("RETRO_NCOUNTER_FILE")
pass_file <- path_env("RETRO_QC_FILE")
gene_qc_file <- path_env("GENE_QC_FILE")
signature_file <- path_env("SIGNATURE_FILE")
pathway_file <- path_env("PATHWAY_FILE")
model_file <- path_env("ENET_MODEL_FILE")
bulk_files <- c(TCGA = path_env("TCGA_BULK_FILE"), CGGA = path_env("CGGA_BULK_FILE"), CPTAC = path_env("CPTAC_BULK_FILE"))
label_files <- c(TCGA = path_env("TCGA_LABEL_FILE"), CGGA = path_env("CGGA_LABEL_FILE"), CPTAC = path_env("CPTAC_LABEL_FILE"))
out_dir <- path_env("OUTPUT_DIR")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
classes <- c("C1", "C2", "C3")
read_expression <- function(file, log_transform = FALSE) {
  x <- fread(file, data.table = FALSE, check.names = FALSE)
  m <- as.matrix(x[, -1, drop = FALSE]); storage.mode(m) <- "numeric"; rownames(m) <- as.character(x[[1]])
  if (log_transform) m <- log2(pmax(m, 0) + 1)
  m
}
read_pathways <- function(file) {
  sheets <- excel_sheets(file)
  setNames(lapply(sheets, function(sheet) {
    x <- suppressMessages(read_excel(file, sheet = sheet, col_names = FALSE))
    unique(na.omit(as.character(x[[1]])))
  }), sheets)
}
ssgsea <- function(expression, sets) {
  if ("ssgseaParam" %in% getNamespaceExports("GSVA")) {
    as.matrix(gsva(ssgseaParam(expression, sets, alpha = 0.25, normalize = FALSE, checkNA = "auto"), verbose = FALSE))
  } else {
    as.matrix(gsva(expression, sets, method = "ssgsea", abs.ranking = TRUE, ssgsea.norm = FALSE, verbose = FALSE))
  }
}
rb <- read_expression(retro_bulk_file)
rn <- read_expression(retro_nc_file, TRUE)
pass <- fread(pass_file, data.table = FALSE)$Sample
rb <- rb[, pass, drop = FALSE]; rn <- rn[, pass, drop = FALSE]
gene_qc <- fread(gene_qc_file, data.table = FALSE)
retained <- gene_qc$Gene[gene_qc$Retained_for_ML]
training_genes <- Reduce(intersect, list(retained, rownames(rb), rownames(rn)))
converter <- bind_rows(lapply(training_genes, function(gene) {
  x <- as.numeric(rb[gene, ]); y <- as.numeric(rn[gene, ]); ok <- is.finite(x) & is.finite(y)
  fit <- lm(y[ok] ~ x[ok])
  data.frame(Gene = gene, Intercept = coef(fit)[1], Slope = coef(fit)[2],
             Median = median(y[ok]), Minimum = min(y[ok]), Maximum = max(y[ok]))
}))
rownames(converter) <- converter$Gene
convert_bulk <- function(bulk) {
  result <- matrix(NA_real_, nrow(converter), ncol(bulk), dimnames = list(converter$Gene, colnames(bulk)))
  for (gene in converter$Gene) {
    parameters <- converter[gene, ]
    value <- if (gene %in% rownames(bulk)) parameters$Intercept + parameters$Slope * bulk[gene, ] else rep(parameters$Median, ncol(bulk))
    value[!is.finite(value)] <- parameters$Median
    result[gene, ] <- pmin(pmax(value, parameters$Minimum), parameters$Maximum)
  }
  result
}
signature <- as.data.frame(read_excel(signature_file))
signature_genes <- intersect(as.character(signature[[1]]), training_genes)
pathways <- lapply(read_pathways(pathway_file), function(g) intersect(g, training_genes))
pathways <- pathways[lengths(pathways) >= 3]
model <- readRDS(model_file)
make_input <- function(pseudo) {
  gene_matrix <- t(pseudo[signature_genes, , drop = FALSE])
  colnames(gene_matrix) <- paste0("GENE__", make.names(colnames(gene_matrix), unique = TRUE))
  pathway_matrix <- t(ssgsea(pseudo, pathways))
  colnames(pathway_matrix) <- paste0("PATH__", make.names(colnames(pathway_matrix), unique = TRUE))
  x <- as.data.frame(cbind(gene_matrix, pathway_matrix), check.names = FALSE)
  missing_features <- setdiff(model$Features, names(x))
  for (feature in missing_features) x[[feature]] <- NA_real_
  list(data = x[, model$Features, drop = FALSE], missing = missing_features)
}
metrics <- function(truth, pred, probability) {
  truth <- factor(truth, levels = classes); pred <- factor(pred, levels = classes); tab <- table(truth, pred)
  accuracy <- sum(diag(tab)) / sum(tab); pe <- sum(rowSums(tab) * colSums(tab)) / sum(tab)^2
  f1 <- sapply(classes, function(k) {tp <- tab[k, k]; fp <- sum(tab[, k]) - tp; fn <- sum(tab[k, ]) - tp; if (2 * tp + fp + fn == 0) NA_real_ else 2 * tp / (2 * tp + fp + fn)})
  total <- sum(tab); den <- sqrt((total^2 - sum(colSums(tab)^2)) * (total^2 - sum(rowSums(tab)^2)))
  aucs <- sapply(classes, function(k) as.numeric(auc(roc(as.integer(truth == k), probability[, k], quiet = TRUE, direction = "<"))))
  data.frame(N = length(truth), Accuracy = accuracy, Kappa = (accuracy - pe) / (1 - pe),
             Macro_F1 = mean(f1, na.rm = TRUE), MCC = (sum(diag(tab)) * total - sum(colSums(tab) * rowSums(tab))) / den,
             Macro_AUC = mean(aucs), AUC_C1 = aucs[1], AUC_C2 = aucs[2], AUC_C3 = aucs[3])
}
prediction_list <- list(); performance_list <- list(); coverage_list <- list()
for (cohort in names(bulk_files)) {
  bulk <- read_expression(bulk_files[[cohort]])
  model_input <- make_input(convert_bulk(bulk))
  probability <- predict(model$Fit, model_input$data, type = "prob")[, classes, drop = FALSE]
  predicted <- classes[max.col(probability, ties.method = "first")]
  labels <- fread(label_files[[cohort]], data.table = FALSE)[, c("Sample", "Cluster")]
  labels <- labels[match(rownames(model_input$data), labels$Sample), , drop = FALSE]
  prediction_list[[cohort]] <- data.frame(Dataset = cohort, Sample = rownames(model_input$data),
                                          Cluster = labels$Cluster, Predicted = predicted, probability, check.names = FALSE)
  performance_list[[cohort]] <- cbind(Dataset = cohort, metrics(labels$Cluster, predicted, probability))
  coverage_list[[cohort]] <- data.frame(Dataset = cohort, Required = length(model$Features),
                                         Available = length(model$Features) - length(model_input$missing),
                                         Missing = length(model_input$missing))
}
predictions <- bind_rows(prediction_list)
performance <- bind_rows(performance_list)
coverage <- bind_rows(coverage_list)
fwrite(predictions, file.path(out_dir, "public_cohort_predictions.csv"))
fwrite(performance, file.path(out_dir, "public_cohort_performance.csv"))
fwrite(coverage, file.path(out_dir, "feature_coverage.csv"))
fwrite(converter, file.path(out_dir, "conversion_parameters.csv"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "session_info.txt"))
