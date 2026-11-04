# Independent prospective validation of the frozen ENET classifier
rm(list = ls()); options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({library(data.table); library(dplyr); library(GSVA); library(pROC); library(readxl)})
path_env <- function(key) {value <- Sys.getenv(key, unset = ""); if (!nzchar(value)) stop("Missing environment variable: ", key); value}
nc_file <- path_env("PRO_NCOUNTER_FILE")
label_file <- path_env("PRO_LABEL_FILE")
ffpe_label_file <- path_env("FFPE_LABEL_FILE")
gene_qc_file <- path_env("GENE_QC_FILE")
signature_file <- path_env("SIGNATURE_FILE")
pathway_file <- path_env("PATHWAY_FILE")
model_file <- path_env("ENET_MODEL_FILE")
out_dir <- path_env("OUTPUT_DIR")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
classes <- c("C1", "C2", "C3")
read_expression <- function(file) {
  x <- fread(file, data.table = FALSE, check.names = FALSE)
  m <- as.matrix(x[, -1, drop = FALSE]); storage.mode(m) <- "numeric"; rownames(m) <- as.character(x[[1]])
  log2(pmax(m, 0) + 1)
}
read_labels <- function(file) {
  extension <- tolower(tools::file_ext(file))
  x <- if (extension %in% c("xlsx", "xls")) as.data.frame(read_excel(file)) else fread(file, data.table = FALSE)
  x[, c("Sample", "Cluster")]
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
nc <- read_expression(nc_file)
labels <- bind_rows(read_labels(label_file), read_labels(ffpe_label_file))
nc <- nc[, labels$Sample, drop = FALSE]
gene_qc <- fread(gene_qc_file, data.table = FALSE)
retained <- gene_qc$Gene[gene_qc$Retained_for_ML]
signature <- as.data.frame(read_excel(signature_file))
signature_genes <- Reduce(intersect, list(as.character(signature[[1]]), retained, rownames(nc)))
pathways <- lapply(read_pathways(pathway_file), function(g) Reduce(intersect, list(g, retained, rownames(nc))))
pathways <- pathways[lengths(pathways) >= 3]
gene_matrix <- t(nc[signature_genes, , drop = FALSE])
colnames(gene_matrix) <- paste0("GENE__", make.names(colnames(gene_matrix), unique = TRUE))
pathway_matrix <- t(ssgsea(nc, pathways))
colnames(pathway_matrix) <- paste0("PATH__", make.names(colnames(pathway_matrix), unique = TRUE))
newdata <- as.data.frame(cbind(gene_matrix, pathway_matrix), check.names = FALSE)
model <- readRDS(model_file)
missing_features <- setdiff(model$Features, names(newdata))
for (feature in missing_features) newdata[[feature]] <- NA_real_
newdata <- newdata[, model$Features, drop = FALSE]
probability <- predict(model$Fit, newdata, type = "prob")[, classes, drop = FALSE]
predicted <- classes[max.col(probability, ties.method = "first")]
truth <- labels$Cluster
predictions <- data.frame(Sample = rownames(newdata), Cluster = truth, Predicted = predicted, probability, check.names = FALSE)
performance <- metrics(predictions$Cluster, predictions$Predicted, as.matrix(predictions[, classes]))
fwrite(predictions, file.path(out_dir, "prospective_predictions.csv"))
fwrite(performance, file.path(out_dir, "prospective_performance.csv"))
fwrite(data.frame(Feature = missing_features), file.path(out_dir, "missing_features.csv"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "session_info.txt"))
