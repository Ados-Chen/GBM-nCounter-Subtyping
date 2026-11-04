# Repeated nested cross-validation and frozen ENET training
rm(list = ls()); options(stringsAsFactors = FALSE); set.seed(20260920)
suppressPackageStartupMessages({
  library(caret); library(data.table); library(dplyr); library(glmnet)
  library(kernlab); library(pROC); library(randomForest); library(readxl); library(xgboost)
})
path_env <- function(key) {value <- Sys.getenv(key, unset = ""); if (!nzchar(value)) stop("Missing environment variable: ", key); value}
input_file <- path_env("ML_INPUT_FILE")
pass_file <- path_env("RETRO_QC_FILE")
out_dir <- path_env("OUTPUT_DIR")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
classes <- c("C1", "C2", "C3")
models <- c("ENET", "RF", "SVM_RFE", "XGBoost")
outer_k <- 5L; inner_k <- 4L; n_repeats <- 20L
dat <- as.data.frame(read_excel(input_file, sheet = "ML_nCounter_Input"), check.names = FALSE)
pass <- fread(pass_file, data.table = FALSE)$Sample
dat <- dat[dat$Sample %in% pass & dat$Cluster %in% classes, , drop = FALSE]
features <- grep("^(GENE__|PATH__)", names(dat), value = TRUE)
x <- dat[, features, drop = FALSE]
x[] <- lapply(x, function(v) suppressWarnings(as.numeric(v)))
y <- factor(dat$Cluster, levels = classes); rownames(x) <- dat$Sample
metrics <- function(truth, pred, prob) {
  truth <- factor(truth, levels = classes); pred <- factor(pred, levels = classes); tab <- table(truth, pred)
  accuracy <- sum(diag(tab)) / sum(tab); pe <- sum(rowSums(tab) * colSums(tab)) / sum(tab)^2
  f1 <- sapply(classes, function(k) {tp <- tab[k, k]; fp <- sum(tab[, k]) - tp; fn <- sum(tab[k, ]) - tp; if (2 * tp + fp + fn == 0) NA_real_ else 2 * tp / (2 * tp + fp + fn)})
  total <- sum(tab); den <- sqrt((total^2 - sum(colSums(tab)^2)) * (total^2 - sum(rowSums(tab)^2)))
  aucs <- sapply(classes, function(k) as.numeric(auc(roc(as.integer(truth == k), prob[, k], quiet = TRUE, direction = "<"))))
  c(Accuracy = accuracy, Kappa = (accuracy - pe) / (1 - pe), Macro_F1 = mean(f1, na.rm = TRUE),
    MCC = (sum(diag(tab)) * total - sum(colSums(tab) * rowSums(tab))) / den,
    Macro_AUC = mean(aucs), AUC_C1 = aucs[1], AUC_C2 = aucs[2], AUC_C3 = aucs[3])
}
summary_fun <- function(data, lev = NULL, model = NULL) metrics(data$obs, data$pred, as.matrix(data[, classes, drop = FALSE]))
fit_model <- function(name, x_train, y_train, seed, final = FALSE) {
  set.seed(seed)
  control_args <- list(method = if (final) "repeatedcv" else "cv", number = inner_k,
                       classProbs = TRUE, summaryFunction = summary_fun,
                       savePredictions = "final", allowParallel = FALSE)
  if (final) control_args$repeats <- 5
  control <- do.call(trainControl, control_args)
  preprocessing <- c("medianImpute", "nzv", "center", "scale")
  selected <- colnames(x_train)
  if (name == "SVM_RFE") {
    rfe_control <- rfeControl(functions = caretFuncs, method = "cv", number = inner_k, allowParallel = FALSE)
    sizes <- unique(pmin(c(10, 20, 40, ncol(x_train)), ncol(x_train)))
    selected <- predictors(rfe(x_train, y_train, sizes = sizes, rfeControl = rfe_control,
                               method = "svmRadial", preProcess = preprocessing, tuneLength = 3))
  }
  xx <- x_train[, selected, drop = FALSE]
  if (name == "ENET") {
    grid <- expand.grid(alpha = c(0, 0.5, 1), lambda = 10^seq(-4, 0, length.out = 5))
    fit <- train(xx, y_train, method = "glmnet", trControl = control, metric = "Macro_AUC",
                 preProcess = preprocessing, tuneGrid = grid)
  } else if (name == "RF") {
    grid <- data.frame(mtry = unique(pmax(1, round(c(sqrt(ncol(xx)), ncol(xx) / 5, ncol(xx) / 3)))))
    fit <- train(xx, y_train, method = "rf", trControl = control, metric = "Macro_AUC",
                 preProcess = preprocessing, tuneGrid = grid, ntree = 1000)
  } else if (name == "SVM_RFE") {
    fit <- train(xx, y_train, method = "svmRadial", trControl = control, metric = "Macro_AUC",
                 preProcess = preprocessing, tuneLength = 5)
  } else {
    grid <- expand.grid(nrounds = c(100, 250), max_depth = c(2, 4), eta = c(0.03, 0.1),
                        gamma = 0, colsample_bytree = 0.8, min_child_weight = c(1, 3), subsample = 0.8)
    fit <- train(xx, y_train, method = "xgbTree", trControl = control, metric = "Macro_AUC",
                 preProcess = preprocessing, tuneGrid = grid, verbose = 0, nthread = 1)
  }
  list(fit = fit, features = selected)
}
oof <- list(); index <- 1L
for (repeat_id in seq_len(n_repeats)) {
  set.seed(1000 + repeat_id); folds <- createFolds(y, k = outer_k, returnTrain = FALSE)
  for (fold_id in seq_along(folds)) {
    test <- folds[[fold_id]]; training <- setdiff(seq_len(nrow(x)), test)
    for (model_name in models) {
      object <- fit_model(model_name, x[training, , drop = FALSE], y[training], repeat_id * 10000 + fold_id)
      probability <- predict(object$fit, x[test, object$features, drop = FALSE], type = "prob")[, classes]
      oof[[index]] <- data.frame(Model = model_name, Repeat = repeat_id, Fold = fold_id,
                                 Sample = rownames(x)[test], Truth = y[test], probability, check.names = FALSE)
      index <- index + 1L
    }
  }
}
oof <- bind_rows(oof)
patient <- oof %>% group_by(Model, Sample, Truth) %>% summarise(across(all_of(classes), mean), .groups = "drop")
patient$Predicted <- classes[max.col(as.matrix(patient[, classes]), ties.method = "first")]
performance <- bind_rows(lapply(models, function(name) {
  d <- patient[patient$Model == name, ]
  data.frame(Model = name, t(metrics(d$Truth, d$Predicted, as.matrix(d[, classes]))), check.names = FALSE)
})) %>% arrange(desc(Macro_AUC), desc(MCC), desc(Macro_F1))
final_enet <- fit_model("ENET", x, y, 90001, final = TRUE)
final_enet$fit$trainingData <- NULL; final_enet$fit$pred <- NULL
final_enet$fit$resample <- NULL; final_enet$fit$resampledCM <- NULL
public_model <- list(Model = "ENET", Fit = final_enet$fit, Features = final_enet$features,
                     Classes = classes, CV = "5-fold outer, 4-fold inner, 20 repeats")
saveRDS(public_model, file.path(out_dir, "ENET_sanitized.rds"), version = 3)
fwrite(oof, file.path(out_dir, "nested_cv_predictions.csv"))
fwrite(patient, file.path(out_dir, "patient_level_predictions.csv"))
fwrite(performance, file.path(out_dir, "model_performance.csv"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "session_info.txt"))
