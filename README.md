[Uploading README.md…]()
# A-Validated-FFPE-Compatible-nCounter-Panel-and-Model-for-Organoid-Confirmed-# GBM nCounter classifier

This repository contains four R scripts and a sanitized ENET model for developing and validating an nCounter-based GBM cluster classifier.

## Files

```text
scripts/
  01_CrossPlatform_QC.R
  02_Train_Compare_and_Freeze_Models.R
  03_Prospective_ENET_Validation.R
  04_Public_Cohort_ENET_Validation.R
models/
  ENET_sanitized.rds
```

## Requirements

Install R and the following packages:

```r
install.packages(c("caret", "data.table", "dplyr", "glmnet", "kernlab",
                   "pROC", "randomForest", "readxl", "xgboost"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("GSVA")
```

## Input format

- Expression matrices: genes in the first column and samples in subsequent columns.
- Label tables: columns named `Sample` and `Cluster`.
- Cluster labels: `C1`, `C2`, or `C3`.
- Sample and gene identifiers must already be identical across the corresponding input files.

## Configure file paths

Input paths are supplied through environment variables. Define the variables required by each script in a local `.Renviron` file or in the operating system environment.

Example:

```text
RETRO_BULK_FILE=<local-path>
RETRO_NCOUNTER_FILE=<local-path>
RETRO_LABEL_FILE=<local-path>
ENET_MODEL_FILE=<local-path>/ENET_sanitized.rds
OUTPUT_DIR=<local-output-directory>
```

The required variable names are listed at the beginning of each script. Do not upload the local `.Renviron` file.

## Run

Run the scripts in numerical order:

```bash
Rscript scripts/01_CrossPlatform_QC.R
Rscript scripts/02_Train_Compare_and_Freeze_Models.R
Rscript scripts/03_Prospective_ENET_Validation.R
Rscript scripts/04_Public_Cohort_ENET_Validation.R
```

Use `01` for cross-platform quality control, `02` for model comparison and ENET fitting, `03` for prospective validation, and `04` for pseudo-nCounter validation in public cohorts.

Only the sanitized model should be included in the public repository. Raw data, clinical data, unsanitized models, local paths, and output tables containing sample identifiers should remain local.
