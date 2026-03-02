# Targeted-phospholipidomics-biomarker-pipeline
This repository demonstrates the analysis pipeline used to identify and validate a targeted plasma phospholipid biomarker panel for pneumoconiosis (PN) using:

- **Discovery cohort**: `pnde_infor.csv`
- **External test cohort**: `pnde_test_infor.csv`

## Pipeline overview

1. Environment setup (install/load packages; record `sessionInfo()`).
2. Data import and cleaning
   - Harmonize exposure duration for HC (set to 0 if missing).
   - Encode smoking status as an ordered factor.
   - Create a binary outcome `Group_Binary` (Case=PN, Control=DE/HC).
3. Split discovery cohort into 70% training / 30% internal validation.
4. Feature screening (limma) using a PN-specific contrast:
   - `PN_Specific = 2*PN - DE - HC`
   - Adjust for age, exposure duration, smoking.
5. Feature reduction (LASSO)
   - 10-fold CV, `lambda.min` selection.
   - Seed reset for reproducibility.
6. Cross-cohort overlap: keep features present in the external cohort.
7. Modeling using shared features + covariates
   - Firth logistic regression (`logistf`)
   - SVM (radial kernel; `e1071`)
   - Random forest (`randomForest`)
8. Evaluation
   - 10-fold CV on training set (mean ± SD metrics)
   - Internal validation set
   - External test set
9. Export metrics to `Integrated_Final_Results_Corrected_02.csv`

## How to run

```r
source("analysis_pipeline_EN.R")
```

Make sure the two CSV files are in your working directory (or edit the file paths in the script).

## Notes on standardization (batch/scale handling)

This workflow applies **separate BoxCox + centering + scaling** for training, internal validation, and external test sets to reduce scale shifts across cohorts. If you prefer a deployment-style transform (fit on training and apply to all), modify the `process_for_model()` function accordingly.
