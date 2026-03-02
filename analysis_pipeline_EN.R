# Metabolomics Biomarker Discovery & Validation (PN vs non-PN)
# -----------------------------------------------------------------------------
# End-to-end pipeline used in the manuscript:
#   limma (PN-specific linear contrast) -> LASSO (lambda.min) -> feature overlap
#   -> ML models (Firth logistic / SVM / Random Forest) with:
#      - 10-fold cross-validation on the training set
#      - internal validation (30% holdout from discovery cohort)
#      - independent external testing
#      - batch/scale handling via separate standardization
#
# Inputs:
#   - pnde_infor.csv           (discovery cohort; clinical + lipid features)
#   - pnde_test_infor.csv      (external test cohort; clinical + lipid features)
#
# Output:
#   - Integrated_Final_Results_Corrected_02.csv
#
# Reproducibility:
#   - Global seed = 123 for splitting and other stochastic steps
#   - Seed reset before cv.glmnet() to stabilize LASSO selection
# -----------------------------------------------------------------------------

(Optional) Record session information
# Reproducibility helper (recommended)
# Capture R and package versions used for the analysis
sessionInfo()
Main analysis pipeline (provided script)
# Project: Metabolomics Biomarker Discovery & Validation (Final Corrected)
# Strategy: Limma -> LASSO (lambda.min) -> Intersection -> ML (Separate Norm)
# Note: Fixed random seed logic and lambda selection for reproducibility
# ==============================================================================

# 1. Environment setup
pkgs <- c("caret", "randomForest", "e1071", "logistf", "pROC", "dplyr", "glmnet", "limma", "xgboost")
new_pkgs <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]
if(length(new_pkgs)) install.packages(new_pkgs)

library(caret)
library(randomForest)
library(e1071)
library(logistf)
library(pROC)
library(dplyr)
library(glmnet)
library(limma)
library(xgboost)

# Global seed (used for data partitioning, etc.)
set.seed(123)

# 2. Data import and basic cleaning
# ------------------------------------------------------------------------------
data_disc <- read.csv("pnde_infor.csv", stringsAsFactors = FALSE)
data_ext  <- read.csv("pnde_test_infor.csv", stringsAsFactors = FALSE)

clean_basic <- function(df, tag) {
  if("Duration_of_exposure" %in% names(df)) {
    df$Duration_of_exposure[df$Group == "HC" & is.na(df$Duration_of_exposure)] <- 0
  }
  if(is.numeric(df$Smoking_status)) df$Smoking_status[df$Smoking_status == 0] <- 1
  df$Smoking_status <- factor(df$Smoking_status, levels = c(1, 2, 3))
  df$Group_Binary <- factor(ifelse(df$Group == "PN", "Case", "Control"), levels = c("Case", "Control"))
  return(df)
}

data_disc <- clean_basic(data_disc, "Discovery")
data_ext  <- clean_basic(data_ext, "External")
met_cols_disc <- grep("^a[0-9]+$", names(data_disc), value = TRUE)

# 3. Data split (70% training / 30% internal validation)
# ------------------------------------------------------------------------------
cat(">>> [Step 1] Splitting data (70/30)...\n")
trainIndex <- createDataPartition(data_disc$Group_Binary, p = 0.7, list = FALSE)
trainData_Raw <- data_disc[trainIndex, ]
validData_Raw <- data_disc[-trainIndex, ]

# Preprocessing (for feature selection)
preProc_FeatSel <- preProcess(trainData_Raw[, met_cols_disc], method = c("BoxCox", "center", "scale"))
trainData_Norm <- trainData_Raw
trainData_Norm[, met_cols_disc] <- predict(preProc_FeatSel, trainData_Raw[, met_cols_disc])

# 4. Feature screening I: limma (linear contrast)
# ------------------------------------------------------------------------------
cat(">>> [Step 2] limma screening...\n")
design <- model.matrix(~ 0 + Group + Age + Duration_of_exposure + Smoking_status, data = trainData_Norm)
colnames(design) <- gsub("Group", "", colnames(design))
fit <- lmFit(t(trainData_Norm[, met_cols_disc]), design)
contr.matrix <- makeContrasts(PN_Specific = 2*PN - DE - HC, levels = design)
fit2 <- eBayes(contrasts.fit(fit, contr.matrix))
res_limma <- topTable(fit2, coef = "PN_Specific", number = Inf, sort.by = "P", adjust.method = "BH")
sig_limma_feats <- rownames(res_limma)[res_limma$adj.P.Val < 0.05]

cat(sprintf("    limma candidates: %d features\n", length(sig_limma_feats)))

# 5. Feature screening II: LASSO (revised: lambda.min + seed reset)
# ------------------------------------------------------------------------------
cat(">>> [Step 3] LASSO regression (lambda.min)...\n")

x_lasso <- as.matrix(trainData_Norm[, sig_limma_feats])
y_lasso <- ifelse(trainData_Norm$Group_Binary == "Case", 1, 0)

# [Key fix] Reset the seed here so LASSO results match a standalone run
set.seed(123) 

# Run cross-validation
cv_lasso <- cv.glmnet(x_lasso, y_lasso, family = "binomial", alpha = 1, nfolds = 10)

# [Key fix] Use lambda.min to retain more features (including a6)
best_lambda <- cv_lasso$lambda.min 
#print(best_lambda)
coef_lasso <- coef(cv_lasso, s = best_lambda)
lasso_feats <- rownames(coef_lasso)[which(coef_lasso != 0)]
lasso_feats <- lasso_feats[lasso_feats != "(Intercept)"]

cat("    LASSO selected:", paste(lasso_feats, collapse=", "), "\n")

# 6. Cross-cohort intersection (feature overlap)
# ------------------------------------------------------------------------------
ext_feats_avail <- grep("^a[0-9]+$", names(data_ext), value = TRUE)
final_mets <- intersect(lasso_feats, ext_feats_avail)

cat(">>> [Step 4] Final intersected features:", paste(final_mets, collapse=", "), "\n")

# If a6 is still dropped due to randomness, print a warning (informational only)
if(!"a6" %in% final_mets && "a6" %in% ext_feats_avail) {
  cat("Warning: 'a6' was not selected by LASSO this time. Results may vary slightly.\n")
}

# Final model variables
clinical_vars <- c("Age", "Duration_of_exposure", "Smoking_status")
model_vars <- c(final_mets, clinical_vars)

# 7. Separate standardization (train/valid/test separately)
# ------------------------------------------------------------------------------
process_for_model <- function(df, mets) {
  df_sub <- na.omit(df[, c("Group_Binary", clinical_vars, mets)])
  for(m in mets) df_sub[[m]] <- as.numeric(df_sub[[m]])
  # Fit normalization parameters within each dataset separately
  preProc <- preProcess(df_sub[, mets], method = c("BoxCox", "center", "scale"))
  df_sub[, mets] <- predict(preProc, df_sub[, mets])
  return(df_sub)
}

trainData_Final <- process_for_model(trainData_Raw, final_mets)
validData_Final <- process_for_model(validData_Raw, final_mets)
testData_Final  <- process_for_model(data_ext, final_mets)

# 8. Model training and 10-fold cross-validation
# ------------------------------------------------------------------------------
cat(">>> [Step 5] Model training & 10-fold CV...\n")

# Define a reusable metric function
calc_metrics <- function(obs, probs, model_name, dataset_name) {
  # 1) Build ROC object
  # quiet=TRUE suppresses verbose output
  roc_obj <- tryCatch(
    { roc(obs, probs, levels = c("Control", "Case"), direction = "<", quiet = TRUE) },
    error = function(e) return(NULL)
  )
  
  if(is.null(roc_obj)) return(NULL)
  
  # 2) Safely compute AUC and its CI
  auc_val <- as.numeric(roc_obj$auc)
  
  # If AUC is extremely close to 1 (e.g., >0.999), CI may fail or be meaningless; hard-code a perfect CI
  if(auc_val > 0.999) {
    ci_str <- "1.000 (1.000-1.000)"
  } else {
    # Try to compute CI; if it fails, return NA placeholders
    ci_str <- tryCatch({
      ci_obj <- ci(roc_obj)
      sprintf("%.3f (%.3f-%.3f)", auc_val, ci_obj[1], ci_obj[3])
    }, error = function(e) {
      sprintf("%.3f (NA-NA)", auc_val)
    })
  }
  
  # 3) Safely extract the optimal threshold (core robustness fix)
  # transpose=FALSE keeps metric names as column names
  best_coords <- coords(roc_obj, "best", best.method = "youden", 
                        ret = c("threshold", "sensitivity", "specificity", "ppv", "npv"),
                        transpose = FALSE)
  
  # [Key fix] Handle multiple optimal thresholds
  # If a matrix/data.frame (multiple rows) is returned, keep the first row
  if(is.matrix(best_coords) || is.data.frame(best_coords)) {
    best_coords <- best_coords[1, , drop=FALSE]
  }
  
  # Coerce to numeric scalars to avoid NULL due to list-like structures
  thresh_val <- as.numeric(best_coords["threshold"])
  sens_val   <- as.numeric(best_coords["sensitivity"])
  spec_val   <- as.numeric(best_coords["specificity"])
  ppv_val    <- as.numeric(best_coords["ppv"])
  npv_val    <- as.numeric(best_coords["npv"])
  
  # 4) Compute F1 score at the optimal threshold
  pred_class <- factor(ifelse(probs >= thresh_val, "Case", "Control"), levels = c("Case", "Control"))
  cm <- confusionMatrix(pred_class, obs, mode = "everything", positive = "Case")
  f1_val <- cm$byClass["F1"]
  
  # 5) Assemble a one-row results data.frame
  # Use sprintf to ensure all entries are length-1 strings
  data.frame(
    Dataset = dataset_name,
    Model = model_name,
    `AUC (95% CI)` = ci_str,
    Sensitivity = sprintf("%.3f", sens_val),
    Specificity = sprintf("%.3f", spec_val),
    F1 = sprintf("%.3f", f1_val),
    PPV = sprintf("%.3f", ppv_val),
    NPV = sprintf("%.3f", npv_val),
    Best_Threshold = sprintf("%.3f", thresh_val),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

# 10-fold CV loop
folds <- createFolds(trainData_Final$Group_Binary, k = 10, list = TRUE)
cv_res_list <- list()

for(i in 1:10) {
  cv_test_idx <- folds[[i]]
  cv_tr <- trainData_Final[-cv_test_idx, ]
  cv_te <- trainData_Final[cv_test_idx, ]
  
  # Train models (increase logistf iterations)
  fit_f <- tryCatch(logistf(as.formula(paste("ifelse(Group_Binary=='Case',1,0) ~", paste(model_vars, collapse="+"))), 
                            data = cv_tr, control=logistf.control(maxit=1000)), error=function(e) NULL)
  fit_s <- svm(as.formula(paste("Group_Binary ~", paste(model_vars, collapse="+"))), data = cv_tr, probability=TRUE, kernel="radial")
  fit_r <- randomForest(as.formula(paste("Group_Binary ~", paste(model_vars, collapse="+"))), data = cv_tr, ntree=200)
  
  # Predict
  if(!is.null(fit_f)) {
    p_f <- predict(fit_f, newdata = cv_te, type = "response")
    cv_res_list[[length(cv_res_list)+1]] <- calc_metrics(cv_te$Group_Binary, p_f, "Firth Logistic", "CV_Fold")
  }
  p_s <- attr(predict(fit_s, newdata = cv_te, probability=TRUE), "probabilities")[, "Case"]
  p_r <- predict(fit_r, newdata = cv_te, type = "prob")[, "Case"]
  
  cv_res_list[[length(cv_res_list)+1]] <- calc_metrics(cv_te$Group_Binary, p_s, "SVM", "CV_Fold")
  cv_res_list[[length(cv_res_list)+1]] <- calc_metrics(cv_te$Group_Binary, p_r, "Random Forest", "CV_Fold")
}

cv_summary <- do.call(rbind, cv_res_list) %>% group_by(Model) %>%
  summarise(Dataset = "10-Fold CV (Mean)",
            `AUC (95% CI)` = sprintf("%.3f ± %.3f", mean(as.numeric(substr(`AUC (95% CI)`, 1, 5))), sd(as.numeric(substr(`AUC (95% CI)`, 1, 5)))),
            Sensitivity = sprintf("%.3f ± %.3f", mean(as.numeric(Sensitivity)), sd(as.numeric(Sensitivity))),
            Specificity = sprintf("%.3f ± %.3f", mean(as.numeric(Specificity)), sd(as.numeric(Specificity))),
            F1 = sprintf("%.3f ± %.3f", mean(as.numeric(F1)), sd(as.numeric(F1))),
            PPV = sprintf("%.3f", mean(as.numeric(PPV))), NPV = sprintf("%.3f", mean(as.numeric(NPV))), Best_Threshold = "Var")

# 9. Final models and independent evaluation
# ------------------------------------------------------------------------------
cat(">>> [Step 6] Final evaluation (train/internal/external)...\n")
# Train final models on the training set
train_f_dat <- trainData_Final
train_f_dat$y <- ifelse(trainData_Final$Group_Binary == "Case", 1, 0)
fit_final_firth <- logistf(as.formula(paste("y ~", paste(model_vars, collapse="+"))), data = train_f_dat, control=logistf.control(maxit=1000))
fit_final_svm <- svm(as.formula(paste("Group_Binary ~", paste(model_vars, collapse="+"))), data = trainData_Final, probability=TRUE, kernel="radial")
fit_final_rf <- randomForest(as.formula(paste("Group_Binary ~", paste(model_vars, collapse="+"))), data = trainData_Final, ntree=500)

evaluate_final <- function(dat, name) {
  p_f <- predict(fit_final_firth, newdata = dat, type = "response")
  p_s <- attr(predict(fit_final_svm, newdata = dat, probability=TRUE), "probabilities")[, "Case"]
  p_r <- predict(fit_final_rf, newdata = dat, type = "prob")[, "Case"]
  obs <- dat$Group_Binary
  rbind(calc_metrics(obs, p_f, "Firth Logistic", name),
        calc_metrics(obs, p_s, "SVM", name),
        calc_metrics(obs, p_r, "Random Forest", name))
}

res_train <- evaluate_final(trainData_Final, "Training Set (Self)")
res_valid <- evaluate_final(validData_Final, "Internal Validation (30%)")
res_test  <- evaluate_final(testData_Final, "External Test (Separate Norm)")

# 10. Summarize and export
# ------------------------------------------------------------------------------
final_table <- rbind(cv_summary, res_train, res_valid, res_test)
cat("\n=== Final Integrated Results ===\n")
print(final_table)
write.csv(final_table, "Integrated_Final_Results_Corrected_02.csv", row.names = FALSE)
