# mlr_nested_cv_selection_2021.R
# Re-runs full MLR selection pipeline (VIF filter, correlation filter, regsubsets, significance trim)
# independently within each of 5 training folds instead of validating a single fixed final model

# Bring in published MLR 2021 model
source("mlr_variable_selection_2021.R") 

library(dplyr)
library(car)
library(leaps)
library(lm.beta)

# Same fold assignment logic: 5 folds with 14-15 sites per fold 
data_2021 <- data_2021 %>%
  arrange(index) %>%
  mutate(fold_number = paste0("Fold", ((index - 1) %% 5) + 1))

table(data_2021$fold_number)  # confirm ~14-15 sites per fold, 5 folds total


run_full_selection_pipeline <- function(train_data, candidate_covariates, vif_threshold = 15, corr_threshold = 0.6) {
  # Step 1: fit on all candidates, get VIF
  model_all <- lm(as.formula(paste("thermal_sensitivity ~", paste(candidate_covariates, collapse = " + "))), data = train_data)
  VIFs <- vif(model_all)
  
  # Step 2: drop VIF > threshold
  covariates_low_VIF <- names(VIFs[VIFs <= vif_threshold])
  
  # Step 3: correlation filtering, using a correlation matrix computed on THIS fold's training data only
  # (not the full-data corr_matrix — using the full-data matrix here would leak information from the held-out fold into the selection process, defeating the point of nesting)
  train_corr_matrix <- cor(train_data[, covariates_low_VIF], method = "spearman", use = "pairwise.complete.obs")
  
  corr_filter_result <- find_correlations(corr_matrix = train_corr_matrix, corr_threshold = corr_threshold, input_covariates = covariates_low_VIF, VIFs = VIFs, data = train_data)
  covariates_corr_filtered <- corr_filter_result$final_covariates
  
  # Step 4: regsubsets best-subset search
  best_subset <- regsubsets(thermal_sensitivity ~ ., data = train_data[, c("thermal_sensitivity", covariates_corr_filtered)], nvmax = length(covariates_corr_filtered), method = "exhaustive")
  summary_best_subset <- summary(best_subset)
  best_size <- which.max(summary_best_subset$adjr2)
  best_covariates_logical <- summary_best_subset$which[best_size, ]
  covariates_selected <- names(best_covariates_logical)[best_covariates_logical][-1]
  
  # Step 5: drop non-significant covariates, matching manuscript methodology
  model_selected <- lm(as.formula(paste("thermal_sensitivity ~", paste(covariates_selected, collapse = " + "))), data = train_data)
  coeffs <- summary(model_selected)$coefficients[-1, , drop = FALSE]
  covariates_final <- rownames(coeffs)[coeffs[, "Pr(>|t|)"] < 0.05]
  
  return(list(
    covariates_after_vif = covariates_low_VIF,
    covariates_after_corr = covariates_corr_filtered,
    covariates_regsubsets = covariates_selected,
    covariates_final_significant = covariates_final
  ))
}

print("MLR published model selection completed. No nested cross-validation selection procedure:")
# Full selection pipeline is run from MLR due to sourcing mlr_variable_selection_2021.R, which defines the find_correlations() function used in the correlation filtering step.
#------------------------------
print("Fold 1")
# --- Fold 1 ---
train_fold1 <- filter(data_2021, fold_number != "Fold1")
test_fold1 <- filter(data_2021, fold_number == "Fold1")
selection_fold1 <- run_full_selection_pipeline(train_fold1, mlr_candidate_covariates_2021)
final_covariates_fold1 <- selection_fold1$covariates_final_significant
formula_fold1 <- if (length(final_covariates_fold1) == 0) "thermal_sensitivity ~ 1" else paste("thermal_sensitivity ~", paste(final_covariates_fold1, collapse = " + "))
model_fold1 <- lm(as.formula(formula_fold1), data = train_fold1)
predictions_fold1 <- predict(model_fold1, newdata = test_fold1)
residuals_fold1 <- test_fold1$thermal_sensitivity - predictions_fold1
rmse_fold1 <- sqrt(mean(residuals_fold1^2))
r2_fold1 <- 1 - sum(residuals_fold1^2) / sum((test_fold1$thermal_sensitivity - mean(data_2021$thermal_sensitivity))^2)
print("Fold 1 complete.")

# --- Fold 2 ---
print("Fold 2")
train_fold2 <- filter(data_2021, fold_number != "Fold2")
test_fold2 <- filter(data_2021, fold_number == "Fold2")
selection_fold2 <- run_full_selection_pipeline(train_fold2, mlr_candidate_covariates_2021)
final_covariates_fold2 <- selection_fold2$covariates_final_significant
formula_fold2 <- if (length(final_covariates_fold2) == 0) "thermal_sensitivity ~ 1" else paste("thermal_sensitivity ~", paste(final_covariates_fold2, collapse = " + "))
model_fold2 <- lm(as.formula(formula_fold2), data = train_fold2)
predictions_fold2 <- predict(model_fold2, newdata = test_fold2)
residuals_fold2 <- test_fold2$thermal_sensitivity - predictions_fold2
rmse_fold2 <- sqrt(mean(residuals_fold2^2))
r2_fold2 <- 1 - sum(residuals_fold2^2) / sum((test_fold2$thermal_sensitivity - mean(data_2021$thermal_sensitivity))^2)
print("Fold 2 complete.")

# --- Fold 3 ---
print("Fold 3")
train_fold3 <- filter(data_2021, fold_number != "Fold3")
test_fold3 <- filter(data_2021, fold_number == "Fold3")
selection_fold3 <- run_full_selection_pipeline(train_fold3, mlr_candidate_covariates_2021)
final_covariates_fold3 <- selection_fold3$covariates_final_significant
formula_fold3 <- if (length(final_covariates_fold3) == 0) "thermal_sensitivity ~ 1" else paste("thermal_sensitivity ~", paste(final_covariates_fold3, collapse = " + "))
model_fold3 <- lm(as.formula(formula_fold3), data = train_fold3)
predictions_fold3 <- predict(model_fold3, newdata = test_fold3)
residuals_fold3 <- test_fold3$thermal_sensitivity - predictions_fold3
rmse_fold3 <- sqrt(mean(residuals_fold3^2))
r2_fold3 <- 1 - sum(residuals_fold3^2) / sum((test_fold3$thermal_sensitivity - mean(data_2021$thermal_sensitivity))^2)
print("Fold 3 complete.")

# --- Fold 4 ---
print("Fold 4")
train_fold4 <- filter(data_2021, fold_number != "Fold4")
test_fold4 <- filter(data_2021, fold_number == "Fold4")
selection_fold4 <- run_full_selection_pipeline(train_fold4, mlr_candidate_covariates_2021)
final_covariates_fold4 <- selection_fold4$covariates_final_significant
formula_fold4 <- if (length(final_covariates_fold4) == 0) "thermal_sensitivity ~ 1" else paste("thermal_sensitivity ~", paste(final_covariates_fold4, collapse = " + "))
model_fold4 <- lm(as.formula(formula_fold4), data = train_fold4)
predictions_fold4 <- predict(model_fold4, newdata = test_fold4)
residuals_fold4 <- test_fold4$thermal_sensitivity - predictions_fold4
rmse_fold4 <- sqrt(mean(residuals_fold4^2))
r2_fold4 <- 1 - sum(residuals_fold4^2) / sum((test_fold4$thermal_sensitivity - mean(data_2021$thermal_sensitivity))^2)
print("Fold 4 complete.")

# --- Fold 5 ---
print("Fold 5")
train_fold5 <- filter(data_2021, fold_number != "Fold5")
test_fold5 <- filter(data_2021, fold_number == "Fold5")
selection_fold5 <- run_full_selection_pipeline(train_fold5, mlr_candidate_covariates_2021)
final_covariates_fold5 <- selection_fold5$covariates_final_significant
formula_fold5 <- if (length(final_covariates_fold5) == 0) "thermal_sensitivity ~ 1" else paste("thermal_sensitivity ~", paste(final_covariates_fold5, collapse = " + "))
model_fold5 <- lm(as.formula(formula_fold5), data = train_fold5)
predictions_fold5 <- predict(model_fold5, newdata = test_fold5)
residuals_fold5 <- test_fold5$thermal_sensitivity - predictions_fold5
rmse_fold5 <- sqrt(mean(residuals_fold5^2))
r2_fold5 <- 1 - sum(residuals_fold5^2) / sum((test_fold5$thermal_sensitivity - mean(data_2021$thermal_sensitivity))^2)
print("Fold 5 complete.")

# Fold print results
cat("=== Fold 1 ===\n")
cat("After VIF:", length(selection_fold1$covariates_after_vif), "| After corr:", length(selection_fold1$covariates_after_corr), "| regsubsets:", length(selection_fold1$covariates_regsubsets), "\n")
cat("Final significant:", paste(final_covariates_fold1, collapse=", "), "\n")
cat("RMSE:", round(rmse_fold1, 4), "| R2:", round(r2_fold1, 4), "\n\n")

cat("=== Fold 2 ===\n")
cat("After VIF:", length(selection_fold2$covariates_after_vif), "| After corr:", length(selection_fold2$covariates_after_corr), "| regsubsets:", length(selection_fold2$covariates_regsubsets), "\n")
cat("Final significant:", paste(final_covariates_fold2, collapse=", "), "\n")
cat("RMSE:", round(rmse_fold2, 4), "| R2:", round(r2_fold2, 4), "\n\n")

cat("=== Fold 3 ===\n")
cat("After VIF:", length(selection_fold3$covariates_after_vif), "| After corr:", length(selection_fold3$covariates_after_corr), "| regsubsets:", length(selection_fold3$covariates_regsubsets), "\n")
cat("Final significant:", paste(final_covariates_fold3, collapse=", "), "\n")
cat("RMSE:", round(rmse_fold3, 4), "| R2:", round(r2_fold3, 4), "\n\n")

cat("=== Fold 4 ===\n")
cat("After VIF:", length(selection_fold4$covariates_after_vif), "| After corr:", length(selection_fold4$covariates_after_corr), "| regsubsets:", length(selection_fold4$covariates_regsubsets), "\n")
cat("Final significant:", paste(final_covariates_fold4, collapse=", "), "\n")
cat("RMSE:", round(rmse_fold4, 4), "| R2:", round(r2_fold4, 4), "\n\n")

cat("=== Fold 5 ===\n")
cat("After VIF:", length(selection_fold5$covariates_after_vif), "| After corr:", length(selection_fold5$covariates_after_corr), "| regsubsets:", length(selection_fold5$covariates_regsubsets), "\n")
cat("Final significant:", paste(final_covariates_fold5, collapse=", "), "\n")
cat("RMSE:", round(rmse_fold5, 4), "| R2:", round(r2_fold5, 4), "\n\n")

# Tables
# --- Nested selection results, one row per fold ---
nested_selection_per_fold_results <- data.frame(
  fold = c("Fold1","Fold2","Fold3","Fold4","Fold5"),
  n_after_vif = c(length(selection_fold1$covariates_after_vif), length(selection_fold2$covariates_after_vif), length(selection_fold3$covariates_after_vif), length(selection_fold4$covariates_after_vif), length(selection_fold5$covariates_after_vif)),
  n_after_corr = c(length(selection_fold1$covariates_after_corr), length(selection_fold2$covariates_after_corr), length(selection_fold3$covariates_after_corr), length(selection_fold4$covariates_after_corr), length(selection_fold5$covariates_after_corr)),
  n_regsubsets = c(length(selection_fold1$covariates_regsubsets), length(selection_fold2$covariates_regsubsets), length(selection_fold3$covariates_regsubsets), length(selection_fold4$covariates_regsubsets), length(selection_fold5$covariates_regsubsets)),
  n_significant = c(length(final_covariates_fold1), length(final_covariates_fold2), length(final_covariates_fold3), length(final_covariates_fold4), length(final_covariates_fold5)),
  covariates_selected = c(paste(final_covariates_fold1, collapse="; "), paste(final_covariates_fold2, collapse="; "), paste(final_covariates_fold3, collapse="; "), paste(final_covariates_fold4, collapse="; "), paste(final_covariates_fold5, collapse="; ")),
  held_out_rmse = c(rmse_fold1, rmse_fold2, rmse_fold3, rmse_fold4, rmse_fold5),
  held_out_r2 = c(r2_fold1, r2_fold2, r2_fold3, r2_fold4, r2_fold5)
)
print(nested_selection_per_fold_results)
write_csv(nested_selection_per_fold_results, "results_2021/mlr/exports/nested_selection_per_fold_results.csv")

# --- Covariate stability, exact name ---
fold_membership <- list(
  Fold1 = final_covariates_fold1,
  Fold2 = final_covariates_fold2,
  Fold3 = final_covariates_fold3,
  Fold4 = final_covariates_fold4,
  Fold5 = final_covariates_fold5
)

unique_covariates <- unique(unlist(fold_membership))

covariate_stability_by_exact_name <- data.frame(
  covariate = unique_covariates,
  n_folds_selected = sapply(unique_covariates, function(cov) sum(sapply(fold_membership, function(f) cov %in% f))),
  folds = sapply(unique_covariates, function(cov) paste(names(fold_membership)[sapply(fold_membership, function(f) cov %in% f)], collapse = ", ")),
  in_published_model = unique_covariates %in% covariates_significant_only
)
print(covariate_stability_by_exact_name)
write_csv(covariate_stability_by_exact_name, "results_2021/mlr/exports/covariate_stability_by_exact_name.csv")

# --- Covariate stability, family level (any scale) ---
covariate_family <- function(x) {
  gsub("_upstream_pct$|_reach_pct$|_buffer_pct$|_upstream$|_reach$|_buffer$", "", x)
}

published_families <- unique(sapply(covariates_significant_only, covariate_family))
all_selected_families <- unique(sapply(unlist(fold_membership), covariate_family))

covariate_stability_by_family <- data.frame(
  family = all_selected_families,
  n_folds_selected_any_scale = sapply(all_selected_families, function(fam) sum(sapply(fold_membership, function(f) any(sapply(f, covariate_family) == fam)))),
  in_published_model = all_selected_families %in% published_families
)
print(covariate_stability_by_family)
write_csv(covariate_stability_by_family, "results_2021/mlr/exports/covariate_stability_by_family.csv")

# --- Nested vs. fixed-model CV performance ---
nested_cv_ssr <- sum(residuals_fold1^2, residuals_fold2^2, residuals_fold3^2, residuals_fold4^2, residuals_fold5^2)
nested_cv_sst <- sum((data_2021$thermal_sensitivity - mean(data_2021$thermal_sensitivity))^2)
nested_cv_r2_pooled <- 1 - nested_cv_ssr / nested_cv_sst
nested_cv_rmse_mean <- mean(c(rmse_fold1, rmse_fold2, rmse_fold3, rmse_fold4, rmse_fold5))

nested_vs_fixed_cv_performance <- data.frame(
  Metric = c("Mean RMSE (avg of 5 folds)", "Pooled R-squared"),
  Fixed_Model_CV = c(0.0894, 0.593),
  Nested_CV = c(round(nested_cv_rmse_mean, 4), round(nested_cv_r2_pooled, 4))
)
print(nested_vs_fixed_cv_performance)
write_csv(nested_vs_fixed_cv_performance, "results_2021/mlr/exports/nested_vs_fixed_cv_performance.csv")