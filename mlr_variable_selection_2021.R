# mlr_variable_selection_2021.R
library(readr)
library(lm.beta)
library(car)
library(leaps)

create_model_summary <- function(model, beta_model, corr_matrix, model_name) {
    model_summary <- summary(model)
    f_P_Value <- pf(model_summary$fstatistic[1], model_summary$fstatistic[2], model_summary$fstatistic[3], lower.tail = FALSE)
    std_coeffs <- beta_model$standardized.coefficients
    std_coeffs <- std_coeffs[!is.na(std_coeffs)]
    coeff_table <- model_summary$coefficients[-1, ]
    number_significant_covariates <- sum(coeff_table[, "Pr(>|t|)"] < 0.05)
    covariate_names <- names(std_coeffs)
    ts_corr <- corr_matrix$r["thermal_sensitivity", covariate_names]   # <- renamed from "thermalSensitivity"
    covariate_summary_df <- data.frame(
        `Model Name` = rep(model_name, length(covariate_names)), 
        `Landscape Covariates` = covariate_names, 
        Estimate = coeff_table[covariate_names, "Estimate"], 
        `Standardized Coefficient` = std_coeffs,
        `Std Error` = coeff_table[covariate_names, "Std. Error"], 
        `t Value` = coeff_table[covariate_names, "t value"], 
        `p Value` = coeff_table[covariate_names, "Pr(>|t|)"],
        Significance = ifelse(coeff_table[covariate_names, "Pr(>|t|)"] < 0.001, "***",
                             ifelse(coeff_table[covariate_names, "Pr(>|t|)"] < 0.01, "**",
                                   ifelse(coeff_table[covariate_names, "Pr(>|t|)"] < 0.05, "*",
                                         ifelse(coeff_table[covariate_names, "Pr(>|t|)"] < 0.1, ".", "")))),
        `Correlation w/ thermal sensitivity` = ts_corr, 
        `R Squared` = rep(model_summary$r.squared, length(covariate_names)), 
        `Adjusted R Squared` = rep(model_summary$adj.r.squared, length(covariate_names)),
        `F statistic` = rep(model_summary$fstatistic[1], length(covariate_names)),  
        `F-p Value` = rep(f_P_Value, length(covariate_names)), 
        `Number of Significant Covariates` = rep(number_significant_covariates, length(covariate_names)), 
        `Total Number of Covariates` = rep(length(covariate_names), length(covariate_names)), 
        stringsAsFactors = FALSE
    )
    return(covariate_summary_df)
}

find_correlations <- function(corr_matrix, corr_threshold, input_covariates, VIFs, data) {
    over_corr_threshold_pairs <- which(abs(corr_matrix) >= corr_threshold & corr_matrix != 1, arr.ind = TRUE)
    if(nrow(over_corr_threshold_pairs) > 0) {
        no_duplicates <- over_corr_threshold_pairs[over_corr_threshold_pairs[,1] < over_corr_threshold_pairs[,2], , drop = FALSE]
        removed_covariates <- c()
        cat(sprintf("Found %d covariate pairs with |correlation| >= %.1f\n", nrow(no_duplicates), corr_threshold))
        for(i in 1:nrow(no_duplicates)) {
            covariate_1 <- rownames(corr_matrix)[no_duplicates[i,1]]
            covariate_2 <- colnames(corr_matrix)[no_duplicates[i,2]]
            if(covariate_1 %in% removed_covariates || covariate_2 %in% removed_covariates) next
            cat(sprintf("\nHigh correlation pair: %s & %s (r = %.3f)\n", covariate_1, covariate_2, corr_matrix[covariate_1, covariate_2]))
            covariate_1_VIF <- VIFs[covariate_1]
            covariate_2_VIF <- VIFs[covariate_2]
            covariate_1_adj_R_squared <- summary(lm(thermal_sensitivity ~ get(covariate_1), data = data))$adj.r.squared   # <- renamed
            covariate_2_adj_R_squared <- summary(lm(thermal_sensitivity ~ get(covariate_2), data = data))$adj.r.squared   # <- renamed
            cat(sprintf("  %s: VIF = %.2f, Adj R² = %.3f\n", covariate_1, covariate_1_VIF, covariate_1_adj_R_squared))
            cat(sprintf("  %s: VIF = %.2f, Adj R² = %.3f\n", covariate_2, covariate_2_VIF, covariate_2_adj_R_squared))
            VIF_diff <- abs(covariate_1_VIF - covariate_2_VIF)
            if(VIF_diff > 2) {
                removed_covariate <- ifelse(covariate_1_VIF > covariate_2_VIF, covariate_1, covariate_2)
            } else {
                removed_covariate <- ifelse(covariate_1_adj_R_squared < covariate_2_adj_R_squared, covariate_1, covariate_2)
            }
            removed_covariates <- c(removed_covariates, removed_covariate)
            cat(sprintf("Removing %s\n", removed_covariate))
        }
        final_covariates <- input_covariates[!input_covariates %in% removed_covariates]
        cat(sprintf("\nRemoved %d variables due to high correlations: %s\n", length(removed_covariates), paste(removed_covariates, collapse = ", ")))
    } else {
        final_covariates <- input_covariates
        cat("No high correlations found\n")
    }
    return(list(
        final_covariates = final_covariates,
        removed_covariates = if(exists("removed_covariates")) removed_covariates else c(),
        numRemoved_covariates = if(exists("removed_covariates")) length(removed_covariates) else 0
    ))
}

# Run MLR variable selection for 2021 data (with solar_exposure removed)
# Load cleaned 2021 data and correlation martix 
data_2021 <- read.csv("clackamas_thermal_sensitivity_covariates_2021.csv", stringsAsFactors = FALSE)
nrow(data_2021)  # should be 72
corr_matrix <- readRDS("results_2021/correlation/RDS/all_corr_matrix.RDS")


# Redefine covariates
covariates_site_specific_2021 <- c("channel_slope", "solar_exposure", "elevation_m", "base_flow_index", "summer_mean_max_air_temp_c", "summer_max_air_temp_c", "annual_precip_mm", "summer_precip_mm", "summer_mean_air_temp_c", "wet_season_precip_mm")
covariates_upstream_2021 <- c("developed_upstream_pct", "lakes_upstream_pct", "agricultural_upstream_pct", "burned_upstream_pct", "road_density_upstream", "high_cascades_upstream_pct", "wetlands_upstream_pct", "veg_cover_upstream_pct", "veg_height_upstream_m", "forest_upstream_pct", "shrub_upstream_pct", "upstream_area_km2")
covariates_reach_2021 <- c("burned_reach_pct", "agricultural_reach_pct", "wetlands_reach_pct", "lakes_reach_pct", "high_cascades_reach_pct", "developed_reach_pct", "road_density_reach", "veg_cover_reach_pct", "veg_height_reach_m")
covariates_buffer_2021 <- c("developed_buffer_pct", "agricultural_buffer_pct", "burned_buffer_pct", "wetlands_buffer_pct", "lakes_buffer_pct", "high_cascades_buffer_pct", "road_density_buffer", "veg_height_buffer_m", "veg_cover_buffer_pct")
all_covariates_2021 <- c(covariates_site_specific_2021, covariates_upstream_2021, covariates_reach_2021, covariates_buffer_2021)
length(all_covariates_2021)  # should be 40

# MLR-specific exclusion: solar_exposure caused a direction mismatch between its 
# standardized beta and Spearman correlation sign across every VIF/correlation threshold tested (manuscript Section 2.3.3)
# Exclude solar_exposure from the MLR candidate pool
mlr_candidate_covariates_2021 <- setdiff(all_covariates_2021, "solar_exposure")
length(mlr_candidate_covariates_2021)  # should be 39

# Model on all 39 candidates, before any collinearity filtering
model_all_candidates <- lm(as.formula(paste("thermal_sensitivity ~", paste(mlr_candidate_covariates_2021, collapse = " + "))), data = data_2021)
beta_model_all_candidates <- lm.beta(model_all_candidates)

model_all_candidates_summary <- create_model_summary(model_all_candidates, beta_model_all_candidates, corr_matrix, "model_all_candidates")
write_csv(model_all_candidates_summary, "results_2021/mlr/exports/model_all_candidates_summary.csv")

VIFs <- vif(model_all_candidates)
cat("All candidate covariates with VIF values:\n")
print(VIFs)

VIFs_df <- data.frame(Covariate = names(VIFs), VIF = VIFs, row.names = NULL)
write_csv(VIFs_df, "results_2021/mlr/exports/model_all_candidates_vifs.csv")

covariates_high_VIF <- names(VIFs[VIFs > 15])
cat("\nCovariates with VIF > 15 (removed):\n")
print(covariates_high_VIF)

covariates_low_VIF <- setdiff(mlr_candidate_covariates_2021, covariates_high_VIF)
length(covariates_low_VIF)

# Model refit using only covariates with VIF <= 15
model_vif_filtered <- lm(as.formula(paste("thermal_sensitivity ~", paste(covariates_low_VIF, collapse = " + "))), data = data_2021)
beta_model_vif_filtered <- lm.beta(model_vif_filtered)

model_vif_filtered_summary <- create_model_summary(model_vif_filtered, beta_model_vif_filtered, corr_matrix, "model_vif_filtered")
write_csv(model_vif_filtered_summary, "results_2021/mlr/exports/model_vif_filtered_summary.csv")

#-------------------------------
# Subset the correlation matrix to just the 15 VIF-surviving covariates
covariates_low_VIF_corr_matrix <- corr_matrix$r[covariates_low_VIF, covariates_low_VIF]


# Run correlation filtering: within any pair |r| >= 0.6, remove the higher-VIF
# covariate, or (if VIFs are within 2 of each other) the one with lower individual
# adjusted R^2 against thermal sensitivity
correlation_filter_result <- find_correlations(corr_matrix = covariates_low_VIF_corr_matrix, corr_threshold = 0.6, input_covariates = covariates_low_VIF, VIFs = VIFs, data = data_2021)

covariates_corr_filtered <- correlation_filter_result$final_covariates
length(covariates_corr_filtered)  # compare against original model3EVs count (11)
print(covariates_corr_filtered)

# Model refit using only covariates with VIF <= 15 AND correlation < 0.6
model_corr_filtered <- lm(as.formula(paste("thermal_sensitivity ~", paste(covariates_corr_filtered, collapse = " + "))), data = data_2021)
beta_model_corr_filtered <- lm.beta(model_corr_filtered)

model_corr_filtered_summary <- create_model_summary(model_corr_filtered, beta_model_corr_filtered, corr_matrix, "model_corr_filtered")
write_csv(model_corr_filtered_summary, "results_2021/mlr/exports/model_corr_filtered_summary.csv")

#-------------------------------
# Exhaustive best-subset search across the 11 correlation-filtered covariates
best_subset <- regsubsets(thermal_sensitivity ~ ., data = data_2021[, c("thermal_sensitivity", covariates_corr_filtered)], nbest = 1, nvmax = length(covariates_corr_filtered), force.in = NULL, force.out = NULL, method = "exhaustive")

summary_best_subset <- summary(best_subset)



as.data.frame(summary_best_subset$outmat)

# Which model size (number of predictors) maximizes adjusted R^2
best_size <- which.max(summary_best_subset$adjr2)
cat("Recommended number of predictors:", best_size, "\n")

# Which covariates are in that best-sized model
best_covariates_logical <- summary_best_subset$which[best_size, ]
covariates_final_selected <- names(best_covariates_logical)[best_covariates_logical][-1]  # drop intercept
print(covariates_final_selected)

# Final model using the regsubsets-selected covariates
model_final_selected <- lm(as.formula(paste("thermal_sensitivity ~", paste(covariates_final_selected, collapse = " + "))), data = data_2021)
beta_model_final_selected <- lm.beta(model_final_selected)

model_final_selected_summary <- create_model_summary(model_final_selected, beta_model_final_selected, corr_matrix, "model_final_selected")
write_csv(model_final_selected_summary, "results_2021/mlr/exports/model_final_selected_summary.csv")

summary(beta_model_final_selected)
summary(beta_model_final_selected)$r.squared
summary(beta_model_final_selected)$adj.r.squared

# model_final_selected = final MLR model 
# Started with 39 candidate covariates (removed solar_exposure),
# 24 covariates removed due to VIF > 15 with 15 covariates surviving VIF filtering
# 4 covariates removed due to |r| >= 0.6 with 11 covariates surviving correlation filtering
# regsubsets selected 6 covariates for the final model
# Drop covariates that were not statistically significant (p >= 0.05) in the
# regsubsets-selected model — this is the actual published model
model_final_selected_coeffs <- summary(model_final_selected)$coefficients[-1, , drop = FALSE]
covariates_significant_only <- rownames(model_final_selected_coeffs)[model_final_selected_coeffs[, "Pr(>|t|)"] < 0.05]

model_published <- lm(as.formula(paste("thermal_sensitivity ~", paste(covariates_significant_only, collapse = " + "))), data = data_2021)
beta_model_published <- lm.beta(model_published)

model_published_summary <- create_model_summary(model_published, beta_model_published, corr_matrix, "model_published")
write_csv(model_published_summary, "results_2021/mlr/exports/model_published_summary.csv")
