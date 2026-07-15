# spearmans_correlation_analysis_2021.R
rm(list = ls())
while (!is.null(dev.list())) dev.off()

library(tidyverse)
library(Hmisc)
library(PerformanceAnalytics)
library(gridExtra)
library(reshape2)
# Load 2021 Clackamas thermal sensitivity metrics and covariates data
data_2021 <- read.csv("clackamas_thermal_sensitivity_covariates_2021.csv", stringsAsFactors = FALSE)
nrow(data_2021)  # should be 72

# Create response variable vectors
all_rvs_2021 <- c("thermal_sensitivity", "mean_stream_temp_c", "mean_air_temp_c", "min_stream_temp_c", "max_stream_temp_c", "min_air_temp_c", "max_air_temp_c", "range_stream_temp_c", "range_air_temp_c")
ts_rvs_2021 <- c("thermal_sensitivity", "mean_stream_temp_c", "mean_air_temp_c")

# Create covariate vectors divided by spatial scale: site-specific, upstream, reach, buffer
covariates_site_specific_2021 <- c("channel_slope", "solar_exposure", "elevation_m", "base_flow_index", "summer_mean_max_air_temp_c", "summer_max_air_temp_c", "annual_precip_mm", "summer_precip_mm", "summer_mean_air_temp_c", "wet_season_precip_mm")
covariates_upstream_2021 <- c("developed_upstream_pct", "lakes_upstream_pct", "agricultural_upstream_pct", "burned_upstream_pct", "road_density_upstream", "high_cascades_upstream_pct", "wetlands_upstream_pct", "veg_cover_upstream_pct", "veg_height_upstream_m", "forest_upstream_pct", "shrub_upstream_pct", "upstream_area_km2")
covariates_reach_2021 <- c("burned_reach_pct", "agricultural_reach_pct", "wetlands_reach_pct", "lakes_reach_pct", "high_cascades_reach_pct", "developed_reach_pct", "road_density_reach", "veg_cover_reach_pct", "veg_height_reach_m")
covariates_buffer_2021 <- c("developed_buffer_pct", "agricultural_buffer_pct", "burned_buffer_pct", "wetlands_buffer_pct", "lakes_buffer_pct", "high_cascades_buffer_pct", "road_density_buffer", "veg_height_buffer_m", "veg_cover_buffer_pct")
all_covariates_2021 <- c(covariates_site_specific_2021, covariates_upstream_2021, covariates_reach_2021, covariates_buffer_2021)

length(all_covariates_2021)  # should be 40


compute_and_save_corr_matrix <- function(data, response_vars, covariate_vars, matrix_name, output_dir = "results_2021/correlation") {
  # Subset data: response variables + this scale's covariates
  subset_data <- data[, c(response_vars, covariate_vars)]
  
  # Compute Spearman's correlation matrix
  corr_matrix <- rcorr(as.matrix(subset_data), type = "spearman")
  
  # Save RDS
  saveRDS(corr_matrix, file.path(output_dir, "RDS", paste0(matrix_name, "_corr_matrix.RDS")))
  
  # Coefficients: move rownames into a column, write CSV
  coeffs_df <- as.data.frame(corr_matrix$r)
  coeffs_df$Variable <- rownames(coeffs_df)
  coeffs_df <- coeffs_df[, c("Variable", setdiff(names(coeffs_df), "Variable"))]
  
  # Round only export copy
  coeffs_export <- coeffs_df
  numeric_cols <- setdiff(names(coeffs_export), "Variable")
  coeffs_export[numeric_cols] <- round(coeffs_export[numeric_cols], 3)
  write_csv(coeffs_export, file.path(output_dir, "exports", paste0(matrix_name, "_coeffs.csv")))
  
  # P-values: same pattern
  pvals_df <- as.data.frame(corr_matrix$P)
  pvals_df$Variable <- rownames(pvals_df)
  pvals_df <- pvals_df[, c("Variable", setdiff(names(pvals_df), "Variable"))]
  
  # Round only export copy
  pvals_export <- pvals_df
  numeric_cols_p <- setdiff(names(pvals_export), "Variable")
  pvals_export[numeric_cols_p] <- lapply(pvals_export[numeric_cols_p], function(col) {
    ifelse(col < 0.001, "<0.001", as.character(round(col, 3)))
  })
  write_csv(pvals_df, file.path(output_dir, "exports", paste0(matrix_name, "_pvalues.csv")))
  
  return(corr_matrix)
}


site_specific_corr_matrix <- compute_and_save_corr_matrix(data = data_2021, response_vars = ts_rvs_2021, covariate_vars = covariates_site_specific_2021, matrix_name = "site_specific")
upstream_corr_matrix <- compute_and_save_corr_matrix(data = data_2021, response_vars = ts_rvs_2021, covariate_vars = covariates_upstream_2021, matrix_name = "upstream")
reach_corr_matrix <- compute_and_save_corr_matrix(data = data_2021, response_vars = ts_rvs_2021, covariate_vars = covariates_reach_2021, matrix_name = "reach")
buffer_corr_matrix <- compute_and_save_corr_matrix(data = data_2021, response_vars = ts_rvs_2021, covariate_vars = covariates_buffer_2021, matrix_name = "buffer")
all_corr_matrix <- compute_and_save_corr_matrix(data = data_2021, response_vars = ts_rvs_2021, covariate_vars = all_covariates_2021, matrix_name = "all")

analyze_corr_matrix <- function(corr_matrix, output_file) {
    # Extract coefficients & p values from correlation matrix
    coefficients <- corr_matrix$r 
    pValues <- corr_matrix$P
    # Reshape data to tall, normal structure
    coefficientsLong <- melt(coefficients, varnames = c("Variable 1", "Variable 2"), value.name = "Correlation Coefficient")
    pValuesLong <- melt(pValues, varnames = c("Variable 1", "Variable 2"), value.name = "p Value")
    # Combine correlation coefficient & p value
    table <- merge(coefficientsLong, pValuesLong, by = c("Variable 1", "Variable 2")) 
    # Remove diagonal (self-correlation) 
    table <- table[table$`Variable 1` != table$`Variable 2`, ]
    # Remove duplicate - not a pair anymore
    table <- table[!duplicated(t(apply(table[,1:2], 1, sort))), ] 
    # Find absolute value of correlation coefficient
    table$`Absolute Value of Correlation` <- abs(table$`Correlation Coefficient`)
    # Categorize correlation strength
    table$`Correlation Strength` <- cut(table$`Absolute Value of Correlation`, breaks = c(0, 0.19, 0.39, 0.59, 0.79, 1.00), labels = c("Very Weak", "Weak", "Moderate", "Strong", "Very Strong"), include.lowest = TRUE)
    # Find direction of correlation
    table$Direction <- ifelse(table$`Correlation Coefficient` > 0, "Positive", "Negative")
    # Sort by highest significance & highest strength
    table <- table[order(table$`p Value`, -table$`Absolute Value of Correlation`), ]
    # Round only export copy
    table_export <- table
    table_export$`Correlation Coefficient` <- round(table_export$`Correlation Coefficient`, 3)
    table_export$`Absolute Value of Correlation` <- round(table_export$`Absolute Value of Correlation`, 3)
    table_export$`p Value` <- ifelse(table_export$`p Value` < 0.001, "<0.001", as.character(round(table_export$`p Value`, 3)))
    # Save table
    write_csv(table_export, output_file)
    return(table)
}

summary_tables <- list(
  site_specific = analyze_corr_matrix(site_specific_corr_matrix, "results_2021/correlation/exports/site_specific_summary.csv"),
  upstream      = analyze_corr_matrix(upstream_corr_matrix, "results_2021/correlation/exports/upstream_summary.csv"),
  reach         = analyze_corr_matrix(reach_corr_matrix, "results_2021/correlation/exports/reach_summary.csv"),
  buffer        = analyze_corr_matrix(buffer_corr_matrix, "results_2021/correlation/exports/buffer_summary.csv"),
  all           = analyze_corr_matrix(all_corr_matrix, "results_2021/correlation/exports/all_summary.csv")
)

ts_correlations <- summary_tables$all[summary_tables$all$`Variable 1` == "thermal_sensitivity" | summary_tables$all$`Variable 2` == "thermal_sensitivity",]

ts_correlations_export <- ts_correlations
ts_correlations_export$`Correlation Coefficient` <- round(ts_correlations_export$`Correlation Coefficient`, 3)
ts_correlations_export$`Absolute Value of Correlation` <- round(ts_correlations_export$`Absolute Value of Correlation`, 3)
ts_correlations_export$`p Value` <- ifelse(ts_correlations_export$`p Value` < 0.001, "<0.001", as.character(round(ts_correlations_export$`p Value`, 3)))

write_csv(ts_correlations_export, "results_2021/correlation/exports/ts_correlations.csv")