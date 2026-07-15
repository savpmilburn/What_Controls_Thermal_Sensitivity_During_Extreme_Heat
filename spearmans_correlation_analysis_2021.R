# spearmans_correlation_analysis_2021.R
rm(list = ls())
while (!is.null(dev.list())) dev.off()

library(tidyverse)
library(Hmisc)
library(car)
library(PerformanceAnalytics)
library(gridExtra)
library(psych)
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
  write_csv(coeffs_df, file.path(output_dir, "exports", paste0(matrix_name, "_coeffs.csv")))
  
  # P-values: same pattern
  pvals_df <- as.data.frame(corr_matrix$P)
  pvals_df$Variable <- rownames(pvals_df)
  pvals_df <- pvals_df[, c("Variable", setdiff(names(pvals_df), "Variable"))]
  write_csv(pvals_df, file.path(output_dir, "exports", paste0(matrix_name, "_pvalues.csv")))
  
  return(corr_matrix)
}


site_specific_corr_matrix <- compute_and_save_corr_matrix(data = data_2021, response_vars = ts_rvs_2021, covariate_vars = covariates_site_specific_2021, matrix_name = "site_specific")
upstream_corr_matrix <- compute_and_save_corr_matrix(data = data_2021, response_vars = ts_rvs_2021, covariate_vars = covariates_upstream_2021, matrix_name = "upstream")
reach_corr_matrix <- compute_and_save_corr_matrix(data = data_2021, response_vars = ts_rvs_2021, covariate_vars = covariates_reach_2021, matrix_name = "reach")
buffer_corr_matrix <- compute_and_save_corr_matrix(data = data_2021, response_vars = ts_rvs_2021, covariate_vars = covariates_buffer_2021, matrix_name = "buffer")
all_corr_matrix <- compute_and_save_corr_matrix(data = data_2021, response_vars = ts_rvs_2021, covariate_vars = all_covariates_2021, matrix_name = "all")