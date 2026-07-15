# clean_data.R
# Read original SortedTSAndEVs2021.csv and produces renamed, 
# documented, reviewer-facing version. Renames columns and drops Max7DADM, GPS_long_, GPS_lat_.
library(dplyr)

# Load original data
old_data <- read.csv("SortedTSAndEVs2021.csv", stringsAsFactors = FALSE)
cat("Original columns (", ncol(old_data), "):\n", sep = "")
print(colnames(old_data))

# Define old column names to new column name mapping
name_map <- c(
    index = "index",
    site = "site", 
    longitude = "x", 
    latitude = "y", 
    stream_name = "Stream_Nam", 
    gridcode = "GRIDCODE",
    thermal_sensitivity = "thermalSensitivity", 
    ts_intercept = "intercept", 
    ts_r_squared = "rSquared", 
    ts_adj_r_squared = "adjRSquared",
    ts_rmse = "rmse",
    mean_stream_temp_c = "meanStreamTemp",
    mean_air_temp_c = "meanAirTemp",
    sd_stream_temp_c = "sdStreamTemp",
    sd_air_temp_c = "sdAirTemp",
    min_stream_temp_c = "minStreamTemp",
    max_stream_temp_c = "maxStreamTemp",
    min_air_temp_c = "minAirTemp",
    max_air_temp_c = "maxAirTemp",
    range_stream_temp_c = "rangeStreamTemp",
    range_air_temp_c = "rangeAirTemp",

    channel_slope = "SLOPE",
    solar_exposure = "Solar", 
    elevation_m = "Elev",
    base_flow_index = "BFI",

    developed_upstream_pct = "h2oDevelop", 
    lakes_upstream_pct = "h2oLakesPe",
    agricultural_upstream_pct = "h2oAgricul",
    burned_upstream_pct = "h2oBurnPer",
    road_density_upstream = "h2oRdDens",
    high_cascades_upstream_pct = "h2oHiCascP",
    wetlands_upstream_pct = "h2oWetland",
    veg_cover_upstream_pct = "h2oVegCov",
    veg_height_upstream_m = "h2oVegHt",
    forest_upstream_pct = "Forest21",
    shrub_upstream_pct = "Shrub21",
    upstream_area_km2 = "h2oKm2",

    burned_reach_pct = "BurnRCA",
    agricultural_reach_pct = "AgricultRC",
    wetlands_reach_pct = "WetlandsRC",
    lakes_reach_pct = "LakesRCA",
    high_cascades_reach_pct = "HiCascRCA",
    developed_reach_pct = "DevelopRCA",
    road_density_reach = "RoadsRCA",
    veg_cover_reach_pct = "VegCover",
    veg_height_reach_m = "VegHeight_",

    developed_buffer_pct = "DevelopBuf",
    agricultural_buffer_pct = "AgBuf",
    burned_buffer_pct = "BurnBuf",
    wetlands_buffer_pct = "WetlandBuf",
    lakes_buffer_pct = "LakesBuf",
    high_cascades_buffer_pct = "HiCascBuf",
    road_density_buffer = "RoadsBuf",
    veg_height_buffer_m = "VegHtBuf",
    veg_cover_buffer_pct = "VegCovBuf",

    summer_mean_max_air_temp_c = "MeanMaxAir",
    summer_max_air_temp_c = "MaxAir_C",
    annual_precip_mm = "Precip_mm",
    summer_precip_mm = "SumPrecip",
    summer_mean_air_temp_c = "MeanAirJJA",
    wet_season_precip_mm = "WetPrecip"
)

# Drop unused columns before renaming
unused_cols <- c("Max7DADM", "GPS_long_", "GPS_lat_")

missing_unused_cols <- setdiff(unused_cols, colnames(old_data))
if (length(missing_unused_cols) > 0) {
    warning("The following unused columns are not present in the original data and will be skipped: ", paste(missing_unused_cols, collapse = ", "))
}

data_trimmed <- old_data %>% select(-any_of(unused_cols))

# Verify before renaming
# Every old name in mapping must exist in trimmed data
old_names_in_map <- unname(name_map)
missing_from_data <- setdiff(old_names_in_map, colnames(data_trimmed))
if (length(missing_from_data) > 0) {
    stop("The following columns in the name mapping are not present in the trimmed data: ", paste(missing_from_data, collapse = ", "))
}

# Every column in trimmed data should be accounted for in mapping
unaccounted_for <- setdiff(colnames(data_trimmed), old_names_in_map)
if (length(unaccounted_for) > 0) {
    warning("The following columns in the trimmed data are not accounted for in the name mapping and will be kept with original names: ", paste(unaccounted_for, collapse = ", "))
}

cat("\nVerification passed: all mapped old-names found in data. \n")
if (length(unaccounted_for) > 0) {
    cat("Warning: some columns in data are not in mapping and will be kept with original names.\n")
}

# Apply renaming mappng 
clean_data <- data_trimmed %>% 
    select(all_of(old_names_in_map)) %>%
    rename(!!!name_map)

cat("\nRenamed columns (", ncol(clean_data), "):\n", sep = "")
print(colnames(clean_data))

stopifnot(nrow(clean_data) == nrow(old_data))
stopifnot(all(clean_data$thermal_sensitivity == old_data$thermalSensitivity))
stopifnot(all(clean_data$high_cascades_upstream_pct == old_data$h2oHiCascP))
stopifnot(all(clean_data$veg_height_upstream_m == old_data$h2oVegHt))

cat("\nAll checks passed. Writing cleaned data to SortedTSAndEVs2021_cleaned.csv\n")

write.csv(clean_data, "clackamas_thermal_sensitivity_covariates_2021.csv", row.names = FALSE)
saveRDS(clean_data, "clackamas_thermal_sensitivity_covariates_2021.RDS")