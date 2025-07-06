# -----------------------------------------------------------------------------
# R SCRIPT FOR ANALYZING NDVI TRENDS IN BURNED VS. UNBURNED AREAS
# -----------------------------------------------------------------------------
#
# OBJECTIVE:  
# This script analyzes the correlation between vegetation trends (from an NDVI
# trend raster) and the occurrence of forest fires (from shapefiles). It
# extracts trend values within fire perimeters, calculates summary statistics,
# compares them to unburned areas, and visualizes the results.
#
# REQUIRED FILES:
# 1. ARIZONA_TREND_full_23_24.tif: A raster where pixel values represent the
#    NDVI trend (e.g., slope of change). Negative values indicate decline.
# 2. Fire_2023_ARIZONA.shp: A shapefile of fire perimeters for 2023.
# 3. Fire_2024_ARIZONA.shp: A shapefile of fire perimeters for 2024.
#
# Author: Paul Arellano, PhD
# Date:   July, 04 2025
#
# -----------------------------------------------------------------------------

# --- 1. SETUP: Install and Load Libraries ---

# Install packages if they are not already installed
# (run these lines once if needed)
# install.packages("terra")
# install.packages("sf")
# install.packages("dplyr")
# install.packages("ggplot2")
# install.packages("tidyr")

# Load the libraries
library(terra)
library(sf)
library(dplyr)
library(ggplot2)
library(tidyr)

cat("✅ Libraries loaded successfully.\n")


# --- 2. CONFIGURATION: Define File Paths ---

# !!! IMPORTANT: UPDATE THESE PATHS TO MATCH YOUR FILE LOCATIONS !!!
path_ndvi_trend <- "C:/Users/pa589/NAU/TREE_STRESS/TreeStress_detection/MODELS/GEE_change_detection/ANOMALIES/Priority_Analysis/NDVI_TREND/ARIZONA_TREND_full_23_24.tif"
path_fire_2023  <- "C:/Users/pa589/NAU/TREE_STRESS/TreeStress_detection/MODELS/GEE_change_detection/ANOMALIES/Priority_Analysis/DATA/FIRES/Fires_2023_ARIZONA.shp"
path_fire_2024  <- "C:/Users/pa589/NAU/TREE_STRESS/TreeStress_detection/MODELS/GEE_change_detection/ANOMALIES/Priority_Analysis/DATA/FIRES/Fires_2024_ARIZONA.shp"
output_dir      <- "C:/Users/pa589/NAU/TREE_STRESS/TreeStress_detection/MODELS/GEE_change_detection/ANOMALIES/Priority_Analysis/Fire_analysis/OUTPUT" # Directory to save results

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


# --- 3. DATA LOADING AND PREPARATION ---

# Load the NDVI trend raster
cat("🔄 Loading NDVI trend raster...\n")
ndvi_trend_raster <- rast(path_ndvi_trend)

# Load the fire perimeter shapefiles
cat("🔄 Loading fire perimeter shapefiles...\n")
fires_2023 <- st_read(path_fire_2023)
fires_2024 <- st_read(path_fire_2024)

# Add a 'year' column to each dataset before merging
fires_2023$year <- 2023
fires_2024$year <- 2024

# Combine the two fire datasets into a single feature collection
all_fires <- rbind(fires_2023, fires_2024)

# --- CRS Check and Reprojection ---
# Ensure the fire polygons have the same Coordinate Reference System (CRS) as the raster
cat("🔄 Checking and aligning Coordinate Reference Systems (CRS)...\n")
if (st_crs(all_fires) != st_crs(ndvi_trend_raster)) {
  cat("CRS mismatch detected. Reprojecting fire polygons...\n")
  all_fires <- st_transform(all_fires, st_crs(ndvi_trend_raster))
}

cat("✅ Data loaded and prepared.\n")
print(ndvi_trend_raster)
print(all_fires)


# --- 4. ANALYSIS: Extract NDVI Trends and Calculate Statistics ---

# --- 4.1. Extract Trends within Burned Areas ---
cat("🔄 Extracting NDVI trend values from within fire perimeters...\n")
# This function returns a list, where each element is the vector of pixel values for a polygon
burned_pixel_values <- terra::extract(ndvi_trend_raster, vect(all_fires))

# Rename the second column (the actual trend values) to "ndvi_trend". The first is "ID".
names(burned_pixel_values)[2] <- "ndvi_trend"


# Combine the extracted values with the fire polygon attributes
# The ID column from the extraction links back to the row number of `all_fires`
burned_data <- burned_pixel_values %>%
  as_tibble() %>%
  # Add fire attributes (like year, name, etc.) using the ID
  mutate(fire_id = ID) %>%
  left_join(st_drop_geometry(all_fires) %>% mutate(fire_id = row_number()))

# --- 4.2. Extract Trends from Unburned Areas (for comparison) ---
cat("🔄 Sampling NDVI trend values from unburned areas for comparison...\n")
# Create a mask of the burned areas
burned_mask <- terra::rasterize(vect(all_fires), ndvi_trend_raster, field = 1)
unburned_raster <- terra::mask(ndvi_trend_raster, burned_mask, inverse = TRUE)

# Take a random sample of pixels from the unburned areas
# We sample a large number to get a representative distribution
unburned_pixel_values <- spatSample(unburned_raster, size = 50000, na.rm = TRUE, as.df = TRUE)
names(unburned_pixel_values) <- "ndvi_trend"

# --- 4.3. Combine Data for Plotting and Summary ---
burned_data$type <- "Burned"
unburned_pixel_values$type <- "Unburned"

# Combine all pixel values into one dataframe for easy plotting
combined_df <- bind_rows(
  burned_data %>% select(ndvi_trend, type),
  unburned_pixel_values %>% select(ndvi_trend, type)
) %>% drop_na()

cat("✅ Trend extraction complete.\n")


# --- 5. STATISTICAL SUMMARY ---
cat("🔄 Calculating summary statistics...\n")

# Define a threshold for what is considered a "significant" trend
# (adjust if your trend values are on a different scale)
trend_threshold <- 0.001

summary_stats <- combined_df %>%
  group_by(type) %>%
  summarise(
    pixel_count = n(),
    mean_trend = mean(ndvi_trend, na.rm = TRUE),
    median_trend = median(ndvi_trend, na.rm = TRUE),
    sd_trend = sd(ndvi_trend, na.rm = TRUE),
    min_trend = min(ndvi_trend, na.rm = TRUE),
    max_trend = max(ndvi_trend, na.rm = TRUE),
    neg_trend_pixels = sum(ndvi_trend < -trend_threshold, na.rm = TRUE),
    pos_trend_pixels = sum(ndvi_trend > trend_threshold, na.rm = TRUE),
    stable_pixels = sum(abs(ndvi_trend) <= trend_threshold, na.rm = TRUE)
  ) %>%
  mutate(
    percent_neg_trend = (neg_trend_pixels / pixel_count) * 100,
    percent_pos_trend = (pos_trend_pixels / pixel_count) * 100,
    percent_stable = (stable_pixels / pixel_count) * 100
  )

cat("--- Summary Statistics ---\n")
print(as.data.frame(summary_stats))

# Save summary statistics to a CSV file
output_csv_path <- file.path(output_dir, "ndvi_trend_fire_summary_stats.csv")
write.csv(summary_stats, output_csv_path, row.names = FALSE)
cat(paste0("✅ Summary statistics saved to: ", output_csv_path, "\n"))


# --- 6. VISUALIZATION ---
cat("🔄 Generating visualizations...\n")

# --- 6.1. Density Plot ---
# Compares the distribution of NDVI trend values
density_plot <- ggplot(combined_df, aes(x = ndvi_trend, fill = type)) +
  geom_density(alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  scale_fill_manual(values = c("Burned" = "#d95f02", "Unburned" = "#1b9e77")) +
  labs(
    title = "Distribution of NDVI Trends in Burned vs. Unburned Areas",
    x = "NDVI Trend (Slope)",
    y = "Density",
    fill = "Area Type"
  ) +
  theme_minimal()

print(density_plot)
ggsave(file.path(output_dir, "ndvi_trend_density_plot.png"), density_plot, width = 8, height = 6)

# --- 6.2. Box Plot ---
# Provides a clear comparison of the central tendency and spread
box_plot <- ggplot(combined_df, aes(x = type, y = ndvi_trend, fill = type)) +
  geom_boxplot(alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_fill_manual(values = c("Burned" = "#d95f02", "Unburned" = "#1b9e77")) +
  labs(
    title = "Comparison of NDVI Trends in Burned vs. Unburned Areas",
    x = "Area Type",
    y = "NDVI Trend (Slope)"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

print(box_plot)
ggsave(file.path(output_dir, "ndvi_trend_boxplot.png"), box_plot, width = 6, height = 6)

cat(paste0("✅ Visualizations saved to: ", output_dir, "\n"))


# --- 7. SUGGESTIONS FOR FURTHER ANALYSIS ---
cat("--- Suggestions for Further Analysis ---\n")
cat("1.  Burn Severity Correlation: If you have burn severity data (e.g., from MTBS or calculated dNBR), correlate the mean/median NDVI trend within each fire perimeter with the severity class. This can validate if a stronger negative NDVI trend corresponds to higher burn severity.\n\n")
cat("2.  Pre-Fire vs. Post-Fire Analysis: The current analysis likely shows the *impact* of the fire. To test if pre-existing vegetation stress made areas more susceptible to burning, you would need to calculate an NDVI trend using data *only from before* the 2023-2024 fire season and see if those trends differ inside vs. outside the fire perimeters.\n\n")
cat("3.  Topographic Influence: Incorporate a Digital Elevation Model (DEM). Analyze if post-fire NDVI trends are correlated with elevation, slope, or aspect. For example, do south-facing slopes (more sun exposure) show slower recovery (less positive trend) after a fire?\n\n")
cat("4.  Land Cover Analysis: Use a land cover map (like NLCD) to stratify your analysis. Do forests show a different trend response to fire compared to grasslands or shrublands?\n\n")
cat("5.  Time-Series Breakpoint Analysis: Instead of a simple linear trend, use methods like BFAST or LandTrendr on the full NDVI time-series for each fire polygon. This can pinpoint the exact date of the disturbance (the fire) and quantify the magnitude of the drop and the subsequent recovery trajectory.\n\n")

cat("--- Script Finished ---\n")
