# UNIVERSAL KRIGING FOR FUNCTIONAL PC's

################################################################################
#
#         Step 1: Load Required R Libraries, setting global parameters
#
################################################################################

# Spatial data handling
library(sf)         # For modern spatial vector data (stations)
library(sp)         # For legacy spatial classes used by gstat
library(terra)      # For reading and working with raster data
library(raster)     # Alternative raster package (used by some functions)

# Geostatistics and kriging
library(gstat)      # Core kriging and variogram modeling
library(automap)    # Automated kriging model selection (uses gstat underneath)

# Data wrangling
library(dplyr)      # Data manipulation
library(tidyr)      # Data reshaping
library(readr)      # Reading CSV and text data

# Cross-validation and evaluation
library(caret)      # Optional: for RMSE, MAE, etc.
library(Metrics)    # Another option for error metrics

# Visualization
library(ggplot2)    # For plotting
library(tmap)       # For thematic maps (raster + vector)
library(viridis)    # Nice color palettes

# Optional: Progress tracking during loops
library(pbapply)    # Adds progress bars to apply-like functions

# Optional: for cool country plots
library(rnaturalearth)
library(rnaturalearthdata)

# Global parameters
mainland_crs <- st_crs(32632)  # UTM zone 32N, WGS84




################################################################################
#
#           Step 2.1: Define the "mainland Italy" geometric object
#
################################################################################

# --- Load and transform country boundaries ---
italy_sf <- ne_countries(scale = "large", country = "Italy", returnclass = "sf") %>%
  st_transform(crs = mainland_crs)
san_marino_sf <- ne_countries(scale = "large", country = "San Marino", returnclass = "sf") %>%
  st_transform(crs = mainland_crs)
combined_sf <- rbind(italy_sf, san_marino_sf)

# --- Define the reference marker point ---
marker_point <- st_sf(geometry = st_sfc(st_point(c(10.1, 41.4)), crs = 4326)) %>%
  st_transform(crs = mainland_crs)

# --- Clip SW region below and to the left of the marker ---
clip_box <- st_polygon(list(rbind(
  c(-10, -10),                  # far SW
  c(10.1, -10),                 # west of marker
  c(10.1, 41.4),                # marker longitude & latitude
  c(-10, 41.4),                 # far west at marker latitude
  c(-10, -10)                   # close polygon
))) %>%
  st_sfc(crs = 4326) %>%
  st_transform(crs = mainland_crs)

mainland_clipped_sf <- st_difference(combined_sf, clip_box)

# --- Switch to interactive mode ---
tmap_mode("view")

# --- Plot interactive map with clipped mainland, stations, and marker ---
tm_shape(mainland_clipped_sf) +
  tm_polygons(col = "lightgrey", border.col = "black") +
  #tm_shape(stations_sf) +
  #tm_dots(col = "red", size = 0.5, shape = 21, border.col = "black", id = "fdsn_code_uni") +
  tm_shape(marker_point) +
  tm_dots(col = "blue", size = 0.5) +
  tm_layout(title = "Mainland Italy + San Marino (SW Removed)")


################################################################################
#
#           Step 2.2: Load the stations data
#
################################################################################

# Load station data
stations_df <- read_csv("ITACAs2s_2_flatfile/DATASET_WITH_SCORES.csv")

# Inspect structure
glimpse(stations_df)

# Convert to sf object (assuming lon/lat columns)
stations_sf <- st_as_sf(stations_df, coords = c("st_longitude", "st_latitude"), crs = 4326)

# Optional: transform to projected CRS (e.g., UTM32N for Italy)
stations_sf <- st_transform(stations_sf, crs = 32632)  # EPSG:32632 = UTM zone 32N

# Extract x and y coordinates to use them as covariates
stations_sf <- stations_sf %>%
  mutate(x = st_coordinates(geometry)[,1],
         y = st_coordinates(geometry)[,2])

# --- Plot mainland Italy + all stations (no marker) ---
tmap_mode("view")

tm_shape(mainland_clipped_sf) +
  tm_polygons(col = "lightgrey", border.col = "black") +
  tm_shape(stations_sf) +
  tm_dots(col = "red", size = 0.5, shape = 21, border.col = "black", id = "fdsn_code_uni") +
  tm_layout(title = "Mainland Italy + Seismic Stations")

################################################################################
#
#                  Step 2.3: Load raster predictors
#
################################################################################

# Load Mori Vs30 map
Vs30_Mori_path <- "INGV_data/Vs30_Mori/Vs30-Mori/VS30MORI50mUTM33.tif"
Vs30_Mori = raster(Vs30_Mori_path)



# Load Brunelli Vs30 map
Vs30_Brunelli_path <- "INGV_data/Vs30_Brunelli_new/vs30ms" 
vs30_Brunelli = raster(Vs30_Brunelli_path)

Vs30_1std_Brunelli_path <- "INGV_data/Vs30_Brunelli_new/vs30ms_1std" 
Vs30_1std_Brunelli = raster(Vs30_1std_Brunelli_path)

Vs30_m1std_Brunelli_path <- "INGV_data/Vs30_Brunelli_new/vs30ms_m1std" 
Vs30_m1std_Brunelli = raster(Vs30_m1std_Brunelli_path)



# Load Brunelli Vs30 map
slope_50_path <- "INGV_data/Slope_50" 
slope_50 = raster(slope_50_path)



# Load Brunelli Vs30 map
lito_Bucci_agreggata_path <- "INGV_data/Carla_lito_Bucci_aggregata/nuova_mappa_lito/new_lit_11_24_matlab.tif"
lito_Bucci_agreggata = raster(lito_Bucci_agreggata_path)





################################################################################
#
#                  Step 3 : computing Functional PC's
#
################################################################################

library(tidyverse)
library(fda)

# 1. Select the functional columns, including dS2S_SA_pga
functional_vars <- stations_sf %>%
  st_drop_geometry() %>%
  dplyr::select(dS2S_SA_pga, starts_with("dS2S_SA_T"))  # Include PGA as period = 0

# 2. Remove rows with NA
functional_data <- na.omit(functional_vars)

# 3. Define grid of argument values (e.g., periods T), including 0 for PGA
periods <- names(functional_vars) %>%
  str_replace("dS2S_SA_pga", "T0_000") %>%  # Treat PGA as T=0.000
  str_extract("T\\d+_\\d+|T\\d+") %>%
  str_replace("T", "") %>%
  str_replace("_", ".") %>%
  as.numeric()

# Order data by increasing period
order_idx <- order(periods)
periods <- periods[order_idx]
data_matrix <- as.matrix(functional_data)[,order_idx]

# 4. Create basis and functional data object
nbasis <- 15
rangeval <- range(periods)
basis_fd <- create.bspline.basis(rangeval = rangeval, nbasis = nbasis)

fd_obj <- Data2fd(argvals = periods, y = t(data_matrix), basisobj = basis_fd)

# 5. Perform Functional PCA
fpca_res <- pca.fd(fd_obj, nharm = 15, centerfns = TRUE)  # Increased nharm from 6 to 8

# 6. Add scores back to spatial dataframe
scores <- fpca_res$scores
colnames(scores) <- paste0("FPC", 1:ncol(scores))
stations_sf[paste0("FPC", 1:ncol(scores))] <- NA
stations_sf[rownames(functional_data), paste0("FPC", 1:ncol(scores))] <- scores

# 7. Scree plot (Improved)
var_explained <- fpca_res$values / sum(fpca_res$values)
cum_var_explained <- cumsum(var_explained)

cat("Percentage of variance explained by each PC:\n")
print(round(100 * var_explained[1:8], 2))
cat("Cumulative percentage of variance explained:\n")
print(round(100 * cum_var_explained[1:8], 2))

plot(var_explained[1:10], type = "b", pch = 16, lwd = 2,
     xlab = "Principal Component", ylab = "Proportion of Variance Explained",
     main = "Functional PC's for acceleration response spectra", xlim = c(1, 10), ylim = c(0, max(var_explained[1:10]) * 1.1))
grid()

# 8. Plot first 8 eigenfunctions (Improved layout and appearance)
layout(matrix(1:4, nrow = 2, ncol = 2, byrow = TRUE))
y_range <- range(sapply(1:4, function(i) eval.fd(periods, fpca_res$harmonics[i])))

for (i in 1:4) {
  plot(fpca_res$harmonics[i], xlab = "", ylab = "", main = paste("FPC", i),
       ylim = y_range, lwd = 2, col = "black")
  grid()
  if (i %% 2 == 1) mtext("Amplitude", side = 2, line = 2.5)
  if (i > 6) mtext("Period (s)", side = 1, line = 2.5)
}
layout(1)






# To reconstruct from FPCs :

# # Example: reconstruct from predicted FPC scores
# K <- length(predicted_scores)
# harmonics <- fpca_res$harmonics[1:K]
# mean_fd <- fpca_res$meanfd
# 
# # Evaluate the mean and harmonics at the original period points
# mean_vals <- eval.fd(periods, mean_fd)
# harmonic_vals <- sapply(1:K, function(k) eval.fd(periods, harmonics[k]))
# 
# # Reconstruct the function as mean + linear combination of eigenfunctions
# reconstructed_vals <- mean_vals + harmonic_vals %*% predicted_scores
# 
# # Optionally: convert to a named vector with period labels
# names(reconstructed_vals) <- paste0("dS2S_SA_T", str_replace_all(sprintf("%.3f", periods), "\\.", "_"))
# 
# # Result: reconstructed_vals is your estimated dS2S curve at each period







################################################################################
#
#                Step N : Training Kriging models
#
################################################################################


# ---- SETUP ----

# 1. Define covariates (adjust duplicates if needed)
covariates <- c("Vs30_Brunelli_simple", "slope50_simple", "Lito_Bucci_simple", "st_elevation", "slope", "x", "y")
left_part_of_formula = "sqrt(Vs30_Brunelli_simple) + sqrt(slope50_simple) + x + y + x*y"
fpc_names <- grep("^FPC\\d+$", names(stations_sf), value = TRUE)
dS2S_names <- names(stations_sf) %>%
  stringr::str_subset("^dS2S_SA_T\\d+_\\d+|^dS2S_SA_pga$")  # Includes PGA

# 2. Filter complete cases
kriging_data <- stations_sf %>%
  filter(if_all(all_of(c(fpc_names, covariates)), ~ !is.na(.))) %>%
  mutate(row_id = row_number())  # Track original rows

# 3. Initialize result dataframe with predicted + reconstructed columns
result <- kriging_data %>%
  st_drop_geometry() %>%
  dplyr::select(row_id, all_of(fpc_names), all_of(dS2S_names)) %>%
  mutate(
    across(all_of(fpc_names), ~NA_real_, .names = "{.col}_predicted"),
    across(all_of(dS2S_names), ~NA_real_, .names = "{.col}_predicted"),
    across(all_of(dS2S_names), ~NA_real_, .names = "{.col}_reconstructed")
  )


# ---- LOO-CV UNIVERSAL KRIGING ----

pb <- txtProgressBar(min = 0, max = nrow(kriging_data), style = 3) 

for (i in 1:nrow(kriging_data)) {
  train_data <- kriging_data[-i, ]
  test_point <- kriging_data[i, ]
  
  predicted_fpcs <- numeric(length(fpc_names))
  
  # Loop over each FPC to predict it via Kriging
  for (j in seq_along(fpc_names)) {
    fpc <- fpc_names[j]
    formula_str <- as.formula(paste(fpc, "~", left_part_of_formula))
    
    krig_model <- gstat::gstat(formula = formula_str, data = train_data)
    
    # Suppress messages and predict
    pred <- suppressMessages(
      suppressWarnings(
        capture.output(
          result_pred <- predict(krig_model, newdata = test_point),
          file = NULL
        )
      )
    )
    
    predicted_fpcs[j] <- result_pred[[1]]
    row_index <- kriging_data$row_id[i]
    result[[paste0(fpc, "_predicted")]][row_index] <- result_pred[[1]]
  }
  
  # ---- RECONSTRUCT dS2S FROM PREDICTED FPCs ----
  
  # Evaluate mean and harmonics
  K <- length(predicted_fpcs)
  mean_vals <- eval.fd(periods, fpca_res$meanfd)
  harmonic_vals <- sapply(1:K, function(k) eval.fd(periods, fpca_res$harmonics[k]))
  
  # Reconstruct functional values
  reconstructed_vals <- mean_vals + harmonic_vals %*% predicted_fpcs
  
  # Store reconstructed dS2S values
  # Store reconstructed values from predicted FPCs
  for (k in seq_along(dS2S_names)) {
    result[[paste0(dS2S_names[k], "_reconstructed")]][row_index] <- reconstructed_vals[k]
  }
  
  # ---- Direct Kriging for each dS2S variable ----
  for (k in seq_along(dS2S_names)) {
    dS2S_var <- dS2S_names[k]
    
    # Train data must exclude test station and contain no NA for target
    valid_train <- train_data %>%
      filter(!is.na(.data[[dS2S_var]]))
    
    if (nrow(valid_train) >= 5) {  # Ensure enough data for kriging
      formula_str <- as.formula(paste(dS2S_var, "~", left_part_of_formula))
      krig_model <- gstat::gstat(formula = formula_str, data = valid_train)
      
      pred_direct <- suppressMessages(
        suppressWarnings(
          capture.output(
            krig_result <- predict(krig_model, newdata = test_point),
            file = NULL
          )
        )
      )
      
      result[[paste0(dS2S_var, "_predicted")]][row_index] <- krig_result[[1]]
    }
  }
  setTxtProgressBar(pb, i)
}
close(pb)


################################################################################
#
#                Step N+1 : Evaluating Kriging models
#
################################################################################



# Helper function to compute metrics
compute_metrics <- function(observed, predicted) {
  errors <- predicted - observed
  data.frame(
    MAE = mean(abs(errors), na.rm = TRUE),
    MSE = mean(errors^2, na.rm = TRUE),
    RMSE = sqrt(mean(errors^2, na.rm = TRUE)),
    MBE = mean(errors, na.rm = TRUE)
  )
}

# Loop over dS2S names and compute both methods
metrics_list <- lapply(dS2S_names, function(name) {
  observed <- result[[name]]
  predicted_direct <- result[[paste0(name, "_predicted")]]
  predicted_reconstructed <- result[[paste0(name, "_reconstructed")]]
  
  metrics_direct <- compute_metrics(observed, predicted_direct)
  metrics_reconstructed <- compute_metrics(observed, predicted_reconstructed)
  
  metrics_direct$period_variable <- name
  metrics_direct$method <- "direct"
  
  metrics_reconstructed$period_variable <- name
  metrics_reconstructed$method <- "reconstructed"
  
  rbind(metrics_direct, metrics_reconstructed)
})

# Combine all into one data frame
evaluation_df <- bind_rows(metrics_list)

# Add human-readable labels
evaluation_df <- evaluation_df %>%
  mutate(period_label = str_replace(period_variable, "dS2S_SA_T", "T = "),
         period_label = str_replace(period_label, "dS2S_SA_pga", "PGA"),
         period_label = str_replace_all(period_label, "_", ".")) %>%
  relocate(method, .after = period_variable) %>%
  relocate(period_label, .after = period_variable)

# View results
print(evaluation_df)



## PLOT THE RESULTS


# Reorder periods numerically for plotting
evaluation_df <- evaluation_df %>%
  mutate(period_numeric = as.numeric(str_extract(period_label, "[\\d.]+")),
         period_numeric = if_else(is.na(period_numeric), 0, period_numeric))  # For "PGA"

# Reshape to long format for ggplot
evaluation_long <- evaluation_df %>%
  pivot_longer(cols = c(MAE, MSE, RMSE, MBE), names_to = "metric", values_to = "value")

# Plot all 4 metrics in separate facets
ggplot(evaluation_long, aes(x = period_numeric, y = value, color = method)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ metric, ncol = 1, scales = "free_y") +
  scale_x_log10(breaks = c(0.01, 0.1, 1, 10), labels = c("0.01", "0.1", "1", "10")) +
  labs(x = "Period (s)", y = "Metric Value", color = "Method",
       title = "Performance Metrics by Period and Method") +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold"))


################################################################################
#
#                Step N+1 : Plot residuals on a map
#
################################################################################


## ERRORS ON MAPS

# Define the period and metric you want to visualize
selected_period <- "dS2S_SA_T0_400"  # <- CHANGE THIS
residual_type <- "abs_error"         # options: "error", "abs_error"

# First, ensure both have `row_id` and drop overlapping columns from stations_sf
stations_with_residuals <- stations_sf %>%
  mutate(row_id = row_number()) %>%
  dplyr::select(-all_of(intersect(names(stations_sf), names(result)))) %>%  # Drop FPCs and dS2S columns that will be joined
  left_join(result, by = "row_id")



# Compute residuals
stations_with_residuals <- stations_with_residuals %>%
  mutate(
    error_direct = .data[[paste0(selected_period, "_predicted")]] - .data[[selected_period]],
    error_reconstructed = .data[[paste0(selected_period, "_reconstructed")]] - .data[[selected_period]],
    abs_error_direct = abs(error_direct),
    abs_error_reconstructed = abs(error_reconstructed)
  )

# Choose which residual type to plot
residual_direct <- paste0(residual_type, "_direct")
residual_reconstructed <- paste0(residual_type, "_reconstructed")

library(leaflet)
library(sf)
library(dplyr)
library(RColorBrewer)
library(leafsync)  # if you want synced maps

# -- Step 1: Transform to EPSG:4326 for leaflet
stations_geo <- stations_with_residuals %>%
  st_transform(crs = 4326)

# -- Step 2: Compute shared range of residuals
combined_residuals <- c(stations_geo$abs_error_direct, stations_geo$abs_error_reconstructed)
pal_shared <- colorBin("YlOrRd", domain = combined_residuals, bins = 7, na.color = "transparent")

# -- Step 3: Map for direct residuals
map_direct <- leaflet(stations_geo) %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addCircleMarkers(
    radius = 5,
    fillColor = ~pal_shared(abs_error_direct),
    stroke = FALSE,  # remove black outlines
    fillOpacity = 0.8,
    label = ~paste0("Direct Residual: ", round(abs_error_direct, 3))
  ) %>%
  addLegend("bottomright", pal = pal_shared, values = ~abs_error_direct,
            title = "Absolute Residual (Direct)", opacity = 1)

# -- Step 4: Map for reconstructed residuals
map_recon <- leaflet(stations_geo) %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addCircleMarkers(
    radius = 5,
    fillColor = ~pal_shared(abs_error_reconstructed),
    stroke = FALSE,  # remove black outlines
    fillOpacity = 0.8,
    label = ~paste0("Reconstructed Residual: ", round(abs_error_reconstructed, 3))
  ) %>%
  addLegend("bottomright", pal = pal_shared, values = ~abs_error_reconstructed,
            title = "Absolute Residual (Reconstructed)", opacity = 1)

# -- Step 5: Show maps side-by-side and synced
leafsync::sync(map_direct, map_recon)





################################################################################
#
#                Step N+1 : Plot reconstructed spectrum
#
################################################################################







# Ensure only one plot is shown
par(mfrow = c(1, 1))

# --- Parameters ---
i <- 34  # Index of the station to plot
station_id <- result$row_id[i]  # Ensure consistent indexing if needed

# --- Extract values ---
observed_vals      <- as.numeric(result[i, dS2S_names])
direct_vals        <- as.numeric(result[i, paste0(dS2S_names, "_predicted")])
reconstructed_vals <- as.numeric(result[i, paste0(dS2S_names, "_reconstructed")])

# --- Plot all curves ---
par(mfrow = c(1, 1))  # Ensure single panel

plot(periods, observed_vals, type = "b", pch = 16, col = "blue",
     xlab = "Period (s)", ylab = expression(delta[S2S]),
     main = paste("δS2S Comparison – Station", i),
     ylim = range(c(observed_vals, direct_vals, reconstructed_vals), na.rm = TRUE))

lines(periods, direct_vals, type = "b", pch = 17, col = "red")
lines(periods, reconstructed_vals, type = "b", pch = 18, col = "darkgreen")

legend("topright", legend = c("Observed", "Direct Kriging", "Reconstructed FPC"),
       col = c("blue", "red", "darkgreen"),
       pch = c(16, 17, 18), lwd = 1)

grid()






# Path to your folder (adjust to your actual path)
save_path <- "results/Pinar_model_selection.RData"

# Save both objects into one file
save(covariates, result, metrics_list, file = save_path)

# To load later:
# load("results/my_objects.RData")




