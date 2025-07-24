library(mvtnorm) # to deal with multivariate normal distributions
library(car) # "Companion to Applied Regression" for regression analysis
library(dplyr)
library( GGally)
library(maps)
library(leaflet) #to use the map
library(corrplot)

#Load data
df_itaca = read.csv("ITACAs2s_2_flatfile/ITACAs2s_SA_2_0.csv", sep=";")
col_names = names(df_itaca)
col_names
colcarnames = col_names[c(1,2,3,4,5,11,12,13,14,15,16,22,23,25)] #names of categorical variables
colcarnames

#NUMERICAL REGIONAL ANALYSIS
#Add region labels (Central, North, South)
df_itaca$region <- with(df_itaca, ifelse(
  # CENTRAL: between 41???44 lat and 10???16 lon
  st_latitude >= 41 & st_latitude <= 44 & st_longitude >= 10 & st_longitude <= 16,
  "Central",
  # NORTH: lat > 44 OR (41???44 lat but west of lon < 10)
  ifelse(st_latitude > 44 | 
           (st_latitude >= 41 & st_latitude <= 44 & st_longitude < 10), 
         "North", 
         # SOUTH: remaining
         "South")
))


#Map visualization-- all data
map <- leaflet(data = df_itaca) %>%
  addTiles() %>%
  setView(lng = 12.5, lat = 41.9, zoom = 6) %>%
  addCircleMarkers(
    lng = ~st_longitude, 
    lat = ~st_latitude,
    color = ~ifelse(region == "Central", "blue",
                    ifelse(region == "North", "green", "pink")),
    radius = 5,
    fill = TRUE,
    fillOpacity = 0.7,
    label = ~region  
  )
map
##
cols_to_use <- c("st_latitude", "st_longitude", "st_elevation", 
                 "slope", "depth_to_bedrock", 
                 "vs30_m_s_from_topography", "dS2S_SA_pga")

# Regional Correlation Matrices(ggpairs)
for (r in unique(df_itaca$region)) {
  cat("Plotting:", r, "\n")
  
  df_region <- df_itaca[df_itaca$region == r, ]
  df_clean <- df_region[, cols_to_use]
  df_clean <- na.omit(df_clean)
  
  if (nrow(df_clean) > 2) {
    p <- ggpairs(df_clean, title = paste("Region:", r, "- 7 Variable Correlation Matrix"))
    print(p)
  } else {
    cat("Not enough data:", r, "\n")
  }
}

#----------------------------------#
#         Raster Data Merge        #
#----------------------------------#
library(sp)
library(sf)
library(raster)
library(terra)

#Utility Functions

reproject_stations = function(stations, raster_map) {
  raster_CRS <- CRS(projection(raster_map))
  stations_new_CRS <- spTransform(stations, raster_CRS)
  return(stations_new_CRS)
}

extract_values = function(raster_map, stations, method, buffer = NULL) {
  values <- extract(raster_map, stations, method = method, buffer = buffer)
  return(values)
}

print_summary_stats = function(values, method_name) {
  cat(paste0("Summary statistics for method '", method_name, "':\n"))
  print(summary(values))
  cat("\n")
}
compare_extractions = function(values1, values2, method1_name, method2_name) {
  common_NA <- which(is.na(values1) & is.na(values2))
  cat(paste("How many NA values in common:", length(common_NA)), "\n")
  cat(paste("Indexes of common NA values:", toString(common_NA)), "\n\n")
  
  valid_indices <- which(!is.na(values1) & !is.na(values2))
  
  if (length(valid_indices) > 0) {
    errors <- values1[valid_indices] - values2[valid_indices]
    mean_abs_error <- mean(abs(errors))
    sd_error <- sd(errors)
    
    cat(paste0("Mean absolute error between '", method1_name, "' and '", method2_name, "': ", round(mean_abs_error, 4)), "\n")
    cat(paste0("Standard deviation of error: ", round(sd_error, 4)), "\n\n")
  } else {
    cat("No valid data points to compare.\n")
  }
}

extract_map = function(stations, path_to_map, method = c("simple", "bilinear"), buffer = NULL, check_NA = TRUE) {
  raster_map <- raster(path_to_map)
  stations_new_CRS <- reproject_stations(stations, raster_map)
  extracted_values <- list()
  
  for (m in method) {
    extracted_values[[m]] <- extract_values(raster_map, stations_new_CRS, method = m, buffer = buffer)
    if (check_NA) {
      cat(paste0("Number of NA values in '", m, "' extraction: ", sum(is.na(extracted_values[[m]]))), "\n")
    }
    print_summary_stats(extracted_values[[m]], m)
  }
  
  if (length(method) == 2) {
    compare_extractions(extracted_values[[method[1]]], extracted_values[[method[2]]], method[1], method[2])
  }
  
  return(extracted_values)
}

base_CRS <- CRS("+proj=longlat +datum=WGS84")
stations <- SpatialPoints(df_itaca[, c("st_longitude", "st_latitude")], proj4string = base_CRS)


#---- Mori
Vs30_Mori_path <- "INGV_data/Vs30_Mori/Vs30-Mori/VS30MORI50mUTM33.tif"
Mori_extract <- extract_map(stations, Vs30_Mori_path)
df_itaca$Vs30_MORI_simple   <- Mori_extract$simple
df_itaca$Vs30_MORI_bilinear <- Mori_extract$bilinear

#---- Brunelli
Vs30ms_Brunelli_path <- "INGV_data/Vs30_Brunelli_new/vs30ms" #exact predictions
Brunelli_Vs30ms_extract <- extract_map(stations, Vs30ms_Brunelli_path)
df_itaca$Vs30_Brunelli_simple   <- Brunelli_Vs30ms_extract$simple
df_itaca$Vs30_Brunelli_bilinear <- Brunelli_Vs30ms_extract$bilinear
#if we want to work with uncertainity 
vs30ms_1std_Brunelli_path <- "INGV_data/Vs30_Brunelli_new/vs30ms_1std" #Vs30 - 1 std??? optimistic
Brunelli_vs30ms_1std_extract <- extract_map(stations, vs30ms_1std_Brunelli_path)
df_itaca$Vs30_Brunelli_1std_simple   <- Brunelli_vs30ms_1std_extract$simple
df_itaca$Vs30_Brunelli_1std_bilinear <- Brunelli_vs30ms_1std_extract$bilinear

vs30ms_m1std_Brunelli_path <- "INGV_data/Vs30_Brunelli_new/vs30ms_m1std" #Vs30 - 1 std??? pesimistic
Brunelli_vs30ms_m1std_extract <- extract_map(stations, vs30ms_m1std_Brunelli_path)
df_itaca$Vs30_Brunelli_m1std_simple   <- Brunelli_vs30ms_m1std_extract$simple
df_itaca$Vs30_Brunelli_m1std_bilinear <- Brunelli_vs30ms_m1std_extract$bilinear

#Merged data frame
write.csv(df_itaca, "ITACAs2s_SA_2_0_with_vs30frommaps.csv", row.names = FALSE)

#Representative period: dS2S_SA_T0_100
# SA(1s) represents the long-period range where site amplification effects are most pronounced and corresponds to the resonance period of mid- to high-rise buildings.
# Therefore, as in Sgobba et al. (2024), SA(1s) is selected as the representative period in this analysis.
cols_pca <- c("st_latitude", "st_longitude", "st_elevation", 
              "slope", "depth_to_bedrock", 
              "Vs30_Brunelli_bilinear", 
              "dS2S_SA_T0_100")  # SA(1s)


for (r in unique(df_itaca$region)) {
  cat("\n===== PCA for Region:", r, "=====\n")
  
  df_region <- df_itaca[df_itaca$region == r, ]
  df_region_clean <- df_region[, cols_pca]
  df_region_clean <- df_region_clean[complete.cases(df_region_clean), ]
  
  if (nrow(df_region_clean) > 5) {
    df_scaled <- scale(df_region_clean)
    pca_region <- princomp(df_scaled, scores = TRUE)
    loadings_region <- pca_region$loadings
    scores <- pca_region$scores
    num_pc <- min(8, ncol(scores))  # Yaln??zca mevcut bile??enler kadar ??izim yap??l??r
    
    # 1. Boxplot of standardized variables and PC scores
    par(mfrow = c(2, 1))
    boxplot(df_scaled, las = 2, col = 'gold', main = paste("Original Variables -", r))
    boxplot(scores, las = 2, col = 'gold', main = paste("Principal Components -", r))
    
    # 2. Significant PCA loadings (|loading| > 0.3)
    par(mar = c(3, 4, 2, 1), mfrow = c(4, 1))
    for (i in 1:min(4, num_pc)) {
      barplot(ifelse(abs(loadings_region[, i]) < 0.3, 0, loadings_region[, i]),
              main = paste("Region:", r, "- PC", i),
              ylim = c(-1, 1),
              col = "grey")
      abline(h = 0)
    }
    
    # 3. Projection of original data using first 1???num_pc PCs
    meanF <- colMeans(df_scaled)
    projection <- matrix(meanF, nrow(df_scaled), ncol(df_scaled), byrow = TRUE)
    par(mfrow = c(2, 5))
    matplot(t(df_scaled), type = 'l', main = paste(r, "- Original Data"), ylim = range(df_scaled))
    matplot(meanF, type = 'l', main = paste(r, "- First 0 PCs"), lwd = 2, ylim = range(df_scaled))
    for (i in 1:num_pc) {
      projection <- projection + scores[, i] %*% t(loadings_region[, i])
      matplot(t(projection), type = 'l', main = paste(r, "- First", i, "PCs"), ylim = range(df_scaled))
      matplot(meanF, type = 'l', lwd = 2, add = TRUE)
    }
    
    # 4. QQ-plot and Shapiro-Wilk test for normality
    results <- data.frame(variable = character(), p_value = numeric())
    for (var in cols_pca) {
      vec <- df_region_clean[[var]]
      sh_test <- shapiro.test(vec)
      qqnorm(vec, main = paste("QQ-Plot for", var, "in", r))
      qqline(vec, col = "red")
      results <- rbind(results, data.frame(variable = var, p_value = sh_test$p.value))
    }
    print(results)
  }
}
  