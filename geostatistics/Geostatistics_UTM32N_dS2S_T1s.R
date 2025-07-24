rm(list=ls())
graphics.off()
library(mvtnorm) # to deal with multivariate normal distributions
library(car) # "Companion to Applied Regression" for regression analysis
library(dplyr)
library( GGally)
library(maps)
library(leaflet) #to use the map
library(corrplot)
library(sp)           ## Data management
library(lattice)      ## Data management
library(geoR)         ## Geostatistics
library(gstat)        ## Geostatistics
library(sf)
library(raster)
# install.packages("rnaturalearth")
# install.packages("rnaturalearthdata")
library(rnaturalearth)
library(rnaturalearthdata)
v.f <- function(x, ...){100-cov.spatial(x, ...)}
v.f.est<-function(x,C0, ...){C0-cov.spatial(x, ...)}
itac = read.csv("FINAL_DATASET.csv") 
#Change of reference sys
coordinates = na.omit(itac[,6:7]) #Taking only lat and long
psf = st_as_sf(coordinates, coords = c("st_longitude", "st_latitude"), crs = 4326)
# Transform in UTM 32N 
new_coordinates = st_transform(psf, 32632) 
newcoords = as.data.frame(st_coordinates(new_coordinates))
colnames(newcoords) = c('x_no_noise', 'y_no_noise')
itaca=itac
itaca = cbind (newcoords, itaca) #New dataset with coordinates
ita = itaca

#generate some noise
set.seed(42)
newcoords = as.data.frame(st_coordinates(new_coordinates))
colnames(newcoords) = c('x', 'y')
# Aggiungi jitter per evitare problemi numerici nel kriging universale
newcoords$x = newcoords$x + rnorm(nrow(newcoords), mean = 0, sd = 1)
newcoords$y = newcoords$y + rnorm(nrow(newcoords), mean = 0, sd = 1)
itaca = cbind (newcoords, itaca) #New dataset with coordinates
ita = itaca

#To use the grouping of manova: livello → nuova classe numerica
mapping <- c(1, 2, 2, 2, 3, 1, 1, 1, 3, 3, 3)
itaca$Lito_Bucci_grouped <- factor(mapping[as.integer(itaca$Lito_Bucci_simple)], levels = c(1, 2, 3))
ita=itaca

attach(itaca)


#Maps: Universal Kriging for dS2S_pga

# ---- Model with both vs30 and slope ----
#dS2S_T_1s ~sqrt(Vs30_Brunelli_simple) + 1/slope50_simple
f = dS2S_SA_T1_000 ~ 1/slope50_simple + sqrt(Vs30_Brunelli_simple)
ita = itaca[!is.na(itaca$slope50_simple) & !is.na(itaca$Vs30_Brunelli_simple) &!is.na(itaca$dS2S_SA_T1_000), ] #take the initial datframe for simplicity (in order to remove only necessary NA). OSS. Don't remove NA of x,y, there are not
coordinates(ita)=c('x','y') #Create a spatial dataframe
proj4string(ita) = CRS("+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs") #set CRS

v.gls<-variogram(f, ita) 
plot(v.gls) #Definetly similar
v.gls.fit <- fit.variogram(v.gls, vgm(0.3, "Sph", 250000, 0.2))
plot(v.gls, v.gls.fit, pch = 19)
g.gls.tr <- gstat(formula = f, data = ita, model = v.gls.fit, nmax = 50, set = list(gls=1))


#s0=s1 known (1st observation)
predict(g.gls.tr, ita[1,])
s0.new=data.frame(x=ita@coords[1,1], y=ita@coords[1,2], Vs30_Brunelli_simple=483.8733, slope50_simple=5.50358343)
coordinates(s0.new)=c('x','y')
s0.new <- as(s0.new, 'SpatialPointsDataFrame')
proj4string(s0.new) = CRS(proj4string(ita)) #set same CRS
predict(g.gls.tr, s0.new)


#map
#vs30
italy = ne_countries(country = "Italy", returnclass = "sf")
italy_utm = st_transform(italy, crs = 32632)  # UTM 32N as vs30 map 
xi = seq(6, 19, by = 0.1) #10km accuracy grid of Italy
yi = seq(36, 48, by = 0.1)
ita.grid = expand.grid(x = xi, y = yi)
dim(ita.grid)
ita.grid = st_as_sf(ita.grid, coords = c("x", "y"), crs = 4326) #latitude longitude WGS84 
ita.grid = st_transform(ita.grid, crs = 32632) #transform in UTM32N as vs30 map
#crs(ita.grid)
vs30_r = raster("INGV_data/Vs30_Brunelli_new/vs30ms") #get vs30 map
#crs(vs30_r)
inside_italy = st_within(ita.grid, italy_utm, sparse = FALSE)[,1] #redundant: uses only values inside italian territory
ita.grid = ita.grid[inside_italy, ]
ita.grid_sp = as(ita.grid, "Spatial") 
vs30_val_italy = extract(vs30_r, ita.grid_sp)
ita.grid$Vs30_Brunelli_simple = vs30_val_italy #assign to the pts of the grid the corresponding vs30
#Not take Sardegna
sard_bbox_ll <- st_as_sfc(st_bbox(c(
  xmin = 8.0, xmax = 10.2,
  ymin = 38.5, ymax = 41.5
), crs = 4326))
sard_bbox_utm <- st_transform(sard_bbox_ll, 32632)
in_sardinia <- st_within(ita.grid, sard_bbox_utm, sparse = FALSE)[,1]
ita.grid$Vs30_Brunelli_simple[in_sardinia] <- NA
#Slope
slope_r = raster("INGV_data/Slope_50") #get slope map
ita.grid_sp = as(ita.grid, "Spatial") 
slope_val_italy = extract(slope_r, ita.grid_sp)
ita.grid$slope50_simple = slope_val_italy #assign to the pts of the grid the corresponding slope
#Not take Sardegna
sard_bbox_ll <- st_as_sfc(st_bbox(c(
  xmin = 8.0, xmax = 10.2,
  ymin = 38.5, ymax = 41.5
), crs = 4326))
sard_bbox_utm <- st_transform(sard_bbox_ll, 32632)
in_sardinia <- st_within(ita.grid, sard_bbox_utm, sparse = FALSE)[,1]
ita.grid$slope50_simple[in_sardinia] <- NA

ita.grid = as(ita.grid, "Spatial")

g.gls.tr <- gstat(formula = f, data = ita, model = v.gls.fit, set = list(gls=1))
lz.ok = predict(g.gls.tr, ita.grid, BLUE = FALSE)
spplot(lz.ok, col.regions = terrain.colors(9))
lz.sf <- st_as_sf(lz.ok)
#map for the DS2S
ggplot(lz.sf) +
  geom_sf(aes(color = var1.pred)) +  
  scale_color_viridis_c(na.value = "grey90") +
  theme_minimal() +
  labs(color = "Predicted Value", title = "Universal Kriging Prediction: dS2S for T=1s")
#map for the var
ggplot(lz.sf) +
  geom_sf(aes(color = var1.var)) +  
  scale_color_viridis_c(na.value = "grey90") +
  theme_minimal() +
  labs(color = "Predicted Value", title = "Universal Kriging Prediction: variability of dS2S for T=1s")


#CV (DON'T RUN if not necessary)
cv=krige.cv(formula = f, locations = ita, model = v.gls.fit)
dim(cv)
sum(!is.na(cv@data[["var1.pred"]])) 
cv_residual = na.omit(cv$observed - cv$var1.pred)
length(cv_residual)
rmse = sqrt(mean(cv_residual^2))
mae  <- mean(abs(cv_residual))
me   <- mean(cv_residual)
loo_mse <- mean(cv_residual^2)
loo_mse #0.2018675






# ---- Model with Lito Bucci ----
#4) dS2S_T_1s ~sqrt(Vs30_Brunelli_simple) + 1/slope50_simple + LitoBucci
itaca$Lito_Bucci_simple=factor(itaca$Lito_Bucci_simple)
levels(itaca$Lito_Bucci_simple)
f = dS2S_SA_T1_000 ~ 1/slope50_simple + sqrt(Vs30_Brunelli_simple) + Lito_Bucci_simple
ita = itaca[!is.na(itaca$slope50_simple) & !is.na(itaca$Vs30_Brunelli_simple) &!is.na(itaca$dS2S_SA_T1_000) & !is.na(itaca$Lito_Bucci_simple), ] 
coordinates(ita)=c('x','y') #Create a spatial dataframe
proj4string(ita) = CRS("+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs") #set CRS

v.gls<-variogram(f, ita) 
plot(v.gls) 
v.gls.fit <- fit.variogram(v.gls, vgm(0.1, "Sph", 250000, 0.2))
plot(v.gls, v.gls.fit, pch = 19)
g.gls.tr <- gstat(formula = f, data = ita, model = v.gls.fit, nmax = 50, set = list(gls=1))


#map
#vs30
italy = ne_countries(country = "Italy", returnclass = "sf")
italy_utm = st_transform(italy, crs = 32632)  # UTM 32N as vs30 map 
xi = seq(6, 19, by = 0.1) #10km accuracy grid of Italy
yi = seq(36, 48, by = 0.1)
ita.grid = expand.grid(x = xi, y = yi)
dim(ita.grid)
ita.grid = st_as_sf(ita.grid, coords = c("x", "y"), crs = 4326) #latitude longitude WGS84 
ita.grid = st_transform(ita.grid, crs = 32632) #transform in UTM32N as vs30 map
#crs(ita.grid)
vs30_r = raster("INGV_data/Vs30_Brunelli_new/vs30ms") #get vs30 map
#crs(vs30_r)
inside_italy = st_within(ita.grid, italy_utm, sparse = FALSE)[,1] #redundant: uses only values inside italian territory
ita.grid = ita.grid[inside_italy, ]
ita.grid_sp = as(ita.grid, "Spatial") 
vs30_val_italy = extract(vs30_r, ita.grid_sp)
ita.grid$Vs30_Brunelli_simple = vs30_val_italy #assign to the pts of the grid the corresponding vs30
#Not take Sardegna
sard_bbox_ll <- st_as_sfc(st_bbox(c(
  xmin = 8.0, xmax = 10.2,
  ymin = 38.5, ymax = 41.5
), crs = 4326))
sard_bbox_utm <- st_transform(sard_bbox_ll, 32632)
in_sardinia <- st_within(ita.grid, sard_bbox_utm, sparse = FALSE)[,1]
ita.grid$Vs30_Brunelli_simple[in_sardinia] <- NA
#Slope
slope_r = raster("INGV_data/Slope_50") #get slope map
ita.grid_sp = as(ita.grid, "Spatial") 
slope_val_italy = extract(slope_r, ita.grid_sp)
ita.grid$slope50_simple = slope_val_italy #assign to the pts of the grid the corresponding slope
#Not take Sardegna
sard_bbox_ll <- st_as_sfc(st_bbox(c(
  xmin = 8.0, xmax = 10.2,
  ymin = 38.5, ymax = 41.5
), crs = 4326))
sard_bbox_utm <- st_transform(sard_bbox_ll, 32632)
in_sardinia <- st_within(ita.grid, sard_bbox_utm, sparse = FALSE)[,1]
ita.grid$slope50_simple[in_sardinia] <- NA
#LitoBucci
Lito_r = raster("INGV_data/Carla_lito_Bucci_aggregata/nuova_mappa_lito/new_lit_11_24_matlab.tif") #get Lito map
ita.grid_sp = as(ita.grid, "Spatial") 
Lito_val_italy = extract(Lito_r, ita.grid_sp)
ita.grid$Lito_Bucci_simple = factor(Lito_val_italy, levels = 1:11) #assign to the pts of the grid the corresponding Litology class
#Not take Sardegna
sard_bbox_ll <- st_as_sfc(st_bbox(c(
  xmin = 8.0, xmax = 10.2,
  ymin = 38.5, ymax = 41.5
), crs = 4326))
sard_bbox_utm <- st_transform(sard_bbox_ll, 32632)
in_sardinia <- st_within(ita.grid, sard_bbox_utm, sparse = FALSE)[,1]
ita.grid$Lito_Bucci_simple[in_sardinia] <- NA

ita.grid = as(ita.grid, "Spatial")

g.gls.tr <- gstat(formula = f, data = ita, model = v.gls.fit, set = list(gls=1))
lz.ok = predict(g.gls.tr, ita.grid, BLUE = FALSE)
spplot(lz.ok, col.regions = terrain.colors(9))
lz.sf <- st_as_sf(lz.ok)
#map for the DS2S
ggplot(lz.sf) +
  geom_sf(aes(color = var1.pred)) +  
  scale_color_viridis_c(na.value = "grey90") +
  theme_minimal() +
  labs(color = "Predicted Value", title = "Universal Kriging Prediction: dS2S for T=1s")
#map for the var
ggplot(lz.sf) +
  geom_sf(aes(color = var1.var)) +  
  scale_color_viridis_c(na.value = "grey90") +
  theme_minimal() +
  labs(color = "Predicted Value", title = "Universal Kriging Prediction: variability of dS2S for T=1s")

#CV (DON'T RUN if not necessary)
cv=krige.cv(formula = f, locations = ita, model = v.gls.fit)
dim(cv)
sum(!is.na(cv@data[["var1.pred"]])) 
cv_residual = na.omit(cv$observed - cv$var1.pred)
length(cv_residual)
rmse = sqrt(mean(cv_residual^2))
mae  <- mean(abs(cv_residual))
me   <- mean(cv_residual)
loo_mse <- mean(cv_residual^2)
loo_mse #0.1967176





# ---- Model with Lito Bucci grouped ----
#5) dS2S_T_1s ~sqrt(Vs30_Brunelli_simple) + 1/slope50_simple + LitoBucciGrouped
itaca$Lito_Bucci_grouped=factor(itaca$Lito_Bucci_grouped)
levels(itaca$Lito_Bucci_grouped)
f = dS2S_SA_T1_000 ~ 1/slope50_simple + sqrt(Vs30_Brunelli_simple) + Lito_Bucci_grouped
ita = itaca[!is.na(itaca$slope50_simple) & !is.na(itaca$Vs30_Brunelli_simple) &!is.na(itaca$dS2S_SA_T1_000) & !is.na(itaca$Lito_Bucci_grouped), ] 
coordinates(ita)=c('x','y') #Create a spatial dataframe
proj4string(ita) = CRS("+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs") #set CRS

v.gls<-variogram(f, ita) 
plot(v.gls) 
v.gls.fit <- fit.variogram(v.gls, vgm(0.1, "Sph", 250000, 0.2))
plot(v.gls, v.gls.fit, pch = 19)
g.gls.tr <- gstat(formula = f, data = ita, model = v.gls.fit, nmax = 50, set = list(gls=1))


#map
#vs30
italy = ne_countries(country = "Italy", returnclass = "sf")
italy_utm = st_transform(italy, crs = 32632)  # UTM 32N as vs30 map 
xi = seq(6, 19, by = 0.1) #10km accuracy grid of Italy
yi = seq(36, 48, by = 0.1)
ita.grid = expand.grid(x = xi, y = yi)
dim(ita.grid)
ita.grid = st_as_sf(ita.grid, coords = c("x", "y"), crs = 4326) #latitude longitude WGS84 
ita.grid = st_transform(ita.grid, crs = 32632) #transform in UTM32N as vs30 map
#crs(ita.grid)
vs30_r = raster("INGV_data/Vs30_Brunelli_new/vs30ms") #get vs30 map
#crs(vs30_r)
inside_italy = st_within(ita.grid, italy_utm, sparse = FALSE)[,1] #redundant: uses only values inside italian territory
ita.grid = ita.grid[inside_italy, ]
ita.grid_sp = as(ita.grid, "Spatial") 
vs30_val_italy = extract(vs30_r, ita.grid_sp)
ita.grid$Vs30_Brunelli_simple = vs30_val_italy #assign to the pts of the grid the corresponding vs30
#Not take Sardegna
sard_bbox_ll <- st_as_sfc(st_bbox(c(
  xmin = 8.0, xmax = 10.2,
  ymin = 38.5, ymax = 41.5
), crs = 4326))
sard_bbox_utm <- st_transform(sard_bbox_ll, 32632)
in_sardinia <- st_within(ita.grid, sard_bbox_utm, sparse = FALSE)[,1]
ita.grid$Vs30_Brunelli_simple[in_sardinia] <- NA
#Slope
slope_r = raster("INGV_data/Slope_50") #get slope map
ita.grid_sp = as(ita.grid, "Spatial") 
slope_val_italy = extract(slope_r, ita.grid_sp)
ita.grid$slope50_simple = slope_val_italy #assign to the pts of the grid the corresponding slope
#Not take Sardegna
sard_bbox_ll <- st_as_sfc(st_bbox(c(
  xmin = 8.0, xmax = 10.2,
  ymin = 38.5, ymax = 41.5
), crs = 4326))
sard_bbox_utm <- st_transform(sard_bbox_ll, 32632)
in_sardinia <- st_within(ita.grid, sard_bbox_utm, sparse = FALSE)[,1]
ita.grid$slope50_simple[in_sardinia] <- NA
#LitoBucci
Lito_r = raster("INGV_data/Carla_lito_Bucci_aggregata/nuova_mappa_lito/new_lit_11_24_matlab.tif") #get Lito map
ita.grid_sp = as(ita.grid, "Spatial") 
Lito_val_italy = extract(Lito_r, ita.grid_sp)
ita.grid$Lito_Bucci_simple = factor(Lito_val_italy, levels = 1:11) #assign to the pts of the grid the corresponding Litology class
mapping <- c(1, 2, 2, 2, 3, 1, 1, 1, 3, 3, 3) # grouping of the levels of Lito_Bucci_simple
ita.grid$Lito_Bucci_grouped <- factor(mapping[as.integer(ita.grid$Lito_Bucci_simple)], levels = c(1, 2, 3))
#Not take Sardegna
sard_bbox_ll <- st_as_sfc(st_bbox(c(
  xmin = 8.0, xmax = 10.2,
  ymin = 38.5, ymax = 41.5
), crs = 4326))
sard_bbox_utm <- st_transform(sard_bbox_ll, 32632)
in_sardinia <- st_within(ita.grid, sard_bbox_utm, sparse = FALSE)[,1]
ita.grid$Lito_Bucci_grouped[in_sardinia] <- NA

ita.grid = as(ita.grid, "Spatial")

g.gls.tr <- gstat(formula = f, data = ita, model = v.gls.fit, set = list(gls=1))
lz.ok = predict(g.gls.tr, ita.grid, BLUE = FALSE)
spplot(lz.ok, col.regions = terrain.colors(9))
lz.sf <- st_as_sf(lz.ok)
#map for the DS2S
ggplot(lz.sf) +
  geom_sf(aes(color = var1.pred)) +  
  scale_color_viridis_c(na.value = "grey90") +
  theme_minimal() +
  labs(color = "Predicted Value", title = "Universal Kriging Prediction: dS2S for T=1s")
#map for the var
ggplot(lz.sf) +
  geom_sf(aes(color = var1.var)) +  
  scale_color_viridis_c(na.value = "grey90") +
  theme_minimal() +
  labs(color = "Predicted Value", title = "Universal Kriging Prediction: variability of dS2S for T=1s")

#CV (DON'T RUN if not necessary)
cv=krige.cv(formula = f, locations = ita, model = v.gls.fit)
dim(cv)
sum(!is.na(cv@data[["var1.pred"]])) 
cv_residual = na.omit(cv$observed - cv$var1.pred)
length(cv_residual)
rmse = sqrt(mean(cv_residual^2))
mae  <- mean(abs(cv_residual))
me   <- mean(cv_residual)
loo_mse <- mean(cv_residual^2)
loo_mse #0.1977965




detach(itaca)
