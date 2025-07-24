library(mvtnorm) # to deal with multivariate normal distributions
library(car) # "Companion to Applied Regression" for regression analysis
library(dplyr)
library( GGally)
library(maps)
library(leaflet) #to use the map
library(corrplot)

#Load data
itac = read.csv("ITACAs2s_2_flatfile/ITACAs2s_SA_2_0.csv", sep=";")
View(itac)
col_names = names(itac)
col_names
colcarnames = col_names[c(1,2,3,4,5,11,12,13,14,15,16,22,23,25)] #names of categorical variables
colcarnames
#map visualization
italy <- map_data("italy")
map <- leaflet() %>%
  addTiles() %>%  # Add default OpenStreetMap tiles
  setView(lng = 12.5, lat = 41.9, zoom = 6)  # Center map on Italy (longitude, latitude, zoom level)
# Add points on the map (using st_longitude and st_latitude)
map <- map %>%
  addCircleMarkers(data = itac, lng = ~st_longitude, lat = ~st_latitude, color = "red", radius = 5, fill = TRUE, fillOpacity = 0.7)
# Show the map
map


#NUMERICAL REGIONAL ANALYSIS
#function to work with rectangular regions
itaca = itac[!is.na(itac$st_latitude) & !is.na(itac$st_longitude), ]
#ALPI OCCIDENTALI: lat 44-45.5, long 6.5-8.5


lat_min <- 41
lat_max <- 44
lon_min <- 10
lon_max <- 16
itaca_region <- itaca[itaca$st_latitude >= lat_min & itaca$st_latitude <= lat_max & itaca$st_longitude >= lon_min & itaca$st_longitude <= lon_max, ]
dim(itaca_region)

#Map visualization
map <- leaflet() %>%
  addTiles() %>%  # Add default OpenStreetMap tiles
  setView(lng = 12.5, lat = 41.9, zoom = 6)  # Center map on Italy (longitude, latitude, zoom level)
# Add points on the map (using st_longitude and st_latitude)
map <- map %>%
  addCircleMarkers(data = itaca_region, lng = ~st_longitude, lat = ~st_latitude, color = "red", radius = 5, fill = TRUE, fillOpacity = 0.7)
# Show the map
map

#CENTRAL ITALY
#Working without categorical variables:  using as representative DS2S PGA and SA(T=5s)
itaca_region = itaca_region [, sapply(itaca_region,is.numeric)]
itaca_region = itaca_region [, -c(5,10)]

itaca_region = itaca_region[!is.na(itaca_region$vs30_from_profile), ]

#itaca_region = na.omit(itaca_region) #not too much value to work without all NA (vs_30_from_profile has quite all values null)
itaca_region_num = itaca_region[,c(1,2,3,4,9,27,58)] #using vs_30_from_topography (because of availability, instead of vs_30_from_profile) and not considering depth_to_bedrock
itaca_region_num = na.omit(itaca_region_num)
dim(itaca_region_num)

ggpairs(itaca_region_num)
#LM doesn't work
li~lm(itaca_region_num$dS2S_SA_pga ~ itaca_region_num$st_elevation + itaca_region_num$slope + itaca_region_num$vs30_m_s_from_topography, data=itaca_region_num)
summary(li)

