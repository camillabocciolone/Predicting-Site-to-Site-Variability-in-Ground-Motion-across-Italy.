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