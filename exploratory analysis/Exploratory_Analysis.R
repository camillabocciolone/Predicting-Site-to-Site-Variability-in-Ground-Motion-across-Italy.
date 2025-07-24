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
dim(itac)
#attach(itac)
col_names = names(itac)
col_names
colcarnames = col_names[c(1,2,3,4,5,12,13,14,15,16,22,23,25)] #names of categorical variables
colcarnames
summary(itac)


#-------map visualization--------
italy <- map_data("italy")
map <- leaflet() %>%
  addTiles() %>%  # Add default OpenStreetMap tiles
  setView(lng = 12.5, lat = 41.9, zoom = 6)  # Center map on Italy (longitude, latitude, zoom level)
# Add points on the map (using st_longitude and st_latitude)
map <- map %>%
  addCircleMarkers(data = itac, lng = ~st_longitude, lat = ~st_latitude, color = "red", radius = 5, fill = TRUE, fillOpacity = 0.7)
# Show the map
map

#---------MISSING VALUES CONSIDERATIONS-------
sum(is.na(itac)) #4573 NA (4,23% of total data) 
na = colSums(is.na(itac)) #vector of all cov with respective sum of NA values
naeq0 = na[na>0] #cols with missing values 
naeq0 #OBS: we observe that the main NA values comes from meas. of VS30, depth_to_bedrock, ...
nadif0 = na [na==0] #cols without missing values
nadif0


#Working without NA (dropping col and rows)
itaca = itac [ , colSums(is.na(itac)) == 0]
dim(itaca)
itaca = itac [ rowSums(is.na(itac)) == 0, ] # dataframe without NA
dim(itaca) #too much reduction

#Working without categorical var
itacanum = itaca [, sapply(itaca,is.numeric)]
itacanum = itacanum [, -c(5,10)] #removing also codes variables: proximity code and vs30_calc_method
dim(itacanum)
attach(itacanum)
c = colnames(itacanum)
c
summary(itacanum)

#new map of stations
map <- leaflet() %>% addTiles() %>%  # Add default OpenStreetMap tiles
  setView(lng = 12.5, lat = 41.9, zoom = 6)  # Center map on Italy (longitude, latitude, zoom level)
# Add points on the map (using st_longitude and st_latitude)
map <- map %>% addCircleMarkers(data = itacanum, lng = ~st_longitude, lat = ~st_latitude, color = "red", radius = 5, fill = TRUE, fillOpacity = 0.7)
map #huge dimensionality reduction

mean = colMeans(itacanum)
var = sapply(itacanum, var)
#sd = sapply(itacanum, sd)
#cov = cov(itacanum)
#cor = cor (itacanum)


#---------FIRST APPROACH TO DATA--------
plot(mean)
x_ticks <- pretty(1:100, n = 16)  #adding more x-axis reference numbers
y_ticks <- pretty(1:700, n = 15) #adding more y's
axis(1, at = x_ticks)
axis(2, at = y_ticks)
colnames(itacanum) #starting from fa_pga (i = 14) they all have about the same mean, as we have thought

plot(var)
x_ticks <- pretty(1:100, n = 16) 
y_ticks <- pretty(1:700, n = 100) 
axis(1, at = x_ticks)
axis(2, at = y_ticks)
#Comment of these plots
mean[mean<100] #cov with magnitude scale of 10s (majority of data) 
mean[mean>100] #cov with magnitude scale of 100s (made by VS30)
var[var>10] #High variability terms (includes all VS30 terms)

#let's see their range
summary(st_elevation) #Highest variance. Measurements from sea level to 7000m, max=9999m quite unlikely (there are outliers) see below
#let's remove the outlier:
which(st_elevation>2000)
st_elevation
plot(st_elevation)
itacanum = itacanum[-which(st_elevation==9999),]
detach(itacanum)
attach(itacanum)
dim(itacanum)
summary(st_elevation)
plot(st_elevation) #Now works better 

summary(slope) #ok
summary(vs30_from_profile) 
summary(depth_to_bedrock) 
summary(vseq_bedrock) #VS30 modified based on the depth of bedrock (consider it in the first phase (then not))
summary(vseq_ntc)
summary(vs30_m_s_from_topography)
#n° recordings seems to be quite ok:
summary(nrec_max)
plot(nrec_max)
summary(nrec_min)
plot(nrec_min)
which(nrec_max>200) #the 125th obs is an outlier in terms of nrec_max >> and nrec_min >>
which(nrec_min>100)



#-----------GGpairs-------------
plot(st_longitude, st_latitude, xlim = c(5,20), ylim = c(35,50)) #a map of Italy 
which(st_elevation>1338) #max elevation point is 58th obs., at a latitude = 42.5267, longitude = 13.3509
ggpairs(itacanum[,c(2,1,3,4)]) #covariates not correlated between them (fine)

ggpairs(itacanum[,c(5,6,7,8,9)]) #as we know they are correlated, thus apport the same information 
ggpairs(itacanum[,c(27,28,40,63,64)]) 

#let's compare all vs30 with different  for linear regression
ggpairs(itacanum[, c(5, 7, 8, 9, 32)]) #DS2_T0_450
ggpairs(itacanum[, c(5, 7, 8, 9, 63)]) #DS2_T10_000: correlation is acceptable (vs30 form profile) <-- in this try it is the most correlated case

#Representative periods pga, T=0.04, T=1s, T=5s
dfitac= read.csv("DATASET_WITH_SCORES.csv") 
df_sub <- dfitac[, c(48, 126, 132)] #includer depth to bedrock 23
df_sub$litoclasses <- as.factor(dfitac$Lito_Bucci_simple)
litoclasses= as.factor(dfitac$Lito_Bucci_simple)
ggpairs(df_sub, aes(color = litoclasses)) +
  scale_color_manual(values = c("red", "blue", "green", "orange", "purple", "brown",
                                "pink", "cyan", "yellow", "darkgreen", "grey"))
palette_colors <- c("red", "blue", "green", "orange", "purple", "brown",
                    "pink", "cyan", "yellow", "darkgreen", "grey")[1:length(levels(litoclasses))]
detach(itacanum)
#some more ggpairs
# covariates - dS2S T=0s
grid()
dfitac= read.csv("DATASET_WITH_SCORES.csv") 
colnames(dfitac) #68,22,14,23; 68,126,132,23; 68,128,132,23
df_sub <- dfitac[, c(68, 126, 132,23)]
colnames(df_sub)=c('dS2S SA', 'Vs30', 'Slope', 'Depth to bedrock')
ggpairs(df_sub)
df_sub <- dfitac[, c(68, 22, 14,23)]
ggpairs(df_sub)
df_sub <- dfitac[, c(68, 127, 132,23)]
colnames(df_sub)=c('dS2S SA', 'Vs30', 'Slope', 'Depth to bedrock')
ggpairs(df_sub)
df_sub <- dfitac[, c(79, 126, 132,23)]
ggpairs(df_sub) 

#transformed covariates - dS2S T=0s
dfitac= read.csv("DATASET_WITH_SCORES.csv") 
dfitac=dfitac[,c(84,23,126,132)]
dfitac=na.omit(dfitac)
df_sub <- data.frame(
  dS2S = dfitac$dS2S_SA_T10_000,
  sqrt_Vs30 = sqrt(dfitac$Vs30_Brunelli_simple),
  inv_slope = 1/dfitac$slope50_simple,
  inv_depth = 1/dfitac$depth_to_bedrock
)
colnames(df_sub)=c('dS2S SA', 'sqrt(Vs30)', '1/Slope', '1/Depth')
plot(df_sub, aes(x = inv_slope, y = dS2S)) 
ggpairs(df_sub)
par(mfrow=c(3,3))



#--------------Normality of covariates-------------------
#Let's see the normal distribution of all covariates 
results <- data.frame(period = character(), p_value = numeric(), stringsAsFactors = FALSE) #creates an empty dataset
#ignore "fa_pga", "fa_01_05_mean", "fa_04_08_mean", "fa_07_11_mean", "fa_pgv", "fa_pga_max", "fa_01_05_max", "fa_04_08_max", "fa_07_11_max", 
#ignore "fa_pgv_max", "fa_pga_min","fa_01_05_min","fa_04_08_min","fa_07_11_min","fa_pgv_min",
cols = c("st_elevation", "slope", "vs30_from_profile", "depth_to_bedrock", "vseq_bedrock", "vseq_ntc", "vs30_m_s_from_topography", 
         "dS2S_SA_pga",              "dS2S_SA_T0_010",           "dS2S_SA_T0_025",           "dS2S_SA_T0_040",
         "dS2S_SA_T0_050",           "dS2S_SA_T0_070",           "dS2S_SA_T0_100",           "dS2S_SA_T0_150",           "dS2S_SA_T0_200",
         "dS2S_SA_T0_250",           "dS2S_SA_T0_300",           "dS2S_SA_T0_350",           "dS2S_SA_T0_400",           "dS2S_SA_T0_450",
         "dS2S_SA_T0_500",           "dS2S_SA_T0_600",           "dS2S_SA_T0_700",           "dS2S_SA_T0_750",           "dS2S_SA_T0_800",
         "dS2S_SA_T0_900",           "dS2S_SA_T1_000",           "dS2S_SA_T1_200",           "dS2S_SA_T1_400",           "dS2S_SA_T1_600",
         "dS2S_SA_T1_800",           "dS2S_SA_T2_000",           "dS2S_SA_T2_500",           "dS2S_SA_T3_000",           "dS2S_SA_T3_500",
         "dS2S_SA_T4_000",           "dS2S_SA_T4_500",           "dS2S_SA_T5_000",           "dS2S_SA_T6_000",           "dS2S_SA_T7_000",
         "dS2S_SA_T8_000",           "dS2S_SA_T9_000",           "dS2S_SA_T10_000",          "dS2S_SA_pgv")
for (col in cols){
  vec = itacanum[[col]]
  sh_test <- shapiro.test(vec)
  qqnorm(vec, main = paste("QQ-Plot for", col))
  qqline(vec, col = "red")
  cat("Results of the Shapiro-Wilk test for", col, ":\n")
  print(sh_test)
  cat("\n")
  results <- rbind(results, data.frame(period = col, p_value = sh_test$p.value))
  readline(prompt = "Press [Enter] to continue...")
}
par(mfrow=c(1,1))




#-------NEW COVARIATES-------
itac = read.csv("FINAL_DATASET.csv") 
newcov=read.table("ampli_f0.txt", header=TRUE, sep=';')
View(newcov)

pos <- match(itac$network_code, newcov$network_code)
sum(is.na(pos))
pos <- match(itac$station_code, newcov$station_code)
sum(is.na(pos))

itac$f01_hz <- NA #Create the column corresponding to the fundamental frequency [Hz] of vibration at the site in the first direction
itac$f01_hz[!is.na(pos)] <- newcov$f01_hv_microtremor_hz[pos[!is.na(pos)]]
sum(is.na(itac[,'f01_hz'])) #642 NA

itac$f01_ampli <- NA #Create the column corresponding to the resonance peak of HSVR curve 1
pos <- match(itac$station_code, newcov$station_code)
itac$f01_ampli[!is.na(pos)] <- newcov$f01_hv_microtremor_ampli[pos[!is.na(pos)]]
sum(is.na(itac[,'f01_ampli'])) #599 NA

itac$f02_hz <- NA # resonance peak of HSVR curve 2
pos <- match(itac$station_code, newcov$station_code)
itac$f02_hz[!is.na(pos)] <- newcov$f02_hv_microtremor_hz[pos[!is.na(pos)]]
sum(is.na(itac[,'f02_hz'])) #866 NA

itac$f02_ampli <- NA
pos <- match(itac$station_code, newcov$station_code)
itac$f02_ampli[!is.na(pos)] <- newcov$f02_hv_microtremor_ampli[pos[!is.na(pos)]]
sum(is.na(itac[,'f02_ampli'])) #861 NA
#So many NA

#write.csv(itac, "DATASET_WITH_NEW_COV.csv", row.names = FALSE) #create datset with f0 covariates

range(na.omit(itac$f01_hz))
range(na.omit(itac$f01_ampli))
range(na.omit(itac$f02_hz))
range(na.omit(itac$f02_ampli))


ggpairs(itac[,c('dS2S_SA_pga','dS2S_SA_T0_500','dS2S_SA_T1_000','f01_hz')]) #highest correlation for periods around 1/1hz = 1s
ggpairs(itac[,c('dS2S_SA_pga','dS2S_SA_T0_500','dS2S_SA_T1_000','f01_ampli')]) #may suggest a correlation (lower wrt f01_hz) with a transformation of the var
ggpairs(itac[,c('dS2S_SA_pga','dS2S_SA_T0_500','dS2S_SA_T1_000','f02_hz')]) #very correlated with period around 1/2hz=0.5s
ggpairs(itac[,c('dS2S_SA_pga','dS2S_SA_T0_500','dS2S_SA_T1_000','f02_ampli')])

Tn = (c(0,0.01,0.025,0.04,0.05,0.07,0.1,0.15,0.2,0.25,0.3,0.35,0.40,0.45,0.5,0.6,0.7,0.75,0.8,0.9,1,1.2,1.4,1.6,1.8,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10))
cor_values_1hz <- sapply(itac[, 43:79], function(x) cor(x, itac$f01_hz, use = "complete.obs"))
plot(Tn, cor_values_1hz, type = "o", pch = 19, col = "steelblue", xlab = "Period Tn [s]", ylab = "Correlation", main = "Correlation dS2S with f01Hz (f0 along 1st dir) across Periods")
cor_values_2hz <- sapply(itac[, 43:79], function(x) cor(x, itac$f02_hz, use = "complete.obs"))
plot(Tn, cor_values_2hz, type = "o", pch = 19, col = "steelblue", xlab = "Period Tn [s]", ylab = "Correlation", main = "Correlation dS2S with f02Hz (f0 along 2nd dir) across Periods")
cor_values_1ampli <- sapply(itac[, 43:79], function(x) cor(x, itac$f01_ampli, use = "complete.obs"))
plot(Tn, cor_values_1ampli, type = "o", pch = 19, col = "steelblue", xlab = "Period Tn [s]", ylab = "Correlation", main = "Correlation dS2S with f01ampli (resonance peak HSVR) ")
cor_values_2ampli <- sapply(itac[, 43:79], function(x) cor(x, itac$f02_ampli, use = "complete.obs"))
plot(Tn, cor_values_2ampli, type = "o", pch = 19, col = "steelblue", xlab = "Period Tn [s]", ylab = "Correlation", main = "Correlation dS2S with f02ampli (resonance peak HSVR) ")




