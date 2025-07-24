library(mvtnorm) # to deal with multivariate normal distributions
library(car) # "Companion to Applied Regression" for regression analysis
library(dplyr)
library( GGally)
library(maps)
library(leaflet) #to use the map
library(corrplot)

#Load data
itac = read.csv("FINAL_DATASET.csv")
View(itac)
#####################SOME COMPARISON WITH SGOBBA########################
#Obs. on vs30_from_profile
(916-sum(is.na(itac$vs30_from_profile)))/916 #Less than 29% of avaailable vs_30_from_profile (Sgobba <25%)

#--------MATPLOT for lithology classes--------------
d0 = itac[,43:79]
d = na.omit(d0)
d0 = na.omit(d0)
colnames(d) = c(
  "pga", "t0_010", "t0_025", "t0_040", "t0_050", "t0_070", "t0_100",
  "t0_150", "t0_200", "t0_250", "t0_300", "t0_350", "t0_400", 
  "t0_450", "t0_500", "t0_600", "t0_700", "t0_750", "t0_800", "t0_900", 
  "t1_000", "t1_200", "t1_400", "t1_600", "t1_800",
  "t2_000", "t2_500", "t3_000", "t3_500", "t4_000", "t4_500", 
  "t5_000", "t6_000", "t7_000", "t8_000", "t9_000", "t10_000")
nomi_colonne <- colnames(d)
converti_periodi <- function(nomi) {
  sapply(nomi, function(x) {
    if (x == "pga") {
      return(0)  # pga = periodo 0
    } else {
      # Rimuove 't' e sostituisce '_' con '.', poi converte in numero
      return(as.numeric(gsub("_", ".", sub("t", "", x))))
    }
  })
}
periodi <- converti_periodi(colnames(d))
y_all <- d[1:879, ]
y_all <- as.matrix(y_all)
y_all_t <- t(y_all)
matplot(periodi, y_all_t, type = "l",lty = 1, col = rgb(0.6, 0.6, 0.6, 0.3),
        xlab = "Period T [s]", ylab = "dS2S SA [cm/s^2]", main = "Spectra comparison based on Bucci lithological classes") #How to use log scale
#Use lithological classes
d00 = itac[as.numeric(rownames(d)),]
li = as.numeric(levels(as.factor(d00$Lito_ISPRA_simple)))
lb = as.numeric(levels(as.factor(d00$Lito_Bucci_simple)))
for(i in lb){
  
  d0i = d00[d00$Lito_Bucci_simple==i, 43:79];
  d0i = na.omit(d0i)
  m = apply(d0i,2,median)
  lines(periodi, as.numeric(m), col = rainbow(11)[i], lwd = 2)
  
}
m = apply(d00[,43:79],2,median)
lines(periodi, as.numeric(m), col = "brown",lty =3, lwd = 2) #line of the average among all sites
Bucci = as.factor(d00$Lito_Bucci_simple)
levels(Bucci) = c("EMR", "ARTR", "ASM", "CCS","FGS","IMR","LIDO","LM","PYD","PYR","SMG")
legend("topright", legend = c(levels(Bucci), "Median","Variance"), col = c(rainbow(11), "brown"), pch = 10, lwd=2, cex= 0.7,title = "Lithology class", ncol=2)
v = apply(d00[,43:79],2,var)
plot(periodi, v, type='b',lty=1,pch=16,col="darkred",xlab='Period T [s]', ylab='Variance of dS2S SA')
grid()
graphics.off()



#---------------plots vs30 vs vs30Mori vs vs30Bucci----------------
#Unique plot
vs30df = na.omit(itac[,c(17,119,121,123,125)])
dim = dim(vs30df)[1]
tc = rgb(c(1, 0, 0, 0, 0), c(0, 1, 0, 0, 0), c(0, 0, 0, 1, 0), alpha = c(0.5, 0.5, 0.5, 0.3, 0.7))

intervals = seq(min(vs30df), max(vs30df),length.out = 10)
vs30df2 = cut(as.vector(vs30df[,2]), breaks = intervals, inlcude.lowest = TRUE)
rel_freq = table(vs30df2)/dim
midpoints = (head(intervals, -1) + tail(intervals, -1)) / 2
barplot(rel_freq, names.arg = round(midpoints, 1), col = tc[2], xlab= "vs30", ylab="relative frequencies", main = "Vs30 comparison" )

intervals = seq(min(vs30df), max(vs30df),length.out = 10)
vs30df2 = cut(as.vector(vs30df[,3]), breaks = intervals, inlcude.lowest = TRUE)
rel_freq = table(vs30df2)/dim
midpoints = (head(intervals, -1) + tail(intervals, -1)) / 2
barplot(rel_freq, names.arg = round(midpoints, 1), col = tc[3], add = TRUE)

intervals = seq(min(vs30df), max(vs30df),length.out = 10)
vs30df2 = cut(as.vector(vs30df[,1]), breaks = intervals, inlcude.lowest = TRUE)
rel_freq = table(vs30df2)/dim
midpoints = (head(intervals, -1) + tail(intervals, -1)) / 2
barplot(rel_freq, names.arg = round(midpoints, 1), col = tc[1], add = TRUE)

legend("topright", legend = c("from profile", "Mori", "Brunelli"), fill = tc[1:3])

#Different plots
par(mfrow=c(3,1))
hist(vs30df[,1], xlim = c(0,2000), ylim=c(0,100), col = "#00A600", xlab = "vs30 from profile", ylab = "absolute frequency", main = "Vs30 Comparison")
hist(vs30df[,2], xlim = c(0,2000), ylim=c(0,100), col = "#E6E600", xlab = "vs30 Mori", ylab = "absolute frequency")
hist(vs30df[,3], xlim = c(0,2000), ylim=c(0,100), col = "#EAB64E", xlab = "vs30 Brunelli", ylab = "absolute frequency")     

#heat map of ≠vs30
par(mfrow=c(1,3))
#vs30 from profile
par(mar = c(4, 4, 2, 2))
map("italy", fill = TRUE, col = "lightgray", bg = "lightblue", mar = c(0, 0, 0, 0))
vs30df = na.omit(itac[,c(7,6,17)])
color_scale = colorRampPalette(c("blue", "green", "yellow", "red"))(100)
col_in = as.numeric(cut(vs30df$vs30_from_profile, breaks = 100))
point_colors = color_scale[col_in]
points(vs30df$st_longitude, vs30df$st_latitude, col = point_colors, pch = 20, cex = 1.5)

min_value <- min(vs30df$vs30_from_profile, na.rm = TRUE)
max_value <- max(vs30df$vs30_from_profile, na.rm = TRUE)
mean_value <- mean(vs30df$vs30_from_profile, na.rm = TRUE)
median_value <- median(vs30df$vs30_from_profile, na.rm = TRUE)

legend("bottomleft", legend = c(sprintf("Min: %.2f", min_value),
                                sprintf("Median: %.2f", median_value),
                                sprintf("Mean: %.2f", mean_value),
                                sprintf("Max: %.2f", max_value)),
       fill = c("blue", "green", "yellow", "red"), title = "Vs30 from profile")
#vs30 Mori
map("italy", fill = TRUE, col = "lightgray", bg = "lightblue", mar = c(0, 0, 0, 0))
vs30df = na.omit(itac[,c(7,6,119)])
vs30df = vs30df[which(vs30df[,3]>min(vs30df)),]
color_scale = colorRampPalette(c("blue", "green", "yellow", "red"))(100)
col_in = as.numeric(cut(vs30df$Vs30_MORI_simple, breaks = 100))
point_colors = color_scale[col_in]
points(vs30df$st_longitude, vs30df$st_latitude, col = point_colors, pch = 20, cex = 1.5)

min_value <- min(vs30df$Vs30_MORI_simple, na.rm = TRUE)
max_value <- max(vs30df$Vs30_MORI_simple, na.rm = TRUE)
mean_value <- mean(vs30df$Vs30_MORI_simple, na.rm = TRUE)
median_value <- median(vs30df$Vs30_MORI_simple, na.rm = TRUE)

legend("bottomleft", legend = c(sprintf("Min: %.2f", min_value),
                                sprintf("Median: %.2f", median_value),
                                sprintf("Mean: %.2f", mean_value),
                                sprintf("Max: %.2f", max_value)),
       fill = c("blue", "green", "yellow", "red"), title = "Vs30 Mori")
#vs30 Brunelli
map("italy", fill = TRUE, col = "lightgray", bg = "lightblue", mar = c(0, 0, 0, 0))
vs30df = na.omit(itac[,c(7,6,121)])
color_scale = colorRampPalette(c("blue", "green", "yellow", "red"))(100)
col_in = as.numeric(cut(vs30df$Vs30_Brunelli_simple, breaks = 100))
point_colors = color_scale[col_in]
points(vs30df$st_longitude, vs30df$st_latitude, col = point_colors, pch = 20, cex = 1.5)

min_value <- min(vs30df$Vs30_Brunelli_simple, na.rm = TRUE)
max_value <- max(vs30df$Vs30_Brunelli_simple, na.rm = TRUE)
mean_value <- mean(vs30df$Vs30_Brunelli_simple, na.rm = TRUE)
median_value <- median(vs30df$Vs30_Brunelli_simple, na.rm = TRUE)

legend("bottomleft", legend = c(sprintf("Min: %.2f", min_value),
                                sprintf("Median: %.2f", median_value),
                                sprintf("Mean: %.2f", mean_value),
                                sprintf("Max: %.2f", max_value)),
       fill = c("blue", "green", "yellow", "red"), title = "Vs30 Brunelli")
#vs30 Brunelli
map("italy", fill = TRUE, col = "lightgray", bg = "lightblue", mar = c(0, 0, 0, 0))
vs30df = na.omit(itac[,c(7,6,121)])
color_scale = colorRampPalette(c("blue", "green", "yellow", "red"))(100)
col_in = as.numeric(cut(vs30df$Vs30_Brunelli_simple, breaks = 100))
point_colors = color_scale[col_in]
points(vs30df$st_longitude, vs30df$st_latitude, col = point_colors, pch = 20, cex = 1.5)

min_value <- min(vs30df$Vs30_Brunelli_simple, na.rm = TRUE)
max_value <- max(vs30df$Vs30_Brunelli_simple, na.rm = TRUE)
mean_value <- mean(vs30df$Vs30_Brunelli_simple, na.rm = TRUE)
median_value <- median(vs30df$Vs30_Brunelli_simple, na.rm = TRUE)

legend("bottomleft", legend = c(sprintf("Min: %.2f", min_value),
                                sprintf("Median: %.2f", median_value),
                                sprintf("Mean: %.2f", mean_value),
                                sprintf("Max: %.2f", max_value)),
       fill = c("blue", "green", "yellow", "red"), title = "Vs30 Brunelli")
#vs30 Brunelli+std
map("italy", fill = TRUE, col = "lightgray", bg = "lightblue", mar = c(0, 0, 0, 0))
vs30df = na.omit(itac[,c(7,6,123)])
color_scale = colorRampPalette(c("blue", "green", "yellow", "red"))(100)
col_in = as.numeric(cut(vs30df$Vs30_Brunelli_1std_simple, breaks = 100))
point_colors = color_scale[col_in]
points(vs30df$st_longitude, vs30df$st_latitude, col = point_colors, pch = 20, cex = 1.5)

min_value <- min(vs30df$Vs30_Brunelli_1std_simple, na.rm = TRUE)
max_value <- max(vs30df$Vs30_Brunelli_1std_simple, na.rm = TRUE)
mean_value <- mean(vs30df$Vs30_Brunelli_1std_simple, na.rm = TRUE)
median_value <- median(vs30df$Vs30_Brunelli_1std_simple, na.rm = TRUE)

legend("bottomleft", legend = c(sprintf("Min: %.2f", min_value),
                                sprintf("Median: %.2f", median_value),
                                sprintf("Mean: %.2f", mean_value),
                                sprintf("Max: %.2f", max_value)),
       fill = c("blue", "green", "yellow", "red"), title = "Vs30 Brunelli+std")
#vs30 Brunelli-std
map("italy", fill = TRUE, col = "lightgray", bg = "lightblue", mar = c(0, 0, 0, 0))
vs30df = na.omit(itac[,c(7,6,125)])
color_scale = colorRampPalette(c("blue", "green", "yellow", "red"))(100)
col_in = as.numeric(cut(vs30df$Vs30_Brunelli_m1std_simple, breaks = 100))
point_colors = color_scale[col_in]
points(vs30df$st_longitude, vs30df$st_latitude, col = point_colors, pch = 20, cex = 1.5)

min_value <- min(vs30df$Vs30_Brunelli_m1std_simple, na.rm = TRUE)
max_value <- max(vs30df$Vs30_Brunelli_m1std_simple, na.rm = TRUE)
mean_value <- mean(vs30df$Vs30_Brunelli_m1std_simple, na.rm = TRUE)
median_value <- median(vs30df$Vs30_Brunelli_m1std_simple, na.rm = TRUE)

legend("bottomleft", legend = c(sprintf("Min: %.2f", min_value),
                                sprintf("Median: %.2f", median_value),
                                sprintf("Mean: %.2f", mean_value),
                                sprintf("Max: %.2f", max_value)),
       fill = c("blue", "green", "yellow", "red"), title = "Vs30 Brunelli-std")
par(mfrow=c(1,2))
#Vs30Mori b
map("italy", fill = TRUE, col = "lightgray", bg = "lightblue", mar = c(0, 0, 0, 0))
vs30df = na.omit(itac[,c(7,6,122)])
color_scale = colorRampPalette(c("blue", "green", "yellow", "red"))(100)
col_in = as.numeric(cut(vs30df$Vs30_Brunelli_bilinear, breaks = 100))
point_colors = color_scale[col_in]
points(vs30df$st_longitude, vs30df$st_latitude, col = point_colors, pch = 20, cex = 1.5)

min_value <- min(vs30df$Vs30_Brunelli_bilinear, na.rm = TRUE)
max_value <- max(vs30df$Vs30_Brunelli_bilinear, na.rm = TRUE)
mean_value <- mean(vs30df$Vs30_Brunelli_bilinear, na.rm = TRUE)
median_value <- median(vs30df$Vs30_Brunelli_bilinear, na.rm = TRUE)

legend("bottomleft", legend = c(sprintf("Min: %.2f", min_value),
                                sprintf("Median: %.2f", median_value),
                                sprintf("Mean: %.2f", mean_value),
                                sprintf("Max: %.2f", max_value)),
       fill = c("blue", "green", "yellow", "red"), title = "Vs30 Mori b")


#PIE CHARTS OF CLASSIFICATION OF THE DATASET
par(mar=c(5, 4, 4, 2))
l = length(levels(as.factor(itac$Lito_ISPRA_simple)))
pie(table(as.factor(itac[,"Lito_ISPRA_simple"])), col=terrain.colors(l), radius = 0.8, cex=0.5, main = "ISPRA Classification")
legend("bottomright", legend = levels(as.factor(itac$Lito_ISPRA_simple)), inset = c(-0.05, 0), fill = terrain.colors(l), cex=0.8, ncol=2, xpd = TRUE)
l = length(levels(as.factor(itac$Lito_Bucci_simple)))
pie(table(as.factor(itac[,"Lito_Bucci_simple"])), col=terrain.colors(l), radius = 1, cex=0.5, main = "BUCCI Classification")
legend("bottomright", legend = levels(as.factor(itac$Lito_Bucci_simple)), inset = c(-0.05, 0), fill = terrain.colors(l), cex=0.8, ncol=2, xpd = TRUE)

#HEAT MAP DS2S
par(mar = c(1,1,1,1))
map("italy", fill = TRUE, col = "lightgray", bg = "lightblue", mar = c(0, 0, 0, 0))
ds2s = na.omit(itac[,c(7,6,43,63,74)]) #915 rows
color_scale = colorRampPalette(c("blue", "green", "yellow", "red"))(100)
col_in = as.numeric(cut(ds2s$dS2S_SA_pga, breaks = 100))
point_colors = color_scale[col_in]
points(ds2s$st_longitude, ds2s$st_latitude, col = point_colors, pch = 20, cex = 1.5)

min_value <- min(ds2s$dS2S_SA_pga, na.rm = TRUE)
max_value <- max(ds2s$dS2S_SA_pga, na.rm = TRUE)
mean_value <- mean(ds2s$dS2S_SA_pga, na.rm = TRUE)
median_value <- median(ds2s$dS2S_SA_pga, na.rm = TRUE)

legend("bottomleft", legend = c(sprintf("Min: %.2f", min_value),
                                sprintf("Median: %.2f", median_value),
                                sprintf("Mean: %.2f", mean_value),
                                sprintf("Max: %.2f", max_value)),
       fill = c("blue", "green", "yellow", "red"), title = "ds2s PGA")


map("italy", fill = TRUE, col = "lightgray", bg = "lightblue", mar = c(0, 0, 0, 0))
ds2s = na.omit(itac[,c(7,6,43,63,74)]) #915 rows
color_scale = colorRampPalette(c("blue", "green", "yellow", "red"))(100)
col_in = as.numeric(cut(ds2s$dS2S_SA_T1_000, breaks = 100))
point_colors = color_scale[col_in]
points(ds2s$st_longitude, ds2s$st_latitude, col = point_colors, pch = 20, cex = 1.5)

min_value <- min(ds2s$dS2S_SA_T1_000, na.rm = TRUE)
max_value <- max(ds2s$dS2S_SA_T1_000, na.rm = TRUE)
mean_value <- mean(ds2s$dS2S_SA_T1_000, na.rm = TRUE)
median_value <- median(ds2s$dS2S_SA_T1_000, na.rm = TRUE)

legend("bottomleft", legend = c(sprintf("Min: %.2f", min_value),
                                sprintf("Median: %.2f", median_value),
                                sprintf("Mean: %.2f", mean_value),
                                sprintf("Max: %.2f", max_value)),
       fill = c("blue", "green", "yellow", "red"), title = "ds2s T=1s")


map("italy", fill = TRUE, col = "lightgray", bg = "lightblue", mar = c(0, 0, 0, 0))
ds2s = na.omit(itac[,c(7,6,43,63,74)]) #915 rows
color_scale = colorRampPalette(c("blue", "green", "yellow", "red"))(100)
col_in = as.numeric(cut(ds2s$dS2S_SA_T5_000, breaks = 100))
point_colors = color_scale[col_in]
points(ds2s$st_longitude, ds2s$st_latitude, col = point_colors, pch = 20, cex = 1.5)

min_value <- min(ds2s$dS2S_SA_T5_000, na.rm = TRUE)
max_value <- max(ds2s$dS2S_SA_T5_000, na.rm = TRUE)
mean_value <- mean(ds2s$dS2S_SA_T5_000, na.rm = TRUE)
median_value <- median(ds2s$dS2S_SA_T5_000, na.rm = TRUE)

legend("bottomleft", legend = c(sprintf("Min: %.2f", min_value),
                                sprintf("Median: %.2f", median_value),
                                sprintf("Mean: %.2f", mean_value),
                                sprintf("Max: %.2f", max_value)),
       fill = c("blue", "green", "yellow", "red"), title = "ds2s T=5s")

legend("bottomleft", legend = c("Low", "", "", "High"),
       fill = c("blue", "green","yellow", "red"), title="ds2s T=5s")



