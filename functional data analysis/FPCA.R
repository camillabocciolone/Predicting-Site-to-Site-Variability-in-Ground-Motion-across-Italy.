rm(list=ls())
graphics.off()
library(fda)

itac = read.csv("FINAL_DATASET.csv") 
itaca = itac[,43:79] #consider only dS2S terms
rowsva <- which(complete.cases(itaca))
NArows <- setdiff(1:916, rowsva) #save NA rows for dS2S to rebuild dataset with scores
itaca = na.omit(itaca) #dim(itaca)=879,37
Tn = (c(0,0.01,0.025,0.04,0.05,0.07,0.1,0.15,0.2,0.25,0.3,0.35,0.40,0.45,0.5,0.6,0.7,0.75,0.8,0.9,1,1.2,1.4,1.6,1.8,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10))
matplot(Tn,t(itaca),type='l',col='grey',lty=1,main='Site-to-site SA',xlab='T [s]',ylab='Spectral Acceleration cm/s^2')
m<- apply(itaca, 2, median)
lines(Tn, as.numeric(m), col = "red", lwd = 2)
#Smoothing
bs <- c(seq(0, 3, by = 0.2), seq(3.5, 10, by = 1)) #Create more knots for low periods
basis <- create.bspline.basis(rangeval = c(0, 10), breaks = bs, norder = 6)
functionalPar <- fdPar(fdobj=basis, Lfdobj=4, lambda=1e-9) #already check lamda 1e-9 is fine 
Xss <- smooth.basis(Tn, t(itaca), functionalPar) #smooth data 

#fpca
fpca <- pca.fd(Xss$fd,nharm=5,centerfns=TRUE) #setting the number of harmonics (PCs to consider) = 5 (actually not specifying only returns 2 cols of scores)
plot(fpca$values[1:27],xlab='j',ylab='Eigenvalues') #consider 27 eigenvalues 
plot(cumsum(fpca$values)[1:27]/sum(fpca$values),xlab='j',ylab='CPV',ylim=c(0.8,1)) #suggests to use 3 FPCs

scrs <- matrix(NA, nrow = 916, ncol = 5)
scrs[rowsva, ] <- fpca$scores[,1:5]
colnames(scrs) <- c('FPC1', 'FPC2', 'FPC3', 'FPC4', 'FPC5')
itac = cbind(scrs,itac)
View(itac)
#write.csv(itac, "DATASET_WITH_SCORES.csv", row.names = FALSE)

#LM on FPCs
itaca = itac[,c(1,2,3,9,10,11,12,20,21,22,23,46:82,122,124,126,132,136)]
attach(itaca)

#lm FPC1
#NA should be removed
itaca=itaca[!is.na(itaca$Vs30_Brunelli_simple),]
itaca=itaca[!is.na(itaca$slope50_simple),]
m=lm(FPC1 ~ sqrt(Vs30_Brunelli_simple) + sqrt(slope50_simple), itaca)
summary(m)

#Diagnostic
plot(m, which=1 ) 
res = m$residuals/summary(m)$sigma
qqnorm(res)
qqline(res, col = 'red')
shapiro.test(res)
#other plots
par(mfrow=c(2,2))
plot(m)
shapiro.test(residuals(m))
par(mfrow=c(1,1))

#Influential points
#Leverages <--
lev=hat(model.matrix(m))
p=m$rank
n=dim(itaca)[1]
pts_lev=lev[which(lev>2*p/n)]
m=lm(FPC1 ~ sqrt(Vs30_Brunelli_simple) + sqrt(slope50_simple),itaca, subset=lev<2*p/n)
summary(m)
#Cook distance
Cdist=cooks.distance(m)
i=which(Cdist>4/(n-p))
which_Cdist=Cdist[i]
id_to_keep=!(1:n%in%i)
m=lm(FPC1 ~ sqrt(Vs30_Brunelli_simple) + sqrt(slope50_simple),itaca[id_to_keep,])
summary(m) 

m=lm(FPC1 ~ sqrt(Vs30_Brunelli_simple) + sqrt(slope50_simple), itaca, subset=lev<2*p/n)
summary(m) 
plot(m, which=1 ) 
res = m$residuals/summary(m)$sigma
qqnorm(res)
qqline(res, col = 'red')
shapiro.test(res)

