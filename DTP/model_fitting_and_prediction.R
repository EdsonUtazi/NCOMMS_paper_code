
#Load required libraries
library(spBayes)
library(raster)
library(parallel)
library(fields)
library(rgdal)


#Load data
vaxdata <- read.csv("filename.csv", header=TRUE)#File containing processed cluster-level vaccination data
							#Each row contains the DHS cluster number, and for each dose, the number of 
							#children surveyed at the cluster and the numbers who received the dose
 
vaxcov <- read.csv("filename.csv", header=TRUE)	#File containing cluster-level covariate information

#Delete clusters where no child was surveyed 
zero.clust <- which(vaxdata$dtp1_Tot==0)  #dtp1_Tot conatins total numbers of children
if (length(zero.clust)>0){   			#in each cluster 
  vaxdata <- vaxdata[-zero.clust,]
  vaxcov <- vaxcov[-zero.clust,]
}


q <- 3  #No of dependent variables, i.e, DTP1-3

#Assign numbers vaccinated and numbers surveyed for each dose to new objects
Numvacc1    <- vaxdata$dtp1_Vax; weights1    <- vaxdata$dtp1_Tot
Numvacc2    <- vaxdata$dtp2_Vax; weights2    <- vaxdata$dtp2_Tot
Numvacc3    <- vaxdata$dtp3_Vax; weights3    <- vaxdata$dtp3_Tot

#Extract and delete data points with anomalies. E.g. where probability of receiving 
#DTP2 is greater than that of DTP1
p1   		<- Numvacc1/weights1; p2 <- Numvacc2/weights2; p3 <- Numvacc3/weights3
dif1 		<- p1-p2; dif2 <- p1-p3; dif3 <- p2-p3
del  		<- unique(c(which(dif1<0),which(dif2<0),which(dif3<0)))
vaxcov 	<- vaxcov[-del,]
Numvacc1 	<- Numvacc1[-del]; Numvacc2 <- Numvacc2[-del]; Numvacc3 <- Numvacc3[-del]
weights1 	<- weights1[-del]; weights2 <- weights2[-del]; weights3 <- weights3[-del]


#Load file with the coordinates (lat-lon) of the survey clusters
#This could be included in any of the data files
f.coords <- read.csv("....")
coords   <- cbind(f.coords$LONGNUM,f.coords$LATNUM)
coords   <- coords[-del,]

#############################################################################################
#For cross-validation. This can be repeated by creating a loop.
#This part can be excluded if not needed. 
ll 	<- length(weights1)
nc 	<- (10/100)*ll         #Take 10% for validation - this was used to estimate the cov prob.
samp.c <- sample(1:nrow(coords), nc, replace=FALSE)
coords.nc 	<- coords[samp.c,]
vaxcov.nc	<- vaxcov[samp.c,]

Numvacc1.nc	<- Numvacc1[samp.c]; weights1.nc	<- weights1[samp.c]
Numvacc2.nc	<- Numvacc2[samp.c]; weights2.nc	<- weights2[samp.c]
Numvacc3.nc	<- Numvacc3[samp.c]; weights3.nc	<- weights3[samp.c]

#Use the rest for model estimation
coords  <- coords[-samp.c,]
vaxcov  <- vaxcov[-samp.c,]
##############################################################################################

set.seed(100)
#Covariates  
VAR1_NAME	   <- vaxcov$VAR1_NAME
VAR2_NAME      <- vaxcov$VAR2_NAME
VAR3_NAME      <- vaxcov$VAR3_NAME
VAR4_NAME      <- vaxcov$VAR4_NAME
VAR5_NAME      <- vaxcov$VAR5_NAME

#Prepare data for initializating the algorithm
xx.1 = xx.2 = xx.3 = cbind(1, VAR1_NAME, VAR2_NAME, VAR3_NAME, VAR4_NAME, VAR5_NAME)   #Note order of covs
x <- mkMvX(list(xx.1,xx.2,xx.3))

Numvacc = weights = numeric(nrow(xx.1)*q) # q = no of dep vars
Numvacc[seq(1,nrow(xx.1)*q, q)] <- Numvacc1; weights[seq(1,nrow(xx.1)*q, q)] <- weights1
Numvacc[seq(2,nrow(xx.1)*q, q)] <- Numvacc2; weights[seq(2,nrow(xx.1)*q, q)] <- weights2 
Numvacc[seq(3,nrow(xx.1)*q, q)] <- Numvacc3; weights[seq(3,nrow(xx.1)*q, q)] <- weights3

#Formula
form <- list(Numvacc1 ~  VAR1_NAME  + VAR2_NAME + VAR3_NAME + VAR4_NAME + VAR5_NAME,  
             Numvacc2 ~  VAR1_NAME  + VAR2_NAME + VAR3_NAME + VAR4_NAME + VAR5_NAME,
             Numvacc3 ~  VAR1_NAME  + VAR2_NAME + VAR3_NAME + VAR4_NAME + VAR5_NAME)

form.2 <- (Numvacc/weights) ~ VAR1_NAME  + VAR2_NAME + VAR3_NAME + VAR4_NAME + VAR5_NAME

#Initial values of some parameters
fit <- glm((Numvacc/weights)~x-1, weights=weights, family="binomial")
beta.starting <- coefficients(fit)
beta.tuning   <- t(chol(vcov(fit)))
A.starting    <- diag(1,q)[lower.tri(diag(1,q), TRUE)]

n.batch      <- 1000         #100*1000 equals run length
batch.length <- 100
n.samples    <- n.batch*batch.length
burn.in      <- 0.1*n.samples + 1

#Sample knot locations using a space-filling design  
n.knots <- 200
bb      <- cover.design(coords, n.knots)  
knots   <- as.matrix(bb$design)

#Model
mod <- spMvGLM(form, family="binomial", weights = cbind(weights1,weights2,weights3), 
              coords=coords, starting=list("beta"=beta.starting, "phi"=rep(0.5,q), "w"=0, "A"=A.starting),
              tuning=list("beta"=beta.tuning, "phi"=rep(0.5,q), "w"=0.5, "nu" = 0.1, "A"=rep(0.1,length(A.starting))),
              priors=list("beta.Normal"=list(rep(rep(0,6),q), rep(rep(1000,6),q)), "phi.Unif"=list(rep(0.166,q), rep(1.11,q)), 
                          "sigma.sq.IG"=c(2, 1), "nu.Unif"=c(0,1), "K.IW"=list(q, diag(1,q))),
              amcmc=list("n.batch"=n.batch, "batch.length"=batch.length, "accept.rate"=0.43),
              knots=knots, cov.model="exponential", verbose=TRUE, n.report=10)  #1.11 - 300 km, 0.166 - 2000 km

#Save model for further analysis
save(mod, file="model_name.rda")

#Read in prediction covariates 
VAR1_NAME       <- raster(".....tif"); VAR2_NAME       <- raster(".....tif")
VAR3_NAME       <- raster(".....tif"); VAR4_NAME       <- raster(".....tif")
VAR5_NAME       <- raster(".....tif")

#Convert covariate rasters to vectors
VAR1_NAME       <- getValues(VAR1_NAME); VAR2_NAME       <- getValues(VAR2_NAME)
VAR3_NAME       <- getValues(VAR3_NAME); VAR4_NAME       <- getValues(VAR4_NAME)
VAR5_NAME       <- getValues(VAR5_NAME)

#Obtain prediction coordinates from one of the raster files
VAR_NAME.p <- raster(".....tif")
Pred_grid2 <- coordinates(VAR_NAME.p)

#Combine grid and covariates
pred.dat <- cbind(Pred_grid2, VAR1_NAME, VAR2_NAME, VAR3_NAME, VAR4_NAME, VAR5_NAME)

#Extract grid cells with missing values and exclude these from the prediction step
ind <- apply(pred.dat, 1, function(x) any(is.na(x)))

miss    <- which(ind==TRUE)
nonmiss <- which(ind==FALSE)

pred.dat.1 <- pred.dat[nonmiss, ]
pred.dat.1 <- data.frame(pred.dat.1)

splitsize <- 10000 ## Change to adjust size of subsets

#create appropriate samples
big 	  <- dim(pred.dat.1)[[1]]
nsplits <- round(big/splitsize)
allspl  <- 1:big
splits  <- cut(allspl, nsplits)

#-------------------------------------------------------------#
#Do this if the prediction locations are to be sampled randomly
#s.samp <- sample(allspl, length(allspl), replace=F)
#splits <- cut(s.samp, nsplits) 
#-------------------------------------------------------------#

#split_dat contains the subsets of prediction data
split_dat <- split(pred.dat.1, splits)

#Define quantiles to be calculated
ff1=function(x) quantile(x,0.025)
ff2=function(x) quantile(x,0.975)

# initialize the cluster
cl <- makeCluster(detectCores()-1) # use the maximum number of processors available
#Load the libraries required for processing
clusterEvalQ(cl, library(spBayes))
  
#Send data/objects needed for analyses to the different cores
#filePathData - directory/folder where files are saved
clusterExport(cl, varlist=c("mod","split_dat","n.samples","burn.in","filePathData","ff1","ff2","q"))

## Main processing Loop
out <- clusterApply(cl, 1:length(split_dat), function(i){
  print(i)
  pred_dat    <- split_dat[[i]] # get the subset of data for prediction
  pred.coords <- pred_dat[,1:2]
  pred.covars.a <- data.matrix(cbind(rep(1,nrow(pred.coords)),pred_dat[,3:ncol(pred_dat)])) #Note constant term added
  pred.covars <- mkMvX(list(pred.covars.a,pred.covars.a,pred.covars.a))
  
  mod.pred    <- spPredict(mod, pred.coords, pred.covars, start=burn.in, n.samples, thin=1, 
                           verbose=FALSE, n.report=100) #change thin
  
  y.hat.median.1 <- apply(mod.pred$p.y.predictive.samples[seq(1,nrow(pred.coords)*q, q),], 1, median)
  y.hat.low.1    <- apply(mod.pred$p.y.predictive.samples[seq(1,nrow(pred.coords)*q, q),], 1, ff1)
  y.hat.up.1     <- apply(mod.pred$p.y.predictive.samples[seq(1,nrow(pred.coords)*q, q),], 1, ff2)
  y.hat.sd.1     <- apply(mod.pred$p.y.predictive.samples[seq(1,nrow(pred.coords)*q, q),], 1, sd)
  y.hat.mean.1     <- apply(mod.pred$p.y.predictive.samples[seq(1,nrow(pred.coords)*q, q),], 1, mean)
  
  y.hat.median.2 <- apply(mod.pred$p.y.predictive.samples[seq(2,nrow(pred.coords)*q, q),], 1, median)
  y.hat.low.2    <- apply(mod.pred$p.y.predictive.samples[seq(2,nrow(pred.coords)*q, q),], 1, ff1)
  y.hat.up.2     <- apply(mod.pred$p.y.predictive.samples[seq(2,nrow(pred.coords)*q, q),], 1, ff2)
  y.hat.sd.2     <- apply(mod.pred$p.y.predictive.samples[seq(2,nrow(pred.coords)*q, q),], 1, sd)
  y.hat.mean.2     <- apply(mod.pred$p.y.predictive.samples[seq(2,nrow(pred.coords)*q, q),], 1, mean)
  
  y.hat.median.3 <- apply(mod.pred$p.y.predictive.samples[seq(3,nrow(pred.coords)*q, q),], 1, median)
  y.hat.low.3    <- apply(mod.pred$p.y.predictive.samples[seq(3,nrow(pred.coords)*q, q),], 1, ff1)
  y.hat.up.3     <- apply(mod.pred$p.y.predictive.samples[seq(3,nrow(pred.coords)*q, q),], 1, ff2)
  y.hat.sd.3     <- apply(mod.pred$p.y.predictive.samples[seq(3,nrow(pred.coords)*q, q),], 1, sd)
  y.hat.mean.3     <- apply(mod.pred$p.y.predictive.samples[seq(3,nrow(pred.coords)*q, q),], 1, mean)
  
  
  prediction  <- cbind(y.hat.median.1, y.hat.low.1, y.hat.up.1, y.hat.sd.1, y.hat.mean.1, 
                       y.hat.median.2, y.hat.low.2, y.hat.up.2, y.hat.sd.2, y.hat.mean.2,
                       y.hat.median.3, y.hat.low.3, y.hat.up.3, y.hat.sd.3, y.hat.mean.3)
  
  
  ## CHANGE filepath HERE
  # write the predictions to separate files
  write.csv(prediction, file=paste0(filePathData,"file_name",i,".csv"), row.names=F, col.names=F)
  # post-processing step needed to merge these files
  rm(mod.pred)
  rm(prediction)
  return(1)
})

stopCluster(cl) # close the cluster

#Combine files after prediction
combine <- data.frame()
for (i in 1:nsplits){
  split.dat <- read.csv(file=paste0(filePathData,"file_name",i,".csv")) 
  combine <- rbind(combine, split.dat)
}

#-----------------------------------------------#
#Do this if prediction data are resampled
#Reorder sampled rows
#combine.ord <- combine[order(s.samp),]
#-----------------------------------------------#

combine.ord <- combine

#Add missing data
pred.out <- matrix(0, nrow=length(miss)+length(nonmiss), ncol=5*q)

pred.out[miss,] <- NA
pred.out[nonmiss,1] <- combine.ord[,1]; pred.out[nonmiss,6] <- combine.ord[,6]; pred.out[nonmiss,11] <- combine.ord[,11] 
pred.out[nonmiss,2] <- combine.ord[,2]; pred.out[nonmiss,7] <- combine.ord[,7]; pred.out[nonmiss,12] <- combine.ord[,12]
pred.out[nonmiss,3] <- combine.ord[,3]; pred.out[nonmiss,8] <- combine.ord[,8]; pred.out[nonmiss,13] <- combine.ord[,13]
pred.out[nonmiss,4] <- combine.ord[,4]; pred.out[nonmiss,9] <- combine.ord[,9]; pred.out[nonmiss,14] <- combine.ord[,14]
pred.out[nonmiss,5] <- combine.ord[,5]; pred.out[nonmiss,10] <- combine.ord[,10]; pred.out[nonmiss,15] <- combine.ord[,15] 

#Convert to rasters
#Sample raster
qw <- raster("filename.tif"))

#Median
median.pred.1 <- raster(qw); values(median.pred.1) <- pred.out[,1]
median.pred.2 <- raster(qw); values(median.pred.2) <- pred.out[,6]
median.pred.3 <- raster(qw); values(median.pred.3) <- pred.out[,11]

#lower CI
low.pred.1 <- raster(qw); values(low.pred.1) <- pred.out[,2]
low.pred.2 <- raster(qw); values(low.pred.2) <- pred.out[,7]
low.pred.3 <- raster(qw); values(low.pred.3) <- pred.out[,12]

#upper CI
up.pred.1 <- raster(qw); values(up.pred.1) <- pred.out[,3]
up.pred.2 <- raster(qw); values(up.pred.2) <- pred.out[,8]
up.pred.3 <- raster(qw); values(up.pred.3) <- pred.out[,13]

#sd
sd.pred.1 <- raster(qw); values(sd.pred.1) <- pred.out[,4]
sd.pred.2 <- raster(qw); values(sd.pred.2) <- pred.out[,9]
sd.pred.3 <- raster(qw); values(sd.pred.3) <- pred.out[,14]

#mean
mean.pred.1 <- raster(qw); values(mean.pred.1) <- pred.out[,5]
mean.pred.2 <- raster(qw); values(mean.pred.2) <- pred.out[,10]
mean.pred.3 <- raster(qw); values(mean.pred.3) <- pred.out[,15]


writeRaster(median.pred.1,"filename.tif", overwrite=TRUE)
writeRaster(low.pred.1,"filename.tif", overwrite=TRUE)
writeRaster(up.pred.1,"filename.tif", overwrite=TRUE)
writeRaster(sd.pred.1,"filename.tif", overwrite=TRUE)
writeRaster(mean.pred.1,"filename.tif", overwrite=TRUE)

writeRaster(median.pred.2,"filename.tif", overwrite=TRUE)
writeRaster(low.pred.2,"filename.tif", overwrite=TRUE)
writeRaster(up.pred.2,"filename.tif", overwrite=TRUE)
writeRaster(sd.pred.2,"filename.tif", overwrite=TRUE)
writeRaster(mean.pred.2,"filename.tif", overwrite=TRUE)

writeRaster(median.pred.3,"filename.tif", overwrite=TRUE)
writeRaster(low.pred.3,"filename.tif", overwrite=TRUE)
writeRaster(up.pred.3,"filename.tif", overwrite=TRUE)
writeRaster(sd.pred.3,"filename.tif", overwrite=TRUE)
writeRaster(mean.pred.3,"filename.tif", overwrite=TRUE)



#Source code for further analysis
source("subcode1.R")

#Source code2 for cross-validation calculations
source("subcode2.R")





