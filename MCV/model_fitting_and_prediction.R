
#Load required libraries

library(spBayes)
library(parallel)
library(raster)
library(rgdal)
library(fields)


#loading the data
vaxdata <- read.csv("filename.csv", header=TRUE) #File containing processed cluster-level vaccination data
								 #Each row contains the DHS cluster number, and for each age group,
								 #the number of children surveyed and the numbers who were vaccinated 

vaxcov <- read.csv("filename.csv",header=TRUE)	 #File containing cluster-level covariate information

#Delete clusters where no child was surveyed
zero.clust <- which(vaxdata$Age9_59_Tot<=1) 	#Age9_59_Tot contains total number of children for the
if (length(zero.clust)>0){				#9-59 month age group
  vaxdata <- vaxdata[-zero.clust,]
  vaxcov <- vaxcov[-zero.clust,]
}

#Extract vaccination data for age group 9-59 months
Numvacc    <- vaxdata$Age9_59_Vax
weights    <- vaxdata$Age9_59_Tot

#Geographical coordinates of clusters
coords     <- cbind(vaxdata$lon,vaxdata$lat)



#############################################################################################
#For cross-validation. This can be repeated by creating a loop.
#This part can be excluded if not needed. 
ll <- length(weights)
nc <- (10/100)*ll         #Take 10% for validation - this was used to estimate the cov prob.
samp.c <- sample(1:nrow(coords), nc, replace=FALSE)
coords.nc 	<- coords[samp.c,]
Numvacc.nc	<- Numvacc[samp.c]
weights.nc	<- weights[samp.c]
vaxcov.nc	<- vaxcov[samp.c,]

#Use the rest for model estimation
coords  <- coords[-samp.c,]
Numvacc <- Numvacc[-samp.c] 
weights <- weights[-samp.c]
vaxcov  <- vaxcov[-samp.c,]
#############################################################################################

set.seed(500)

#Covariates  
VAR1_NAME	   <- vaxcov$VAR1_NAME
VAR2_NAME      <- vaxcov$VAR2_NAME
VAR3_NAME      <- vaxcov$VAR3_NAME
VAR4_NAME      <- vaxcov$VAR4_NAME

#Formula
form 	 <- Numvacc ~  VAR1_NAME  + VAR2_NAME + VAR3_NAME + VAR4_NAME
form.2 <- (Numvacc/weights) ~ VAR1_NAME  + VAR2_NAME + VAR3_NAME + VAR4_NAME  

#Initial values of some parameters
fit 			<- glm(form.2, weights=weights, family="binomial")
beta.starting 	<- coefficients(fit)
beta.tuning 	<- t(chol(vcov(fit)))

n.batch 		<- 1000         ##100*1000 equals run length
batch.length 	<- 100
n.samples 		<- n.batch*batch.length
burn.in 		<- 0.1*n.samples + 1

#Sample knot locations using a space-filling design 
n.knots 	<- 100   #Can increase to 100 if need be
bb       	<- cover.design(coords, n.knots)  
knots 	<- as.matrix(bb$design)


#COMPARE DIFFERENT COV MODELS LATER!
mod <- spGLM(form, family="binomial", weights = weights, 
            coords=coords, starting=list("beta"=beta.starting, "phi"=1,"sigma.sq"=1, "w"=0, "nu"=1),
            tuning=list("beta"=beta.tuning, "phi"=0.5, "sigma.sq"=0.5, "w"=0.5, "nu" = 0.1),
            priors=list("beta.Normal"=list(c(rep(0,5)),c(rep(1000,5))), "phi.Unif"=c(0.475, 3.3), "sigma.sq.IG"=c(2, 1), "nu.Unif"=c(0.1,1)),
            amcmc=list("n.batch"=n.batch, "batch.length"=batch.length, "accept.rate"=0.43),
            knots=knots,cov.model="exponential", verbose=TRUE, n.report=10) # 0.475 - 700 km,  3.3 - 100 km


#Save model for analysis
save(mod, file="model_name.rda")

#Read in prediction covariates 
VAR1_NAME       <- raster(".....tif"); VAR2_NAME       <- raster(".....tif")
VAR3_NAME       <- raster(".....tif"); VAR4_NAME       <- raster(".....tif")

#Convert covariate rasters to vectors
VAR1_NAME       <- getValues(VAR1_NAME); VAR2_NAME       <- getValues(VAR2_NAME)
VAR3_NAME       <- getValues(VAR3_NAME); VAR4_NAME       <- getValues(VAR4_NAME)

#Obtain prediction coordinates from one of the raster files
VAR_NAME.p <- raster(".....tif")
Pred_grid2 <- coordinates(VAR_NAME.p)

#Combine grid and covariates
pred.dat <- cbind(Pred_grid2, VAR1_NAME, VAR2_NAME, VAR3_NAME, VAR4_NAME)

#Extract grid cells with missing values and exclude these from the prediction step
ind <- apply(pred.dat, 1, function(x) any(is.na(x)))

miss    <- which(ind==TRUE)
nonmiss <- which(ind==FALSE)

pred.dat.1 <- pred.dat[nonmiss, ]
pred.dat.1 <- data.frame(pred.dat.1)

splitsize <- 20000 ## CHANGE HERE to adjust size of subsets

#create appropriate samples
big <- dim(pred.dat.1)[[1]]
nsplits <- round(big/splitsize)
allspl <- 1:big
splits <- cut(allspl, nsplits)

#-------------------------------------------------------------#
#Do this if the prediction locations are to be sampled randomly
#s.samp <- sample(allspl, length(allspl), replace=F)
#splits <- cut(s.samp, nsplits) ## CHECK NUMBER OF PARTITIONS
#-------------------------------------------------------------#

# the split_dat object contains the subsets of prediction data
split_dat <- split(pred.dat.1, splits)

#Define quantiles to be calculated
ff1=function(x) quantile(x,0.025)
ff2=function(x) quantile(x,0.975)

# initialize the cluster
cl <- makeCluster(detectCores()-1) # use the maximum number of processors available
# tell the cluster to load the libraries required for processing

clusterEvalQ(cl,library("spBayes"))

# this would need to be spBayes, for example
# send data/objects needed for analyses to the different cores
clusterExport(cl, varlist=c("mod","split_dat","n.samples","burn.in","filePathData","ff1","ff2"))

## Main processing Loop
out <- clusterApply(cl, 1:length(split_dat), function(i){
  print(i)
  pred_dat    <- split_dat[[i]] # get the subset of data for prediction
  pred.coords <- pred_dat[,1:2]
  #pred.covars <- pred_dat[,3:ncol(pred_dat)]
  pred.covars <- data.matrix(cbind(rep(1,nrow(pred.coords)),pred_dat[,3:ncol(pred_dat)])) #Note const added
  mod.pred    <- spPredict(mod, pred.coords, pred.covars, start=burn.in, n.samples, thin=5, 
                           verbose=FALSE, n.report=100) #change thin
  y.hat.median <- apply(mod.pred$p.y.predictive.samples, 1, median)
  y.hat.low    <- apply(mod.pred$p.y.predictive.samples, 1, ff1)
  y.hat.up     <- apply(mod.pred$p.y.predictive.samples, 1, ff2)
  y.hat.sd     <- apply(mod.pred$p.y.predictive.samples, 1, sd)
  y.hat.mean     <- apply(mod.pred$p.y.predictive.samples, 1, mean)
  prediction  <- cbind(y.hat.median, y.hat.low, y.hat.up, y.hat.sd, y.hat.mean)
  
  
  ## CHANGE filepath HERE
  # write the predictions to separate files
  write.csv(prediction, file=paste0(filePathData,"file_name",i,".csv"), row.names=F, col.names=F)
  # post-processing step needed to merge these files
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
#Do this if prrediction data are resampled
#Reorder sampled rows
#combine.ord <- combine[order(s.samp),]
#-----------------------------------------------#

#Convert to raster
combine.ord <- combine

#Add missing data
pred.out <- matrix(0, nrow=length(miss)+length(nonmiss), ncol=5)

pred.out[miss,] <- NA
pred.out[nonmiss,1] <- combine.ord[,1]
pred.out[nonmiss,2] <- combine.ord[,2]
pred.out[nonmiss,3] <- combine.ord[,3]
pred.out[nonmiss,4] <- combine.ord[,4]
pred.out[nonmiss,5] <- combine.ord[,5]

#Convert to rasters
#Sample raster
qw <- raster("filename.tif"))

#Median
median.pred 	 	<- raster(qw)
values(median.pred) 	<- pred.out[,1]

#lower CI
low.pred 			<- raster(qw)
values(low.pred) 		<- pred.out[,2]

#upper CI
up.pred 			<- raster(qw)
values(up.pred) 		<- pred.out[,3]

#sd
sd.pred 			<- raster(qw)
values(sd.pred) 		<- pred.out[,4]

#Mean
mean.pred 			<- raster(qw)
values(mean.pred) 	<- pred.out[,5]

writeRaster(median.pred, "filename.tif", overwrite=TRUE)
writeRaster(low.pred,"filename.tif", overwrite=TRUE)
writeRaster(up.pred,"filename.tif", overwrite=TRUE)
writeRaster(sd.pred,"filename.tif", overwrite=TRUE)
writeRaster(mean.pred,"filename.tif", overwrite=TRUE)  

#Source code for further analysis
source("subcode1.R")

#Source code2 for cross-validation calculations
source("subcode2.R")





