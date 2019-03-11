
#Load these libraries
library(Metrics)
library(plyr)
library(xtable)
library(ggplot2)
library(reshape2)
library(MASS)
library(bestglm)
require(dplyr) 
library(car)
library(gtools)

#Set working directory
setwd("......")

#Read in processed data sets
covs <- read.csv("covariate_data.csv", header = TRUE) #Covariate data
vaxdat <- read.csv("Measles_vax_data.csv", header=TRUE) #Vaccination data

#Delete Class_rainfed_MaJORITY - contains many missing values
#covs <- covs[,-10]

names(covs)
names(vaxdat)

#Get the coords of the clusters
vaxcoord <- covs[,3:4]
head(vaxcoord)

#Vaccination variables
TotChild <- vaxdat$Age9_59_Tot # Total number of children surveyed in each cluster
TotVax   <- vaxdat$Age9_59_Vax # Number of children vaccinated

#Extract columns with covariates
covars <- covs[,6:ncol(covs)]

############################################################
#Histograms

par(mfrow=c(2,5))
for (i in 1:10){
hist(covars[,i], main = names(covars)[i])
}

windows()
par(mfrow=c(2,5))
for (i in 11:20){
hist(covars[,i], main = names(covars)[i])
}

windows()
par(mfrow=c(2,5))
for (i in 21:24){
hist(covars[,i], main = names(covars)[i])
}

#par(mfrow=c(1,3))
#for (i in 31:33){
#hist(covars[,i], main = names(covars)[i])
#}

##############################################################


#Take logs of heavily skewed covariates
llog <- c(1,4,7,8,9,11:15,17,18,20,21,22,24)
for (i in 1:length(llog)){
covars[,llog[i]] <- log(covars[,llog[i]] + 0.5)
names(covars)[llog[i]] <- paste0("l_",names(covars)[llog[i]])
}

##############################################################
#Linearity checks
pp = TotVax/TotChild
pp[pp==0] <- 0.001
pp[pp==1] <- 0.99
par(mfrow=c(2,5))
for (i in 1:10){
plot(covars[,i],logit(pp), main = names(covars)[i], ylab="logit(Probability)", xlab="")
abline(lm(logit(pp)~covars[,i]))
}

windows()
par(mfrow=c(2,5))
for (i in 11:20){
plot(covars[,i],logit(pp+0.0001), main = names(covars)[i], ylab="logit(Probability)", xlab="")
abline(lm(logit(pp+0.0001)~covars[,i]))
}

windows()
par(mfrow=c(2,5))
for (i in 21:24){
plot(covars[,i],logit(pp+0.0001), main = names(covars)[i], ylab="logit(Probability)", xlab="")
abline(lm(logit(pp+0.0001)~covars[,i]))
}
###############################################################


#-------------------------Covariate selection starts here--------------------#

#Single covariate models first
Data <- cbind(TotChild, TotVax, covars)

#Delete clusters where TotChild is zero
zero.clust <- which(Data$TotChild==0)
if (length(zero.clust)>0) Data <- Data[-zero.clust,]

Data.1 <- Data[,-1]  #Delete TotChild for model
covlist <- names(Data.1)
covlist=covlist[-grep("TotVax", covlist)]


##Single Covariate models
AICs=dev=nulldev=n.par=r2=pr2=iter=numeric()
model=response=character()
n.iter=5
propsub=0.8
resp="cbind(TotVax, TotChild-TotVax)"
resp1 = "TotVax"
weights = "TotChild"

for(i in 1:n.iter){
  subind=sample(1:nrow(Data), propsub*nrow(Data), replace=F)
  subdat=Data[subind,]
  preddat=Data[-subind,]
  ## the null model, for comparison
  nullmod=glm(formula(paste(resp, "~1")), data=Data, family = binomial(logit))
  for(j in covlist){
    form=formula(paste(resp, "~", j))
    pmod=glm(form, data=subdat, family = binomial(logit))
    fullmod=glm(form, data=Data, family = binomial(logit))
    ## pr2 checks the predictive power of the model against a 'new' subset of the data
    pr2=c(pr2,cor(pmod$fitted, subdat[[resp1]]/subdat[[weights]])^2) 
    ## AIC for the model on the full data
    AICs=c(AICs, AIC(fullmod))
    dev=c(dev, deviance(fullmod))
    n.par=c(n.par,1)
    iter=c(iter,i)
    r2=c(r2,cor(fullmod$fitted, Data[[resp1]]/Data[[weights]])^2) 
    model=c(model,j)
    response=c(response, resp1)
    #print(paste(i, j, sep="X"))
    }
}
op=data.frame(model, response, dev, AICs, r2, pr2, n.par )
singles=ddply(op, .(model), summarise,
              model=model[1],
              response=response[1],
              dev=dev[1],
              pr2=mean(pr2),
              AICs=AICs[1],
              r2=r2[1])

### Add deviance reduction
singles$devred=singles$dev/deviance(nullmod)

## order by AICs (but could be devred or pr2)
#singles=singles[order(singles$AICs),]

## order by pr2
singles=singles[order(singles$AIC, decreasing=TRUE),]
singles

keepsingles=singles   #can decide which single covariates to keep at this stage

#Ranks
singles <- within(singles, ranks <- 1:nrow(singles))

#Covariates and ranks
cov.rank <- data.frame(covariate=as.character(singles$model), ranks = singles$ranks)

#Output file 1
#write.csv(singles, file="Single.Cov.Models.csv")
#####End of single covariate models

#################################################################################
##Detection of multicollinearity through correlations between covariates and VIF analysis

#Correlations
#Determine correlations between the covariates and extract highly correlated pairs
#for screening and elimination. 
#Flag covariate pairs with correlations >= 0.8   #change to a higher value if necessary

corrs <- cor(covars[,-1])
bigcors <- matrix(0,1, 2) 
for (i in 2: nrow(corrs)){
  for (j in 1:(i-1)){
  if (abs(corrs[i,j])>=0.8) bigcors <- rbind(bigcors, c(i,j))
  }
}

bigcors <- bigcors[-1,]
bigcors.dat <- matrix(0, nrow(bigcors), 3)
for (i in 1:nrow(bigcors)){
  bigcors.dat[i,] <- c(rownames(corrs)[bigcors[i,1]],colnames(corrs)[bigcors[i,2]], 
                       round(corrs[bigcors[i,1],bigcors[i,2]],3))
}

#Pairs of covariates with high correlations
bigcors.dat

#Select between pairs of covariates using their ranks
all.covs <-rownames(corrs)
for (i in 1:nrow(bigcors.dat)){
 #print(i)
 name.cov <- bigcors.dat[i,1:2]
 if (name.cov[1]%in%all.covs && name.cov[2]%in%all.covs){
	r1 <- which(cov.rank$covariate==name.cov[1])
	r2 <- which(cov.rank$covariate==name.cov[2])
	if (r1>r2) all.covs <- all.covs[-which(all.covs==name.cov[2])]
	if (r2>r1) all.covs <- all.covs[-which(all.covs==name.cov[1])]
	}
 }

#Check that all remaining covariates are not highly correlated
corrs  <- cor(covars[,all.covs])
corrs

######################################################
###Covariates to be deleted from "covars" at stage 1 (|rho|>0.8))
#del <- c(9,16,19,23)
#These are:
#5 - Chirps, 14 - Distance_highway_MEAN, 10 - GPW, 19 - WC_preci, 3- Evapotrans
#23 - Elevation, 12 - Distance_urban, 25 - Nighttime_lights, 26 - Protected_areas
#covars.1 <- covars[,-del]
#sel.cov.names <- names(covars.1)
######################################################

covars.1      <- covars[,all.covs]
sel.cov.names <- names(covars.1)

#Selected covariates from singles 
singles1 <- singles[singles$model %in% sel.cov.names,]

#VIF analysis using the remaining covariates
covnames <- as.character(singles1$model)

Data <- cbind(TotChild, TotVax, covars.1)
head(Data)

#Delete clusters where TotChild is zero
zero.clust <- which(Data$TotChild==0)
if (length(zero.clust)>0) Data <- Data[-zero.clust,]

form    <- paste("cbind(TotVax, TotChild-TotVax)","~", paste(covnames, collapse=" + "))
mod.vif <- glm(form, data=Data, family = binomial(logit))
summary(mod.vif)
vif(mod.vif)
#####################################################################################



#Covariate selection in a binomial regression model
keepsingles <- subset(singles1, pr2>0) #Use all slected covariates or a subset based on 
                                       #pr2 (predictive R-square) values

covnames <- as.character(keepsingles$model)

#P-value approach - akin to backward elimination
#Fit all models
covnames <- as.character(keepsingles$model)
form <- paste("cbind(TotVax, TotChild-TotVax)","~", paste(covnames, collapse=" + "))
fit  <- glm(form, data = Data, family = binomial(logit))
summary(fit)
p.values <- summary(fit)$coefficients[,4]     #p-values
p.values <- p.values[-1]
max.p.no <- as.numeric(which(p.values==max(p.values)) )

max.p <- as.numeric(p.values[as.numeric(max.p.no)])
while (max.p > 0.05){
#	print(covnames)
	covnames <- covnames[-max.p.no]
	form <- paste("cbind(TotVax, TotChild-TotVax)","~", paste(covnames, collapse=" + "))
	fit  <- glm(form, data = Data, family = binomial(logit))
	p.values <- summary(fit)$coefficients[,4]     #p-values
	p.values <- p.values[-1]
	max.p.no <- as.numeric(which(p.values==max(p.values))) 
	max.p <- as.numeric(p.values[as.numeric(max.p.no)])
	print(max.p)
}


#Additional checks with selected covs
form <- paste("cbind(TotVax, TotChild-TotVax)","~", paste(covnames, collapse=" + "))
fit  <- glm(form, data = Data, family = binomial(logit))
#summary(fit)

#R^2
p.est <- as.numeric(fit$fitted)
p.obs <- Data[["TotVax"]]/Data[["TotChild"]]
p.obs.bar <- mean(p.obs)
r2 <- (cor(p.est, p.obs))^2

#P-values
p.values <- summary(fit)$coefficients[,4]; p.values <- p.values[-1]

#Clear r console using ctrl + L

#Print r-squared, Print names of selected covariates, Print VIF values and p-values
cat("r-squared for age group 0-59 is:", r2, sep="  ", "\n")
cat("Selected covariates for age group 0-59 are:", covnames, sep="  ", "\n")
cat("VIF values of selected covariates for 0-59 months are:", vif(fit), sep="  ", "\n")
cat("P-values of selected covariates for 0-59 months are:", p.values, sep="  ", "\n")


##Fill in manually
#Names of selected covariates for all age groups
covnames <- c("l_Distance_GUF_MEAN", "NPP_MEAN", "Aridity_MEAN", 
"l_Goat_density_MEAN", "l_Distance_rail_MEAN", "l_Nighttime_lights_MEAN")

Data <- cbind(TotChild, TotVax, covars.1)

#Export selected covariates for all age groups
cov.exp <- subset(Data, select=covnames)
head(cov.exp)
write.csv(cov.exp, "covars_selected_log_measles.csv")
#--------------------------------------End of covariate selection---------------------------------#





