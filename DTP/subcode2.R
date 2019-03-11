#For calculating cross-validation statistics - 95% coverage rate, VMSE and percentage bias 


#Predict at validation locations
  pred.coords   <- coords.nc
  pred.covars.a <- data.matrix(cbind(rep(1,nrow(pred.coords)),vaxcov.nc)) #Note const added
  pred.covars   <- mkMvX(list(pred.covars.a,pred.covars.a,pred.covars.a))
  mod.pred      <- spPredict(mod, pred.coords, pred.covars, start=burn.in, n.samples, thin=1, 
                           verbose=FALSE, n.report=100) #change thin

#Compute estimates of numbers of vaccinated children for the calculation
#of the coverage rates

hh  <- mod.pred$p.y.predictive.samples
hh1 <- hh[seq(1,nrow(hh),q),]
hh2 <- hh[seq(2,nrow(hh),q),]
hh3 <- hh[seq(3,nrow(hh),q),]
  
rr1=rr2=rr3=matrix(0,nrow(hh1),ncol(hh1))
for (i in 1:nrow(hh1)){
	for (j in 1:ncol(hh1)){
	rr1[i,j] <- rbinom(1,weights1.nc[i],hh1[i,j])
	rr2[i,j] <- rbinom(1,weights2.nc[i],hh2[i,j])
	rr3[i,j] <- rbinom(1,weights3.nc[i],hh3[i,j])
	}
}

  y.hat.mean1  <- apply(rr1, 1, mean);  y.hat.sd1    <- apply(rr1, 1, sd)
  y.hat.low1   <- apply(rr1, 1, function(x) quantile(x,0.025)); y.hat.up1 <- apply(rr1, 1, function(x) quantile(x,0.975))   #Try later
  y.obs1 <- Numvacc1.nc
 
  y.hat.mean2  <- apply(rr2, 1, mean);  y.hat.sd2    <- apply(rr2, 1, sd)
  y.hat.low2   <- apply(rr2, 1, function(x) quantile(x,0.025)); y.hat.up2 <- apply(rr2, 1, function(x) quantile(x,0.975))   #Try later
  y.obs2 <- Numvacc2.nc
  
  y.hat.mean3  <- apply(rr3, 1, mean);  y.hat.sd3    <- apply(rr3, 1, sd)
  y.hat.low3   <- apply(rr3, 1, function(x) quantile(x,0.025)); y.hat.up3 <- apply(rr3, 1, function(x) quantile(x,0.975))   #Try later
  y.obs3 <- Numvacc3.nc
  
  
count1 = count2 = count3 = 0
for(i in 1:nc){
if ((y.obs1[i] >= y.hat.low1[i]) && (y.obs1[i] <= y.hat.up1[i])) count1 <- count1 + 1
if ((y.obs2[i] >= y.hat.low2[i]) && (y.obs2[i] <= y.hat.up2[i])) count2 <- count2 + 1
if ((y.obs3[i] >= y.hat.low3[i]) && (y.obs3[i] <= y.hat.up3[i])) count3 <- count3 + 1
}


#Coverage rate
cov.rate1 <- (count1/nc)*100
cov.rate2 <- (count2/nc)*100
cov.rate3 <- (count3/nc)*100

#Validation mean square error
p.hat.mean1  <- apply(hh1, 1, mean); p.obs1 <- Numvacc1.nc/weights1.nc
vmse1 <- sum((p.hat.mean1-p.obs1)^2)/nc

p.hat.mean2  <- apply(hh2, 1, mean); p.obs2 <- Numvacc1.nc/weights1.nc
vmse2 <- sum((p.hat.mean2-p.obs2)^2)/nc

p.hat.mean3  <- apply(hh3, 1, mean); p.obs3 <- Numvacc1.nc/weights1.nc
vmse3 <- sum((p.hat.mean3-p.obs3)^2)/nc

#Percentage bias
bias.1 <- (sum(p.hat.mean1-p.obs1)/sum(p.obs1))*100
bias.2 <- (sum(p.hat.mean2-p.obs2)/sum(p.obs2))*100
bias.3 <- (sum(p.hat.mean3-p.obs3)/sum(p.obs3))*100









