#For calculating cross-validation statistics - 95% coverage rate, VMSE and percentage bias 

#Predict at validation locations
pred.coords <- coords.nc
pred.covars <- data.matrix(cbind(rep(1,nrow(pred.coords)),vaxcov.nc)) #Note const added
mod.pred    <- spPredict(mod, pred.coords, pred.covars, start=burn.in, n.samples, thin=1, 
                         verbose=FALSE, n.report=100) #change thin

#Compute estimates of numbers of vaccinated children for the calculation
#of the coverage rates
hh <- mod.pred$p.y.predictive.samples
rr <- matrix(0,nrow(hh),ncol(hh))
for (i in 1:nrow(hh)){
  for (j in 1:ncol(hh)){
    rr[i,j] <- rbinom(1,weights.nc[i],hh[i,j])
  }
}

y.hat.mean  <- apply(rr, 1, mean)
y.hat.sd    <- apply(rr, 1, sd)
y.hat.low   <- apply(rr, 1, function(x) quantile(x,0.025))
y.hat.up    <- apply(rr, 1, function(x) quantile(x,0.975))

y.obs <- Numvacc.nc
count <- 0
for(i in 1:nc){
  if ((y.obs[i] >= y.hat.low[i]) && (y.obs[i] <= y.hat.up[i])) count <- count + 1
}

#Coverage rate
cov.rate <- (count/nc)*100

#Validation mean square error
p.hat.mean  <- apply(mod.pred$p.y.predictive.samples, 1, mean)
p.obs       <- Numvacc.nc/weights.nc
vmse <- sum((p.hat.mean-p.obs)^2)/nc

#Percentage bias
bias <- (sum(p.hat.mean-p.obs)/sum(p.obs))*100





