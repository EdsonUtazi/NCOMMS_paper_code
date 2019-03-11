#For generating parameter estimates and calculating the R-squared values of the 
#fitted models 

ncovs 	<- 4    #No of covariates
burn.in 	<- 0.5*n.samples + 1
sub.samps 	<- burn.in:n.samples
print(summary(window(mod$p.beta.theta.samples, start=burn.in)))

#Parameter estimates
beta.hat <- mod$p.beta.theta.samples[sub.samps,]
param.out <- t(apply(beta.hat, 2, function (x) c(mean(x), sd(x), quantile(x, probs = c(0.025, 0.25, 0.5,0.75,0.975)))))
colnames(param.out) <- c("Mean", "SD", "2.5%", "25%", "50%", "75%", "97.5%")

#Calculate p-hat
w.hat <- mod$p.w.samples[,sub.samps]
X <- mod$X

linear.pred <- matrix(0, nrow=nrow(X), ncol=length(sub.samps))
for(i in 1:length(sub.samps)){
  linear.pred[,i] <- X%*%matrix(beta.hat[i,1:(ncovs+1)],nrow=(ncovs+1), ncol=1)
}
 
p.hat  <- 1/(1+exp(-(linear.pred+w.hat)))  
p.hat.1 <- apply(p.hat, 1, function(x) c(mean(x), sd(x)))

#Calculate R-squared 
p.est <- p.hat.1[1,]
p.obs <- Numvacc/weights
p.obs.bar <- mean(p.obs)
r2 <- (cor(p.est, p.obs))^2
#bias <- (sum(p.est-p.obs)/sum(p.obs))*100





