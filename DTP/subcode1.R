#For generating parameter estimates and calculating the R-squared values of the 
#fitted models

ncovs 	<- 5    #No of covariates
burn.in 	<- 0.5*n.samples + 1
sub.samps 	<- burn.in:n.samples
print(summary(window(mod$p.beta.theta.samples, start=burn.in)))

beta.hat <- mod$p.beta.theta.samples[sub.samps,]

#Parameter estimates
param.out <- t(apply(beta.hat, 2, function (x) c(mean(x), sd(x), quantile(x, probs = c(0.025, 0.25, 0.5,0.75,0.975)))))
colnames(param.out) <- c("Mean", "SD", "2.5%", "25%", "50%", "75%", "97.5%")

#Calculate p-hat
beta.hat.1 <- t(mod$p.beta.theta.samples[sub.samps,1:((ncovs+1)*3)])
w.hat <- mod$p.w.samples[,sub.samps]
p.hat  <- 1/(1+exp(-(x%*%beta.hat.1+w.hat))) 
p.hat.a <- apply(p.hat, 1, function(x) c(mean(x), sd(x)))

p.hat.1 <-  p.hat.a[1,seq(1,ncol(p.hat.a),q)]
p.hat.2 <-  p.hat.a[1,seq(2,ncol(p.hat.a),q)]
p.hat.3 <-  p.hat.a[1,seq(3,ncol(p.hat.a),q)]

#Calculate R-squared 
p.est1 <- p.hat.1
p.obs1 <- Numvacc1/weights1
p.obs.bar1 <- mean(p.obs1)
r2.1 <- (cor(p.est1,p.obs1))^2
#bias.1 <- (sum(p.est1-p.obs1)/sum(p.obs1))*100

p.est2 <- p.hat.2
p.obs2 <- Numvacc2/weights2
p.obs.bar2 <- mean(p.obs2)
r2.2 <- (cor(p.est2,p.obs2))^2
#bias.2 <- (sum(p.est2-p.obs2)/sum(p.obs2))*100

p.est3 <- p.hat.3
p.obs3 <- Numvacc3/weights3
p.obs.bar3 <- mean(p.obs3)
r2.3 <- (cor(p.est3,p.obs3))^2
#bias.3 <- (sum(p.est3-p.obs3)/sum(p.obs3))*100


