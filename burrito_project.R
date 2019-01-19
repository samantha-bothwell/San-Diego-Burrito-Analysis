##burrito project

rm(list=ls())
setwd("D:/School/Spring 2018/Stat 440/Labs/Lab 3")
load("burritodata.Rda")

#make data frame for top ingredients
count=rep(NA,33)
for(i in 1:33){
  count[i]=sum(burrito[,i+17])
}
list1=rbind(count,ingredients)
list1=t(list1)
list1=as.data.frame(list1)
list1$count=as.numeric(as.character(list1$count))
list1=list1[order(-list1$count),]
list1[1:15,]

#remove values 113,135
which(is.na(burrito$Cost))
burrito=burrito[complete.cases(burrito$Cost),]

## main ingredients: Beef, Pork, Chicken, Fish, Shrimp
burrito.meat=burrito[c("Location","Cost","Beef","Pork","Chicken","Fish","Shrimp")]

#declare function for model
BLMM <- function(y,x,z,S=1000,mu=0,tau=10,a1=0.01, a2=0.01, b1=0.01, b2=0.01){
  
  # create the variables now to avoid repeating this calculation each step
  n <- nrow(x)
  p <- ncol(x)
  if(is.null(p)) p <- 1
  g <- length(unique(z))
  X <- cbind(1,x,model.matrix(~z-1))  # add intercept and random effect indicators
  Xy <- t(X)%*%y
  XX <- t(X)%*%X
  
  # Places to store results
  beta_results <- matrix(NA, S, ncol(X))
  sig2inv_results <- rep(NA,S)
  kappa2inv_results <- rep(NA,S)
  
  # starting values
  sig2inv <- 1/var(y)
  kappa2inv <- 1/var(aggregate(y,by=list(z),mean)$x)
  
  # MCMC!!!!
  for(s in 1:S){
    
    # update beta including intercept
    tau2inv <- c(0,rep(1/tau^2,p), rep(kappa2inv,g))
    v <- chol2inv(chol((sig2inv*XX+diag(tau2inv))))
    beta <- t(chol(v)) %*% rnorm(ncol(X)) + v %*% ( sig2inv*Xy +tau2inv*mu  )
    
    #update sig2inv
    sig2inv <- rgamma(1, a1 + length(y)/2, a2 + sum((y-X%*%beta)^2)/2 )
    
    #update sig2inv
    kappa2inv <- rgamma(1, b1 + g/2, b2 + sum((beta[-c(1:(1+p))])^2)/2 )
    beta[2:3]
    
    # store results
    beta_results[s,] <- beta
    sig2inv_results[s] <- sig2inv
    kappa2inv_results[s] <- kappa2inv
    
  }
  
  return(list(alpha=beta_results[,1], beta=beta_results[,2:(1+p)], gamma=beta_results[,-c(1:(1+p))], 
              sigma=sqrt(1/sig2inv_results), kappa=sqrt(1/kappa2inv_results)))
}

#make x,y,z for meat data
y.meat=burrito.meat$Cost
x.meat=cbind(burrito.meat$Beef,burrito.meat$Pork,burrito.meat$Chicken,
             burrito.meat$Fish,burrito.meat$Shrimp)
z.meat=burrito.meat$Location

meat=BLMM(y=y.meat,x=x.meat,z=z.meat)

#which restaurants have lowest cost
table(apply(meat$gamma,1,which.min))/1000

#results for meat data
meat.results=data.frame(
  mean=c(mean(meat$alpha),mean(meat$beta[,1]),mean(meat$beta[,2]),mean(meat$beta[,3]),
         mean(meat$beta[,4]),mean(meat$beta[,5]),mean(meat$sigma),mean(meat$kappa)),
  sd=c(sd(meat$alpha),sd(meat$beta[,1]),sd(meat$beta[,2]),sd(meat$beta[,3]),sd(meat$beta[,4]),
       sd(meat$beta[,5]),sd(meat$sigma),sd(meat$kappa)),
  quantile_lower = c(quantile(meat$alpha, 0.025), quantile(meat$beta[,1], 0.025), 
                     quantile(meat$beta[,2], 0.025), quantile(meat$beta[,3], 0.025), 
                     quantile(meat$beta[,4], 0.025), quantile(meat$beta[,5], 0.025), 
                     quantile(meat$sigma, 0.025),quantile(meat$kappa,0.025)),
  quantile_upper = c(quantile(meat$alpha, 0.975), quantile(meat$beta[,1], 0.975), 
                     quantile(meat$beta[,2], 0.975), quantile(meat$beta[,3], 0.975), 
                     quantile(meat$beta[,4], 0.975), quantile(meat$beta[,5], 0.975), 
                     quantile(meat$sigma, 0.975),quantile(meat$kappa,0.025)) 
)
row.names(meat.results)=c("alpha","Beef","Pork","Chicken","Fish","Shrimp","sigma","kappa")
round(meat.results,2)

#just looking at significance of beta_i's
meat.plot=meat.results[c(2:6),]

library(plotrix)

#CI graph for beta_i
par(mfrow=c(1,1))
plotCI(x=meat.plot$mean,y=NULL,uiw=meat.plot$quantile_upper-meat.plot$mean,
       liw = meat.plot$mean-meat.plot$quantile_lower,xaxt="n",
       xlab="",ylab="",main="CI for Meat")
axis(1,at=seq(1,5,by=1),labels=c("Beef","Pork","Chicken","Fish","Shrimp"))
abline(h=0,col="red")

#traceplots
par(mfrow=c(2,4))
plot(meat$alpha,type = "l",las=1,main = "alpha")
plot(meat$sigma,type = "l",las=1,main = "sigma")
plot(meat$kappa,type = "l",las=1,main = "kappa")
plot(meat$beta[,1],type = "l",las=1,main = "beta1 (Beef)")
plot(meat$beta[,2],type = "l",las=1,main = "beta2 (Pork)")
plot(meat$beta[,3],type = "l",las=1,main = "beta3 (Chicken)")
plot(meat$beta[,4],type = "l",las=1,main = "beta4 (Fish)")
plot(meat$beta[,4],type = "l",las=1,main = "beta5 (Shrimp)")


#looking at top 10 other ingredients

burrito.other=burrito[c("Location","Cost","Cheese","Pico","Guac","Fries","Sour_cream",
                        "Sauce","Rice","Beans","Onion","Cilantro")]

y.other=burrito.other$Cost
x.other=cbind(burrito.other$Cheese,burrito.other$Pico,burrito.other$Guac,burrito.other$Fries,
              burrito.other$Sour_cream,burrito.other$Sauce,burrito.other$Rice,burrito.other$Beans,
              burrito.other$Onion,burrito.other$Cilantro)
z.other=burrito.other$Location

other=BLMM(y=y.other,x=x.other,z=z.other)

other.results=data.frame(
  mean=c(mean(other$alpha),mean(other$beta[,1]),mean(other$beta[,2]),mean(other$beta[,3]),
         mean(other$beta[,4]),mean(other$beta[,5]),mean(other$beta[,6]),mean(other$beta[,7]),
         mean(other$beta[,8]),mean(other$beta[,9]),mean(other$beta[,10]),mean(other$sigma),
         mean(other$kappa)),
  sd=c(sd(other$alpha),sd(other$beta[,1]),sd(other$beta[,2]),sd(other$beta[,3]),sd(other$beta[,4]),
       sd(other$beta[,5]),sd(other$beta[,6]),sd(other$beta[,7]),sd(other$beta[,8]),sd(other$beta[,9]),
       sd(other$beta[,10]),sd(other$sigma),sd(other$kappa)),
  quantile_lower = c(quantile(other$alpha, 0.025), quantile(other$beta[,1], 0.025), 
                     quantile(other$beta[,2], 0.025), quantile(other$beta[,3], 0.025), 
                     quantile(other$beta[,4], 0.025), quantile(other$beta[,5], 0.025),
                     quantile(other$beta[,6], 0.025),quantile(other$beta[,7], 0.025),
                     quantile(other$beta[,8], 0.025),quantile(other$beta[,9], 0.025),
                     quantile(other$beta[,10], 0.025),quantile(other$sigma, 0.025),
                     quantile(other$kappa, 0.025)),
  quantile_upper = c(quantile(other$alpha, 0.975), quantile(other$beta[,1], 0.975), 
                     quantile(other$beta[,2], 0.975), quantile(other$beta[,3], 0.975), 
                     quantile(other$beta[,4], 0.975), quantile(other$beta[,5], 0.975),
                     quantile(other$beta[,6], 0.975), quantile(other$beta[,7], 0.975),
                     quantile(other$beta[,8], 0.975), quantile(other$beta[,9], 0.975), 
                     quantile(other$beta[,10], 0.975), quantile(other$sigma, 0.975),
                     quantile(other$kappa,0.025))
)
row.names(other.results)=c("alpha","Cheese","Pico","Guac","Fries","Sour_cream",
                           "Sauce","Rice","Beans","Onion","Cilantro","sigma","kappa")
round(other.results,2)
other.plot=other.results[c(2:11),]

#CI for other ingredients
plotCI(x=other.plot$mean,y=NULL,uiw=other.plot$quantile_upper-other.plot$mean,
       liw = other.plot$mean-other.plot$quantile_lower,xaxt="n",ylab="",xlab="",
       main="CI for other ingredients")
axis(1,at=seq(1,10,by=1),labels=c("Cheese","Pico","Guac","Fries","Sour_cream",
                                  "Sauce","Rice","Beans","Onion","Cilantro"))
abline(h=0,col="red")


