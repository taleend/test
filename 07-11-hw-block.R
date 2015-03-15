source('clarkFunctions.R')

#----produce data-------------
n <- 50
mu <- 4
sd <- 1
y <- rnorm(n,mu,sd)    #produce data
lo <- c(-5,0.001)
hi <- c(5,4)

#----simulation---------------
library(mvtnorm)
ng <- 5000
gibbs <- matrix(NA,ng,2)
g <-c(0,0.1)  #initial values for mean and variance
Covmatrix <- matrix(c(.015,0,0,.01),2,2) # prior covariance
accept <- 0

for (j in 1:ng)
{
 prop <- rmvnorm(1,g,Covmatrix)
 while (prop[2]<0)  prop <- rmvnorm(1,g,Covmatrix) #discard the negative value for variance
 pnow <- sum(dnorm(y,g[1],g[2],log=T))         #log Pr current value
 pnew <- sum(dnorm(y,prop[1],prop[2],log=T))      #log Pr proposal
 atmp <- acceptMH(pnow,pnew,g,prop,BLOCK=T)
 g <- atmp$x
 gibbs[j,1] <- g[1]
 gibbs[j,2] <- g[2]
 accept <- accept+atmp$accept
}

#----Plots-----------------------
processPars(gibbs[,1],xtrue=mu,CPLOT=T,DPLOT=T,burnin=500)
x11()
processPars(gibbs[,2],xtrue=sd,CPLOT=T,DPLOT=T,burnin=500)
x11()
plot(gibbs[,1],gibbs[,2])
