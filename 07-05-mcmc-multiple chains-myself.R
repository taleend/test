#--------MCMC------------------------------------------
#--------withdraw data from normal distribution--------
n=50
mu=1
sd=1
td <- rnorm(n,mu,sd)


#-------a metropolis sampler---------------
#-------produce data-----------------------
n <- 50
mu <- 1
sd <- 1
y <- rnorm(n,mu,sd)

ng <- 1000          #number of iterations
mgibbs <- matrix(data=0,nrow=3,ncol=ng) #vector to hold samples
mug <- 1            #initial value for the mean
v <- .2             #variance for the proposal distribution
accept <- rep(0,3)

for (i in 1:3)
{
 for (g in 1:ng)
 {
  mprop <- rnorm(1,mug,v)             #proposal
  pnow <- sum(dnorm(y,mug,sd,log=T))  #log Pr current value
  pnew <- sum(dnorm(y,mprop,sd,log=T))#log Pr proposal

  atmp <- acceptMH(pnow,pnew,mug,mprop,BLOCK=T)
  mug <- atmp$x
  accept[i] <- accept[i] + atmp$accept
  mgibbs[i,g] <- mug                    #store value
 }
x11()
processPars(mgibbs[i,],xtrue=mu,CPLOT=T,DPLOT=T)
}

#--------------autocorrelation--------------
rho <- ar(mgibbs[1,],order.max=10,aic=F)$ar
rho <- rho[1:min(c(1:10)[rho<.01])]
#processPars(mgibbs[1,],xtrue=mu,CPLOT=T,DPLOT=T)
#x11()
#processPars(mgibbs[2,],xtrue=mu,CPLOT=T,DPLOT=T)
#x11()
#processPars(mgibbs[3,],xtrue=mu,CPLOT=T,DPLOT=T)
