#--------MCMC------------------------------------------
#--------withdraw data from normal distribution--------
n=50
mu=1
sd=1
td <- rnorm(n,mu,sd)

#-------a function --------
# accept for M, M-H
# if BLOCK, then accept as a block, otherwise, accept individually
acceptMH <- function(p0,p1,x0,x1,BLOCK=F) 
{
 nz <- length(x0)  # number to acept
 if(BLOCK)nz <- 1
 a <- exp(p1-p0)  # acceptance PR

 z <- runif(nz,0,1)
 keep <- which(z<a, arr.ind=T)

if(BLOCK & length(keep) > 0)
x0 <- x1
if(!BLOCK)
x0[keep] <- X1[keep]
accept <- length(keep)

list(x=x0,accept=accept)
}

#-------a metropolis sampler---------------
#-------produce data-----------------------
n <- 50
mu <- 1
sd <- 1
y <- rnorm(n,mu,sd)

#-------sampler---------------------------
#-------make nc chains--------------------
ng <- 1000          #number of iterations
nc <- 5
mgibbs <- matrix(0,nc,ng) #vector to hold samples
mug <- 1            #initial value for the mean
v <- .2             #variance for the proposal distribution
accept <- rep(0,nc)         #counter for number of acceptances

for (h in 1:nc)
{
for (g in 1:ng)
{
 mprop <- rnorm(1,mug,v)             #proposal
 pnow <- sum(dnorm(y,mug,sd,log=T))  #log Pr current value
 pnew <- sum(dnorm(y,mprop,sd,log=T))#log Pr proposal

 atmp <- acceptMH(pnow,pnew,mug,mprop,BLOCK=T)
 mug <- atmp$x
 accept[h] <- accept[h] + atmp$accept
 mgibbs[h,g] <- mug                    #store value
}
}

#-------------check if the average values of each chain are similar------
mugib <- rep(0,nc)
sdgib <- rep(0,nc)
for (l in 1:nc)
{
mugib[l] <- mean(mgibbs[l,])
sdgib[l] <- sd(mgibbs[l,])
rho <- ar(mgibbs[l,],order.max=10,aic=F)$ar
rho <- rho[1:min(c(1:10)[rho<.01])]
x11()
processPars(mgibbs[l,],xtrue=mu,CPLOT=T,DPLOT=T)
}

#--produce graphs that all chains are ploted together to see similarity--
plot(mgibbs[1,],type='l')
for (k in 2:nc)
{
lines(mgibbs[k,],col=k)
}



#--------------autocorrelation--------------
rho <- ar(mgibbs,order.max=10,aic=F)$ar
rho <- rho[1:min(c(1:10)[rho<.01])]
processPars(mgibbs,xtrue=mu,CPLOT=T,DPLOT=T)
