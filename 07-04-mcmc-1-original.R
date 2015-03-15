source('clarkFunctions.R')
n <- 50
mu <- 4
sd <- 1
y <- rnorm(n,mu,sd)    #produce data


ng <- 1000             #no. of iterations
mgibbs <- rep(0,ng)    #vector to hold samples
mug <- 1               #initial value for the mean
v <- .2                #variance for the proposal distribution
accept <- 0            #counter for number of acceptances
for(g in 1:ng)
{
mprop <- rnorm(1,mug,v)                   #proposal
pnow <- sum(dnorm(y,mug,sd,log=T))        #log Pr current value
pnew <- sum(dnorm(y,mprop,sd,log=T))      #log Pr proposal
atmp <- acceptMH(pnow,pnew,mug,mprop,BLOCK=T)
mug <- atmp$x
accept <- accept + atmp$accept
mgibbs[g] <- mug                          #store value
}