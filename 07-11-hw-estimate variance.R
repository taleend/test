source('clarkFunctions.R')
#----Produce data----------------------
n <- 50
mu <- 4
sd <- 1
y <- rnorm(n,mu,sd)    #produce data

#----Simulation------------------------
V <- c(0.01, 0.1, 0.5, 1, 2, 5)
ng <- 2000                           #no. of iterations
vgibbs <- matrix(NA,length(V),ng)    #vector to hold samples
accept <- rep(0,length(V))
vg <- 0                              # initial value for v

for (j in 1:length(V))
{
for (g in 1:ng)
{
 vprop <- tnorm(1,0,4,vg,V[j])            #proposal
 pnow <- sum(dnorm(y,mu,vg,log=T))+sum(dinvGamma(y,1,1,log=T))   #log Pr current value
 pnew <- sum(dnorm(y,mu,vprop,log=T))+sum(dinvGamma(y,1,1,log=T)) #log Pr proposal
 atmp <- acceptMH(pnow,pnew,vg,vprop,BLOCK=T)
 vg <- atmp$x
 accept[j] <- accept[j] + atmp$accept
 vgibbs[j,g] <-vg
}
}

for (j in 1:length(V))
{
x11()
processPars(vgibbs[j,],xtrue=sd,CPLOT=T,DPLOT=T)
}