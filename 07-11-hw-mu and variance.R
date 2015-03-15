source('clarkFunctions.R')

#----produce data-------------
n <- 50
mu <- 4
sd <- 1
y <- rnorm(n,mu,sd)    #produce data

#----simulation---------------
ng <- 2000
mgibbs <- rep(NA,ng)
vgibbs <- mgibbs
accept <- rep(0,2)
Vm <- 0.2
Vs <- 0.15

mg <- 0
vg <- 0.1

for (i in 1:ng)
{
 mprop <- rnorm(1,mg,Vm)                   #proposal
 pnow <- sum(dnorm(y,mg,vg,log=T))         #log Pr current value
 pnew <- sum(dnorm(y,mprop,vg,log=T))      #log Pr proposal
 atmp <- acceptMH(pnow,pnew,mg,mprop,BLOCK=T)
 mg <- atmp$x
 mgibbs[i] <- mg
 accept[1] <- accept[1]+atmp$accept   

 vprop <- rnorm(1,vg,Vs)              #proposal
 pnow <- sum(dnorm(y,mg,vg,log=T))    #log Pr current value
 pnew <- sum(dnorm(y,mg,vprop,log=T)) #log Pr proposal
 atmp <- acceptMH(pnow,pnew,vg,vprop,BLOCK=T)
 vg <- atmp$x
 vgibbs[i] <-vg
 accept[2] <- accept[2]+atmp$accept  
}

#----Plots--------------------------------
processPars(mgibbs,xtrue=mu,CPLOT=T,DPLOT=T,burnin=200)
x11() 
processPars(vgibbs,xtrue=sd,CPLOT=T,DPLOT=T,burnin=200)
x11() 
plot(mgibbs,vgibbs)
