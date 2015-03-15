#------produce data------------------
mu <- 1 #simulate data
sigma <- 1
y <- rnorm(100,mu,sigma)
#----
ng<-5000
sgibbs <- rep(0,ng)
sg <- .1
jump <- .1
accept <- 0
s1 <- .1
s2 <- .1

for(g in 1:ng)
{
sprop <- rlnorm(1,log(sg),jump) #propose from the lognormal

#Pr for current and proposed values, corrected for the asymmetry

pnow <- sum(dnorm(y,mu,sg,log=T)) + #likelihood

dgamma(sg,s1,s2,log=T) - #prior

dlnorm(sg,log(sprop),jump,log=T) #jump correction

pnew <- sum(dnorm(y,mu,sprop,log=T)) +
dgamma(sprop,s1,s2,log=T) -
dlnorm(sprop,log(sg),jump,log=T)

atmp <- acceptMH(pnow,pnew,sg,sprop,BLOCK=T)
sg <- atmp$x
sgibbs[g] <- sg
accept <- accept + atmp$accept
}