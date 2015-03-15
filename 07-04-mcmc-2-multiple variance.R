n <- 50
mu <- 1
sd <- 1
y <- rnorm(n,mu,sd)    #produce data


ng <- 10000                      #no. of iterations
mug <- 1                         #initial value for the mean
v <- c(0.01,0.05,.2,0.8,2)       #variance for the proposal distribution
nc <- length(v)
mgibbs <- matrix(0,nc,ng)        #vector to hold samples
accept <- rep(0,nc)              #counter for number of acceptances

mugib <- rep(0,nc)               #average for each row
sdgib <- rep(0,nc)               #sd for each row

for (h in 1:nc)
{
for(g in 1:ng)
{
 mprop <- rnorm(1,mug,v[h])                   #proposal
 pnow <- sum(dnorm(y,mug,sd,log=T))        #log Pr current value
 pnew <- sum(dnorm(y,mprop,sd,log=T))      #log Pr proposal
 atmp <- acceptMH(pnow,pnew,mug,mprop,BLOCK=T)
 mug <- atmp$x
 accept[h] <- accept[h] + atmp$accept
 mgibbs[h,g] <- mug                          #store value
}
x11()
processPars(mgibbs[h,],xtrue=mu,CPLOT=T,DPLOT=T)
}

