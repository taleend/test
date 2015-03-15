source('clarkFunctions.R')

#----Produce data----------------------
n <- 12
m <- rnbinom(n,size=1,mu=20) + 1
theta <- rbeta(n,3,10)
y <- rbinom(n,m,theta)

#----simulation-------------------------

ng <- 1000
v <- 0.029
tgibbs <- matrix(NA,length(v),ng)
accept <- rep(0,length(v))
tg <- 0.001

for (i in 1:length(v))
{
for (j in 1:ng)
{
 tprop <- tnorm(1,0,1,tg,v[i])
 pnow <- sum(dbinom(y,m,tg,log=T))+sum(dbeta(tg,1,1,log=T))
 pnew <- sum(dbinom(y,m,tprop,log=T))+sum(dbeta(tprop,1,1,log=T))
 atmp <- acceptMH(pnow,pnew,tg,tprop,BLOCK=T)
 tg <-atmp$x
 accept[i] <- accept[i] + atmp$accept
 tgibbs[i,j] <- tg
}
x11()
processPars(tgibbs[i,],xtrue=0.23,CPLOT=T,DPLOT=T)
}
