source('clarkFunctions.R')
library('mvtnorm')
#------function for the equation------------------
fx <- function(x,pars){pars[1]*x*(1-x/pars[2])}


#-------------data simulation------------------------
nt <- 100                         #no. time intervals
nm <- 20                          #number of missing values
wm <- sort(sample(c(2:nt),nm))    #which are missing
r <- .1                           #parameter values
k <- 100
sig <- 10                         #process err variance
x0 <- 1                           #1st value
x <- rep(x0,nt)
dt <- runif(nt,.4,2)              #random time increments
time <- cumsum(dt); dt <- dt[-nt] #observation times
tau <- 5^2      #obs error variance


for(t in 1:(nt-1))
{
x[t+1] <- tnorm(1,.1,1000,x[t] + fx(x[t],c(r,k))*dt[t],sqrt(sig*dt[t]))
}
y <- tnorm(nt,0,Inf,x,sqrt(tau))     #simulate observations
y[wm] <- NA
plot(time,y)
lines(time,x)

#---------modeling---------------------------------
mus <- 10 #mean for process error variance
s1 <- 10 #weak prior
s2 <- mus*(s1 - 1)
mut <- 25 #mean for obs error variance
v1 <- 2*nt #wt comparable to data
v2 <- mut*(v1 - 1)
priorB <- matrix(c(.2,100),ncol=1)
priorVB <- diag(c(1,10))
loB <- c(0,0)
hiB <- c(1,1000)


#tmp <- initialStatesSS(y)
#xg <- tmp$x
library(stats)
xg <- y    
xg[wm] <- predict(smooth.spline(time[-wm],y[-wm]),time)$y[wm]
notMiss <- c(1:length(y))
notMiss <- notMiss[-wm]
propSd <- sd(diff(y),na.rm=T)

propx <- propSd

sg <- mus #initialize
tg <- mut
rg <- .1
kg <- 150
bg <- matrix(c(rg,kg),ncol=1)
propVar <- diag(c(.01,.1))
ap <- 0 #acceptance rates for (r,K)
ax <- 0 #acceptance rates for x
ng <- 5000
pgibbs <- matrix(0,ng,2); colnames(pgibbs) <- c('r','K')
vgibbs <- matrix(0,ng,2); colnames(vgibbs) <- c('sigma','tau')
xgibbs <- matrix(0,ng,nt)
cvar <- matrix(c(.0001,.03,.03,100),2,2)
for(g in 1:ng){
#parameters
tmp <- updateSSparsMet(bg,priorB,priorVB,
loB,hiB,xg,sg,dt,propVar)
bg <- tmp$b
ap <- ap + tmp$aa

tmp <- updateSSstatesMet(xg,y,bg,sg,tg,lo=0,hi=1000,
notMiss=notMiss,
dt=dt,propSd=propx)
xg <- tmp$x
ax <- ax + tmp$aa
if(g %in% c(50,100,200,500,1000)){ #tune proposals
if(ax/g < .2)propx <- propx/2
if(ax/g > .4)propx <- propx*2
print(propx)
}
#sample variances
sg <- updateSigmaIG(xg[-1],xg[-nt]+fx(xg[-nt],bg),s1,s2)
tg <- updateSigmaIG(xg[notMiss],y[notMiss],v1,v2)
pgibbs[g,] <- bg
vgibbs[g,] <- c(sg,tg)
xgibbs[g,] <- xg
}

processPars(cbind(pgibbs ,vgibbs),xtrue=c(r,k,sig,tau),CPLOT=T)
x11()
processStates(xgibbs,y)
lines(time,x,col=2)
plot(time[wm],y[wm],col=4)