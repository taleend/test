source('clarkFunctions.R')

tmp <- inData('dataAlces.txt',ynames='density',
tname='year',na.rm=F)
y <- as.vector(tmp$y)
time <- tmp$t
nt <- length(time)
dt <- diff(time)

library(stats)
n <- length(y)
wm   <- which(is.na(y))
notMiss <- c(1:n)
notMiss <- notMiss[-wm]
xg <- y
xg[wm] <- predict(smooth.spline(time[-wm],y[-wm]),time)$y[wm]
propSd <- sd(diff(y),na.rm=T)
propx <- propSd/nt

priorB <- matrix(c(0.2,100),ncol=1)
priorVB <- diag(c(100,100))
loB <- c(0,0)
hiB <- c(1,1000)

#plot(time,y)
#points(time[wm],xg[wm],col=2)


mus <- 4^2 #mean for process error variance
s1 <- 1 #weak prior
s2 <- 1
mut <- 10 #mean for obs error variance
v1 <- 1 #wt comparable to data
v2 <- 1



sg <- mus #initialize
tg <- mut
rg <- .1
kg <- 150
bg <- matrix(c(rg,kg),ncol=1)
propVar <- diag(c(.01,1))
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
x11()
processPars(cbind(pgibbs ,vgibbs),CPLOT=T)
x11()
processStates(xgibbs,y)