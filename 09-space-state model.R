source('clarkFunctions.R')
library(stats)
nt <- 100 #no. time intervals
time <- c(1:nt) #sample times
nm <- 30 #number of missing values
wm <- sort(sample(c(2:nt),nm)) #which are missing
sig <- 10 #process err variance
tau <- 20 #obs error variance
x0 <- 10 #1st value
x <- rep(x0,nt)
for(t in 1:(nt-1))x[t+1] <- x[t] + rnorm(1,0,sqrt(sig))
y <- rnorm(nt,x,sqrt(tau)) #simulate observations
y[wm] <- NA
plot(time,y)
lines(time,x)

xg <- y
xg[wm] <- predict(smooth.spline(time[-wm],y[-wm]),time)$y[wm]

mus <- 2 #mean for process error variance
s1 <- 2 #weak prior
s2 <- mus*(s1 - 1)

mut <- 20 #mean for obs error variance
v1 <- 2*nt #wt comparable to data
v2 <- mut*(v1 - 1)
sg <- mus #initialize variances
tg <- mut
ng <- 1000
vgibbs <- matrix(0,ng,2)
xgibbs <- matrix(0,ng,nt)
for(g in 1:ng){
xg <- updateStateSpaceRW(xg,y,wm,tg,sg) #sample x's
sg <- updateSigmaIG(xg[-1],xg[-nt],s1,s2)
tg <- updateSigmaIG(y[-wm],xg[-wm],v1,v2)
vgibbs[g,] <- c(sg,tg)
xgibbs[g,] <- xg
}

processPars(vgibbs,xtrue=c(sig,tau),CPLOT=T,DPLOT=T)

processStates(xgibbs,y)