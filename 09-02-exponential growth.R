source('clarkFunctions.R')

nt <- 100 #no. time intervals
time <- c(1:nt) #sample times
nm <- 30 #number of missing values
wm <- sort(sample(c(2:nt),nm)) #which are missing
sig <- 10 #process err variance
tau <- 20 #obs error variance
x0 <- 10 #1st value
x <- rep(x0,nt)
beta <- 0.1
for(t in 1:(nt-1))x[t+1] <- x[t] + beta + rnorm(1,0,sqrt(sig))

y <- rnorm(nt,x,sqrt(tau)) #simulate observations

y[wm] <- NA
plot(time,y)
lines(time,x)

#-------simulation------------
#xg <- initialStatesSS(y)$x
library(stats)
xg <- y    
xg[wm] <- predict(smooth.spline(time[-wm],y[-wm]),time)$y[wm]

mus <- 100 #mean for process error variance
s1 <- 10 #weak prior
s2 <- mus*(s1 - 1)

mut <- 10 #mean for obs error variance
v1 <- 2*nt #wt comparable to data
v2 <- mut*(v1 - 1)

priorB <- 0 #non-informative prior
priorIVB <- 10 #prior precision
sg <- mus #initialize
tg <- mut
bg <- 0
ng <- 1000
bgibbs <- rep(0,ng)
vgibbs <- matrix(0,ng,2)
xgibbs <- matrix(0,ng,nt)

for(g in 1:ng){
bg <- updateSSB1(xg,sg,priorB,priorIVB)
xg <- updateSSRW(xg,y,wm,tg,sg)
sg <- updateSigmaIG(xg[-1],xg[-nt] + bg,s1,s2)
tg <- updateSigmaIG(y[-wm],xg[-wm],v1,v2)
vgibbs[g,] <- c(sg,tg)
xgibbs[g,] <- xg
bgibbs[g] <- bg
}

processPars(cbind(bgibbs,vgibbs),xtrue=c(beta,sig,tau),CPLOT=T,DPLOT=T)
x11()
processStates(xgibbs,y)
lines(x,col=2)