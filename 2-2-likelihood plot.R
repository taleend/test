source('clarkFunctions.R')

#----observation data------------------------------
t1 <- read.table("arc1ready.txt",header=T)
yob <- t1$V3       #--observed data of soil respiration
y <- log(yob)      #--logged soil respiration data for estimation
x <- cbind(rep(1,length(y)),t1$V4,t1$V5)
t <- x[,2]
w <- x[,3]
tw <- (x[,2] - mean(x[,2]))*(x[,3] - mean(x[,3]))
x <- cbind(x,tw)


#----read parameter estimates-------------
nx <- length(tw)
parset <- read.table("2-4bs.txt",header=T)

np <- length(parset[,1])
yres <- matrix(NA,np,nx)

L <- rep(0,np)
sdre <- rep(0,np)
R2 <- rep(0,np)

#---plot prediction vs. observation, get residual---------
layout(matrix(1:18,3,6))
for (i in 1:np)
{
k <- parset[i,2]
b0 <- parset[i,4]
b1 <- parset[i,5]
b2 <- parset[i,6]
b3 <- parset[i,7]
R0 <- w/(w+k)
Q10 <- exp(b0+b1*t+b2*w+b3*tw)
yp <- R0*Q10

plot(yob,type='l',col=2,xlab='time order',ylab='soil respiration') # plot the data in time order
title(paste('k =',k,sep=' '))
lines(yp,col=4)   

yres[i,] <- yob-yp
sdre <- sd(yres[i,])
L[i] <- -sum(dnorm(yres[i,],0,sdre,log=T)) # the smaller, the better
#----get R2-------------------
lmy <- lm(yob~0+yp)     #linear models
slmy <- summary(lmy)
R2[i] <- slmy$r.squared
}

parL <- cbind(parset,L,R2)



#----selected parset------

parL <- parL[7:15,]

np <- length(parL[,1])

layout(matrix(1:9,3,3))
for (i in 1:np)
{
k <- parL[i,2]
b0 <- parL[i,4]
b1 <- parL[i,5]
b2 <- parL[i,6]
b3 <- parL[i,7]
R0 <- w/(w+k)
Q1 <- exp(b0+b1*t+b2*w+b3*tw)
Q2 <- exp(b0+b1*t+b2*w)
yd <- R0*(Q1-Q2)
hist(yd,main=NULL)
title(paste('k=',k,sep=' '))
}