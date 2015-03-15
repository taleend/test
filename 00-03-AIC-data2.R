source('clarkFunctions.R')

getZ <- function(km){log(-w) - log(-km - w)}

T1 <- read.table("arc1.txt",header=F)
T1 <- na.omit(T1)
t1 <- T1[which(T1$V5<0),] #--exclude positive values for column 5(water table data)
yob <- t1$V3       #--observed data of soil respiration
y <- log(yob)      #--logged soil respiration data for estimation
x <- cbind(rep(1,length(y)),t1$V4,t1$V5)
t <- x[,2]
w <- x[,3]
#tw <- (x[,2] - mean(x[,2]))*(x[,3] - mean(x[,3]))
#x <- cbind(x[,1],t,tw-24*w)


parset <- read.table('00-data-arc1-2model.txt',header=T)
parset <- cbind(parset[,2],parset[,4:6])

np <- length(parset[,1])
parset <- parset[1:np-1,]
np <- np-1
R2 <- rep(0,np)            #---R square

L <- rep(0,np)

layout(matrix(1:18,3,6))

for (i in 1:np)
{
k <- parset[i,1]
be0 <- parset[i,2]
be1 <- parset[i,3]
be2 <- parset[i,4]
R0 <- w/(w+k)
Q10 <- exp(be0+be1*t+be2*w)
yp <- R0*Q10

yres <- yob-yp
yrsd <- sd(yres)
L[i] <- sum(dnorm(yres,0,yrsd,log=T)) #loged likelihood

lmy <- lm(yob~0+yp)     #linear models
slmy <- summary(lmy)
R2[i] <- slmy$r.squared

plot(yob,type='l',col=2) # plot the data in time order
title(paste('k =',k,sep=' '))
lines(yp,col=4)    
}

ML <- L[np]
D <- -2*(L[1:np]-ML)
cbind(parset[,1],D)
x11()
plot(parset[1:np-1,1],D[1:np-1],xlab='k values',ylab='D')