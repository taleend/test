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
tw <- (x[,2] - mean(x[,2]))*(x[,3] - mean(x[,3]))
x <- cbind(x[,1],t,tw-24*w)


parset <- read.table('00-data-arc1-1model-combine.txt',header=T)

np <- length(parset[,1])
R2 <- rep(0,np)            #---R square

L <- rep(0,np)

layout(matrix(1:24,4,6))

for (i in 1:np)
{
k <- parset[i,1]
v <- parset[i,2]
be1 <- parset[i,3]
be2 <- parset[i,4]
R0 <- w*v/(w+k)
Q10 <- be1*exp(be2*(w+10))
yp <- R0*Q10^(0.1*t-2.4)

yres <- yob-yp
yrsd <- sd(yres)
L[i] <- sum(dnorm(yres,0,yrsd,log=T)) #loged likelihood

lmy <- lm(yob~0+yp)     #linear models
slmy <- summary(lmy)
R2[i] <- slmy$r.squared

#plot(yob-yp)   #plot residuals
plot(yob,type='l',col=2) # plot the data in time order
lines(yp,col=4)          
#plot(yob,yp,ylim=c(0,25)) # scatter plot
title(paste('k =',k,sep=' '))
#lines(c(0,30),c(0,30),col=2)
}

ML <- L[np]
D <- -2*(L[1:np]-ML)
cbind(parset[,1],D)
plot(parset[1:np-1,1],D[1:np-1],xlab='k values',ylab='D')