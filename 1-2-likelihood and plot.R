source('clarkFunctions.R')
#----simulate data and function getZ---------------------------
getZ <- function(km){log(-w) - log(-km - w)}

V <- 13.58
be1 <- 3.63
be2 <- 0.053

b0 <- signif((log(V) - 2.4*log(be1) - 24*be2),5)
b1 <- signif(0.1*log(be1)+be2,5)
b2 <- -2.4*be2
b3 <- 0.1*be2

n    <- 1000
beta <- matrix(c(b0,b1,b2,b3),4,1) #--values of parameters from MLE
km   <- -7.2                                 #--value of km from MLE
sd   <- .5
sg   <- sd^2

loX <- c(1,3,-40)    #int, temp, h2o, interaction
hiX <- c(1,25,-0.001)      

x <- simX(n,loX,hiX)     #simulate x and y
t <- x[,2]
w <- x[,3]
x <- cbind(x,(x[,2] - mean(x[,2]))*(x[,3] - mean(x[,3])))
z <- getZ(km)
y <- rnorm(n,z + x%*%beta,sd)
#-----------------------------------------------

#----prediction vs. observation plots-------------
parset <- read.table("1-simpar-dif k.txt",header=T)
n <- length(parset[,1])

parset <- parset[1:n-1,]

n <- length(parset[,1])

yp <- matrix(NA,n,1000)     #--for prediction of y
resi <- matrix(NA,n,1000)   #--for residual (y-yp)

layout(matrix(1:6,2,3))
for (i in 1:n)
{
inik <- parset[i,1]
k <- parset[i,3]
b0 <- parset[i,4]
b1 <- parset[i,5]
b2<- parset[i,6]
b3<- parset[i,7]
beta2 <- matrix(c(b0,b1,b2,b3),4,1)

z2 <- getZ(k)
yp[i,] <- z2+x%*%beta2
resi[i,] <- y-yp[i,]
#----prediction vs. observation--------------
plot(y,yp[i,],ylab='prediction',xlab='simulated data')
title(paste('initial k =',inik,sep=' '))
lines(c(-25,25),c(-25,25),col=2) #--for predition vs. observation
#----residual plots--------------------------
#plot(y,y-yp[i,],ylab='residual',xlab='observation')


#abline(h=0,col=2)  #--for residual plots
}

#-----likelihood and R2-----------
L <- rep(0,n)
sdre <- rep(0,n)
R2 <- rep(0,n)

for (i in 1:n)
{
#----get logged likelihood----------
sdre[i] <- sd(resi[i,])
#L[i] <- prod(dnorm(resi[i,],0,sdre[i]))  # some of them go too small
L[i] <- -sum(dnorm(resi[i,],0,sdre[i],log=T)) # the smaller, the better

#----get R2-------------------
lmy <- lm(y~0+yp[i,])     #linear models
slmy <- summary(lmy)
R2[i] <- slmy$r.squared

}
cbind(parset,L)