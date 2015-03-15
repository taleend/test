source('clarkFunctions.r')
tmp <- inData('dataAphidTurns.txt',xnames='angle',
ynames='counts')
x <- tmp$x
y <- tmp$y
plot(x,y,type='s',xlab='angle (radians)')

#yy <- c(8,4,7,10,5,4,5,4,9,7,4,12,9,11,10,9,9,18,12,14,9,6,8,8,10,8,7,6,3,7,6,9,5,8,8,7)
#yyf <- yy/sum(yy)
#barplot(yyf)

n <- sum(y)
m <- length(x)

fitturn <- function(par)
{
b0 <- par[1]
b1 <- par[2]
b2 <- par[3]
mu <- (b0 + b1*cos(b2 + x))/2/pi/b0 #mean function
-sum(dpois(y,2*pi/m*mu,log=T)) # - ln likelihood
}

p0 <- c(5,2,-1) #initial values for (b0,b1,b2)
fit <- nlminb(p0,fitturn,lower=c(2,1,-2),upper=c(10,5,1))

b0 <- fit$par[1]
b1 <- fit$par[2]
b2 <- fit$par[3]
tseq <- seq(-pi,pi,length=100)
lines(tseq, (b0 + b1*cos(b2 + tseq)),col=2) #unnormalized function


#calculate frequency and devided by diff(x)
yf <- y/sum(y)/0.175
plot(x,yf,type='s',xlab='angle (radians)')
lines(tseq, (b0 + b1*cos(b2 + tseq))/2/pi/b0,col=2)

C <- 70
f <- 1/2/pi
cf <- C*f

rturn <- function(n,cf)
{ #generate random deviates
xx <- numeric(0)
nx <- 0
while(nx < n)
{
q <- runif(n-nx,-pi,pi)
a <- (b0 + b1*cos(b2 + q))/cf #unnormalized function
z <- runif(n-nx,0,1)
xxz <- q[z < a]
xx <- c(xx,xxz)
nx <- length(xx)
}
xx
}

#make it to 4 cf values. xz1 is a 4*1000 matrix
n <- 1000
C1 <- c(5,70,220,500)
cf1 <-C1*f
xz1 <- matrix(data=NA,nrow=4,ncol=1000)
layout(matrix(1:4,2,2))

for (i in 1:4)
{
xz1[i,] <- rturn(n,cf1[i])
hist(xz1[i,],nclass=m,pr=T)
lines(tseq,(b0+b1*cos(b2+tseq))/2/pi/b0,col=2)
}


