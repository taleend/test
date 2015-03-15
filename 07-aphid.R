source('clarkFunctions.R')
tmp <- inData('dataAphidTurns.txt', xnames='angle', ynames='counts')

x <- tmp$x
y <- tmp$y
plot(x,y,type='s', xlab='angle(radians)')

n <- sum(y)
m <- length(x)

fitturn <- function(par){
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
lines(tseq, (b0 + b1*cos(b2 + tseq))) #unnormalized function

C <- 70
f <- 1/2/pi #the uniform density
cf <- C*f #the basis for rejection

rturn <- function(n,cf){ #generate random deviates
x <- numeric(0)
nx <- 0
while(nx < n){
q <- runif(n-nx,-pi,pi)
a <- (b0 + b1*cos(b2 + q))/cf #unnormalized function
z <- runif(n-nx,0,1)
xz <- q[z < a]
x <- c(x,xz)
nx <- length(x)
}
x
}
n <- 1000
xz <- rturn(n,cf)
hist(xz,nclass=m,probability=T)
lines(tseq, (b0 + b1*cos(b2 + tseq))/2/pi/b0) #normalized function

