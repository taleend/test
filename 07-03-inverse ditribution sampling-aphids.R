#------------begin with norm distribution--------------------
mu <- 15
sd <- 2
tseq <- seq(1,30, length=100)
layout(matrix(1:2,2,1))
plot(tseq,dnorm(tseq,mu,sd),type='l')
plot(tseq,pnorm(tseq,mu,sd),type='l')

#------draw vertical/horizontal lines on the cdf curve
m <- 20
z <- runif(m,0,1)
theta <- qnorm(z,mu,sd)
for (i in 1:m)
{
 lines(c(theta[i],theta[i]),c(0,z[i]),col=i)  #vertical lines
 lines(c(theta[i],0),c(z[i],z[i]),col=i)  #horizontal lines
}

z <- runif(1000,0,1)
theta <- hist(qnorm(z,mu,sd),nclass=100)
lines(theta$mids,theta$density,type='s')

--------------truncate function---------------
lo <- 0.5
hi <- 2
z <- runif(n, pnorm(lo,mu,sd),pnorm(hi,mu,sd))
p <- qnorm(z,mu,sd)

rtrunc <- function(n,lo,hi,p1,p2,FN)
{
pf <- paste('p',FN,sep='')
qf <- paste('q',FN,sep='')
pf1 <- match.fun(pf)
qf1 <- match.fun(qf)

z <- runif(n,pf1(lo,p1,p2),pf1(hi,p1,p2))
qf1(z,p1,p2)
}

rv <-rtrunc(1000,.3,.7,10,2,'beta')


#-------------aphid data-------------------------------------

#------------get values of parameters---------------
source('clarkFunctions.r')
tmp <- inData('dataAphidTurns.txt',xnames='angle',
ynames='counts')
x <- tmp$x
y <- tmp$y
plot(x,y,type='s',xlab='angle (radians)')

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


#-------------inverse distribution--------
m <- 50
lo <- -pi
hi <- pi
theta <- seq(lo, hi, length=m)

p <- (b0+b1*cos(b2+theta))/2/pi/b0 #mean function

par(mfrow=c(2,1))
P <- cumsum(p)
plot(theta,P,type='l')

z <- runif(1000,P[1],P[length(P)])
tvar <- theta[findInterval(z,P)]
th <- hist(tvar,plot=F,nclass=20)
plot(th$mids,th$density,type='s')


#my method for cumsum
#pp <- rep(0,m)
#pp[1]=p[1]
#for (i in 2:m)
#{
# pp[i]=pp[i-1]+p[i]
#}





