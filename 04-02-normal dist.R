x11()
layout(matrix(1:2,2,1))
layout.show(2)

#----an observation y and its likelihood----
sigma <- 1
y <- rnorm(1,0,sqrt(sigma))  #--simulate 1 observation
mseq <- seq(-5,5,length=100) 
plot(mseq,dnorm(mseq,y,sqrt(sigma)),type='l',xlab='Parameter mu', ylim=c(0,1))
abline(v=y,lty=2)

#----a prior density of the mean mu, a normal density too, with mu0 and tau---
mu0 <- -1
tau <- .5
lines(mseq,dnorm(mseq,mu0,sqrt(tau)),col=2)
abline(v=mu0,lty=3,col=2)

#----posterior formula----------
v <- y/sigma + mu0/tau
V <- 1/(1/sigma + 1/tau)
lines(mseq, dnorm(mseq,V*v,sqrt(V)),col=3)
legend(3,0.9,pch='-',col=c(1,2,3),c('likelihood','prior', 'posterior'))
abline(v=V*v,col=3,lty=3)

n <- 15
sigma <- 2
mu <- 1
y <- rnorm(n,mu,sqrt(sigma))
mu0 <- 0
tau <- 1
V <- 1/(n/sigma + 1/tau)
v <- sum(y)/sigma + mu0/tau
prior <- dnorm(mseq,mu0,sqrt(tau))
like <- dnorm(mseq,mean(y),sqrt(sigma/n))
post <- dnorm(mseq,V*v,sqrt(V))
plot(mseq,post,type='l',lty=3,col=3)
lines(mseq,prior,lty=2,col=2)
lines(mseq,like)

#----estimate variance, with mean known---------
source('clarkFunctions.r')
a <- 1 #prior parameter values
b <- 1
sseq <- seq(.001,5,length=1000) #sequence of values
S <- sum((y - mu)^2) #sum of squares
prior <- dinvGamma(sseq,a,b)
like <- (S/2)^(n-1)*exp(-S/2/sseq)/gamma(n-1)/sseq^n
post <- dinvGamma(sseq,a + n/2,b + S/2)
x11()
plot(sseq,post,type='l',ylim=(c(0,0.65)),lty=3)
lines(sseq,prior,lty=2,col=2)
lines(sseq,like,col=3)

