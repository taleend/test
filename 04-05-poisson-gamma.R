par(mfrow=c(2,1))

#----produce data-----------------
n <- 10
lamda <- 3.5
y <- rpois(n,lamda)#simulated observations

#----simulation-------------------
Y <- sum(y)
mu <- 2 #prior mean
b <- 1
a <- mu*b
lseq <- seq(0,10,length=1000) #sequence of values
prior <- dgamma(lseq,a,b)
like <- n^(S+1)*lseq^Y*exp(-n*lseq)/gamma(Y+1)
post <- dgamma(lseq,a + Y,b + n)
plot(lseq,post,type='l')
lines(lseq,prior,lty=2)
lines(lseq,like,lty=3)

b <- 100
a <- mu*b
prior <- dgamma(lseq,a,b)
post <- dgamma(lseq,a + Y,b + n)
plot(lseq,post,type='l')
lines(lseq,prior,lty=2)
lines(lseq,like,lty=3)

