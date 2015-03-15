par(mfrow=c(2,1))
n <- 10
lamda <- 17.8
y <- rpois(n,lamda) #simulated observations
Y <- sum(y)
mu <- 10 #prior mean
b <- 2
a <- mu*b
lseq <- seq(0,30,length=1000) #sequence of values
prior <- dgamma(lseq,a,b)
like <- dgamma(lseq,1 + Y,n)
post <- dgamma(lseq,a + Y,b + n)
plot(lseq,post,type='l')
lines(lseq,prior,lty=2)
lines(lseq,like,lty=3)


a1 <- a + Y
b1 <- b + n
pr <- b1/(1 + b1)

postMean <- a1/b1
yseq <- c(0:30)
plot(yseq,dpois(yseq,postMean),type='s')
lines(yseq,dnbinom(yseq,size=a1,prob=pr),type='s',col=2)
abline(v=qnbinom(c(.025,.975),size=a1,prob=pr),lty=2,
col=2)
