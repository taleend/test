#----confidence interval of introductory stats
mu <- 0
y <- rnorm(1,mu,1)
mseq <- seq(-8,8,length=100)
like <- dnorm(mseq,y,1)
plot(mseq,like,type='l')
abline(v=y,lty=2)

alpha <- .05
ci <- qnorm(c(alpha/2,1-alpha/2),y)
abline(v=ci,lty=3)
lines(mseq,dnorm(mseq,ci[1]),col=2)
lines(mseq,dnorm(mseq,ci[2]),col=3)