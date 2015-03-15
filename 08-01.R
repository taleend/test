n <- 12
m <- rnbinom(n,size=1,mu=20) + 1
theta <- rbeta(n,3,10)
y <- rbinom(n,m,theta)

a <- 1
b <- 1
tseq <- seq(0,1,length=1000)
plot(tseq,dbeta(tseq,a + sum(y),b + sum(m - y)),type='l',lwd=2)
lines(tseq,dbeta(tseq,a,b),lty=2)

for(j in 1:n)lines(tseq,dbeta(tseq,a + y[j],b + m[j] - y[j]),col=j)

