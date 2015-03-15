mu <- 0.34
sd <- 2
n <-60
m <- 10000
y <- rep(NA,n)
s <- rep(NA,m)


ssd <- rep(NA,100)

for (i in 1:m)
{
y=rnorm(n,mu,sd)
s[i]=sd(y)*(n-1)/n
}

hist(s,breaks=seq(1,3,0.1),ylim=c(0,3),pr=T)
ss <- seq(1,3,length=100)
lines(ss,dnorm(ss,mean(s),sd(s)),col=2)