cv1 <- matrix(data=c(1,.5,.5,1),nrow=2,ncol=2)
mu <- matrix(c(0,0),1,2)
thetanow <- matrix(c(-1,1),2,1)
result <- rep(NA,1000)
for (i in 1:1000)
{
thetanew <- t(rmvnorm(1,me1,cv1))
pnow <- dmvnorm(thetanow,mu,cv1)
pnew <- dmvnorm(thetanew,mu,cv1)
a <- pnew/pnow
if (a>1)
{
result[i] <- pnew
thetanow <- thetanew
}
}