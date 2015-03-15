n <- 10
a <- rexp(n,.1) #simulated ages
ma <- mean(a)
mle <- 1/ma
rseq <- seq(.5,2,length=100)*mle
normlik <- n*(log(ma) + log(rseq)) + n*(1 - rseq*ma) #normed likelihood
par(mfrow=c(3,1))
plot(rseq,normlik,type='l')
abline(v=mle,lty=2)
D <- -2*normlik         #Deviance
plot(rseq,D,type='l')
abline(v=mle,lty=2)
abline(h=c(qchisq(0.05,1),qchisq(0.95,1)),lty=2,col=2)    #value of D at P = 0.05
P <- 1 - pchisq(D,1)    #P value
plot(rseq,P,type='l')
abline(v=mle,lty=2)
abline(h=.05,lty=2)