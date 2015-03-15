library(mvtnorm) #needed for MVNormal
par(mfrow=c(2,1))
n <- 100
sigma <- 5 #variance is known
xy <- rmvnorm(n,c(0,0),diag(sigma,2)) #simulated (x,y)
r <- sqrt(xy[,1]^2 + xy[,2]^2) #transform to distance
a <- 10 #prior parameter values
b <- 40
sseq <- seq(.001,30,length=1000) #sequence of values
S <- sum(r^2) #sum of squares
prior <- dinvGamma(sseq,a,b)
like <- dinvGamma(sseq,n-1,S/2)
post <- dinvGamma(sseq,a + n,b + S/2)
plot(sseq,post,type='l')
lines(sseq,prior,lty=2)
lines(sseq,like,lty=3)


spost <- (b + S/2)/(a + n + 1) #most probable value
kfit <- sseq/spost*exp(-.5*sseq^2/spost) #fitted model
hist(r,probability=T,xlim=c(0,20),ylim=c(0,0.3))
lines(sseq,kfit)