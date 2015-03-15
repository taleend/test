source('clarkFunctions.R')

#----Produce data----------------------
n <- 12
m <- rnbinom(n,size=1,mu=20) + 1
theta <- rbeta(n,3,10)
y <- rbinom(n,m,theta)

#----simulation-------------------------
a <- 1
b <- 1

a1 <- .1
a2 <- .1
b1 <- .1
b2 <- .1
ag <- bg <- a
ng <- 100000
tgibbs <- matrix(NA,ng,n)
pgibbs <- matrix(NA,ng,2)

for(g in 1:ng)
{

tg <- rbeta(n,ag + y,bg + m - y)
as <- tnorm(1,0,10,ag,.4)
bs <- tnorm(1,0,10,bg,1.2)

pnow <- sum(dbeta(tg,ag + y,bg + m - y,log=T)) +
dgamma(ag,a1,a2,log=T) + dgamma(bg,b1,b2,log=T)

pnew <- sum(dbeta(tg,as + y,bs + m - y,log=T)) +
dgamma(as,a1,a2,log=T) + dgamma(bs,b1,b2,log=T)

tmp <- acceptMH(pnow,pnew,c(ag,bg),c(as,bs),BLOCK=T)$x

ag <- tmp[1]
bg <- tmp[2]
tgibbs[g,] <- tg
pgibbs[g,] <- tmp

}

x11()
processPars(pgibbs,CPLOT=T)
x11
plotObsPred(theta,apply(tgibbs,2,mean),apply(tgibbs,2,sd))


ab <- apply(pgibbs,2,mean)
x11()
tseq <- seq(0,1,length=1000)
par(mfrow=c(3,4))

for(j in 1:n){
d0 <- dbeta(tseq,a + sum(y),b + sum(m - y)) # posterior from option A
di <- dbeta(tseq,a + y[j],b + m[j] - y[j])  # posterior from option B
plot(tseq,d0,type='l',lwd=2.5)
lines(tseq,di)

tmp <- density(tgibbs[,j])
lines(tmp$x,tmp$y,col=2)
lines(tseq,dbeta(tseq,a,b),lty=2,col=3) #Option C: prior for a and b
lines(tseq,dbeta(tseq,mean(pgibbs[,1]),mean(pgibbs[,2])),
col=4,lty=2)                            #Option C: posterior for a and b
lines(tseq,dbeta(tseq,3,10),col=5,lty=2) # the original dist for theta
title(paste('m =',m[j],sep=' '))
}