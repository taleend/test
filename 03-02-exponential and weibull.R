#----compare geometric and exponential distribution--

nt <- 100
tDisc <- c(1:nt)
tCont <- seq(0,nt,length=1000)
lambda <- .1
d <- dgeom(tDisc-1,lambda)
c <- dexp(tCont,lambda)
plot(tDisc-.5,d,type='s')
lines(tCont,c,col=2)


#--exponential and weibull distribution----
n <- 100
lambda <- .05
c <- 3

fexp <- rexp(n,lambda)          #sample exponential
fwei <- rweibull(n,c,1/lambda)  #sample Weibull
tseq <- seq(0,200,by=1)         #sequence for t

hist(fexp,probability=T,breaks=tseq,ylim=c(0,.2))
wcurve <- hist(fwei,breaks=tseq,plot=F)
lines(wcurve$mids,wcurve$density,type='s',col='red')
#----add means to the two distributions-----
abline(v=1/lambda,lty=2)
abline(v=1/lambda*gamma(1/c+1),col='red')

#---plot risk over time-----
abline(h=lambda,lty=2)
lines(tseq,c*lambda^c*tseq^(c-1),col=4,lty=2)

#--plot density--------
plot(tseq,dexp(tseq,lambda),type='l',ylim=c(0,.15))
lines(tseq,dweibull(tseq,c,1/lambda),col='red')

#--Example 3.4 MLE for continuous data-fire intervals--

fires <- read.table('dataFires.txt',header=T)
dim(fires)
#int <- fires$interval
#num <- fires$number
#plot(int,num,type='s')

years <- c(0:max(fires[,'interval']))
nfire <- years*0
nfire[fires[,'interval']] <- fires[,'number']
plot(years,nfire,type='s')

#--estimate the parameters for Weibull dist.---

likeWeibull <- function(param)
{
-sum(dweibull(y,param[1],1/param[2],log=T))
}
par0 <- c(.1,1)
out <- nlminb(par0,likeWeibull,lower=c(.01,.01),upper=c(10,3))

lw <- out$par[2]
cw <- out$par[1]
hazard <- cw*lw^cw*tseq^(cw-1)
fweib <- dweibull(tseq,cw,1/lw) #density of events
plot(years,nfire/sum(nfire),type='s',ylim=c(0,0.2))
lines(tseq,hazard,lty=2,col='red')
lines(tseq,fweib,col=3)
abline(h=1/mean(y),lty=2,col=4)
lines(tseq,dexp(tseq,lambda),col=5,lty=2)