par(mfrow=c(2,2))
xnames <- c('diam','trt')
tmp <- inData('dataTreeFACE.txt',xnames,'cones',na.rm=T,INTERCEPT=F)
x <- tmp$x
y <- tmp$y
wlo <- which(x[,'trt'] == 0) # low CO2 treatment
xlo <- x[wlo,'diam']
xhi <- x[-wlo,'diam'] # high CO2
ylo <- y[wlo]
yhi <- y[-wlo]

a <- 1 #priors
b <- 1
X2lo <- sum(xlo^2)
Ylo <- sum(ylo)
bseq <- seq(0,.01,length=1000)

prior <- dgamma(bseq,a,b)
post1 <- dgamma(bseq,a + Ylo,b + X2lo)
like1 <- dgamma(bseq,1 + Ylo,X2lo)
plot(bseq,post1,type='l',lty=3)
lines(bseq,prior,lty=2)
lines(bseq,like1)

X2hi <- sum(xhi^2)
Yhi <- sum(yhi)
post2 <- dgamma(bseq,a + Yhi,b + X2hi)
like2 <- dgamma(bseq,1 + Yhi,X2hi)
lines(bseq,like2,col='red')
lines(bseq,post2,col='red',lty=3)

mu_lo <- (a + Ylo)/(b + X2lo)
mu_hi <- (a + Yhi)/(b + X2hi)
sd_lo <- sqrt( (a + Ylo - 1)/((b + X2lo)^2) )
sd_hi <- sqrt( (a + Yhi - 1)/((b + X2hi)^2) )

CI_lo <- qgamma(c(0.5,.025,.975),a + Ylo,b + X2lo)
CI_hi <- qgamma(c(0.5,.025,.975),a + Yhi,b + X2hi)
fit <- cbind(c(mu_lo,mu_hi),c(sd_lo,sd_hi),rbind(CI_lo,CI_hi))
colnames(fit) <- c('mu','sd','0.5','0.025','0.975')
rownames(fit) <- c('ambient','elevated')
signif(fit,4)