par(mfrow=c(2,2))
xnames <- c('diam','trt')
tmp <- inData('dataTreeFACE.txt',xnames,'cones',
na.rm=T,INTERCEPT=F)
x <- tmp$x
y <- tmp$y
wlo <- which(x[,'trt'] == 0) # low CO2 treatment
xlo <- x[wlo,'diam']
ylo <- y[wlo]
a <- 1 #priors
b <- 1
b1 <- b + sum(xlo^2)
a1 <- a + sum(ylo)
p1 <- b1/(b1 + xlo^2)
cmu <- qnbinom(.5,size=a1,prob=p1)
clo <- qnbinom(.025,size=a1,prob=p1)
chi <- qnbinom(.975,size=a1,prob=p1)
plot(ylo,cmu,ylim=c(0,10))
for(i in 1:length(xlo))lines(c(ylo[i],ylo[i]),c(clo[i],chi[i]))
abline(0,1,col=3)
abline(h=mean(ylo),col=2)
