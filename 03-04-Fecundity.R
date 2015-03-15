source(clarkFunctions.R)
ynames <- 'cones'
xnames <- c('diam','trt')
tmp <- inData('dataFACE.txt', xnames, ynames, tname = 'tstat',na.rm = T,INTERCEPT = F)
x <- tmp$x
y <- tmp$y          #response is no. cones
mstat <- tmp$t      #maturation index
mindex <- which(mstat == 1)        #index for mature
crange <- range(y)                 #range of counts
cvalues <- c(0:(crange[2]+1))
wlo <- which(x[,'trt'] == 0 & mstat == 1) #index for mature tree...
whi <- which(x[,'trt'] == 1 & mstat == 1) #...with lo and hi CO2
slo <- summary(y[wlo])
shi <- summary(y[whi])
rbind(slo,shi)
hist(y[wlo],probability=T,breaks=c(0:crange[2]))#hist of y for low CO2
hi <- hist(y[whi],breaks=cvalues,plot=F)
lines(cvalues[-length(cvalues)],hi$density,type='s',col='red')

#--plot diameter aganst y (number of cones)
par(mfrow=c(3,1))
plot(x[,'diam'],y) #all by diameter
#diameter by CO2 trt
plot(x[,'diam'],y,ylim=crange,col=(x[,'trt']+1))
#log scale (will not show zeros)
plot(x[,'diam'],y,col=(x[,'trt']+1),log='xy')


#----MLE for the low trt---
Y <- sum(y[wlo])
X2 <- sum(x[wlo,'diam']^2)
bmle <- Y/X2
bse <- sqrt(Y)/X2
c(bmle,bse)
#----likelihood----------
L1 <- sum(dpois(y[wlo],bmle*x[wlo,'diam']^2,log=T),na.rm=T)
#----plot--------
par(mfrow=c(3,1))
bseq <- seq(0,.01,length=1000)
plot(bseq,dnorm(bseq,bmle,bse),type='l')
ci95 <- c(bmle - 1.96*bse,bmle + 1.96*bse)
abline(v=ci95,lty=2)

#----for mature trees-----------

Y <- sum(y[mindex])
X2 <- sum(x[mindex,'diam']^2)
bmle <- Y/X2
bse <- sqrt(Y)/X2
c(bmle,bse)
#--likelihood
L1 <- sum(dpois(y[mindex],bmle*x[mindex,'diam']^2,log=T))
plot(bseq,dnorm(bseq,bmle,bse),type='l')
ci95 <- c(bmle - 1.96*bse,bmle + 1.96*bse)
abline(v=ci95,lty=2)


#-----for all tree and years----------

tmp <- inData('dataFACE.txt', xnames, ynames, tname = 't',iname = 'i', na.rm = T, INTERCEPT = F)
year <- tmp$t #indicator for year
tree <- tmp$i #indicator for tree
mattable <- treebytime(tree,year,mstat)
conetable <- treebytime(tree,year,y) #tree by yr table
yrmean <- apply(conetable*mattable,2,mean,na.rm=T)
treemean <- apply(conetable*mattable,1,mean,na.rm=T)

par(mfrow=c(2,1))
#----jitter: add a small amount of noise to the data----
plot(jitter(year),jitter(y))
lines(yrmean,col='red')
#----hist of treemean-----------------
hist(treemean,breaks=c(0:40),main='Individual means')
nyr <- ncol(conetable) #no. years
ntree <- nrow(conetable) #no. trees
yrvec <- c(1:nyr)
diamtable <- treebytime(tree,year,x[,'diam'])
bmle2 <- apply(conetable*mattable,2,mean,na.rm=T)/
apply(diamtable^2*mattable,2,mean,na.rm=T)
L2 <- sum(dpois(y[mindex],bmle2[year[mindex]]*
x[mindex,'diam']^2,log=T))
signif(bmle2,3)
#----plot again
par(mfrow=c(2,2))
pred.y1 <- bmle*x[,'diam']^2
pred.y2 <- bmle2[year]*x[,'diam']^2
plot(y,pred.y1);abline(0,1,lty=2)
plot(y,pred.y2);abline(0,1,lty=2)
#CO2 effect
bmleCO2 <- c(sum(y[wlo])/sum(x[wlo,'diam']^2),
(sum(y[whi])/sum(x[whi,'diam']^2)))
L3 <- sum(dpois(y,bmleCO2[x[,'trt']+1]*x[,'diam']^2,log=T))
pred.y3 <- bmleCO2[x[,'trt'] + 1]*x[,'diam']^2
plot(y,pred.y3);abline(0,1,lty=2)
