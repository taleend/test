source('clarkFunctions.R')
filename <- 'dataTreeFACE.txt'
xnames <- c('diam','trt','nfert')
tmp <- inData(filename,xnames,'dnow',na.rm = T,INTERCEPT = T)
x <- tmp$x
y <- tmp$y
n <- nrow(x)
p <- ncol(x)
par(mfrow=c(2,2))
for(j in 2:p){
plot(x[,j],y)
title(colnames(x)[j])
}

#----determine the rank--------------------
qr(x)$rank
cov(x[,'diam'],y)/var(x[,'diam'])


#----a function for any x and y--------
linearModel <- function(x,y){ #objects for a linear regression
n <- length(y)
b <- solve(crossprod(x))%*%crossprod(x,y) #coefficients
py <- x%*%b #predictions
plot(y,py,xlab='Observed',ylab='Predicted')
abline(0,1,lty=2)
s <- (crossprod(y - py)/(n - ncol(x)))[1] #variance estimate
varb <- s*solve(crossprod(x)) #parameter covar
se <- sqrt(diag(varb)) #standard errors

ci <- matrix(NA,ncol(x),2);colnames(ci) <- c('.025','.975')
for(j in 1:ncol(x))ci[j,] <- c(b[j] - 1.96*se[j],b[j] + 1.96*se[j])
pmat <- cbind(b,se,ci)
colnames(pmat)[1] <- 'estimate'
list(coeff = pmat, shat = s, parVar = varb)
}
#----end of the function---------------

linearModel(x,y)
cor(x[,-1])


#----interaction----
dXtrt <- x[,'diam']*x[,'trt']
dXnfert <- x[,'diam']*x[,'nfert']
x <- cbind(x,dXtrt,dXnfert)
p <- ncol(x) #more variables in x


#----A plot of the normal distribution described by the estimated...
#----...mean and variance for intercept and diameter slope----
library(cluster)
bhat <- tmp$coeff[1:2,1]
varb <- tmp$parVar[1:2,1:2]
ci <- tmp$coeff[1:2,c('.025','.975')]
plot(ellipsoidPoints(varb,d2 = qchisq(0.95, df = 2),loc= bhat),
type='l')
lines(ellipsoidPoints(varb,d2 = qchisq(0.68, df = 2),loc= bhat))
abline(v=ci[1,],lty=2,col=2)
abline(h=ci[2,],lty=2,col=3)