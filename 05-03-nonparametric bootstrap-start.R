setwd("C:/Users/WEn/Documents/01-work/2-Classes/02-2012 spring/Med/Rwk")
source('clarkFunctions.R')
filename <- ('dataFACE.txt')
ynames <- 'cones'
xnames <- c('diam','trt','nfert')
tmp <- inData(filename, xnames, ynames,na.rm = T,INTERCEPT = T)
x <- tmp$x
y <- tmp$y
n <- nrow(x)
p <- ncol(x)

loCO2 <- which(x[,'trt'] == 0)
xlo <- x[loCO2,'diam']
ylo <- y[loCO2]
nlo <- length(ylo)
Y <- sum(ylo)
X2 <- sum(xlo^2)
bmle <- Y/X2
bse <- sqrt(Y)/X2
c(bmle,bse)

par(mfrow=c(4,1))
bseq <- seq(0,.01,length=1000)
plot(bseq,dnorm(bseq,bmle,bse),type='l')

#----bootstrap-----------------
nboot <- 2000
bvals <- rep(0,nboot)
for (b in 1:nboot)
{
 bindex <- sample(c(1:nlo),replace=T)
 ytmp <- ylo[bindex]
 xtmp <- xlo[bindex]
 bvals[b] <- sum(ytmp)/sum(xtmp^2) #---MLE
}
#----work on bval----------------
sd(bvals)
hist(bvals,pr=T,nclass=20,xlim=c(0,0.01))
lines(density(bvals))
cil <- quantile(bvals,c(0.025,0.975))
abline(v=cil,lty=2)

