source('clarkFunctions.R')
ynames <- 'cones'
xnames <- 'diam'
tmp <- inData("dataTreeFACE.txt",xnames,ynames,na.rm=T,INTERCEPT=T)
x <- tmp$x            #design matrix
diam <- x[,'diam']
y <- tmp$y #response
y[y > 1] <- 1         #binary response
n <- nrow(x)
plot(diam,jitter(y))

par0 <- c(-5,.2)
lo <- c(-10,-10)
hi <- c(10,10)

likeBernLogit <- function(pars)
{
g <- x %*% pars
theta <- inv.logit(g)
s <- y*log(theta ) + (1 - y)*log(1 - theta )
-sum(s)
}

out <- nlminb(par0,likeBernLogit,lower=lo,upper=hi)
b1 <- out$par
L1 <- out$objective


#----a null model-------------------
theta <- sum(y)/n
b0 <- logit(theta)
L0 <- -sum(dbinom(y,1,theta,log=T))
abline(h=theta,lty=2,col=2)
abline(h=c(0,1),col=3)

dseq <- seq(min(diam),max(diam),length=100) #sequence diameters
y1 <- b1[1] + b1[2]*dseq #predicted logit theta
lines(dseq,inv.logit(y1),col=2) #invert and plot

