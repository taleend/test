library(mvtnorm)
n <- 100 #no. wells
J <- 20 #no. sources
alpha <- 20
p <- 2 #no. predictors + intercept
beta <- matrix(c(1,2),p,1)
sigma <- .6^2
tau <- .1^2
x <- matrix(1,J,p) #design matrix
x[,2] <- runif(J,0,5)
f <- rnorm(J,x%*%beta,sqrt(sigma)) #sources
ix <- runif(n,0,100) #locations of wells & sources
iy <- runif(n,0,100)
jx <- runif(J,0,100)
jy <- runif(J,0,100)
dij <- distmat(jx,jy,ix,iy) #distance matrix
k <- getKernel(alpha) #kernel matrix
S <- getKernSigma(k,sigma,tau) #covariance matrix
y <- t(myrmvnorm(1,t(k%*%x%*%beta),S))

b <- matrix(0,p,1)
B <- diag(100,p)
BInv <- solve(B)
amax <- 100
smax <- .2
mtau <- .1^2
t1 <- n
t2 <- mtau*(t1 - 1)
bg <- b #initial values
sg <- tg <- .1
ag <- 10

vas <- diag(c(.1,.001^2)) #proposal for alpha, sigma
ng <- 5000
bgibbs <- matrix(0,ng,p)
pgibbs <- matrix(0,ng,3)
colnames(pgibbs) <- c('a','s2','t2')
fhat <- matrix(0,ng,J)
g0 <- 1
for(g in 1:ng){
bg <- updateKernB()
tg <- updateKernT()
tmp <- updateKernAS()
ag <- tmp[1]
sg <- tmp[2]

fhat[g,] <- rnorm(J,x%*%bg,sqrt(sg)) #sources
if(g %in% c(100,300,500)){
vas <- var(pgibbs[g0:g,c('a','s2')])
g0 <- g
}
bgibbs[g,] <- bg
pgibbs[g,] <- c(ag,sg,tg)
}
processPars(bgibbs,beta,CPLOT=T)
processPars(pgibbs,c(alpha,sigma,tau),CPLOT=T)
fci <- predVsObs(f,fhat)
fplot <- 2*(fci[1,] - min(fci[1,]))/(max(fci[1,]) - min(fci[1,]))
yplot <- 2*(y - min(y))/(max(y) - min(y))
symbols(ix,iy,circles=yplot,inches=F)
symbols(jx,jy,circles=fplot,inches=F,
add=T,bg=1)
dii <- distmat(ix,iy,ix,iy) #distance matrix
kMat <- getKernel(mean(pgibbs[,'a']))
sigMat <- getKernSigma(kMat,mean(pgibbs[,'s2']),mean(pgibbs[,'t2']))
dist2 <- cov2Dist(sigMat)
plot(dii,dist2,ylim=c(0,1))
ns <- 100
dseq <- seq(0,100,length=ns)
khat <- exp(-(matrix(dseq,ng,ns,byrow=T)/
matrix(pgibbs[,'a'],ng,ns))^2)
khat <- apply(khat,2,quantile,c(.5,.025,.975))
lines(dseq,khat[1,])
for(j in 2:3)lines(dseq,khat[j,],lty=2)
