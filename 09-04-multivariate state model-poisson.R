source('clarkFunctions.R')

#---data producing--------------------
ns <- 5 - 1 #no. species
n <- 8 #no. sites
p <- 3 #no. columns in X
nt <- 50 #no. time incremetns
rangeB <- c(-.8,.8) #covariate parameters
rangeA <- c(-.5,.5) #species interaction parameters
tmp <- simMVNData(n,p,ns,rangeB) #simulate MVN data
x <- tmp$x
b <- tmp$b #parameters for covariates
sigma <- tmp$s
a <- sigma*0 + runif(length(sigma),rangeA[1],rangeA[2]) #spp pars

a[abs(a) < .3] <- 0 #make some equal zero

y1 <- tmp$y
specnames <- tmp$specnames
covars <- colnames(x)[-1]
sitenames <- paste('site',c(1:n),sep='-')
timenames <- paste('t',c(1:nt),sep='-')
y <- matrix(0,nt*n,ns) #time by species-site
colnames(y) <- specnames
rownames(y) <- outer(timenames,sitenames,paste,sep='_')
xtime <- matrix(rep(x,each=nt),nt*n,ncol(x))
colnames(xtime) <- colnames(x)
rownames(xtime) <- rownames(y)

tindex <- c(0:(n-1))*nt
y[(tindex+1),] <- as.vector(y1)

for(t in 2:nt)
{
ti <- t + tindex
xtime[ti-1,-1] <- rnorm(n*(p-1),xtime[ti-1,-1],.5)
mu <- xtime[ti-1,]%*%b + y[ti-1,]%*%a
y[ti,] <- myrmvnorm(n,mu,sigma)
}


#----model simulation------------------
#Poisson likelihood
likelihood <- 'dpois'
yobs <- matrix(rpois(length(y),exp(y)),nrow(y),ncol(y))
tmp <- initialStatesSS('dpois',yobs,mus,mw)
yg <- tmp$x


prior.W <- diag(1,ns) #prior covariance matrix, for Wishart
prior.WDF <- p + n*ns*nt
pB <- matrix(0,p,ns)

#initial values
sg <- prior.W
bg <- b*0
ag <- a*0
ng <- 3000
bgibbs <- matrix(0,ng,p*ns)
colnames(bgibbs) <- as.vector(outer(colnames(x),colnames(y),
paste,sep='_'))

agibbs <- matrix(0,ng,ns*ns)
colnames(agibbs) <- as.vector(outer(colnames(y),colnames(y),
paste,sep='_'))
sgibbs <- matrix(0,ng,ns*ns)
colnames(sgibbs) <- outer(colnames(y),colnames(y),paste,sep='_')
ypred <- y*0
ypred2 <- ypred
accept <- 0
t1 <- tindex+1
t2 <- tindex+nt

for(g in 1:ng)
{
bg <- bUpdateMVNorm(xtime[-t2,],yg[-t1,] - yg[-t2,]%*%ag,
bg,sg)
ag <- bUpdateMVNorm(yg[-t2,],yg[-t1,] - xtime[-t2,]%*%bg,
ag,sg)
S <- crossprod(yg[-t1,]-xtime[-t2,]%*%bg-yg[-t2,]%*%ag)
sinv <- rwish(ns + prior.WDF, solve(S + prior.W*prior.WDF))
sg <- solve(sinv)
tmp <- updateSSstatesMVN(n,nt,tindex,xtime,
yg,yobs,bg,ag,obsError)
yg <- tmp$y
accept <- tmp$accept
if(likelihood == 'dnorm')obsError <- sigmaUpdate(yg,yobs,s1,s2)
bgibbs[g,] <- as.vector(bg)
agibbs[g,] <- as.vector(ag)
sgibbs[g,] <- as.vector(sg)
ypred <- ypred + yg
ypred2 <- ypred2 + yg^2
print(g)
}

processPars(sgibbs,xtrue=as.vector(sigma),CPLOT=T)
processPars(bgibbs,xtrue=as.vector(b),sigOnly=F,CPLOT=T,burnin=100)
processPars(agibbs,xtrue=as.vector(a),sigOnly=F,CPLOT=T,burnin=100)

par(mfrow=c(2,2))
predVsObs(as.vector(b),bgibbs)
predVsObs(as.vector(a),agibbs)
predVsObs(as.vector(sigma),sgibbs)
ymean <- ypred/ng
ysd <- sqrt(ypred2/ng - ymean^2)
ylo <- ymean - 1.96*ysd
yhi <- ymean + 1.96*ysd
plot(yobs,ymean)
for(s in 1:ns){
points(yobs[,s],ymean[,s],col=s)
for(i in 1:(nt*n)){
lines(c(yobs[i,s],yobs[i,s]),c(ylo[i,s],yhi[i,s]),col=s)
}
}
abline(0,1,lty=2)
 






