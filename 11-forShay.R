
rm(list=ls())

source('clarkFunctions.r')

y2w <- function(y){
	
	for(j in 1:nsize)w[ ,j] <- apply(y[,sizecols[j,]],1,sum)
	w
}

w2q <- function(a,w){
	
	for(j in 1:nsize)q[,j] <- a[j,1] + a[j,2]*w[,j]
	q
}

ysampMVNPois <- function(x,y,z,a,b,s,tg,E){  #sample y's for MVN 1st stage, Pois 2nd stage

  r     <- ncol(y)
  n     <- nrow(y)
  propy <- matrix(rnorm(length(y),y,.003),n,r,byrow=F)
  w     <- y2w(y)
  propw <- y2w(propy)
  muq   <- w2q(a,w)
  mupq  <- w2q(a,propw)
  
  pnow <- rowSums(dnorm(q,muq,sqrt(tg),log=T))        #acoustic data
  pnew <- rowSums(dnorm(q,mupq,sqrt(tg),log=T))
  xb   <- x%*%b
  
  for(i in 1:n){
    pnow[i] <- pnow[i] + mydmvnorm(y[i,],xb[i,],s,log=T) #predictors
    pnew[i] <- pnew[i] + mydmvnorm(propy[i,],xb[i,],s,log=T)
  }
  
  pnow <- pnow + rowSums(dpois(z,E*exp(y),log=T))       #count data
  pnew <- pnew + rowSums(dpois(z,E*exp(propy),log=T))

  aa <- exp(pnew - pnow)
  zz <- runif(n,0,1)
  
  y[zz < aa] <- propy[zz < aa]
  accept <- length(which(zz < aa))

  list(y = y, a = accept)
}


n        <- 200
nspec    <- 2
nsize    <- 3
sizenames <- paste('size',c(1:nsize),sep='')
specnames <- as.vector(outer(sizenames,specnames,paste,sep='-'))


#predictors and count data

#abundance
sizecols <- numeric(0)                      #columns for each size class
for(j in 1:nsize)sizecols <- rbind(sizecols,seq(j,j+nsize,by=nsize))

#acoustic data
w <- q <- matrix(0,n,nsize)
colnames(w) <- colnames(q) <- sizenames
w <- y2w(y)

int   <- runif(nsize,-.1,.1)   #additive, multiplicative bias in acoustic
slope <- runif(nsize,.5,1)
a     <- cbind(int,slope)
rownames(a) <- sizenames
tau   <- .1                    #variance
muq   <- w2q(a,w)
q     <- matrix(rnorm(n*nsize,muq,sqrt(tau)),n,nsize) #acoustic response

colnames(q) <- sizenames

priorA   <- a       #informative prior acoutic parameters
priorIVA <- diag(n*nspec,2)

mut <- tg <- .1
t1  <- n*nspec
t2  <- mut*(t1 - 1)

ag <- a*0


ng     <- 1000
agibbs <- matrix(0,ng,length(ag))
tgibbs <- rep(0,ng)
  
  
  
  wg <- y2w(yg)
  for(j in 1:nsize)ag[j,] <- bUpdateNorm(cbind(rep(1,n),wg[,j]),q[,j],tg,priorA[j,],priorIVA)

  mu <- w2q(ag,wg)
  tg <- sigmaUpdate(as.vector(q),as.vector(mu),t1,t2)

  agibbs[g,] <- as.vector(ag)
  tgibbs[g]  <- tg

accepty <- accept/n/ng


predVsObs(as.vector(b),bgibbs)

predVsObs(as.vector(a),agibbs)

predVsObs(as.vector(sigma),sgibbs)

ymean <- ypred/ng
ysd   <- sqrt(ypred2/ng - ymean^2)