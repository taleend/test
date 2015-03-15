
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
nsize    <- 3r        <- nspec*nsize     #no. spp times 3 size classesk        <- 3               #no. predictorseffort_i <- runif(n,.5,2)   #effort at location idetect_s <- runif(r,.5,1)   #detection for spp seffort   <- matrix(effort_i,n,r)*matrix(detect_s,n,r,byrow=T)sigma    <- diag(.1,r)specnames <- paste('s',c(1:nspec),sep='')
sizenames <- paste('size',c(1:nsize),sep='')
specnames <- as.vector(outer(sizenames,specnames,paste,sep='-'))sampnames <- paste('i',c(1:n),sep='')
b <- bmultiProp(matrix(0,k,r) + 1,diag(.5,k*r))$cloX <- c(1,-1,0)       #range of covariates, length = khiX <- c(1,1,1)

#predictors and count datax <- simX(n,loX,hiX)     #simulate x and yy <- simY(x,b,'mvnorm',r,sigma=sigma)       #log abundancez <- matrix(rpois(n*r,effort*exp(y)),n,r)   #counts

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
rownames(x) <- rownames(y) <- sampnamescolnames(x) <- paste('v',c(1:k),sep='')colnames(y) <- specnames
colnames(q) <- sizenamesprior.W   <- diag(.1,r)  #prior covariance matrixprior.WDF <- r + 2

priorA   <- a       #informative prior acoutic parameters
priorIVA <- diag(n*nspec,2)

mut <- tg <- .1
t1  <- n*nspec
t2  <- mut*(t1 - 1)
bg <- b*0
ag <- a*0sg <- diag(1,r)yg <- log(z + .1)


ng     <- 1000bgibbs <- matrix(0,ng,length(bg))
agibbs <- matrix(0,ng,length(ag))colnames(bgibbs) <- as.vector(outer(colnames(x),colnames(y),paste,sep='-'))sgibbs <- matrix(0,ng,length(sg))
tgibbs <- rep(0,ng)colnames(sgibbs) <- outer(colnames(y),colnames(y),paste,sep='-')ypred  <- z*0names(ypred) <- as.vector(outer(rownames(x),colnames(y),                                paste,sep='-'))ypred2 <- ypredaccept <- 0for(g in 1:ng){  bg    <- bUpdateMVNorm(x,yg,bg,sg)
    sinv  <- rwish(n + prior.WDF,solve(crossprod(yg - x%*%bg) + prior.W*prior.WDF ))  sg    <- solve(sinv)
    tmp   <- ysampMVNPois(x,yg,z,ag,bg,sg,tg,effort)  yg    <- tmp$y  accept <- accept + tmp$a
  
  wg <- y2w(yg)
  for(j in 1:nsize)ag[j,] <- bUpdateNorm(cbind(rep(1,n),wg[,j]),q[,j],tg,priorA[j,],priorIVA)

  mu <- w2q(ag,wg)
  tg <- sigmaUpdate(as.vector(q),as.vector(mu),t1,t2)

  agibbs[g,] <- as.vector(ag)  bgibbs[g,] <- as.vector(bg)  sgibbs[g,] <- as.vector(sg)
  tgibbs[g]  <- tg  ypred  <- ypred + yg  ypred2 <- ypred2 + yg^2}

accepty <- accept/n/ngprocessPars(bgibbs[,1:9],as.vector(b)[1:9],CPLOT=T)
processPars(agibbs,as.vector(a),CPLOT=T)

predVsObs(as.vector(b),bgibbs)

predVsObs(as.vector(a),agibbs)

predVsObs(as.vector(sigma),sgibbs)

ymean <- ypred/ng
ysd   <- sqrt(ypred2/ng - ymean^2)plotObsPred(as.vector(y),ymean,ysd)covS <- matrix(apply(sgibbs,2,mean),r,r) #mean cov estimaterownames(covS) <- specnamescolnames(covS) <- specnamesvarS <- diag(covS)           sdS  <- sqrt(varS)                           #sd estimatescorS <- covS/matrix(sdS,r,r,byrow=T)/             matrix(sdS,r,r,byrow=F)       #cor estimates
