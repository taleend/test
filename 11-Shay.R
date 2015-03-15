source('clarkFunctions.R')

n        <- 100
r        <- 4               #no. spp
k        <- 3               # num of predictors
effort_i <- runif(n,.5,2)   #effort at location i
detect_s <- runif(r,.5,1)   #detection for spp s
effort   <- matrix(effort_i,n,r)*matrix(detect_s,n,r,byrow=T)
sigma    <- diag(.1,r)
specnames <- paste('s',c(1:r),sep='')
sampnames <- paste('i',c(1:n),sep='')

b <- bmultiProp(priorB*0 + 1,diag(.5,k*r))$c

loX <- c(1,-1,0)       #range of covariates, length = k
hiX <- c(1,1,1)

LIKE <- 'mvnorm'          #'norm', 'pois', 'binom', or 'multinom'

x <- simX(n,loX,hiX)     #simulate x and y
y <- simY(x,b,LIKE,r,sigma=sigma)
z <- matrix(rpois(n*r,effort*exp(y)),n,r)

rownames(x) <- rownames(y) <- sampnames
colnames(x) <- paste('v',c(1:k),sep='')
colnames(y) <- specnames

prior.W   <- diag(.1,r)  #prior covariance matrix
prior.WDF <- k + round(n,0)

ng     <- 2000
bg <- b*0
sg <- diag(1,r)
yg <- log(z + .1)

bgibbs <- matrix(0,ng,length(bg))
colnames(bgibbs) <- as.vector(outer(colnames(x),colnames(y),paste,sep='-'))
sgibbs <- matrix(0,ng,length(sg))
colnames(sgibbs) <- outer(colnames(y),colnames(y),paste,sep='-')

ypred  <- z*0
names(ypred) <- as.vector(outer(rownames(x),colnames(y),
                                paste,sep='-'))
ypred2 <- ypred
accept <- 0


for(g in 1:ng){

  bg    <- bMVNorm(x,yg,bg,sg)
  sinv  <- rwish(n + prior.WDF,crossprod(yg - x%*%bg) + prior.W*prior.WDF )
  sg    <- solve(sinv)
  tmp   <- ysampMVNPois(x,yg,z,bg,sg,effort)
  yg    <- tmp$y
  accept <- accept + tmp$a

  bgibbs[g,] <- as.vector(bg)
  sgibbs[g,] <- as.vector(sg)

  ypred  <- ypred + y
  ypred2 <- ypred2 + y^2

}

processPars(bgibbs,as.vector(b),CPLOT=T)
processPars(bgibbs,as.vector(b),DPLOT=T)

plotObsPred(as.vector(b),apply(bgibbs,2,mean),apply(bgibbs,2,sd))

plotObsPred(as.vector(b),apply(bgibbs,2,mean),apply(bgibbs,2,sd))

To examine the residual correlation I execute the following:
covS <- matrix(apply(sgibbs,2,mean),r,r) #mean cov estimate
rownames(covS) <- specnames
colnames(covS) <- specnames
varS <- diag(covS)           
sdS  <- sqrt(varS)                           #sd estimates
corS <- covS/matrix(sdS,r,r,byrow=T)/
             matrix(sdS,r,r,byrow=F)       #cor estimates
