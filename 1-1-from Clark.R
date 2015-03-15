source('clarkFunctions.R')

getZ <- function(km){log(w) - log(km + w)}

updateK <- function(){

  ks <- tnorm(1,0,1,kg,rexp(1,100))

  zn <- getZ(kg)
  zs <- getZ(ks)

  pnow <- sum(dnorm(y,zn + x%*%bg,sqrt(sg),log=T)) + 
          dnorm(kg,priorK,priorKsd,log=T)
  pnew <- sum(dnorm(y,zs + x%*%bg,sqrt(sg),log=T)) + 
          dnorm(ks,priorK,priorKsd,log=T)

  tmp <- acceptMH(pnow,pnew,kg,ks)
  list(kg = tmp$x, accept = tmp$accept)
}


#Here is a simulated data set and Gibbs sampler:

n    <- 1000
beta <- matrix(c(-.2,.1,.3,-.2),4,1)
km   <- .2
sd   <- .1
sg   <- sd^2

loX <- c(1,-10,0)    #int, temp, h20, interaction
hiX <- c(1,20,1)      

x <- simX(n,loX,hiX)     #simulate x and y
w <- x[,3]
x <- cbind(x,(x[,2] - mean(x[,2]))*(x[,3] - mean(x[,3])))
z <- log(w) - log(km + w)
y <- rnorm(n,z + x%*%beta,sd)

k <- 4

priorB   <- beta*0                    #prior means
priorVB  <- diag(1000,k)           #prior covariance
priorIVB <- solve(priorVB)         #prior inv covariance
priorK  <- .1
priorKsd <- 1

s1      <- 1                      #prior values for variance
s2      <- 1

kg <- 5
z <- getZ(kg)
bg <- beta*0

ng <- 10000
ac <- 0

kgibbs <- matrix(0,ng,2)
bgibbs <- matrix(0,ng,k)
for(g in 1:ng){

  bg  <- bUpdateNorm(x,y - z,sg,priorB,priorIVB)
  sg  <- sigmaUpdate(x%*%bg,y - z,s1,s2)
  tmp <- updateK()
  kg  <- tmp$kg
  ac  <- ac + tmp$accept

  bgibbs[g,] <- bg
  kgibbs[g,] <- c(kg,sg)

}

processPars(kgibbs,xtrue=c(km,sd^2),CPLOT=T)
processPars(bgibbs,xtrue=beta,CPLOT=T)
