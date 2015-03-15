source('clarkFunctions.R')

getZ <- function(km){log(-w) - log(-km - w)}
updateK <- function()
{
  ks <- -tnorm(1,0,50,-kg,rexp(1,2))

  zn <- getZ(kg)
  zs <- getZ(ks)

  pnow <- sum(dnorm(y,zn + x%*%bg,sqrt(sg),log=T)) + 
          dnorm(kg,priorK,priorKsd,log=T)
  pnew <- sum(dnorm(y,zs + x%*%bg,sqrt(sg),log=T)) + 
          dnorm(ks,priorK,priorKsd,log=T)

  tmp <- acceptMH(pnow,pnew,kg,ks)
  list(kg = tmp$x, accept = tmp$accept)
}

#----conversion of parameters: from V, beta1, beta2, to b0-b3
V <- 13.58
be1 <- 3.63
be2 <- 0.053

b0 <- log(V) - 2.4*log(be1) - 24*be2
b1 <- 0.1*log(be1)+be2
b2 <- -2.4*be2
b3 <- 0.1*be2


#----Here is a simulated data set and Gibbs sampler:
n    <- 1000
beta <- matrix(c(b0,b1,b2,b3),4,1) #--values of parameters from MLE
km   <- -7.2                                 #--value of km from MLE
sd   <- .5
sg   <- sd^2

loX <- c(1,3,-40)    #int, temp, h2o, interaction
hiX <- c(1,25,-0.001)      

x <- simX(n,loX,hiX)     #simulate x and y
t <- x[,2]
w <- x[,3]
x <- cbind(x,(x[,2] - mean(x[,2]))*(x[,3] - mean(x[,3])))
z <- getZ(km)
y <- rnorm(n,z + x%*%beta,sd)

#----parameter estimation------------------
k <- 4
priorB   <- matrix(0,4,1)          #prior means
priorVB  <- diag(1000,k)           #prior covariance
priorIVB <- solve(priorVB)         #prior inv covariance
priorK  <- .1
priorKsd <- 1

s1      <- 1                      #prior values for variance
s2      <- 1

kgi <- c(-5,-8,-10,-12,-16,-20) # different intial values for kg
nk <- length(kgi)
mkgi <- matrix(0,nk,7)   # matrix for storing results from different inital values for kg

for (j in 1:nk)
{
  kg <- kgi[j]
  z <- getZ(kg)
  bg <- matrix(0,4,1)

  ng <- 10000
  ac <- 0

  kgibbs <- matrix(0,ng,2)
  bgibbs <- matrix(0,ng,k)

  for(g in 1:ng)
  {
    bg  <- bUpdateNorm(x,y - z,sg,priorB,priorIVB) #--proposal for b0, b1, b2, b3
    sg  <- sigmaUpdate(x%*%bg,y - z,s1,s2)         #--proposal for variance
    tmp <- updateK()
    kg  <- tmp$kg
    ac  <- ac + tmp$accept

    bgibbs[g,] <- bg
    kgibbs[g,] <- c(kg,sg)
   }
  x11()
  processPars(kgibbs[,1],CPLOT=T)

  mkgi[j,1] <- kgi[j]
  mkgi[j,2] <- ac/ng
  mkgi[j,3] <-processPars(kgibbs[,1],burnin=3000)$summary[1,1]
  mkgi[j,4:7] <-processPars(bgibbs,burnin=3000)$summary[1:4,1]
}

true <- c(0,0,km,t(beta))
mkgi <- rbind(mkgi,true)
colnames(mkgi) <- c('initial_k','accept_rate','estiamted_k','b0','b1','b2','b3')

write.table(mkgi,file="1-simpar-dif k",row.names=F)
# read a table: read.table("1-simpar-dif k.txt",header=T)


#---converstion to parameters for un-logged model
