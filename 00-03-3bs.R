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


#----Here is a simulated data set and Gibbs sampler:
n    <- 1000
beta <- matrix(c(-1.967,.189,.006),3,1) #--values of parameters from MLE
km   <- -6.1                                 #--value of km from MLE
sd   <- .1
sg   <- sd^2

loX <- c(1,3,-40)    #int, temp, h2o, interaction
hiX <- c(1,25,-0.001)      

x <- simX(n,loX,hiX)     #simulate x and y
t <- x[,2]
w <- x[,3]
tw <- (x[,2] - mean(x[,2]))*(x[,3] - mean(x[,3]))
x <- cbind(x[,1],t,tw-24*w)
z <- log(-w) - log(-km - w)
y <- rnorm(n,z + x%*%beta,sd)

#----estimation------------------
k <- length(beta)
priorB   <- beta*0                    #prior means
priorVB  <- diag(1000,k)           #prior covariance
priorIVB <- solve(priorVB)         #prior inv covariance
priorK  <- .1
priorKsd <- 1

s1      <- 1                      #prior values for variance
s2      <- 1

kgi <- -20 # different intial values for kg
nk <- length(kgi)
mk <- matrix(0,nk,6)   # matrix for storing averaged results from different inital values for kg

for (j in 1:nk)
{
kg <- kgi[j]
z <- getZ(kg)
bg <- matrix(0,4,1)

ng <- 300000
ac <- 0

kgibbs <- matrix(0,ng,2)
bgibbs <- matrix(0,ng,k)

for(g in 1:ng)
{
  bg  <- bUpdateNorm(x,y - z,sg,priorB,priorIVB) #--proposal for b0, b1, b2
  sg  <- sigmaUpdate(x%*%bg,y - z,s1,s2)       #--proposal for variance
  tmp <- updateK()
  kg  <- tmp$kg
  ac  <- ac + tmp$accept

  bgibbs[g,] <- bg
  kgibbs[g,] <- c(kg,sg)
}
x11()
processPars(kgibbs[,1],CPLOT=T)
title(paste('k =',kgi[j],sep=' '))

x11()
processPars(bgibbs,CPLOT=T)
title(paste('k =',kgi[j],sep=' '))

mk[j,1] <- kgi[j]
mk[j,2] <- ac/ng
mk[j,3] <-processPars(kgibbs[,1],burnin=3000)$summary[1,1]
mk[j,4:6] <-processPars(bgibbs,burnin=3000)$summary[1:3,1]
}

true <- c(0,0,km,t(beta))
mk <- rbind(mk,true)
colnames(mk) <- c('initial_k','accept','estiamted_k','b0','b1','b2')

tmk <- matrix(0,(nk+1),4)#matrix for storing values parameter values for the original model
colnames(tmk) <-c('k','v','be1','be2')
tmk[,1] <- mk[,3]
tmk[,4] <- mk[,6]*10
tmk[,3] <- exp(10*mk[,5]-10*tmk[,4])
tmk[,2] <- exp(mk[,4]+ 2.4*log(tmk[,3])+24*tmk[,4])



bset <- mk[,4:6]
beta <- matrix(c(-1.967,.189,.006),3,1)
ML <- -sum(getZ(km)+x%*%beta)
L <- rep(0,nk)
D <- L
for (i in 1:nk)
{
beta <-matrix(bset[i,],3,1)
L[i] <- -sum(getZ(mk[i,3])+x%*%beta)
D[i] <- -2*(L[i]-ML)
}












# plot prediction vs. observation
x11()
yob <- exp(y)
layout(matrix(1:6,2,3))
for (i in 1:(nk+1))
{
K <- tmk[i,1]
V <- tmk[i,2]
be1 <- tmk[i,3]
be2 <- tmk[i,4]

R0 <- w*V/(w+K)
Q10 <- be1*exp(be2*(w+10))
yp <- R0*Q10^(0.1*t-2.4)
plot(yob,yp)
lines(c(0,50),c(0,50),col=2)
}

