source('clarkFunctions.R')
T1 <- read.table("arc1.txt",header=F)
T1 <- na.omit(T1)
#----exclude positive values for col5(water table data)---
t1 <- T1[which(T1$V5<0),]
y <- t1$V3
x <- cbind(rep(1,length(y)),t1$V4,t1$V5)
t <- x[,2]
w <- x[,3]
x <- cbind(x,(x[,2] - mean(x[,2]))*(x[,3] - mean(x[,3])))

getZ <- function(km){log(-w) - log(-km -w)}

updateK <- function()
{
  ks <- -tnorm(1,0,50,-kg,rexp(1,10))
  zn <- getZ(kg)
  zs <- getZ(ks)

  pnow <- sum(dnorm(y,zn + x%*%bg,sqrt(sg),log=T)) + 
          dnorm(kg,priorK,priorKsd,log=T)
  pnew <- sum(dnorm(y,zs + x%*%bg,sqrt(sg),log=T)) + 
          dnorm(ks,priorK,priorKsd,log=T)

  tmp <- acceptMH(pnow,pnew,kg,ks)
  list(kg = tmp$x, accept = tmp$accept)
}

k <- 4
priorB   <- matrix(0,k,1)          #prior means
priorVB  <- diag(1000,k)           #prior covariance
priorIVB <- solve(priorVB)         #prior inv covariance
priorK  <- .1
priorKsd <- 1
sg <- 0.01

s1      <- 1                      #prior values for variance
s2      <- 1

kg <- -1
z <- getZ(kg)
bg <- matrix(0,4,1)

ng <- 5000
ac <- 0

kgibbs <- matrix(0,ng,2)
bgibbs <- matrix(0,ng,k)


for(g in 1:ng)
{
  bg  <- bUpdateNorm(x,y - z,sg,priorB,priorIVB)
  sg  <- sigmaUpdate(x%*%bg,y - z,s1,s2)
  tmp <- updateK()
  kg  <- tmp$kg
  ac  <- ac + tmp$accept

  bgibbs[g,] <- bg
  kgibbs[g,] <- c(kg,sg)
}

processPars(kgibbs)
processPars(bgibbs)


K <- mean(kgibbs[,1])
b0 <- mean(bgibbs[,1])
b1 <- mean(bgibbs[,2])
b2 <- mean(bgibbs[,3])
b3 <- mean(bgibbs[,4])

lny <- log(y)
lnyp <- log(w/(w+K))+b0+b1*t+b2*w+b3*t*w
plot(lny,lnyp)
lines(c(-5,25),c(-5,25),col=2)

yp <- exp(lnyp)
plot(y,yp,ylim=c(0,50000))
lines(c(-5,30),c(-5,30),col=2)


