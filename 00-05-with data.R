source('clarkFunctions.R')

T1 <- read.table("arc1.txt",header=F)
T1 <- na.omit(T1)
t1 <- T1[which(T1$V5<0),] #--exclude positive values for column 5(water table data)
yob <- t1$V3       #--observed data of soil respiration
y <- log(yob)      #--logged soil respiration data for estimation
x <- cbind(rep(1,length(y)),t1$V4,t1$V5)
t <- x[,2]
w <- x[,3]
x <- cbind(x,(x[,2] - mean(x[,2]))*(x[,3] - mean(x[,3])))   #--design matrix

getZ <- function(km){log(-w) - log(-km -w)} #--the minus mark is used because water table data is negative

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

k <- 4
priorB   <- matrix(0,k,1)          #prior means
priorVB  <- diag(1000,k)           #prior covariance
priorIVB <- solve(priorVB)         #prior inv covariance
priorK  <- .1
priorKsd <- 1
sg <- 0.01

s1      <- 1                      #prior values for variance
s2      <- 1

kgi <- (-20):(-1)      #--different intial values for kg
nk <- length(kgi)
mkgi <- matrix(0,nk,7) #--matrix for storing averaged results from different inital values for kg
accept <- rep(0,nk)

kgibbs <- rep(0,2*ng*nk)  
dim(kgibbs) <- c(ng,2,nk)  #--store all simulation for k
bgibbs <- rep(0,ng*k*nk)
dim(bgibbs) <- c(ng,k,nk)  #--store all simulation for b's

for (j in 1:nk)
{
kg <- kgi[j]
z <- getZ(kg)
bg <- matrix(0,4,1)

ng <- 10000
ac <- 0

for(g in 1:ng)
{
  bg  <- bUpdateNorm(x,y - z,sg,priorB,priorIVB)
  sg  <- sigmaUpdate(x%*%bg,y - z,s1,s2)
  tmp <- updateK()
  kg  <- tmp$kg
  ac  <- ac + tmp$accept

  bgibbs[g,,j] <- bg
  kgibbs[g,,j] <- c(kg,sg)
}

#accept[j] <- ac
mkgi[j,1] <- kgi[j]
mkgi[j,2] <-processPars(kgibbs[,1,j],burnin=3000)$summary[1,1]
mkgi[j,3] <- ac/ng
mkgi[j,4:7] <-processPars(bgibbs[,,j],burnin=3000)$summary[1:4,1]
}

true <- c(0,km,0,t(beta))
mkgi <- rbind(mkgi,true)
colnames(mkgi) <- c('initial_k','estiamted_k','accept_rate','b0','b1','b2','b3')

tmk <- matrix(0,(nk+1),4)
colnames(tmk) <-c('k','v','be1','be2')
tmk[,1] <- mkgi[,2]
tmk[,4] <- mkgi[,6]/-2.4
tmk[,3] <- exp(10*mkgi[,5]-10*tmk[,4])
tmk[,2] <- exp(mkgi[,4]+ 2.4*log(tmk[,3])+24*tmk[,4])

layout(matrix(1:20,4,5))
for (i in 2:(nk+1))
{
K <- tmk[i,1]
V <- tmk[i,2]
b1 <- tmk[i,3]
b2 <- tmk[i,4]

R0 <- w*V/(w+K)
Q10 <- b1*exp(b2*(w+10))
yp <- R0*Q10^(0.1*t-2.4)
plot(yob,yp,ylim=c(0,25))
lines(c(0,50),c(0,50),col=2)
}

tmks <- rbind(tmk[7:14,],tmk[21,])
x11()
layout(matrix(1:9,3,))
for (i in 1:9)
{
K <- tmks[i,1]
V <- tmks[i,2]
b1 <- tmks[i,3]
b2 <- tmks[i,4]

R0 <- w*V/(w+K)
Q10 <- b1*exp(b2*10)
yp <- R0*Q10^(0.1*t-2.4)
plot(yob,yp,ylim=c(0,25))
lines(c(0,50),c(0,50),col=2)
}

x11()
layout(matrix(1:9,3,))
for (i in 1:9)
{
K <- tmks[i,1]
V <- tmks[i,2]
b1 <- tmks[i,3]
b2 <- tmks[i,4]

R0 <- w*V/(w+K)
Q10 <- b1*exp(b2*10)
yp <- R0*Q10^(0.1*t-2.4)
plot(yob,yp-yob)
}

x11()
layout(matrix(1:9,3,3))
dif <- matrix(0,9,length(y))
for (i in 1:9)
{
K <- tmks[i,1]
V <- tmks[i,2]
b1 <- tmks[i,3]
b2 <- tmks[i,4]

R0 <- w*V/(w+K)
Q10 <- b1*exp(b2*(w+10))
Q10N <- b1*exp(b2*10)
yp <- R0*Q10^(0.1*t-2.4)
ypN <- R0*Q10N^(0.1*t-2.4)
dif[i,] <- yp-ypN
hist(dif[i,])
#plot(yob,yp-yob,col='forestgreen')
}