source('clarkFunctions.R')

T1 <- read.table("arc1.txt",header=F)
T1 <- na.omit(T1)
t1 <- T1[which(T1$V5<0),] #--exclude positive values for column 5(water table data)
write.table(t1,file="arc1ready.txt")
yob <- t1$V3       #--observed data of soil respiration
y <- log(yob)      #--logged soil respiration data for estimation
x <- cbind(rep(1,length(y)),t1$V4,t1$V5)
t <- x[,2]
w <- x[,3]
tw <- (x[,2] - mean(x[,2]))*(x[,3] - mean(x[,3]))
x <- cbind(x,tw)

#----functions------------------------------------------------------
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

#kgi <- seq(-20,-4,1)
kgi <- (-20):(-3)      #--different intial values for kg
#kgi <- c(-15,-10,-5)
nk <- length(kgi)
mkgi <- matrix(0,nk,7) #--matrix for storing averaged results from different inital values for kg
colnames(mkgi) <- c('initial_k','estiamted_k','accept_rate','b0','b1','b2','b3')
accept <- rep(0,nk)

ng <- 5000

kgibbs <- rep(0,2*ng*nk)  
dim(kgibbs) <- c(ng,2,nk)  #--store all simulation for k
bgibbs <- rep(0,ng*k*nk)
dim(bgibbs) <- c(ng,k,nk)  #--store all simulation for b's

for (j in 1:nk)
{
kg <- kgi[j]
z <- getZ(kg)
bg <- matrix(0,4,1)

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

accept[j] <- ac
mkgi[j,1] <- kgi[j]
mkgi[j,2] <-processPars(kgibbs[,1,j],burnin=2000)$summary[1,1]
mkgi[j,3] <- ac/ng
mkgi[j,4:7] <-processPars(bgibbs[,,j],burnin=2000)$summary[1:4,1]
}

#true <- c(0,km,0,t(beta))
#mkgi <- rbind(mkgi,true)

write.table(mkgi,file="2-4bs.txt")


