source('clarkFunctions.r')
b <- read.csv('dataFecundity.csv',header=T)

yr      <- b[,'year']       #duration
beginyr <- min(yr)
endyr   <- max(yr)
yrvec   <- c(beginyr:endyr)
tall    <- sort(unique(yr))
ti      <- match(yr,tall)

plot  <- b[,'plot']         #plots within groups
group <- b[,'Group']
plotnames  <- sort(unique(plot))
groupnames <- sort(unique(group))
ngroup     <- length(groupnames)
tab  <- table(plot,group)
wi   <- which(tab > 0,arr.ind=T)

gi <- match(group,groupnames)
nt <- length(tall)
n  <- nrow(b)

gmeans <- rep(0,ngroup)   #group and yr effects
tmeans <- rep(0,nt)

#The above code organizes the data for efficient Gibbs sampling.  
#Here are cover values used to determine fecundity:
C <- b[,'cover.t0']
cmat <- matrix(0,n,ngroup)
cmat[cbind(c(1:n),group)] <- C
cMean <- apply(cmat,2,mean)
cMean <- cMean[group]
pc    <- .9
cg    <- pc*C + (1 - pc)*cMean

#Here are predictors and response:

xnames <- c('TmeanSpr1','pptSum1')
p      <- length(xnames) + 1
x      <- matrix(1,n,p)
x[,c(2:p)] <- matrix(unlist(b[,xnames]),n,(p-1))
colnames(x) <- c('intercept',xnames)
y      <- b[,'nrec.t1']
lambda <- y + .1            #initial value of lambda

pairs(cbind(y,x[,-1],C))

#Here are prior parameter values and initial values:

priorb   <- matrix(0,p,1)
priorVb  <- diag(20,p)
priorIVB <- solve(priorVb)
loB      <- priorb*0 - 100
hiB      <- priorb*0 + 100

vg0 <- 1   #variance group effects
sg1 <- ngroup*4
sg2 <- vg0*(sg1 - 1)

vt0 <- 1   
st1 <- nt*4
st2 <- vt0*(st1 - 1)

vg <- vg0             #initial variances
vt <- vt0

priorSize   <- 1
priorSdSize <- .1

sg <- priorSize

XX  <- crossprod(x)
IXX <- solve(XX)

bg <- IXX%*%crossprod(x,log(lambda) - log(cg))  # intial values
teffect <- tmeans[ti]
geffect <- gmeans[gi]

#In the Gibbs sampler I first update the group effects using this function:

updateGroup <- function(){  #metropolis

  z  <- log(lambda) - log(cg) - x%*%bg - teffect #group means
  zm <- matrix(0,n,ngroup)
  zm[cbind(c(1:n),gi)] <- z

  z1 <- zm
  z1[zm != 0] <- 1
  
  zmu <- apply(zm,2,sum)
  nz  <- apply(z1,2,sum) 

  gp <- rnorm(ngroup,zmu/nz,rexp(ngroup,10)) #proposal
  gp <- gp - mean(gp)                        #center

  m1 <- cg*exp(x%*%bg + teffect + geffect)
  m2 <- cg*exp(x%*%bg + teffect + gp[gi])

  pnow <- sum(dnbinom(y,size=sg,mu=m1,log=T)) + 
          sum(dnorm(gmeans,0,sqrt(vg),log=T))
  pnew <- sum(dnbinom(y,size=sg,mu=m2,log=T)) + 
          sum(dnorm(gp,0,sqrt(vg),log=T))
  
  acceptMH(pnow,pnew,gmeans,gp)$x
}

#Here is a sample for year effects:
updateTime <- function(){

  z <- log(lambda) - log(cg) - x%*%bg - geffect
  zm <- matrix(0,n,nt)
  zm[cbind(c(1:n),ti)] <- z
  z1 <- zm
  z1[zm != 0] <- 1
  
  zmu <- apply(zm,2,sum)
  nz  <- apply(z1,2,sum)

  tp <- rnorm(nt,zmu/nz,rexp(nt,100))
  tp <- tp - mean(tp)

  m1 <- cg*exp(x%*%bg + teffect + geffect)
  m2 <- cg*exp(x%*%bg + tp[ti]  + geffect)

  pnow <- sum(dnbinom(y,size=sg,mu=m1,log=T)) + 
          sum(dnorm(tmeans,0,sqrt(vt),log=T))
  pnew <- sum(dnbinom(y,size=sg,mu=m2,log=T)) + 
          sum(dnorm(tp,0,sqrt(vt),log=T))
  
  acceptMH(pnow,pnew,tmeans,tp)$x
}

#Here is a sample for predictor coefficients:

updateB <- function(){  #metropolis to update fixed effects

  vall <- vt + vg
  V    <- IXX*vall

  bp   <- t(myrmvnorm(1,bg,V))
  sp   <- tnorm(1,0,20,sg,.01)

  ci <- pc*C + (1 - pc)*cMean

  m1 <- ci*exp(x%*%bg + teffect + geffect)
  m2 <- ci*exp(x%*%bp + teffect + geffect)

  pnow <- sum(dnbinom(y,size=sg,mu=m1,log=T)) + 
              mydmvnorm(t(bg),priorb,priorVb,log=T) +
              dnorm(sg,priorSize,priorSdSize,log=T)
  pnew <- sum(dnbinom(y,size=sp,mu=m2,log=T)) + 
              mydmvnorm(t(bp),priorb,priorVb,log=T) +
              dnorm(sp,priorSize,priorSdSize,log=T)

  tmp <- acceptMH(pnow,pnew,c(bg,sg),c(bp,sp),BLOCK=T)$x
  bg[1:p,1] <- tmp[1:p]
  sg        <- tmp[p+1]

  list(bg = bg, sg = sg)
}

#Here is an update for the fraction of seed originating within the plot:

updatePC <- function(){

  pd <- tnorm(1,.9,1,pc,.01)

  ci <- pc*C + (1 - pc)*cMean
  cd <- pd*C + (1 - pd)*cMean

  m1 <- ci*exp(x%*%bg + teffect + geffect)
  m2 <- cd*exp(x%*%bg + teffect + geffect)

  pnow <- sum(dnbinom(y,size=sg,mu=m1,log=T)) 
  pnew <- sum(dnbinom(y,size=sg,mu=m2,log=T)) 

  acceptMH(pnow,pnew,pc,pd)$x
}

#Here is an update for the variances:

updateS <- function(s1,s2,s){

  u1 <- s1 + length(s)/2
  u2 <- s2 + .5*sum(s^2)
  1/rgamma(1,u1,u2)
}

#Here is a Gibbs sampler:
timeEffect <- F             #include random yr effects?

ng <- 5000

tgibbs <- matrix(0,ng,nt)
ggibbs <- matrix(0,ng,ngroup)
bgibbs <- matrix(0,ng,p)
vgibbs <- matrix(0,ng,3)
yhat   <- matrix(0,ng,n)


colnames(bgibbs) <- colnames(x)
colnames(vgibbs) <- c('groupVar','timeVar','nbPar')

for(g in 1:ng){

  gmeans <- updateGroup()
  geffect <- gmeans[gi]

  if(timeEffect){
    tmeans <- updateTime()
    teffect <- tmeans[ti]
  }

  tmp <- updateB()
  bg  <- tmp$bg
  sg  <- tmp$sg

  pc <- updatePC()

  cg <- C*pc + cMean*(1 - pc)

  lambda   <- cg*exp(x%*%bg + teffect + geffect)
  yhat[g,] <- rnbinom(n,size=sg,mu=lambda)

  vg <- updateS(sg1,sg2,gmeans)
  vt <- updateS(st1,st2,tmeans)

  if(g %in% c(200,500,1000)){
    bcov <- cov(bgibbs[1:(g-1),])
  }

  tgibbs[g,] <- tmeans
  ggibbs[g,] <- gmeans
  bgibbs[g,] <- bg
  vgibbs[g,] <- c(vg,vt,sg)
  print(bg)
}

processPars(bgibbs,c(1:p)*0,CPLOT=T)

processPars(vgibbs,rep(1,3),CPLOT=T)

processPars(ggibbs,rep(0,ngroup),CPLOT=T)

processPars(tgibbs[,1:9],rep(0,9),CPLOT=T)

predVsObs(y,yhat)

