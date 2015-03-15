load('dataPathogen.RData')
died <- apply(data[,c(1,3,5)],1,sum)     # no. died, survived
survived <- apply(data[,c(2,4,6)],1,sum)
nvec <- data[1,]
nd <- sum(nvec[1:4])                     #those with detection
ns <- sum(nvec[5:6])                     #without detection
v <- .02                                 #variance for proposals
lo <- c(.01,.9,.1,.1)                    #low and high truncation values
hi <- c(.99,1,1,.9)
priorm <- c(.5,.95,.7,.6)                #prior mean
priorV <- diag(c(1,.2,.5,.5))            #prior covariance
rnames <- c('theta','phi','s0','s1')
ng <- 10000
pvec <- tnorm(4,lo,hi,priorm,sqrt(diag(priorV))) #initial from prior
pgibbs <- matrix(0,ng,4)
colnames(pgibbs) <- rnames
a <- 0
for(g in 1:ng){
prop <- tnorm(4,lo,hi,pvec,v)             #proposal
pnew <- sum(nvec[1:4]*log(qprob(prop))) +
dmvnorm(prop,priorm,priorV,log=T)
pnow <- sum(nvec[1:4]*log(qprob(pvec))) +
dmvnorm(pvec,priorm,priorV,log=T)
atmp <- acceptMH(pnow,pnew,pvec,prop,BLOCK=T)
pvec <- atmp$x
a <- a + atmp$accept
pgibbs[g,] <- pvec
}
