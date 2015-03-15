m <- 50 #individuals in sample
p <- 0.6 #frequency of allele a
f <- -.2 #inbreeding coefficient
pr <- pmake(c(p ,f)) #values for (Paa,Pab,Pbb)
y <- rmultinom(1,m,pr) #a random vectors of size m
minf(p)
npop <- 20 #no. populations
y <- rmultinom(npop,m,pr) #n random vectors of size m
g1 <- 1 #beta parameters for p
g2 <- 1
priorF <- 0 #mean and sd for f
priorFSD <- 1
fg <- f #initial value of f
propp <- tnorm(1,.02,.98,p,.002) #propose p
fl <- minf(propp) #lower f limit for pa
propf <- tnorm(1,fl,1,fg,.01) #propose f > fl
ng <- 5000
nchain <- 10
pgibbs <- matrix(0,ng,nchain)
fgibbs <- pgibbs
accept <- 0
for(j in 1:nchain){
pg <- rbeta(1,1,1) #draw initial values from prior
lg <- minf(pg)
fg <- tnorm(1,lg,1,priorF,priorFSD)
for(g in 1:ng){
pf <- update_pf()
pg <- pf$pg
fg <- pf$fg
pgibbs[g,j] <- pg
fgibbs[g,j] <- fg
}
}
par(mfrow=c(2,1))
plot(fgibbs[,1],type='l',ylim=c(-1,1))
for(j in 2:nchain)lines(fgibbs[,j])
abline(h=f)
plot(pgibbs[,1],type='l',ylim=c(0,1))
for(j in 2:nchain)lines(pgibbs[,j])
abline(h=p)
par(mfrow=c(1,1))
plot(pgibbs,fgibbs,xlim=c(0,1),ylim=c(-1,1))
lines(pseq,minf(pseq)
processPars(as.vector(fgibbs[-c(1:1000),]),xtrue=f,DPLOT=T)
