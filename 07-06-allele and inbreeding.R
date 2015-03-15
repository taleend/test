source('clarkFunctions.R')
#---calculate fl, the lower limit of f, pa >=0.5-------------------

pseq <- seq(0.49,1,length=100)
fl <- -(1-pseq)/pseq
plot(pseq,fl,type='l')

#----code for pmake-----
#function(pars)
#{  
#  p   <- pars[1]                   #frequency of a
#  f   <- pars[2]                   #inbreeding coefficient 
#  paa <- p^2 + f*p*(1 - p)        #frequency of p.aa
#  pab <- 2*p*(1 - p)*(1 - f)      #frequency of p.ab
#  pbb <- (1 - p)^2 + f*p*(1 - p)  #frequency of p.bb
#  c(paa,pab,pbb)
#}

#------calcualte lower limit of f, using function minf---------------

pseq <- seq(0.01,0.99,length=100)
plot(pseq,minf(pseq),type='l',xlab='pa',ylab='lower limit of f',col=2)


#--------generate data--------------------------
m <- 50
p <- 0.6
f <- -.2
pr <- pmake(c(p,f))
npop <- 20                  #no. populations
y <- rmultinom(npop,m,pr)   #n random vectors of size m

#----propose intial values for various parameters----
g1 <- 1 #beta parameters for p
g2 <- 1
priorF <- 0 #mean and sd for f
priorFSD <- 1

fg <- f #initial value of f

#----Propose p and f values based on current values pg and fg--

propp <- tnorm(1,.02,.98,p,.002) #propose p
fl <- minf(propp) #lower f limit for pa
propf <- tnorm(1,fl,1,fg,.01) #propose f > fl

#----one chain with a nested inner loop for a single chain

ng <- 5000
nchain <- 10
pgibbs <- matrix(0,ng,nchain)
fgibbs <- pgibbs
accept <- 0

for(j in 1:nchain)
{
pg <- rbeta(1,1,1) #draw initial values from prior
lg <- minf(pg)     # lg for what?
fg <- tnorm(1,lg,1,priorF,priorFSD)

for(g in 1:ng)
{
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
lines(pseq,minf(pseq))

processPars(as.vector(fgibbs[-c(1:1000),]),xtrue=f,DPLOT=T)


