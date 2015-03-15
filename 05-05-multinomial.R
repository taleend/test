source('clarkFunctions.r')
n <- 1000 #sample size
theta <- .6 #Pr{I}
phi <- .8 #Pr{D|I}
s0 <- .9 #Pr{S|I=0}
s1 <- .4 #Pr(S|I=1}
p <- qprob(c(theta,phi,s0,s1)) #Pr{D,S}
nvec <- as.vector(rmultinom(1,n,p)) #simulate data
mle <- nlminb(rep(.5,4),like_multinom, #MLEs
lower=c(.1,.1,.1,.1),upper=c(.8,1,1,1))$par

