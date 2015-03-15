source('clarkFunctions.R')

getZ <- function(km){log(-w) - log(-km - w)}

n    <- 1000
beta <- matrix(c(-1.967,.189,.006),3,1) #--values of parameters from MLE
km   <- -6.1                                 #--value of km from MLE
sd   <- .1
loX <- c(1,3,-40)    #int, temp, h2o, interaction
hiX <- c(1,25,-0.001)      

x <- simX(n,loX,hiX)     #simulate x and y
t <- x[,2]
w <- x[,3]
tw <- (x[,2] - mean(x[,2]))*(x[,3] - mean(x[,3]))
x <- cbind(x[,1],t,tw-24*w)
z <- log(-w) - log(-km - w)
y <- rnorm(n,z + x%*%beta,sd)


parset <- read.table('00-simulationPars.txt')
colnames(parset) <- c('k','b0','b1','b2')
np <- length(parset[,1])

L <- rep(0,np)

for (i in 1:np)
{
k <- parset[i,1]
b0 <- parset[i,2]
b1 <- parset[i,3]
b2 <- parset[i,4]
yp <- getZ(k)+b0+b1*t+b2*(tw-24*w)
L[i] <- sum(dnorm(y-yp,mean=0,sd=0.1))
}

D <- rep(0,np-1)
ML <- L[1]
for (i in 1:np-1)
D[i] <- -2*(L[i+1]-ML)
cbind(parset[2:19,1],D)
