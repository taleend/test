n <- 100 #no. observations
k <- 4 #no. covariates
b <- matrix(c(-.1,2,.5,.3),ncol=1) #parameters
sigma <- .1 #variance parameter
loX <- c(1,0,-10,-3) #range of covariates, length = k
hiX <- c(1,1,10,0)
LIKE <- 'norm' #'norm', 'pois', 'binom', 'multinom', 'mvnorm',
'mvnorm-Pois'
x <- simX(n,loX,hiX) #simulate x and y
y <- simY(x,b,LIKE,r,sigma)
r <- 1
priorB <- b*0 #prior means
priorVB <- diag(1000,k*r) #prior covariance
priorIV <- solve(priorVB) #prior inv covariance
s1 <- 1 #prior values for variance
s2 <- 1
ng <- 1000
tmp <- gibbsLoop(LIKE,ng,x,y,b=priorB,sigma=1)
bg <- tmp$bgibbs
sg <- tmp$sgibbs
processPars(bg,xtrue=b,CPLOT=T)
processPars(sg,xtrue=sigma,CPLOT=T)
plot(y,tmp$ymean)
abline(0,1)
ylo <- tmp$ymean - 1.96*tmp$yse
yhi <- tmp$ymean + 1.96*tmp$yse
for(i in 1:n){
lines(c(y[i],y[i]),c(ylo[i],yhi[i]))
}
