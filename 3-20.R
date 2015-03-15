filename <- ('dataFACE.txt')
ynames <- 'cones'
xnames <- c('diam','trt')
tmp <- inData(filename, xnames, ynames, na.rm = T, INTERCEPT = T)
x <- tmp$x
y <- tmp$y
n <- length(y)
k <- ncol(x)
r <- 1
b <- matrix(0,k,1)
LIKE <- 'pois'
priorB <- b*0 #prior means
priorVB <- diag(1000,k*r) #prior covariance
priorIV <- solve(priorVB) #prior inv covariance
ng <- 10000
tmp <- gibbsLoop(LIKE,ng,x,y,b,0)
processPars(tmp$bg,CPLOT=T)
plotObsPred(y,tmp$ymean,tmp$yse)
