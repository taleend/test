n <- 10
p <- 3
b <- matrix(runif(p,-1,1),p,1)
x <- matrix(runif(n*p,-1,1),n,p)
x[,1] <-1
y1 <- x%*%b

sigma <- .1

y <- rnorm(n,y1,sigma)

bhat <- solve(crossprod(x))%*%crossprod(x,y)

shat <- crossprod(y-x%*%bhat)/(n-p)

coeffs <- cbind(c(b,sigma),c(bhat,shat))
colnames(coeffs) <- c('true','estimated')
rownames(coeffs) <- c(paste('b',0:(p-1),sep='-'),'sigma')