n <- 100
b0 <- 10
b1 <- 2
beta <- matrix(c(b0,b1),2,1)
sd <- 4

x1 <- runif(n,0,20)
x <- cbind(rep(1,n),x1)
y <- rnorm(n,x%*%beta,sd)

plot(x[,2],y)
abline(b0,b1)

priorB <- as.vector(c(0,0))
priorIV <- solve(diag(1000,2))

s1 <- 0.1
s2 <- 0.1

b.update <- function()
{
  V <- solve(crossprod(x)/sg + priorIV)
  v <- crossprod(x,y)/sg + priorIV%*%priorB
  t(rmvnorm(1,V%*%v,V))
}

sg <-10^2

library(mvtnorm)

bg <- b.update()
v.update <- function()
{
  u1 <- s1+n/2
  u2 <- s2 + .5*crossprod(y-x%*%bg)
1/rgamma(1,u1,u2)
}

sg <-1
bg <-b.update()

sg <- v.update()

ng <-1000
bgibbs <- matrix(NA,ng,2)
colnames(bgibbs) <- c('b0','b1')
sgibbs <- rep(0,ng)

for (g in 1:ng)
{
 sg <- v.update()
 bg <- b.update()

 bgibbs[g,] <- bg
 sgibbs[g] <- sg
}

processPars(cbind(bgibbs,sgibbs), xtrue=c(beta,sd^2), CPLOT=T)
processPars(cbind(bgibbs,sgibbs), xtrue=c(beta,sd^2),DPLOT=T,burnin=100)

summary(lm(y~x[,2]))

lm(formula = y~x[,2])