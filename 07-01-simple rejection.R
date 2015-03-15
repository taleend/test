#----simple example for generating simulation for f(x)=(x-0.5)^2/2.4,0<x<5----
#----the maximal f(x) is (5,8.4375), we choose g(x)=0.5

c <-1.875
n <- 20000
x0 <- runif(n,0,2)
f <- (x0-0.5)^2/2.4
g <- 0.5
cc <- f/c/g

f1 <-rep(NA,n)

for (i in 1:n)
{
yy <- runif(1,0,1)
if (yy<cc[i])
f1[i]=x0[i]
}

f2 <-f1[!is.na(f1)]

#-----------------
ff1 <- rep(NA,n)

for (i in 1:n)
{
if (cc[i]<1)
ff1[i]=x0[i]
}

ff2 <-ff1[!is.na(ff1)]