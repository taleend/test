x11()
yob <-t1$V3
layout(matrix(1:9,3,3))
for (i in 1:9)
{
K <- tmk[i,1]
V <- tmk[i,2]
b1 <- tmk[i,3]
b2 <- tmk[i,4]

R0 <- w*V/(w+K)
Q10 <- b1*exp(b2*(w+10))
yp <- R0*Q10^(0.1*t-2.4)
plot(yob,yp,ylim=c(0,25))
lines(c(0,50),c(0,50),col=2)
}

layout(matrix(1:12, 3, 4))
for (j in 1:(nk+1))
{
 K <- mkgi[j,2]
 b0 <- mkgi[j,4]
 b1 <- mkgi[j,5]
 b2 <- mkgi[j,6]
 b3 <- mkgi[j,7]
 yp <- getZ(K)+b0+b1*t+b2*w+b3*(t-mean(t))*(w-mean(w))

 plot(y,yp)
 lines(c(-15,25),c(-15,25),col=2)
}

x11()
layout(matrix(1:8, 2, 4))
for (j in 4:(nk+1))
{
 K <- mkgi[j,2]
 b0 <- mkgi[j,4]
 b1 <- mkgi[j,5]
 b2 <- mkgi[j,6]
 b3 <- mkgi[j,7]
 yp <- getZ(K)+b0+b1*t+b2*w+b3*(t-mean(t))*(w-mean(w))

 plot(y,yp)
 lines(c(-15,25),c(-15,25),col=2)
}


