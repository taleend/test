#---get data--------------
T1 <- read.table("arc1.txt",header=F)
T1 <- na.omit(T1)
#----exclude positive values for col5-------------
t1 <- T1[which(T1$V5<0),]
y <- t1$V3
t <- t1$V4
w <- t1$V5

lny <- log(y)


#---guofang's model
V <-13.58
K <- -7.2
b1 <- 3.63
b2 <- 0.0533

R0 <- w*V/(w+K)
Q10 <- b1*exp(b2*(w+10))

yp <- R0*Q10^(0.1*t-2.4)
x11()
plot(y,yp)
lines(c(0,25),c(0,25),col=2)

lnyp <- log(w/(w+K))+log(V)+(0.1*t-2.4)*(log(b1)+b2*w+10*b2)
plot(lny,lnyp)
lines(c(-5,25),c(-5,25),col=2)

#----my 1st----------
V <-157.748
K <- -5.839
b1 <- 15.142
b2 <- 0.076

R0 <- w*V/(w+K)
Q10 <- b1*exp(b2*(w+10))
yp <- R0*Q10^(0.1*t-2.4)
x11()
plot(y,yp,col=2)

lnyp <- log(w/(w+K))+log(V)+(0.1*t-2.4)*(log(b1)+b2*w+10*b2) 

K <- -5.839
bb0 <- -3.291
bb1 <- 0.348
bb2 <- -0.183
bb3 <- -0.006
lnypb <- log(w/(w+K)) + bb0 + bb1*t + bb2*w +bb3*t*w
plot(lny,lnyp)
lines(c(-5,25),c(-5,25),col=2)


#----my second-------
V <-157.748
K <- -5.839
b1 <- 59.145
b2 <- -0.06

R0 <- w*V/(w+K)
Q10 <- b1*exp(b2*(w+10))
yp <- R0*Q10^(0.1*t-2.4)
x11()
plot(y,yp,col=3)