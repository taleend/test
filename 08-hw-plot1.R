x1 <- matrix(0,512,3)
y1 <- matrix(0,512,3)

for (i in 1:3)
{
 x1[,i] <- density(ff[,i])$x
 y1[,i] <- density(ff[,i])$y
}

x2 <- matrix(0,512,3)
y2 <- matrix(0,512,3)

for (i in 1:3)
{
 x2[,i] <- density(ft[,i])$x
 y2[,i] <- density(ft[,i])$y
}


x3 <- matrix(0,512,3)
y3 <- matrix(0,512,3)

for (i in 1:3)
{
 x3[,i] <- density(tf[,i])$x
 y3[,i] <- density(tf[,i])$y
}

x4 <- matrix(0,512,3)
y4 <- matrix(0,512,3)

for (i in 1:3)
{
 x4[,i] <- density(tt[,i])$x
 y4[,i] <- density(tt[,i])$y
}

x11()
layout(matrix(1:3,3,1))
plot(x1[,1],y1[,1],type='l',xlim=c(-5,6),ylim=c(0,0.5),ylab='density')
title("INTERCEPT")
legend(3.5,.478,c('Group Off Year Off','Group Off Year On',
		'Group On Year Off','Group On Year On'),col=c(1,2,3,4),pch='-')

lines(x2[,1],y2[,1],col=2)
lines(x3[,1],y3[,1],col=3)
lines(x4[,1],y4[,1],col=4)

plot(x1[,2],y1[,2],type='l',xlim=c(-.5,1),ylim=c(0,5),ylab='density')
title("TMEANSPR1")
lines(x2[,2],y2[,2],col=2)
lines(x3[,2],y3[,2],col=3)
lines(x4[,2],y4[,2],col=4)

plot(x1[,3],y1[,3],type='l',xlim=c(-0.01,0.04),ylim=c(0,150),ylab='density')
title("PPTSUM1")
lines(x2[,3],y2[,3],col=2)
lines(x3[,3],y3[,3],col=3)
lines(x4[,3],y4[,3],col=4)

