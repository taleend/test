P1 <- t(matrix(c(.35,.35,.1,.1,.1,.15,.55,.1,.1,.1,.15,.15,.1,.2,.4,.15,.15,.4,.1,.2,.15,.15,.2,.4,.1),5,5))
P2 <- t(matrix(c(.49985,.49985,.0001,.0001,.0001,.29985,.69985,.0001,.0001,.0001,.00015,.00015,.1999,.2999,.4999,.00015,.00015,.4999,.1999,.2999,.00015,.00015,.2999,.4999,.1999),5,5))

# calculate P1(n), where n=2^5
P125 <- P1
for (i in 1:5)
P125 <- P125%*%P125

# calculate P1(n), where n=2^10
P11 <- matrix(0,50,5)
P125 <- P1
for (i in 1:10)
{
P125 <- P125%*%P125
P11[(5*i-4):(5*i),] <- P125
}

# calculate P2(n), where n=2^5
P225 <- P2
for (i in 1:5)
P225 <- P225%*%P225

# calculate P2(n), where n=2^50
P21 <- matrix(0,250,5)
P225 <- P2
for (i in 1:50)
{
P225 <- P225%*%P225
P21[(5*i-4):(5*i),] <- P225
}

# calculate P2(n), where n=2^60
# the value of max(P21) goes very high when i reaches 59
P21 <- matrix(0,300,5)
P225 <- P2
for (i in 1:57)
{
P225 <- P225%*%P225
P21[(5*i-4):(5*i),] <- P225
}

plot(P21[1:280,1],type='l',ylim=c(0,max(P21[1:280])),lty=2)
for (i in 2:5)
lines(P21[1:280,i],col=i,lty=2)

sumpr <- matrix(0,300,1)
for (i in 1:300)
sumpr[i] <- sum(P21[i,])


# try Metropolis algorithm Example 7
M1 <- t(matrix(c(.4,.15,.3,.15,.1,.466667,.3,.133333,.15,.225,.575,.05,.3,.4,.2,.1),4,4))
# calculate M1(n), where n=2^50=1.12E15
Ct <- matrix(0,200,4)
M3 <- M1
for (i in 1:50)
{
M3 <- M3%*%M3
Ct[(4*i-3):(4*i),] <- M3
}
