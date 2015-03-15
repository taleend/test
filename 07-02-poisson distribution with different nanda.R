xseq <- 0:30
plot(xseq,exp(-1)/factorial(xseq),type='l')
for (i in 2:10)
lines(xseq,exp(-i)*i^xseq/factorial(xseq),col=i,type='l')
#legend(15,0.3,c('lambda=1','lambda=2','lambda=3','lambda=4','lambda=5'),col=c(1,2,3,4,5),pch='-')
