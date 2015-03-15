dIG <- function(x,a,b)
{
p <- b^a/factorial(a-1)*x^(-a-1)*exp(-b/x)
}

par(mfrow=c(2,1))
tseq <- seq(0,5,0.05)
plot(tseq,dinvGamma(tseq,2,1),type='l',col=2)
plot(tseq,dIG(tseq,2,1),type='l',col=2,lty=2)