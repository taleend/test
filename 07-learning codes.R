qq <- runif(5,-pi,pi)
aa <- (b0 + b1*cos(b2+qq))/cf
zz <- runif(5,0,1)

zx <- qq[zz<aa]

plot(qq)
points(aa,col=2)
points(zz,col=3,pch=5)
points(zx,col=4,pch=4)

gr <- 1:7
abline(v=gr,col=7)

points(zz-aa,col=6,pch='-')

re <- 0
x1 <- sort(x)
x2 <- sort(xz)
for (i in 1:length(x2))
{
for (j in 1:length(x1))
	{
	if (x2[i]==x1[j])
	re=re+1
	}
}

#an example from r-intro-cn about functions

ricker <- function(nzero, r, K=1, time=100, from=0, to=time) 
{
N <- numeric(time+1)
N[1] <- nzero
for (i in 1:time) N[i+1] <- N[i]*exp(r*(1 - N[i]/K))
Time <- 0:time
plot(Time, N, type="l", xlim=c(from, to))
N
}

layout(matrix(1:3,3,1))
ricker(0.1,1);title("r=1")
ricker(0.1,2);title("r=2")
ricker(0.1,3);title("r=3")

#-----Markov chain Monte Carlo simulation---------------

