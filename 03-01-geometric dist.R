m <- 1000 #no. experiments
first3 <- rep(NA,m) #for 1st 3 in each exp.
for(j in 1:m){ #repeat for m exp's.
count <- 0 #counter
trial <- 0 #no. on a die
while(trial != 3){ #repeat until 3
count <- count + 1
trial <- sample(c(1:6),1) #roll a die
}
first3[j] <- count
}

hist(first3,breaks=c(1:max(first3)),right=F,pr=T)

which(rmultinom(1,1,rep(1/6,6)) ==1)


#-------find the optimal estimate using optimize-------
likeGeom <- function(theta)
-dgeom(4,theta)

optimize(likeGeom,lower=0,upper=1)

#-----how theta changes according to the y values----
thetaseq <- seq(0.001,0.999,length=100)
plot(thetaseq,dgeom(1,thetaseq),type='l',ylim=c(0,0.3))
y <- 2:10
for (i in 1:9)
lines(thetaseq,dgeom((y[i]-1),thetaseq),col=i)

#-----R sentences for dgeom----
likeGeom <- function(y)
{
 q <-1/mean(y)
 x <- sum(dgeom(y,q,log=T))
 list(theta=q,loglik=x)
}



