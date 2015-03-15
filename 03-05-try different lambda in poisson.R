source(clarkFunctions.R)
ynames <- 'cones'
xnames <- c('diam','trt')
tmp <- inData('dataFACE.txt', xnames, ynames, tname = 'tstat',na.rm = T,INTERCEPT = F)
x <- tmp$x
y <- tmp$y          #response is no. cones
mstat <- tmp$t      #maturation index
mindex <- which(mstat == 1)        #index for mature
crange <- range(y)                 #range of counts
cvalues <- c(0:(crange[2]+1))
wlo <- which(x[,'trt'] == 0 & mstat == 1) #index for mature tree under low trt

y <- y[wlo]
x <- x[wlo]
par(mfrow=c(2,1))
plot(x,y)

y1 <- sort(y)
p <- matrix(0,4,length(y))

nx <- seq(5,30,length=4)
for (i in 1:4)
{
lambda <- 0.007*nx[i]^2 #--calulate lambda for poisson dist using beta=0.007
p[i,] <- exp(-lambda)*lambda^y1/factorial(y1)
}

plot(y1,p[1,],ylim=c(0,.9),type='l',xlim=c(0,20))
for (i in 2:4)
lines(y1,p[i,],col=i)
