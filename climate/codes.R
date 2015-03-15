ind <- read.csv('index by yr-difference.csv',header=T)

x11()
layout(matrix(1:4, 2, 2))
plot(ind$yr,-ind$pdsi.d,ylim=c(-5,7))
abline(h=0)
plot(ind$yr,-ind$phdi.d,ylim=c(-5,7))
abline(h=0)
plot(ind$yr,-ind$zndx.d,ylim=c(-5,7))
abline(h=0)
plot(ind$yr,-ind$pmdi.d,ylim=c(-5,7))
abline(h=0)

x11()
layout(matrix(1:4, 2, 2))
boxplot(-pdsi.d~yr,data=ind,ylim=c(-5,7),ylab="adjusted pdsi")
abline(h=0,lty=2)
boxplot(-phdi.d~yr,data=ind,ylim=c(-5,7),ylab="adjusted phdi")
abline(h=0,lty=2)
boxplot(-zndx.d~yr,data=ind,ylim=c(-5,7),ylab="adjusted zndx")
abline(h=0,lty=2)
boxplot(-pmdi.d~yr,data=ind,ylim=c(-5,7),ylab="adjusted pmdi")
abline(h=0,lty=2)


pdm <- with(ind,tapply(-pdsi.d,yr,mean))
phm <- with(ind,tapply(-phdi.d,yr,mean))
znm <- with(ind,tapply(-zndx.d,yr,mean))
pmm <- with(ind,tapply(-pmdi.d,yr,mean))

x11()
yr <- c("00","01","02","03","04","05","06","07","08","09","10","11")
plot(pdm,xlab='year',ylab="adjusted drought index",type="o",ylim=c(-3,2),axes=F)
axis(1,1:12,labels=yr)
axis(2)
abline(h=0,lty=2)

lines(phm,col=2)
points(phm,col=2)

lines(znm,col=3)
points(znm,col=3)

lines(pmm,col=4)
points(pmm,col=4)

legend(8,-1,c("pdsi","phdi","zndx","pmdi"),lty=c(1,1,1),lwd=c(2.5,2.5,2.5),col=c(1,2,3,4),bty="n")



#no adjusted

x11()
layout(matrix(1:4, 2, 2))
plot(ind$yr,ind$pdavg,ylim=c(-6,6))
abline(h=0)
plot(ind$yr,ind$phavg,ylim=c(-6,6))
abline(h=0)
plot(ind$yr,ind$znavg,ylim=c(-6,6))
abline(h=0)
plot(ind$yr,ind$pmavg,ylim=c(-6,6))
abline(h=0)

x11()
layout(matrix(1:4, 2, 2))
boxplot(pdavg~yr,data=ind,ylim=c(-6,6),ylab="pdsi")
abline(h=0,lty=2)
boxplot(phavg~yr,data=ind,ylim=c(-6,6),ylab="phdi")
abline(h=0,lty=2)
boxplot(znavg~yr,data=ind,ylim=c(-6,6),ylab="zndx")
abline(h=0,lty=2)
boxplot(pmavg~yr,data=ind,ylim=c(-6,6),ylab="pmdi")
abline(h=0,lty=2)


pdm <- with(ind,tapply(pdavg,yr,mean))
phm <- with(ind,tapply(phavg,yr,mean))
znm <- with(ind,tapply(znavg,yr,mean))
pmm <- with(ind,tapply(pmavg,yr,mean))

x11()
yr <- c("00","01","02","03","04","05","06","07","08","09","10","11")
plot(pdm,xlab='year',ylab="drought index",type="o",ylim=c(-3,2),axes=F)
axis(1,1:12,labels=yr)
axis(2)
abline(h=0,lty=2)

lines(phm,col=2)
points(phm,col=2)

lines(znm,col=3)
points(znm,col=3)

lines(pmm,col=4)
points(pmm,col=4)

legend(8,-1,c("pdsi","phdi","zndx","pmdi"),lty=c(1,1,1),lwd=c(2.5,2.5,2.5),col=c(1,2,3,4),bty="n")

s09 <- read.csv('2009summer.csv',header=T)
s11 <- read.csv('2011summer.csv',header=T)

boxplot(cbind(s09$pdsi,s09$phdi,s09$zndx,s09$pmdi,s11$pdsi,s11$phdi,s11$zndx,s11$pmdi),ylim=c(-8,8),axes=F)
axis(1,1:8,labels=c("09pd","09ph","09zn","09pm","11pd","11ph","11zn","11pm"))
axis(2,-8:8,labels=-8:8)
abline(h=0,lty=2)

#2000-2011 summer boxplots

s00 <- read.csv('2000summer.csv',header=T)
xlabel <- unique(s00$ym)
boxplot(pdsi~ym,data=s00,ylab="pdsi",ylim=c(-8,8))#axes=F)
#axis(1,1:60,labels=xlabel)
#axis(2,-8:8,labels=-8:8)
abline(h=0,lty=2)

boxplot(phdi~ym,data=s00,ylab="phdi",ylim=c(-8,8))
abline(h=0,lty=2)

boxplot(zndx~ym,data=s00,ylab="zndx",ylim=c(-8,8))
abline(h=0,lty=2)

boxplot(pmdi~ym,data=s00,ylab="pmdi",ylim=c(-8,8))
abline(h=0,lty=2)

