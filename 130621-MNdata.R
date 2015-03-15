setwd("C:/Users/WEn/Documents/01-work/2-Classes/02-2012 spring/Med/Rwk")

M <- read.csv("MNdata.csv",header=T)

m <- matrix(c(1,1,2,1),2,2)
layout(m, widths=c(2, 1),
heights=c(1, 2))
layout.show(2)


plot(M$wavelength,M$sample.1,type="l",xlab="Wave Length",ylab="Intensity",lwd=2)
lines(M$wavelength,(M$sample.2+0.025),type="l",col=2,lwd=2)
lines(M$wavelength,M$sample.3,type="l",col=3,lwd=2)
lines(M$wavelength,(M$sample.4+0.025),type="l",col=4,lwd=2)
lines(M$wavelength,M$sample.5,type="l",col=5,lwd=2)
lines(M$wavelength,M$sample.6,type="l",col=6,lwd=2)
lines(M$wavelength,M$sample.7,type="l",col=7,lwd=2)
lines(M$wavelength,M$sample.8,type="l",col=8,lwd=2)
lines(M$wavelength,M$sample.9,type="l",col=9,lwd=2)
lines(M$wavelength,M$sample.10,type="l",col=10,lwd=2)
lines(M$wavelength,M$sample.11,type="l",col=11,lwd=2)
lines(M$wavelength,M$sample.12,type="l",col=12,lwd=2)
lines(M$wavelength,M$sample.13,type="l",col=13,lwd=2)
lines(M$wavelength,M$sample.14,type="l",col=14,lwd=2)

#plot(M$wavelength,M$sample.1,type="l",xlab="Wave Length",ylab="Intensity")

legend("topright",c("sample 1","sample 2","sample 3","sample 4","sample 5","sample 6","sample 7","sample 8","sample 9",
"sample 10","sample 11","sample 12","sample 13","sample 14"),col=1:14,lty=rep(1,14),bty="n")