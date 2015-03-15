clim <- read.table("C:/Users/CNR/Documents/1-research/1-pinemap/Climate data/02_date.csv",header=T,sep=",")
#c2 <- clim$sp02[1:200]

stnum <- unique(clim$state) # get state number

st1 <- subset(clim,state==1)
plot(st1$ob,st1$pdsi)

c41_6 <-subset(clim,state==41,clim_div==6)

plot(clim$ob, clim$zndx,type='l')
lines(clim$ob,clim$pmdi,col=2)

c1_8 <- read.table("C:/Users/CNR/Documents/1-research/1-pinemap/Climate data/02-date-state1.csv",header=T,sep=",")
dt <- c1_8$ob
zn18 <- c1_8$zndx
pm18 <- c1_8$pmdi
plot(zn18,type='l')
lines(pm18,col=2)