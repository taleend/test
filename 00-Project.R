T1 <- read.table("arc1.txt",header=F)
T1 <- na.omit(T1)
T2 <- read.table("arc3.txt",header=F)
T2 <- na.omit(T2)
T3 <- read.table("arc4.txt",header=F)
T3 <- na.omit(T3)

# exclude positive values for col5
t1 <- T1[which(T1$V5<0),]
t2 <- T2[which(T2$V5<0),]
t3 <- T3[which(T3$V5<0),]

# combine t1, t2, t3 into one dataset t
t <- rbind(t1,t2,t3)

range(t1$V3)
range(t1$V4)
range(t1$V5)

range(t2$V3)
range(t2$V4)
range(t2$V5)

range(t3$V3)
range(t3$V4)
range(t3$V5)



