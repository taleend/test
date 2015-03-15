source('clarkFunctions.R')
filename <- ('dataFACE.txt')
ynames <- 'cones'
xnames <- c('diam','trt','nfert')
tmp <- inData(filename, xnames, ynames,na.rm = T,INTERCEPT = T)
x <- tmp$x
y <- tmp$y
n <- nrow(x)
p <- ncol(x)

tmt <- as.factor(x[,'trt'])
diam <- x[,'diam']
fitb1 <- glm(y ~ tmt*diam ,family=poisson(log))
d.bins <- c(0:30) #bins for histogram
anova.1 <- anova(fitb1,test = "Chisq")
pmat <- cbind(1,contrasts(tmt))%*%matrix(coef(fitb1),2)#coeffs
colnames(pmat) <- c('intercept','slope')

#----prediction for mean response with s.e.
pred.d <- predict.glm(fitb1, data.frame(diam=d.bins,
tmt=factor(rep(0,length(d.bins)),
levels=levels(tmt))),
type="response",se.fit=T)
pred.e <- predict.glm(fitb1,data.frame(diam=d.bins,
tmt=factor(rep(1,length(d.bins)),
levels=levels(tmt))),
type="response",se.fit=T)
par(mfrow=c(3,1))
plot(diam,y)
plot(d.bins,pred.d$fit,type='l')
lines(d.bins,(pred.d$fit + 1.96*pred.d$se),lty=2)
lines(d.bins,(pred.d$fit - 1.96*pred.d$se),lty=2)
lines(d.bins,pred.e$fit,col='red')


dcone <- c(0:40) #vector for poisson probabilities
dval <- seq(5,25,by=10) #several diameters to examine
plot(d.bins,exp(pmat[1,1] + pmat[1,2]*d.bins),type='l',ylim=c(0,20))
for(i in 1:length(dval))
{
pval <- exp(pmat[1,1] + pmat[1,2]*dval[i])
bdist <- dpois(dcone,pval)
lines((20*bdist+dval[i]),dcone,type='s',col=(i+1)) #horizontal scaling
abline(v=dval[i],col=(i+4))
}

#----Probit regression---------------------
icones <- y
icones[y > 0] <- 1
fitb2 <- glm(icones ~ tmt*diam,family=binomial(probit))
anova.2 <- anova(fitb2,test="Chisq")
#recover separate intercepts and slopes using contrasts matrix
cmat <- cbind(1,contrasts(tmt)) %*% matrix(coef(fitb2),2)
dseq <- c(0:30)
#predictions for each treatment
pred.a <- predict.glm(fitb2,data.frame(diam=d.bins,
tmt=factor(rep(0,length(d.bins)),
levels=levels(tmt))),
type="response",se.fit=T)
pred.c <- predict.glm(fitb2,data.frame(diam=d.bins,
tmt=factor(rep(1,length(d.bins)),
levels=levels(tmt))),
type="response",se.fit=T)
plot(dseq,pred.a$fit,type='l',ylim=c(0,1))
lines(d.bins,(pred.a$fit + pred.a$se),lty=2)
lines(d.bins,(pred.a$fit - pred.a$se),lty=2)
lines(d.bins,pred.c$fit,col='red')
