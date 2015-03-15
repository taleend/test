

library(Matrix)

############################### set up scale for a map
mapSetup<- function(xlim,ylim,scale){  #scale is m per inch

  px   <- diff(xlim)/scale  py   <- diff(ylim)/scale  par(pin=c(px,py))

}

############################### map species
mapSpecies <- function(specs,x,y,z,mapx=range(x),mapy=range(y),
                       scale=0,add=F){
	
	specNames <- sort(unique(specs))
	ns        <- length(specNames)
	cc        <- as.vector(sample(c(1:100),ns))
	cvec      <- 1
	if(ns > 1)cvec <- match(specs,specNames)
	
	if(scale > 0)mapsetup(mapx,mapy,scale)
	symbols(x,y,circles=z/10,inches=F,xlim=mapx,ylim=mapy,fg=cc[cvec],add=add)
}
	
############################## empirical distribution functionmy.ecdf <- function(x){  n  <- length(x)  xx <- sort(unique(x))                 # find unique values of x  fx <- cumsum(tabulate(match(x,xx))/n) # number of observations for each x  list(x = xx,fx = fx)                  # return the unique values and edcf}
############################# sample plots from a map
samplePlots <- function(mapx,mapy,wide,mapscale,PLOTIT = T){
  yt      <- seq(mapy[1],mapy[2],by=wide)      #y grid locations  yl      <- length(yt)  xt      <- seq(mapx[1],mapx[2],by=wide)      #x grid locations
  xl      <- length(xt)  mapgrid <- cbind(rep(xt,each=yl),rep(yt,xl)) #x and y locations  loc     <- mapgrid[sample(yl*xl,nplot,replace=F),] #samples from grid  sl <- .5*wide                             #plot edges  xbound <- cbind((loc[,1]-sl),(loc[,1]+sl))  ybound <- cbind((loc[,2]-sl),(loc[,2]+sl))  xindex <- c(1,2,2,1,1)  yindex <- c(1,1,2,2,1)
  
  specname  <- sort(unique(treedata[,'species']))
  nspec     <- length(specname) # plot.data  <- numeric(0)                   #list of observations  tableSpec <- matrix(0,nrow=nspec,ncol=nplot)  rownames(tableSpec) <- specname  if(PLOTIT)plot(-1000,0,xlim=mapx,ylim=mapy)  for(i in 1:nplot){  # extract trees on sample plot i and store them in table.spec    xt <- !is.na(cut(treedata[,'x'],breaks=xbound[i,],exclude=NA))    yt <- !is.na(cut(treedata[,'y'],breaks=ybound[i,],exclude=NA))
    
    tmp <- treedata[xt & yt,]
    if(nrow(tmp) > 0){
    	
      tableSpec[,i] <- table(tmp[,'species'])      if(PLOTIT)symbols(tmp[,'x'],tmp[,'y'],circles=tmp[,'dbh']/20,inches=F,add=T)
    }  # draw a box around each plot
    if(PLOTIT){      xvec <- xbound[i,xindex]      yvec <- ybound[i,yindex]      lines(xvec,yvec)
    }  }

  tableSpec
}
######################################

inData <- function(filename, xnames = NULL, ynames = NULL, tname = NULL, 
                   iname = NULL, na.rm = F, INTERCEPT = F){  
                   	
  #read in data file, return design matrix x, response y
  #xnames, ynames, tname, iname are column headings in filename
  #time indicator t
  #individual indicator i

  data <- read.table(filename,header=T)
  
  if(is.atomic(xnames)){
    x    <- data[,xnames]
    if(!is.matrix(x))x <- as.matrix(x)
    if(INTERCEPT){
      intercept <- rep(1,nrow(data))
  	   x <- cbind(intercept,x)
  	   colnames(x) <- c('intercept',xnames)
    }
  }
  if(is.atomic(ynames)){
    y <- data[,ynames]
    if(!is.matrix(y))y <- as.matrix(y) 
    y  <- matrix(y,nrow(data),length(ynames))
    colnames(y) <- ynames
  }
  
  wf <- c(1:nrow(data))
  
  if(na.rm){
  	 wf <- which(is.finite(apply(x,1,sum)) & is.finite(apply(y,1,sum)))
  	 x  <- x[wf,]
  	 y  <- y[wf,]
  }
  
  z  <- list(x = x, y = y)
  
  if(is.atomic(tname))z$t <- data[wf,tname]
  if(is.atomic(iname))z$i <- data[wf,iname]
  z
}
####################################################
simX <- function(n,loX,hiX){                #generate design matrix

  k <- length(loX)
  x <- matrix(1,n,k)
  for(j in 1:k)x[,j] <- runif(n,loX[j],hiX[j])
  x
}
####################################################
simY <- function(x,b,LIKE,r = 1,size=rep(1,nrow(x)),sigma = 0,Effort = 1){     #simulate response

  u <- x%*%b

  if(LIKE == 'norm')return( rnorm(length(u),u,sqrt(sigma)) )
  if(LIKE == 'pois'){
  	 u <- exp(u)*Effort
  	 return( rpois(nrow(x),u) )
  } 
  if(LIKE == 'binom'){
  	 u <- inv.logit(u)
  	 return( rbinom(nrow(x),1,u) )
  }
  if(LIKE == 'multinom'){
     zs <- apply(exp(u),1,sum)
     z1 <- 1/(1 + zs)
     zm <- exp(u)/ (1 + zs)
     u  <- cbind(zm,z1)
     return( myrmultinom(size,u) )
  }
  if(LIKE == 'mvnorm'){
    u <- myrmvnorm(n,u,sigma)
    return( u )
  }
  if(LIKE == 'mvnorm-multinom'){
    u <- myrmvnorm(n,u,sigma)
    zs <- apply(exp(u),1,sum)
    z1 <- 1/(1 + zs)
    zm <- exp(u)/ (1 + zs)
    u  <- cbind(zm,z1)
    return( myrmultinom(size,u) )
  }
  numeric(0)
}
#########################################
simMVData <- function(LIKE,k,r){
	
   xnames <- c('intercept',paste('x',c(1:(k-1)),sep='-') )
   p      <- length(xnames)
   x      <- matrix(runif(n*k,-4,4),n,k)
   x[,1]  <- 1
   colnames(x) <- xnames
   
   b           <- matrix(runif(k*ns,-2,1),k,r)
   specnames   <- paste('S',c(1:r),sep='-')
   colnames(b) <- specnames
   rownames(b) <- xnames
   sigma <- diag(runif(r,0,.1/r),r)
   rownames(sigma) <- specnames
   colnames(sigma) <- specnames

   mu <- x %*% b
   y     <- myrmvnorm(n,mu,sigma)
   sigma <- crossprod(y - x%*%b)/(n-k)
   y     <- myrmvnorm(n,mu,sigma)
   colnames(y) <- specnames
   
   if(LIKE == 'mvnorm-Pois')    z <- matrix(rpois(length(y),exp(y)),n,r,byrow=F)
   if(LIKE == 'mvnorm-multinom'){
   	
   	  p1 <- matrix(rbeta(n*r,.2,.2),n,r)
   	  p1[cbind(c(1:n),sample(c(1:r),n,replace=T))] <- 10
   	  p1 <- p1/matrix(apply(p1,1,sum),n,r)
   	  zz <- myrmultinom(10,p1)
   	  
   	  b  <- bInitMNomLogit(x,zz,size)
   	
   	  sigma <- sigma[1:(r-1),1:(r-1)]
     y  <- myrmvnorm(n, x %*% b,sigma)
   	  zs <- apply(exp(y),1,sum)
     z1 <- 1/(1 + zs)
     zm <- exp(y)/ (1 + zs)
     y  <- cbind(zm,z1)
   	  z  <- myrmultinom(size,y)
   	  plot(y,jitter(z))
   	}
   
   list(x = x, b = b, y = y, z = z, s = sigma, specnames = specnames)
}

#########################################
myrmultinom <- function(size,p){  

  #n multinomial r.v. for a n by ncol(p) matrix of probs
  #each row of p is a probability vector

  n     <- nrow(p)
  J     <- ncol(p)
  y     <- matrix(0,n,J)
  sizej <- size
  if(length(size) == 1)sizej <- rep(size,n)
  sumj  <- rep(0,n)
  dpj   <- rep(1,n)
  pj    <- p
  wj    <- c(1:n)

  for(j in 1:(J-1)){
    a     <- round(pj[wj,1],10)
    y[wj,j] <- rbinom(length(wj),sizej[wj],a)
    sumj  <- sumj + y[,j]
    sizej <- size - sumj
    dpj   <- dpj - p[,j]
    pj    <- p[,c((j+1):J)]/dpj
    wj    <- which(sumj < size,arr.ind=T) 
  }

  y[,J] <- size - apply(y,1,sum)
  y
}


##############################################################

#like_geom <- function(theta) -sum(dgeom((y-1),theta,log=T)) 

get_lf <- function(bounds,FUN){  #gets the likelihood function

  x <- seq(bounds[1],bounds[2],length=1000)
  q <- x*0
  for(i in 1:1000){q[i] <- FUN(x[i])}
  list(x = x, likelihood = q)
}


treebytime <- function(iindex,tindex,var){ #matrix for var with individuals by time
  vmat <- matrix(NA,max(iindex),max(tindex))
  vmat[cbind(iindex,tindex)] <- var
  vmat
}

like_conedata <- function(b){
  g <- b*x^2
  -sum(dpois(y,g,log=T))
}

like_coneyear <- function(b){
  g     <- b[year]*diamMature^2
  -sum(dpois(coneMature,g,log=T))
}

mleGeom <- function(n){ 
	#mle and lnL for geometric
	
  q <- 1/mean(n)               #ML estimate
  x <- sum(dgeom(n,q,log=T))   #lnL
  list(theta = q, loglik = x)
}


like_weib <- function(param){  #Weibull likelihood  -sum(dweibull(y,param[1],1/param[2],log=T))}

like_exp_binom <- function(rho){  theta <- exp(-rho*20)  -dbinom(6,30,theta,log=T)}

likeBernLogit <- function(pars){
  
  g      <- x %*% pars
  theta  <- inv.logit(g)
  s <- y*log(theta) + (1 - y)*log(1 - theta)
  -sum(s)
}


like_norm <- function(par){
  -sum(dnorm(y,par[1],sqrt(par[2]),log=T))
}

surv_prodlim <- function(surv,life){ #product limit estimator of survival
  ctable <- table(surv,life)         #survival by sample date
  yrsum  <- apply(ctable,2,sum)     
  nt     <- n - cumsum(yrsum)        #no. at risk by date
  nt     <- c(n,nt[-nyr])            #1-yr shift
  dt     <- ctable[1,]              #no. died by date
  lamdat <- dt/nt

  s <- rep(1,nyr)                    #survival function
  for(j in 2:nyr)s[j] <- s[j-1]*(1 - lamdat[j])
  list(lamda = lamdat, survprod = s)
}

surv_weib <- function(b){

  lamda <- b[1]
  c     <- b[2]
  beta  <- b[3:length(b)]
  xb    <- x %*% beta
  h0    <- log(c) + c*log(lamda) + (c-1)*log(life[notc]) #only uncensored
  l0    <- -(lamda*life)^c                                 #all individuals
  h     <- h0 + xb[notc]
  s     <- l0*exp(xb)
  s[notc] <- s[notc] + h
  -sum(s)                      #returns -lnL
}

inv.logit <- function(f)(1/(1 + exp(-f)))

logit <- function(f) log(f/(1 - f))

dinvGamma <- function(x,a,b,log=FALSE){
	
	p <- a*log(b) - lgamma(a) - (a+1)*log(x) - b/x
	if(log)return(p)
	exp(p)
}


qprob <- function(p){  #evaluate likelihood for pathogen model: multinomial    theta <- p[1]    phi   <- p[2]    s0    <- p[3]    s1    <- p[4]    q00 <- (1 - theta)*(1 - s0) + theta*(1 - phi)*(1 - s1)    q01 <- (1 - theta)*s0 + theta*(1 - phi)*s1    q10 <- theta*phi*(1 - s1)    q11 <- theta*phi*s1    c(q00,q01,q10,q11)}

like_multinom <- function(pars){   q <- qprob(pars)  -sum(nvec*log(q))}

like_pois_probit <- function(par){  a0 <- par[1]  a1 <- par[2]  b  <- par[3] #diameters for zeros  xd     <- xtmp[ytmp == 0]  qd     <- b*xd^2  theta1 <- pnorm(xd,a0,a1)               #probit  pz     <- sum(log((1-theta1) + theta1*dpois(0,qd))) #trees having cones  xd     <- xtmp[ytmp > 0]  yd     <- ytmp[ytmp > 0]  qd     <- b*xd^2  theta2 <- pnorm(xd,a0,a1)  pk     <- sum(log(theta2) + dpois(yd,qd,log=T))  -pz - pk}

#elephant movement example

gk <- function(locnow,q,ss){   dx   <- distance[locnow,]         pvec <- q*exp(-(dx^2)/ss)   pvec/sum(pvec)}

move <- function(locnow,q,ss){   gt  <- gk(locnow,q,ss)      #prob moving to other patches   wk  <- rmultinom(1,1,gt)    #movement vector   c(1:npatch) [wk == 1]       #extract new patch location}update_z <- function(){  #update unknown locations

    for(t in tindex){   #only simulate hidden states      pk   <- gk(z[t-1],qual[(t-1),],sg)*(1-pg[z[t]])      if(t < nt)pk <- pk *gk(z[t+1],qual[(t+1),],sg)      pk   <- pk/sum(pk)      pv   <- rmultinom(1,1,pk)      z[t] <- pseq[pv == 1]    }    z[tindex]
}

update_p <- function(){  #detection probabilities    p1 <- sp1 + uk            p2 <- sp2 + tabulate(z,nbins=npatch) - uk    rbeta(npatch,p1,p2)
}

update_sigma <- function(){  #movement parameter

    props <- tnorm(1,0,Inf,sg,2)   #proposal for sigma
    psum1 <- 0    psum2 <- 0    for(t in 2:(nt-1)){     psum1 <- psum1 + log(gk(z[t],qual[t,],sg)[z[t+1]])     psum2 <- psum2 + log(gk(z[t],qual[t,],props)[z[t+1]])    }    a <- exp(psum2 - psum1)    z <- runif(1,0,1)    if(z < a) sg <- props
    sg
}

#simulation

fitturn <- function(par){ #fit aphid turn data  b0 <- par[1]  b1 <- par[2]  b2 <- par[3]  mu <- (b0 + b1*cos(b2 + x))/2/pi/b0  #mean function  -sum(dpois(y,2*pi/m*mu,log=T))       # - ln likelihood}

rturn <- function(n,cf){  #generate random deviates  x  <- numeric(0)  nx <- 0  while(nx < n){    q  <- runif(n-nx,-pi,pi)    a  <- (b0 + b1*cos(b2 + q))/cf                #unnormalized    z  <- runif(n-nx,0,1)
    xz <- q[z < a]    x  <- c(x,xz)    nx <- length(x)  }  x}

rtrunc <- function(n,lo,hi,p1,p2,FN){    #truncated using inv dist sampling   #truncated 2 parameter distribution	
  if(FN == 'invGamma'){
  	z1 <- 1 - pgamma(1/lo,p1,p2)
  	z2 <- 1 - pgamma(1/hi,p1,p2)
   z  <- 1 - runif(n,z1,z2)
   return(1/qgamma(z,p1,p2))
  }
  pf <- paste('p',FN,sep='')  qf <- paste('q',FN,sep='')  pf1 <- match.fun(pf)  qf1 <- match.fun(qf)
  
  z1 <- pf1(lo,p1,p2)
  z2 <- pf1(hi,p1,p2)
  z  <- runif(n,z1,z2)
    qf1(z,p1,p2)}


tnorm <- function(n,lo,hi,mu,sig){   

  #normal truncated lo and hi

  if(length(lo) == 1 & length(mu) > 1)lo <- rep(lo,length(mu))
  if(length(hi) == 1 & length(mu) > 1)hi <- rep(hi,length(mu))

  q1 <- pnorm(lo,mu,sig)
  q2 <- pnorm(hi,mu,sig)

  z <- runif(n,q1,q2)
  z <- qnorm(z,mu,sig)
  z[z == Inf]  <- lo[z == Inf]
  z[z == -Inf] <- hi[z == -Inf]
  z
}


acceptMH <- function(p0,p1,x0,x1,BLOCK=F){   #accept for M, M-H	# if BLOCK, then accept as a block,	# otherwise, accept individually  nz          <- length(x0)  #no. to accept  if(BLOCK)nz <- 1    a    <- exp(p1 - p0)       #acceptance PR  z    <- runif(nz,0,1)  keep <- which(z < a,arr.ind=T)    if(BLOCK & length(keep) > 0)x0 <- x1  if(!BLOCK)                  x0[keep] <- x1[keep]             ac <- length(keep)          list(x = x0, accept = ac)}

pmake <- function(pars){    p   <- pars[1]                   #frequency of a  f   <- pars[2]                   #inbreeding coefficient   paa <- p^2 + f*p*(1 - p)	  #frequency of p.aa  pab <- 2*p*(1 - p)*(1 - f)	  #frequency of p.ab  pbb <- (1 - p)^2 + f*p*(1 - p)  #frequency of p.bb  c(paa,pab,pbb)}

minf <- function(p){   #minimum inbreeding coefficient  lof <- rep(-1,length(p))  lop <- -(1 - p)/p  hip <- -p/(1 - p)  lof[p > (1 - p)] <- lop[p > (1 - p)]  lof[p < (1 - p)] <- hip[p < (1 - p)]  lof}pf.update <- function(p,f){  #log posterior  multlik(c(p,f)) + dbeta(p,g1,g2,log=T) +            dnorm(f,priorF,priorFSD,log=T)}
multlik <- function(pars){  #multinom likelihood   p    <- pmake(pars)  pmat <- matrix(rep(p,npop),3,npop)  sum(y*log(pmat))}


tnorm.mvt <- function(avec,muvec,smat,lo,hi){   

  # truncated multvariate normal
  # muvec is the vector of means
  # smat is the covariance matrix 

  for(k in 1:length(muvec)){
    skk    <- smat[-k,-k]
    piece1 <- smat[-k,k] %*% solve(skk)
    muk <- muvec[k] + piece1 %*% (avec[-k] - muvec[-k])
    sgk <- as.numeric(smat[k,k] - piece1 %*% smat[k,-k])
    avec[k] <- tnorm(1,lo[k],hi[k],muk,sqrt(sgk))
  }
  avec
}

predVsObs <- function(o,p){ 
	
  #o  - length n vector of obs or true values
  #p - ng by n matrix of estimates
  
  n <- length(o)
  y <- apply(p,2,quantile,c(.5,.025,.975))

  plot(o,y[1,])
  for(j in 1:n)lines(c(o[j],o[j]),y[2:3,j])
  abline(0,1,lty=2)
  y
}

update_pf <- function(){

  propp  <- tnorm(1,.02,.98,pg,.002) #propose pa  fl     <- minf(propp)              #lower f limit for pa  propf  <- tnorm(1,fl,1,fg,.01)     #propose f > fl  pnew   <- pf.update(propp,propf)  pnow   <- pf.update(pg,fg)  atmp   <- acceptMH(pnow,pnew,c(pg,fg),c(propp,propf),BLOCK=T)  pg     <- atmp$x[1]  fg     <- atmp$x[2]  ac     <- atmp$accept

  list(pg = pg, fg = fg, accept = ac)
}

##sample from posteriors for a linear regression

b.update <- function(){  V <- solve(crossprod(x)/sg + priorIV)  v <- crossprod(x,y)/sg + priorIV %*% priorB  t(rmvnorm(1,V%*%v,V))}

v.update <- function(){  u1 <- s1 + n/2  u2 <- s2 + .5*crossprod(y - x%*%bg)  1/rgamma(1,u1,u2)}
#######################################################
plotObsPred <- function(obs,yMean,ySE){
	
	plot(obs,yMean)
	ylo <- yMean - 1.96*ySE
	yhi <- yMean + 1.96*ySE
	for(i in 1:length(obs))lines(c(obs[i],obs[i]),c(ylo[i],yhi[i]))
	abline(0,1,lty=2)

}
#######################################################
deviance <- function(y,x,b,s=0,LIKE){
	    
    if(LIKE == 'norm')  dv <- dnorm(y,x%*%b,sqrt(s),log=T) 
    if(LIKE == 'pois')  dv <- dpois(y,exp(x%*%b),log=T)
    if(LIKE == 'binom') dv <- dbinom(y,1,inv.logit(x%*%b),log=T)
    if(LIKE == 'mvnorm')dv <- dmvnorm(y,x%*%b,s,log=T)
    if(LIKE == 'multinom')dv <- multinomLike(y,x,b)
    -2*dv
}
#######################################################
multinomLike <- function(y,x,b){  #log likelihood multinomial logit
	
  tiny <- 1e-20
  huge <- 1 - tiny

	  z <- x%*%b
     zs    <- apply(exp(z),1,sum)
     z1    <- 1/(1 + zs)
     zm    <- exp(z)/ (1 + zs)
     z2  <- cbind(zm,z1)
     z2[z2 < tiny] <- tiny
     z2[z2 > huge] <- huge
     y*log(z2)
}
#######################################################

gibbsLoop <- function(LIKE,ng,x,y,b,sigma = 0,z = numeric(0)){

  kx     <- length(b)
  
  bgibbs <- matrix(NA,ng,kx)
  colnames(bgibbs) <- paste('b',c(1:kx),sep='-')
  sgibbs <- matrix(NA,ng,max(1,length(sigma)))
  
  r <- 1                       #responses
  if(is.matrix(y))r <- ncol(y)
  
  if(sigma[1] == 0)sgibbs <- numeric(0) #variances
  
  sg <- sigma
  
  if(is.matrix(sigma)){                  #Wishart prior
    sg        <- prior.W
    colnames(sgibbs) <- outer(rownames(sigma),rownames(sigma),paste,sep='_') 
    colnames(bgibbs) <- as.vector(outer(colnames(x),colnames(b),paste,sep='_'))
  }
  
  pBVar <- solve(crossprod(x))

  pred <- pred2 <- rep(0,nrow(x)*r)

  bg <- b
  yg <- y
  
  y1 <- yg
  
  if(LIKE == 'pois') bg <- pBVar%*%crossprod(x,log(y + .1))
  if(LIKE == 'binom'){
  	yl <- y
  	yl[yl == 0] <- .1
  	yl[yl == 1] <- .9
  	bg <- pBVar%*%crossprod(x,logit(yl))
  }
  
  if(LIKE == "mvnorm-multinom"){
  	y1   <- prob2Logit(yg)
  }
  if(LIKE == 'multinom')pBVar <- diag(.01,k*(r-1))
  
  dev <- 0

  for(g in 1:ng){

    bg <- bUpdateGibbs(x,y1,bg,LIKE,sigma=sg,pBVar)
    
    if(sigma[1] > 0 & length(sigma) == 1){
    	 sg <- sigmaUpdate(x,y1,bg)
    }
    
    if(length(grep('mvnorm',LIKE)) > 0){
    	 sinv <- wishsamp(x,y1,bg)
    	 sg   <- solve(sinv)
    }
    
    if(LIKE == "mvnorm-multinom"){
    	y1 <- ysampMvnormMultinom(x,y1,z,bg,sg)$y
    	yg <- logit2Prob(y1)
    }
    
    dev <- dev + sum(deviance(y1,x,bg,sg,LIKE))

    py    <- as.vector(simY(x,bg,LIKE,r,size,sg))
    pred  <- pred + py
    pred2 <- pred2 + py^2

    bgibbs[g,] <- bg
    if(length(sgibbs) > 0)sgibbs[g,] <- sg
    
    if(g %in% c(200,500)){
    	pBVar <- cov(bgibbs[10:g,])
    }
  }

  ymean <- pred/ng
  yse   <- sqrt(pred2/ng - ymean^2)
  
  if(r > 1){
  	ymean <- matrix(ymean,n,r)
  	yse   <- matrix(yse,n,r)
  }
  bmean <- apply(bgibbs,2,mean)
  bmean <- matrix(bmean,nrow(bg),ncol(bg))
  smean <- numeric(0)
  if(length(sg) > 1) smean <- matrix(apply(sgibbs,2,mean),nrow(sg),ncol(sg))
  if(length(sg) == 1)smean <- mean(sgibbs)
  
  meanDev <- dev/ng
  
  pd  <- meanDev - sum(deviance(y1,x,bmean,smean,LIKE))
  dic <- 2*pd + meanDev

  list(bgibbs = bgibbs,sgibbs = sgibbs, ymean = ymean, yse = yse, dic = dic)
}
####################################################
bUpdateGibbs <- function(x,y,b,LIKE,sigma = 0,pBVar=0){

  if(LIKE == 'norm')    return( bUpdateNorm(x,y,sigma) )
  if(LIKE == 'multinom')return( bUpdateMNom(x,y,b,pBVar) )
  if(LIKE %in% c('mvnorm','mvnorm-multinom'))return( bupdateMVNorm(x,y,b,sigma) )

  b <- matrix(b,length(b),1)
  c <- t(myrmvnorm(1,t(b),pBVar))

  znow <- x%*%b
  znew <- x%*%c

  if(LIKE == 'pois'){
     pnow <- dpois(y,exp(znow),log=T)
     pnew <- dpois(y,exp(znew),log=T)
  }
  if(LIKE == 'binom'){
     pnow <- dbinom(y,1,inv.logit(znow),log=T)
     pnew <- dbinom(y,1,inv.logit(znew),log=T)
  }

  pnow <- sum(pnow) + mydmvnorm(t(b),priorB,priorVB,log=T)
  pnew <- sum(pnew) + mydmvnorm(t(c),priorB,priorVB,log=T)

  a <- exp(pnew - pnow)
  z <- runif(1,0,1)
  if(z < a)b <- c
  b
}
####################################################
ysampMvnormPois <- function(x,y,b,sigma,Effort){
	
  r <- ncol(y)
  k <- nrow(b)

  propy <- matrix(rnorm(length(y),y,.02),n,r,byrow=F)
  
  pnow <- rep(0,n)
  pnew <- rep(0,n)
  
  for(i in 1:n){
    pnow[i] <- mydmvnorm(y[i,],(x %*% bgM)[i,], sigma,log=T)
    pnew[i] <- mydmvnorm(propy[i,],(x %*% bgM)[i,], sigma,log=T)
  }
  
  pnow <- pnow + rowSums(dpois(z, Effort*exp(y),log=T))
  pnew <- pnew + rowSums(dpois(z, Effort*exp(propy),log=T))

  a <- exp(sum(pnew) - sum(pnow))
  zz <- runif(1,0,1)
  accept <- 0
  if(zz < a){
     y <- propy
     accept <- 1
  }
  list(y = y, a = accept)
}
#########################
ysampMvnormMultinom <- function(x,y,z,b,sigma){  #sample y on MVlogit scale

  r  <- ncol(y)
  k  <- nrow(b)
  
  propy <- matrix(rnorm(length(y),y,.01),n,r,byrow=F)
  
  pnow <- rep(0,n)
  pnew <- rep(0,n)
  
  for(i in 1:n){
    pnow[i] <- mydmvnorm(y[i,],(x %*% b)[i,], sigma,log=T)
    pnew[i] <- mydmvnorm(propy[i,],(x %*% b)[i,], sigma,log=T)
   }

  zs    <- apply(exp(y),1,sum)
  z1    <- 1/(1 + zs)
  zm    <- exp(y)/ (1 + zs)
  znow  <- cbind(zm,z1)

  zs    <- apply(exp(propy),1,sum)
  z1    <- 1/(1 + zs)
  zm    <- exp(propy)/ (1 + zs)
  znew  <- cbind(zm,z1)

  pnow <- apply(pnow + z*log(znow),1,sum)
  pnew <- apply(pnew + z*log(znew),1,sum)

  a  <- exp(pnew - pnow)
  zz <- runif(length(a),0,1)
  y[zz < a,] <- propy[zz < a,]
  accept <- length(zz[zz < a])

  list(y = y, a = accept)
}

###################################
rwish <- function (v, S) {
    if (!is.matrix(S)) 
        S <- matrix(S)
    if (nrow(S) != ncol(S)) {
        stop(message = "S not square in rwish().\n")
    }
    if (v < nrow(S)) {
        stop(message = "v is less than the dimension of S in rwish().\n")
    }
    p  <- nrow(S)
    CC <- chol(S)
    Z  <- matrix(0, p, p)
    diag(Z) <- sqrt(rchisq(p, v:(v - p + 1)))
    Z[diag(Z < 100)] <- 100
    if (p > 1) {
        pseq <- 1:(p - 1)
        Z[rep(p * pseq, pseq) + unlist(lapply(pseq, seq))] <- rnorm(p * 
            (p - 1)/2)
    }
    crossprod(Z %*% CC)
}

y2zMVlogit <- function(y){     #multivar logit to fractions
  zs   <- apply(exp(y),1,sum)
  z1   <- 1/(1 + zs)
  zm   <- exp(y)/ (1 + zs)
  cbind(zm,z1)
}

z2yMVlogit <- function(z){     #fractions to multivar logit
 # log(z[,-(ns+1)]*(1 + z))       
  
  log(z[,-(ns+1)]/
  (1 - apply(z[,-(ns+1)],1,sum)))
  
}

###################################
bUpdateMVNorm <- function(x,y,b,sigma){  # update b's for mvnorm

  r      <- ncol(b)
  k      <- nrow(b)
  
  bigv   <- solve(crossprod(x))      #multivariate
  smallv <- crossprod(x,y)
  mu     <- bigv%*%smallv
  vaa    <- kronecker(sigma,bigv)   
  bg     <- matrix(myrmvnorm(1,matrix(mu,k*r,1),sigma=vaa),k,r,byrow=F)
  
  colnames(bg) <- colnames(b)
  rownames(bg) <- rownames(b)
  bg
}

###################################
bInitMNomLogit <- function(x,y,size){   #initialize coeffs for multinom logit
	
	if(length(size) == 1)y  <- y/size
	if(length(size) == nrow(x))y <- y/matrix(size,nrow(x),ncol(y))
	bj <- matrix(NA,ncol(x),ncol(y)-1)
	rownames(bj) <- colnames(x)
	colnames(bj) <- colnames(y)[-ncol(y)]
	
	for(j in 1:(ncol(y)-1)){
	  yj <- y[,j]
	  y0 <- yj + .01
	  y0[y0 > 1] <- .99
	  y0 <- logit(y0)
	  bj[,j] <- solve(crossprod(x))%*%crossprod(x,y0)
	 }
	 bj
}
	
####################################################
bUpdateMNom <- function(x,y,b,pBVar){

  bvec <- as.vector(b)

  tmp  <- bmultiProp(b,pBVar)
  c    <- tmp$c
  cvec <- tmp$cvec
  
  pnow <- multinomLike(y,x,b)
  pnew <- multinomLike(y,x,c)

  pnow <- sum(pnow) + mydmvnorm(bvec,priorB,priorVB,log=T)
  pnew <- sum(pnew) + mydmvnorm(cvec,priorB,priorVB,log=T)

  a <- exp(pnew - pnow)
  z <- runif(1,0,1)
  if(z < a)b <- c
  b
}

####################################################
bmultiProp <- function(b = matrix(0,k,r-1),pBVar=diag(.1,k*(r-1))){  
	
    bvec <- as.vector(b)
    cvec <- myrmvnorm(1,t(bvec),pBVar)
    c    <- matrix(cvec,nrow(b),ncol(b))
    
    if('lob' %in% ls()){                     #if lob and hib available, use tnorm.mvt
    	cvec <- tnorm.mvt(t(bvec),pBVar,lob,hib)
      c    <- matrix(cvec,nrow(b),ncol(b))
    }

  list(c = c, cvec = cvec)
}

####################################################
bUpdateNorm <- function(x,y,sigma,priorB,priorIV){

  V <- solve(crossprod(x)/sigma + priorIV)
  v <- crossprod(x,y)/sigma + priorIV %*% priorB
  t(myrmvnorm(1,t(V%*%v),V))
}

####################################################
sigmaUpdate <- function(x,mu,s1,s2){

  u1 <- s1 + length(y)/2
  u2 <- s2 + .5*sum( (y - mu)^2 )
  1/rgamma(1,u1,u2)
}
##########################################################

processStates <- function(xchains,y){
	
	xci <- apply(xchains,2,quantile,c(.5,.025,.975))
	
	yr <- range(xci)
	
	plot(xci[1,],type='l',lwd=2,ylim=yr)
	for(j in 2:3)lines(xci[j,],lty=2)
	points(y)
}
##########################################################
processPars <- function(xgb,xtrue=numeric(0),CPLOT=F,DPLOT=F,
                        sigOnly = F,burnin=1,xlimits = NULL){  

  #xg      - matrix of gibbs chains
  #xtrue   - true values (simulated data)
  #CPLOT   - if T, plot chains
  #DPLOT   - if T, plot density
  #burnin  - analyze chains > burnin
  #xlimits - xlimits for plot
  #sigOnly - plot only parameters that 95% CI does not include 0
  
  if(!is.matrix(xgb))xgb <- matrix(xgb,ncol=1)
  if(is.null(colnames(xgb)))colnames(xgb) <- paste('V',c(1:ncol(xgb)),sep='-')
  
  if(sigOnly){
    wi   <- grep('intercept',colnames(xgb))      #extract covariates for plotting
    btmp <- xgb
    if(length(wi) > 0){
    	btmp <- xgb[,-wi]
      if(length(xtrue) > 0)xtrue <- xtrue[-wi]
    }

    wq   <- apply(btmp,2,quantile,c(.025,.975))  #extract parameters != 0
    wq   <- which(wq[1,] < 0 & wq[2,] > 0)
    if(length(wq) > 0){
      xgb  <- btmp[,-wq]
      if(length(xtrue) > 0)xtrue <- xtrue[-wq]
    }
   }

  if(!is.matrix(xgb))xgb <- as.matrix(xgb)
  if(burnin > 1){
  	     if(burnin > (nrow(xgb) + 100))stop("burnin too large")
  	     xgb <- xgb[-c(1:burnin),]
  }
  if(!is.matrix(xgb))xgb <- as.matrix(xgb)
  nc <- ncol(xgb)
  nf <- round(sqrt(nc),0)

  out <- t(rbind(apply(xgb,2,mean),apply(xgb,2,quantile,c(.025,.975))))
  if(!is.null(colnames(xgb)))rownames(out) <- colnames(xgb)
  colnames(out) <- c('mean','0.025','0.975')
  if(length(xtrue) > 0){
    out <- cbind(out,xtrue)
    colnames(out) <- c('mean','0.025','0.975','true value')
  }

  armat <- matrix(0,nc,10)  #for AR model
  
 # for(j in 1:nc)armat[j,] <- ar(xgb[,j],aic=F,order.max = 10)$ar[1:10]

 # if(!is.null(colnames(xgb)))rownames(armat) <- colnames(xgb)
#  colnames(armat) <- paste('AR',c(1:10),sep='-')

  if(CPLOT | DPLOT)par(mfrow=c((nf+1),nf))
  if(CPLOT & DPLOT)par(mfrow=c((nf+1),nc))

  if(CPLOT){
      for(j in 1:nc){
       plot(xgb[,j],type='l')
       abline(h=out[j,],lty=2)
       if(length(xtrue) > 0)abline(h=xtrue[j],col='red')
       title(colnames(xgb)[j])
     }
  }
  xlims <- xlimits
  if(DPLOT){
      for(j in 1:nc){
        xj <- density(xgb[,j])
        if(is.null(xlimits))xlims <- range(xj$x)
        plot(xj$x,xj$y,type='l',xlim=xlims)
        abline(v=out[j,],lty=2)
        if(length(xtrue) > 0)abline(v=xtrue[j],col='red')
        title(colnames(xgb)[j])
     }
  }
  list(summary = signif(out,4))

}


gibbsStates <- function(time,x,y,xtrue=numeric(0)){  nc <- ncol(x)  par(mfrow=c(1,1))  plot(time,y)  xci <- apply(x,2,quantile,c(.5,.025,.975))  lines(time,xci[1,])  lines(time,xci[2,],lty=2)  lines(time,xci[3,],lty=2)  if(length(xtrue) > 0)lines(time,xtrue,col='red')}
#############state space models

updateSSB1 <- function(states,sg,priorB,priorIVB){  #intercept SS model
	
	V <- 1/( (length(y) - 1)/sg + priorIVB)
	v <- (states[length(states)] - states[1])/sg + priorB*priorIVB
	rnorm(1,V*v,sqrt(V))
}
##################################
updateSSRW <- function(states,y,missing,tg,sg){        #state-space random walk 
	#update continuous states, random walk
	#missing times, obs y, obs error tg, process error sg

  for(t in 1:nt){

    VI <- 0
    v  <- 0

    if(!t %in% missing){          #observations
      v  <- y[t]/tg
      VI <- 1/tg
    }

    if(t < nt){              #t+1 term excluded for last 
      v  <- v + states[t+1]/sg 
      VI <- VI + 1/sg
    }

   if(t > 1){                #t-1 term excluded for 1st 
      v  <- v + states[t-1]/sg
      VI <- VI + 1/sg
   }

   V     <- 1/VI
   states[t] <- rnorm(1,V*v,sqrt(V))
  }
  states
}

updateSigmaIG <- function(y,mu,s1,s2){
	
	u1 <- s1 + length(y)/2
	u2 <- s2 + .5*sum( (y - mu)^2 )
	1/rgamma(1,u1,u2)
}

SSupdateS <- function(FUN,...){        #process error
  
  FUN <- match.fun(FUN)
  x2 <- xg[-1]
  x1 <- xg[-nt] +  FUN(...)
  u1 <- s1 + (nt-1)/2
  u2 <- s2 + .5*sum( ((x2 - x1)^2) )
  1/rgamma(1,u1,u2)
}

SSupdateT <- function(){

  u1 <- v1 + (nt-nm)/2
  u2 <- v2 + .5*sum( (y[-wm] - xg[-wm])^2)
  1/rgamma(1,u1,u2)
}

SSRWupdateB <- function(){  #regression parameters in rw model  V <- 1/( (nt - 1)/sg + 1/vb )  v <- (x[nt] - x[1])/sg + b0/vb  rnorm(1,V*v,sqrt(V))}
SSLGupdateB <- function(){  z <- xg[-1] - xg[-nt]  X <- cbind(rep(1,(nt-1)),exp(xg[-nt]))  V <- solve(crossprod(X)	+priorIVb)  v <- crossprod(X,z) + priorIVb%*%priorb  rmvnorm(1,V%*%v,V)  }
SSMetUpdateX <- function(){     #Metropolis for state-space  for(t in 1:nt){    xprop <- rnorm(1,xg[t],.1)    pnow <- 0    pnew <- 0    if(!t %in% wm){          #observations      pnow  <- pnow + dnorm(y[t],xg[t],sqrt(tg),log=T)      pnew  <- pnew + dnorm(y[t],xprop,sqrt(tg),log=T)    }    if(t < nt){              #t+1 term excluded for last       pnow  <- pnow + dnorm(xg[t+1],fx(xg[t]),sqrt(sg),log=T)      pnew  <- pnew + dnorm(xg[t+1],fx(xprop),sqrt(sg),log=T)    }   if(t > 1){                #t-1 term excluded for 1st       pnow  <- pnow + dnorm(xg[t],fx(xg[t-1]),sqrt(sg),log=T)      pnew  <- pnew + dnorm(xprop,fx(xg[t-1]),sqrt(sg),log=T)   }   xg[t] <- acceptMH(pnow,pnew,xg[t],xprop,BLOCK=T)$x   }  xg}

updateSSparsMet <- function(b,priorB,priorVB,
                            lo=rep(-Inf,length(b)),hi=rep(Inf,length(b)),
                            x,sigma,dt=1,propVar){  k <- length(b)
  if(k == 1)pb <- tnorm(1,lo,hi,b,propVar)  if(k > 1) pb <- tnorm.mvt(b,b,propVar,lo,hi) 
  
  xnow <- xg[-nt] + fx(xg[-nt],b)*dt     #predicted x
  xnew <- xg[-nt] + fx(xg[-nt],pb)*dt
	
  pnow <- sum(dnorm(xg[-1],xnow,sqrt(sigma*dt),log=T)) +
          dmvnorm(t(b),priorB,priorVB,log=T)
  pnew <- sum(dnorm(xg[-1],xnew,sqrt(sigma*dt),log=T)) +
	      dmvnorm(t(pb),priorB,priorVB,log=T)
  atmp <- acceptMH(pnow,pnew,b,pb,BLOCK=T)
  b   <- atmp$x

  list(b = b, aa = atmp$accept)

}

initialStatesSS <- function(likelihood,y,priorObsEr,priorObsErWt,bias=1){

  library(stats)
  if(!is.matrix(y))y <- matrix(y,ncol=1)
  r <- ncol(y)

  n    <- length(y)
  time <- c(1:n)
  wm   <- which(is.na(y))
  notMiss <- c(1:n)
  
  x <- y*bias

  if(length(wm) > 0){
   notMiss <- notMiss[-wm]
   x[wm]   <- predict(smooth.spline(time[-wm],y[-wm]),time)$y[wm]
  }

  propSd <- sd(diff(y),na.rm=T)
  
  if(likelihood == 'dpois')x <- log(y + .1)
  
  s1s2 <- c(0,0)
  if(likelihood == 'dnorm'){
   	 s1 <- priorObsErWt
   	 s2 <- priorObsEr*(s1 - 1)
   	 s1s2 <- c(s1,s2)
  }

  list(x = x, propSd = propSd,notMiss = notMiss, miss = wm, s1s2 = s1s2)
}

updateSSstatesMet <- function(x,y,b,sigma,tau,lo=-Inf,hi=Inf,notMiss=c(1:length(x)),
                              dt=1,propSd=.1){   #propose/accept x as a block
  nt   <- length(x)  xp   <- tnorm(nt,lo,hi,x,rexp(nt,1/propSd) )       #proposed x
    xnew <- xp[-nt] + fx(xp[-nt],b)*dt      #mean proposed  xnow <- x[-nt]  + fx(x[-nt],b)*dt       #mean current  pnow <- sum(dnorm(x[-1],xnow,sqrt(sigma*dt),log=T)) +	       sum(dnorm(y[notMiss],x[notMiss],sqrt(tau),log=T))
	         pnew <- sum(dnorm(xp[-1],xnew,sqrt(sigma*dt),log=T)) +	       sum(dnorm(y[notMiss],xp[notMiss],sqrt(tau),log=T))
	       
  tmp <- acceptMH(pnow,pnew,x[-1],xp[-1],BLOCK=T)  x[-1]   <- tmp$x  list(x = x, aa = tmp$accept)}SSLGupdateS <- function(){    xp <- xg[-nt] + fx(xg[-nt],c(rg,kg))*dt    #mean current  u1 <- s1 + (nt-1)/2  u2 <- s2 + .5*sum( ((xg[-1] - xp)^2)/dt )  1/rgamma(1,u1,u2)}myrmvnorm <- function (n, mu = rep(0, nrow(sigma)), sigma = diag(length(mean))){
	
    ev <- eigen(sigma, sym = TRUE)$values
    if (!all(ev >= -sqrt(.Machine$double.eps) * abs(ev[1])))
        warning("sigma is numerically not positive definite")
    sigsvd <- svd(sigma)
    retval <- t(sigsvd$v %*% (t(sigsvd$u) * sqrt(sigsvd$d)))
    retval <- matrix(rnorm(n * ncol(sigma)), nrow = n) %*% retval
    if(min(dim(mu)) == 1 & n > 1)mu <- matrix(mu,n,ncol(sigma),byrow=T)
    if(nrow(retval) == 1 & ncol(mu) == 1)mu <- t(mu)
    retval + mu
}

mydmvnorm <- function(x,mean,sigma,log=FALSE){

  #mv normal density

    if (is.vector(x))x <- matrix(x, ncol = length(x))
    if (is.vector(mean))mean <- matrix(mean, ncol = length(x))

    distval <- mahalanobis(x, mean, sigma)
    logdet <- sum(log(eigen(sigma, symmetric = TRUE, only.values = TRUE)$values))
    logretval <- -(ncol(x) * log(2 * pi) + logdet + distval)/2
    if(log)return(logretval)
    exp(logretval)
}

ZIPparUpdate <- function(){

  cprop <- t(myrmvnorm(1,c(ag,bg),parcov))
  ap <- cprop[1:nv]
  bp <- cprop[(nv+1):(nv+nx)]

  if(DIST == 'Poisson'){
    pnow  <- -ziPoisNegLik(c(ag,bg))
    pnew  <- -ziPoisNegLik(c(ap,bp))
  }
  if(DIST == 'LN'){
    pnow  <- -ziLNNegLik(c(ag,bg,sg))
    pnew  <- -ziLNNegLik(c(ap,bp,sg))
  }

  pnow  <- pnow + mydmvnorm(t(ag),aprior,aVar,log=T) +
                  mydmvnorm(t(bg),bprior,bVar,log=T)
  pnew  <- pnew + mydmvnorm(t(ap),aprior,aVar,log=T) +
                  mydmvnorm(t(bp),bprior,bVar,log=T)
  a <- exp(pnew - pnow)
  z <- runif(1,0,1)
  if(z < a){
    ag <- ap
    bg <- bp
  }

  ss <- 0
  if(DIST == 'LN'){
    ss <- sg
    sp <- tnorm(1,0,10,sg,.005)
    pnow  <- -ziLNNegLik(c(ag,bg,sg))
    pnew  <- -ziLNNegLik(c(ag,bg,sp))
    a <- exp(pnew - pnow)
    z <- runif(1,0,1)
    if(z < a)ss <- sp
  }
     
  list(ag = ag, bg = bg,sg = ss)
}

ziPoisNegLik <- function(par){

  a <- matrix(par[1:nv],nv,1)
  b <- matrix(par[(nv+1):(nv+nx)],nx,1)

  VA   <- vb%*%a
  LB   <- xb%*%b

  theta <- inv.logit(VA)
  lamda <- exp(LB)

  pnow  <- rep(0,n)
  pnow[yb == 0] <- log(1 - theta + theta*dpois(0,A*lamda))[yb == 0]
  pnow[yb > 0]  <- log(theta[yb > 0]) + dpois(yb[yb > 0],(A*lamda)[yb > 0],log=T)
  -sum(pnow)
}
#########################################
ziLNNegLik <- function(par){

  a <- matrix(par[1:nv],nv,1)
  b <- matrix(par[(nv+1):(nv+nx)],nx,1)
  s <- par[length(par)]

  VA   <- vb%*%a
  LB   <- xb%*%b

  theta <- inv.logit(VA)

  pnow  <- rep(0,n)
  pnow[yb == 0] <- log(1 - theta)[yb == 0]
  pnow[yb > 0]  <- log(theta[yb > 0]) + dnorm(log(yb[yb > 0]),LB[yb > 0],s,log=T)
  -sum(pnow)
}

####################################################
bgsampMVN <- function(){  #sample reg pars for MVN 
	
  V   <- solve(crossprod(x)) 
  v   <- crossprod(x,y)
  mu  <- V%*%v
  vaa <- kronecker(sg,V)   
  
 # vaa <- nearPD(vaa)$mat
  matrix( myrmvnorm(1,as.vector(mu),vaa) ,p,ns,byrow=F)

}
######################################
bgTimeMVN <- function(){  #sample reg pars for MVN 
	
  yy  <- yg[-(tindex+1),] - yg[-(tindex+nt),]%*%ag
  xx  <- xtime[-(tindex+nt),]
  
  V   <- solve(crossprod(xx)) 
  v   <- crossprod(xx,yy)
  mu  <- V%*%v
  vaa <- kronecker(sg,V)   
  
  matrix( rmvnorm(1,as.vector(mu),vaa) ,p,ns,byrow=F)

}
######################################
randUpdateMVNorm <- function(x,y,sigma){  #sample random effects pars for MVN 
	  
  V   <- solve(crossprod(x))
  v   <- crossprod(x,y)
  mu  <- V%*%v
  vaa <- kronecker(sigma,V)   
  
  matrix( rmvnorm(1,as.vector(mu),vaa) ,ns,ns,byrow=F)

}
######################################
wishTimesamp <- function(){   #sample from Inv Wishart
	
	yy  <- yg[-(tindex+1),] 
   xx  <- xtime[-(tindex+nt),]%*%bg + yg[-(tindex+nt),]%*%ag

   scp  <- crossprod((yy - xx))
   vmat <- solve(scp + prior.W*prior.WDF)
   v2   <- ns + prior.WDF
   stmp <- myrmvnorm(v2,matrix(0,v2,ns),vmat)
   crossprod(stmp)

}

######################################
wishsamp <- function(x,yy,b){   #sample from Inv Wishart

   r    <- ncol(b)
   scp  <- crossprod((yy - x %*% b))
   vmat <- solve(scp + prior.W*prior.WDF)
   v2   <- ns + prior.WDF
   stmp <- myrmvnorm(v2,matrix(0,v2,r),vmat)
   crossprod(stmp)
}
####################################################
ysampMVNPois <- function(x,y,z,b,s,E){  #sample y's for MVN 1st stage, Pois 2nd stage

  r <- ncol(y)
  propy <- matrix(rnorm(length(y),y,.01),n,r,byrow=F)
  
  pnow <- rep(0,n)
  pnew <- rep(0,n)
  
  for(i in 1:n){
    pnow[i] <- mydmvnorm(y[i,],(x %*% b)[i,],s,log=T)
    pnew[i] <- mydmvnorm(propy[i,],(x %*% b)[i,],s,log=T)
  }
  
  pnow <- pnow + rowSums(dpois(z,E*exp(y),log=T))
  pnew <- pnew + rowSums(dpois(z,E*exp(propy),log=T))

  a <- exp(sum(pnew) - sum(pnow))
  zz <- runif(1,0,1)
  accept <- 0
  if(zz < a){
     y <- propy
     accept <- 1
  }
  list(y = y, a = accept)
  
}
#########################

sampleEffort <- function(){  #if Poisson effort is unknown
	
	u1 <- a1 + apply(z,1,sum)
	u2 <- a2 + apply(exp(y),1,sum)
	a <- rgamma(n,u1,u2)
	matrix(a,n,ns)
}

#########################
simMVNData <- function(n,p,ns,rangeB){
	
   covars <- c('intercept',paste('x',c(1:(p-1)),sep='-') )
   p      <- length(covars)
   x      <- matrix(runif(n*p,-1,1),n,p)
   x[,1]  <- 1
   colnames(x) <- covars
   
   b           <- matrix(runif(p*ns,rangeB[1],rangeB[2]),p,ns)
   specnames   <- paste('S',c(1:ns),sep='-')
   colnames(b) <- specnames
   rownames(b) <- covars
   sigma <- diag(runif(ns,0,.1),ns)
   rownames(sigma) <- specnames
   colnames(sigma) <- specnames

   mu    <- x %*% b
   y     <- myrmvnorm(n,mu,sigma)/10
   sigma <- cov(y)
   sigma <- as.matrix(nearPD(sigma)$mat)
   y     <- myrmvnorm(n,mu,sigma)
   colnames(y) <- specnames
   
   list(x = x, b = b, y = y, s = sigma, specnames = specnames)
}

logit2Prob <- function(y){   #multivar logit to fractions
	
  if(!is.matrix(y))y <- matrix(y,nrow=1)
  
  zs   <- apply(exp(y),1,sum)
  z1   <- 1/(1 + zs)
  zm   <- exp(y)/ (1 + zs)
  cbind(zm,z1)
  
}

prob2Logit <- function(y){     #fractions to multivar logit      
  
  log(y[,-r]/(1 - apply(y[,-r],1,sum)))
  
}

dataModel <- function(FUN,...){
	
	 FUN <- match.fun(FUN)
	 FUN(...)
}

################################
updateSSstatesMVN <- function(n,nt,tindex,x,y,yobs,b,a,obsError){
	
   ns <- ncol(y)
   mu <- xtime%*%b + yg%*%a
    
   accept <- 0
   
   osd <- sqrt(obsError) 
	
	for(t in 1:nt){
				
		ti <- t + tindex
		propy <- matrix(rnorm(n*ns,y[ti,],.1),n,ns)
		
		if(likelihood == 'dnorm'){
		  dnow <- dataModel(likelihood,y[ti,],yobs[ti,],osd,log=T)
		  dnew <- dataModel(likelihood,propy,yobs[ti,],osd,log=T)
		}
		if(likelihood == 'dpois'){
		  dnow <- dataModel(likelihood,yobs[ti,],exp(y[ti,]),log=T)
		  dnew <- dataModel(likelihood,yobs[ti,],exp(propy),log=T)
		}		
		
		pnow <- apply(dnow,1,sum)
		pnew <- apply(dnew,1,sum)
		
		if(t > 1){
			pnow <- pnow + diag(-(y[ti,] - mu[ti-1,])%*%sinv%*%t(y[ti,] - mu[ti-1,])/2)
			pnew <- pnew + diag(-(propy - mu[ti-1,])%*%sinv%*%t(propy - mu[ti-1,])/2)
		}	
		if(t < nt){
			mu1 <- mu[ti,]
			mu2 <- xtime[ti,]%*%b + propy%*%a
			pnow <- pnow + diag(-(y[ti+1,] - mu1)%*%sinv%*%t(y[ti+1,] - mu1)/2)
			pnew <- pnew + diag(-(y[ti+1,] - mu2)%*%sinv%*%t(y[ti+1,] - mu2)/2)
		}	
		aa <- exp(pnew - pnow)
		za <- runif(n,0,1)
		wp <- which(za < aa)
		if(length(wp) > 0){
			y[ti[wp],] <- propy[wp,]
			accept <- accept + length(wp)
		}
   }
   list(y = y, accept = accept/nrow(y))
}

obsErrorMVNlogit <- function(){
	
	u1 <- s1 + .5*n*ns*nt
	u2 <- s2 + .5*sum( (yg - yobs)^2 )
	1/rgamma(1,u1,u2)
}		
		
		
		
		

predYMVN <- function(x,bg,sg){ 
	#predict multivariate y

  vaa <- kronecker(diag(1,n),sg)
  matrix(myrmvnorm(1,as.vector(x%*%bg),vaa),n,ns)
}

#####################
distmat <- function(xt,yt,xs,ys){
    xd <- outer(xt,xs,function(xt,xs) (xt - xs)^2)
    yd <- outer(yt,ys,function(yt,ys) (yt - ys)^2)
    t(sqrt(xd + yd)) 
}
#####################
getSetup <- function(mapx,mapy,grid){   #setup map for movement

  nx   <- diff(mapx)/grid + 1
  ny   <- diff(mapy)/grid + 1
  xs <- seq(mapx[1],mapx[2],length=nx)
  ys <- seq(mapy[1],mapy[2],length=ny)

  pj   <- as.matrix(expand.grid(x=seq(0,100,length=nx),y=seq(0,100,length=ny)))
  dj   <- distmat(pj[,1],pj[,2],pj[,1],pj[,2])

  list(xs = xs, ys =  ys, pj = pj, dj = dj)
}
#################################
getNeighbors <- function(j){   #find neighbor cells for movement model

  dz     <- matrix(distj[j,],nrow=length(j),byrow=F)
  nb     <- t(apply(dz,1,order,decreasing=F))[,1:nm]
  dindex <- cbind(rep(j,each=nm),as.vector(t(nb)))
  distNb <- matrix(distj[dindex],nrow=length(j),byrow=T)

  list(nb = nb, di = distNb)
}
##########################

MoveUpdateB <- function(){  #update reg pars for movement model

   bp   <- t(myrmvnorm(1,t(bg),parcov))

   pnow <- 0
   pnew <- 0

   for(t in 1:(nt-1)){

     znext <- which(hood[t,] == z[t+1])
     xt    <- cbind(xvar[hood[t,]] - xvar[z[t]],distn[t,]) #subtract mean to get gradient
     e     <- exp(xt%*%bg)
     th    <- (e/sum(e))[znext]
     pnow  <- pnow + log(th)

     e     <- exp(xt%*%bp)
     th    <- (e/sum(e))[znext]
     pnew  <- pnew + log(th)

  }

  pnow <- pnow + mydmvnorm(t(bg),bprior,bVar,log=T)
  pnew <- pnew + mydmvnorm(t(bp),bprior,bVar,log=T)

  a <- exp(pnew - pnow)
  z <- runif(1,0,1)
  if(z < a)bg <- bp
  bg
}

getKernel    <- function(a)     exp(-(dij/a)^2)               #kernel for distance dij                
getKernSigma <- function(k,s,t) crossprod(t(k))*s + diag(t,n) #covariance


updateKernAS <- function(){  #update alpha and sigma
	
  asnow <- c(ag,sg)
  asnew <- tnorm.mvt(asnow,asnow,vas,c(0,0),c(amax,smax))
	
  know <- getKernel(ag)	
  knew <- getKernel(asnew[1])
  snow <- getKernSigma(know,sg,tg)
  snew <- getKernSigma(knew,asnew[2],tg)
	
  pnow <- dmvnorm(t(y),know%*%x%*%bg,snow,log=T)
  pnew <- dmvnorm(t(y),knew%*%x%*%bg,snew,log=T)
	
  a <- exp(pnew - pnow)
  z <- runif(1,0,1)
  if(z < a)asnow <- asnew
  asnow

}

updateKernB <- function(){
	
	s  <- getKernSigma(k,sg,tg)
	v1 <- t(x)%*%t(k)%*%s
	V  <- solve(v1%*%k%*%x)
	v  <- v1%*%y
	t(myrmvnorm(1,V%*%v,V))
}

updateKernT <- function(){
	
	e  <- matrix(rnorm(J,0,sqrt(sg)),J,1)
	rs <- y - k%*%(x%*%bg + e)
	u1 <- t1 + n/2
	u2 <- t2 + .5*crossprod(rs)
	1/rgamma(1,u1,u2)
}

cov2Dist <- function(sigma){ #distance induced by covariance
	
	n <- nrow(sigma)
	matrix(diag(sigma),n,n) + matrix(diag(sigma),n,n,byrow=T) - 2*sigma
}

####################

makeMap <- function(xlim,ylim,nx,ny,d,WEIGHT=F,ROWNORMAL=F){  
	#adjacency matrix W
	#rowsum matrix Dw
	#expanded map grid
	
  mapx <- seq(xlim[1],xlim[2],length=nx)
  mapy <- seq(ylim[1],ylim[2],length=ny)
  n    <- nx*ny

  mapgrid <- expand.grid(mapx,mapy)
  plot(mapgrid[,1],mapgrid[,2])

  mIndex <- matrix(c(1:n),nx,ny,byrow=T)
  
  dd <- as.matrix(dist(mapgrid,diag=T))
  diag(dd) <- d + 1
  
  wh <- which(dd <= d,arr.ind=T)
    
  hood <- table(wh[,1],by=wh[,2])
  if(WEIGHT)hood <- hood*1/dd
  if(ROWNORMAL)hood <- hood/matrix(apply(hood,1,sum),nrow(hood),nrow(hood))
  
  W   <- hood
  Dw  <- diag(rowSums(W))

  tmp   <- eigen(Dw)
  dsi   <- solve(tmp$vectors%*%diag(sqrt(tmp$values))%*%solve(tmp$vectors))
  rho   <- range(eigen(dsi%*%W%*%dsi)$values)

  list(mapgrid = mapgrid, W = W, mIndex = mIndex, Dw = Dw, rho = rho)
}

updateBetaCov <- function(x,y,sigInv,priorB,priorIVB){

  x1 <- t(x)%*% sigInv
  V  <- chol2inv(chol(x1%*%x + priorIVB))
  v <- x1%*%y + priorIVB%*%priorB
  t(rmvnorm(1,V%*%v,V))
}

updateTauCAR <- function(mu,y,rho,tau,rlims,tlims,propVar){

  n  <- length(mu)
  lo <- c(rlims[1],tlims[1])
  hi <- c(rlims[2],tlims[2])
  prop <- tnorm.mvt(c(rho,tau),c(rho,tau),propVar,lo,hi)

  si1 <- (D - rho*W)/tau
  si2 <- (D - prop[1]*W)/prop[2]

  ld1 <- -sum(log(eigen(si1, symmetric = T, 
          only.values = T)$values))
  ld2 <- -sum(log(eigen(si2, symmetric = T, 
          only.values = T)$values))

  snow <- t(mu)%*%si1%*%mu
  snew <- t(mu)%*%si2%*%mu

  pnow <- -(n*log(2*pi) + ld1 + snow)/2
  pnew <- -(n*log(2*pi) + ld2 + snew)/2

  tmp <- acceptMH(pnow,pnew,c(rho,tau),prop,BLOCK=T)
  r <- tmp$x[1]
  t <- tmp$x[2]

  c(r,t)
}


