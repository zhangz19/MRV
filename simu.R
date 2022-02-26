

# MRV: median regression with varying coefficients
# simulation program
# Zhen Zhang (zhangquake1@outook.com), 2022-2-16
# note: run the code with Rscript simu 0 first to generate data 1-nsim, then run simu 1-nsim for job array

require(splines)
require(R.matlab)

repl <- as.numeric(commandArgs(TRUE))   #repl is the job ID of job-array runs
# repl <- 0   # test run

nsim <- 500  #number of simulated data sets
design <- 2  #1=setting in Yin et al(2014), 2=simulation designs in the manuscript
simuDataOnly <- ifelse(repl==0, TRUE, FALSE)

# set up simulation scenarios
# scenarios <- do.call(expand.grid, list(errType=1:3, cen_pt=c(30,50), semi_varying=c(TRUE, FALSE), n=1e3))
scenarios <- data.frame(errType=1, cen_pt=c(30), semi_varying=c(TRUE), n=1e3)
dir0 <- getwd()

for(runId in 1:nrow(scenarios)){
  print(scenarios[runId,])
  errType <- scenarios$errType[runId] #1-3 error types in the manuscript
  cen_pt <- scenarios$cen_pt[runId]  #censoring percentage
  semi_varying <- scenarios$semi_varying[runId]  
  #if semi_varying==TRUE: Z_g for gamma (constant effect);  if FALSE: Z=0 (no gamma term, all time-varying [beta])
  n <- scenarios$n[runId]  # number of observations
  
  foldNam <- paste0('./run', runId)
  if(!dir.exists(foldNam)) dir.create(foldNam)
  setwd(foldNam)  #temporarily go to the subfolder
  
  # M <- 2
  # nboot <- 0 #400
  
  ptm <- proc.time()[3]
  
  starts <- 1; ends <- nsim
  if(!simuDataOnly) starts <- ends <- repl
  
  for(id in starts:ends){
    
    set.seed(id*10)
    cat(id,' '); if(!id%%20) cat('\n') 
    
    nknot <- 20 #number of knots for both methods
    cen_depend <- FALSE  #if TRUE: censoring depends on error and covariates
    
    if(design==1){
      # Yin et al(2014): simu 1
      #n <- 200  
      W <- runif(n,-1,1)
      Z <- cbind(1, runif(n), rnorm(n))
      p <- ncol(betas <- cbind(0.5, W^2, cos(3*W))) 
      eps <- rnorm(n,0,.1)
      mu <- matrix(NA,n,3); for(j in 1:p) mu[,j] <- Z[,j]*betas[,j]
      Ts <- exp(rowSums(mu) + eps)
      K <- 14; C <- runif(n,0,K)
      V <- apply(rbind(Ts,C),2,min)  #Ts; #
      Delta <- as.integer(Ts<=C)
      cat('\ncensoring rate = ', mean(1-Delta),'\n', sep='')
      myknots <- c(-.85,-.45,-.05, .45, .85)
      B <- bs(c(W,-0.5,0,0.5), knots = myknots, Boundary.knots=c(-1.2,1.2), degree = 3, intercept = TRUE)
    }
    
    if(design==2){
      gammas <- ifelse(semi_varying, 1, 0)
      
      # CAUTION: set seed inside datGen
      datGen <- function(K){
        set.seed(id*10)
        W <- runif(n,0,1)
        #Z <- cbind(1, rnorm(n), runif(n)); 
        Z <- cbind(1, rnorm(n)) 
        Z_g <- runif(n) #as.numeric(rbinom(n,1,.5))
        p <- ncol(betas <- cbind(W^3, sin(2*pi*W)))   #, 1
        gammas <- ifelse(semi_varying, 1, 0)
        
        cen_depend <- FALSE  #if TRUE: censoring depends on error and covariates
        
        if(errType==1){
          eps <- runif(n,-1,1)
          # if(cen_pt==30) K <- 6  #30% censoring
          # if(cen_pt==50) K <- 3  #50% censoring
        }
        
        if(errType==2){
          p0 <- 0.6; mu0 <- -0.25; lambda0 <- 0.5; sig0 <- -mu0/qnorm(0.5/p0)
          #p0*(mu0^2+sig0^2)+(1-p0)*2*lambda0^2-(p0*mu0 + (1-p0)*lambda0)^2
          tmp <- rbinom(n,1,p0)
          eps <- rnorm(n,mu0,sig0)*tmp + rexp(n,rate=1/lambda0)*(1-tmp)
          #var(eps)
          #hist(eps)
        }
        
        if(errType==3){
          p0 <- 0.6; tmp <- rbinom(n,1,p0)
          eps <- runif(n,-1,1)*tmp + rt(n,6)*(1-tmp)
          cen_depend <- TRUE
          #var(eps)
          #hist(eps) 
        }
        
        mu <- matrix(NA,n,p); for(j in 1:p) mu[,j] <- Z[,j]*betas[,j]
        Ts <- exp(rowSums(mu) + Z_g*gammas + eps)
        if(!cen_depend) 
          C <- runif(n,0,K)
        if(cen_depend){
          # K <- 3.5; if(cen_pt==50) K <- 1.7  
          C <- runif(n,0,K*(0.5*Z[,2]^2+exp(W+Z_g)))
        }
        V <- apply(rbind(Ts,C),2,min) #Ts; #
        #V <- Ts    # no censoring
        Delta <- as.integer(Ts<=C)
        # cat('\ncensoring rate = ', mean(1-Delta),'\n', sep='')
        return(list(censorRate=mean(1-Delta), V=V, Delta=Delta, W=W, Z=Z, Z_g=Z_g, betas=betas, mu=mu, eps=eps))
      }
      
      K <- uniroot(function(K) return(datGen(K)$censorRate-cen_pt/100), interval=c(0,100))$root
      # K
      tmp <- datGen(K)
      V <- tmp$V; Delta <- tmp$Delta; W <- tmp$W; Z <- tmp$Z; Z_g <- tmp$Z_g;  betas <- tmp$betas; mu=tmp$mu; eps=tmp$eps; 
      
      
      # if(nknot==9) myknots <- seq(.1,.9,by=.1) #9-knot design.  c(.2,.4,.6,.8)#c(.1,.3,.5,.7,.9)# #c(.2,.4,.6,.8) 
      # if(nknot==19) myknots <- seq(.05,.95,by=.05)  #19-knot design
      myknots <- seq(1/nknot, 1-1/nknot, length=nknot)
      len0 <- length(u0 <- c(.1,.3,.5,.7,.9)) #specific locations to calculate the bias and RMSE
      B <- bs(c(W,u0), knots = myknots,  degree = 3, intercept = TRUE, Boundary.knots=c(0,1)) #ncol(B)=nknot+3(cubic)+1(intercept)
      
      # # for plotting functional curves
      # B1 <- bs(u1 <- seq(.01,.99,len=50), knots = myknots,  degree = 3, intercept = TRUE, Boundary.knots=c(0,1))
      # if(simuDataOnly) writeMat('B1.mat',B1=B1, u1=u1)
      
    }
    
    if(simuDataOnly){ # save results for MATLAB codes
      
      #which(rowSums(B)==0)
      B0 <- B[(length(W)+1):nrow(B),]
      B <- B[1:length(W), ]
      labs <- rep(1,n)
      
      # # this part is now done in MATLAB without saving the file
      #    if(simuDataOnly){ # save bootstrap samples
      #      nboot2 <- 500
      #      indmat <- matrix(0,n,nboot2)
      #      for(b in 1:nboot2){
      #        if(!b%%100) print(b)
      #        inds <- sample(n,n,replace=T)
      #        tab <- table(factor(inds,levels=1:n))
      #        indmat[,b] <- tab
      #      }
      #      writeMat(paste('indmat_simu',id,'.mat',sep=''),indmat=indmat)
      #    }
      
      # run with semi varying:
      if(semi_varying) writeMat(paste0('sdat',id,'.mat'), V=V, Delta=Delta, X=Z, Z=Z_g, U=W, B=B, 
                                Trues=list(betas=betas, labs=labs, gammas=gammas, mu=mu, eps=eps, B0=B0))
      
      # run with all varying:
      #if(!semi_varying) writeMat(paste('simudat',id,'.mat',sep=''), V=V, Delta=Delta, X=cbind(Z, Z_g), Z=0, U=W, B=B, Trues=list(betas=cbind(betas,gammas), labs=labs, gammas=0, mu=mu, eps=eps, B0=B0))
      
      # totally discard Z
      if(!semi_varying) writeMat(paste0('sdat',id,'.mat'), V=V, Delta=Delta, X=Z, Z=0, U=W, B=B, 
                                 Trues=list(betas=betas, labs=labs, gammas=0, mu=mu, eps=eps, B0=B0))
      
    }
    
    # ==================================  run Yin et al (2014)'s method
    if(!simuDataOnly){
      
      p <- ncol(betas <- cbind(W^3, sin(2*pi*W)))  #, 1
      #Z <- cbind(Z, Z_g)
      
      #len <- length(ws <- seq(-.9, 1, by=0.1)) # 20 equal intervals: knots to evaluate for the kernel method
      # len <- length(ws <- seq(.1,.9,by=.1))
      len <- length(ws <- myknots)
      #hn <- 0.08 #0.06 #W1 <- density(W, bw=hn, kernel='gaussian')
      nbw <- length(hns <- c(.04, .07, .1, .13))  #difference choices of bandwidth for the kernel method
      
      tmp <- readMat(paste0('sdat',id,'.mat'))  #read the generated data set from previous step
      W <- tmp$U;  V <- tmp$V;  Delta <- tmp$Delta;  Z <- tmp$X
      dat0 <- dat <- as.data.frame(cbind(W, V, Delta, Z)) 
      if(semi_varying) p <- ncol(dat <- cbind(dat, Z_g))-3  #first 3 columns are not covariates. 
      
      r <- 1; P <- (r+1)*p
      #xhat <- array(NA, dim=c(P, len, 1+nboot))
      xhat <- array(NA, dim=c(P, len, nbw))
      
      for(loop_hn in 1:nbw){
        cat(loop_hn,' ')
        hn <- hns[loop_hn]
        # source('kernelMethod.R')
        
        kernelMethod <- function(dat){
          # The kernel methods in Yin et al (2014)'s paper
          # dat is as.data.frame(cbind(W, V, Delta, Z)), each variable is n-vector
          
          len <- length(ws)
          ncind <- which(dat$Delta==1) 
          mat <- matrix(NA, P, len)
          for(k in 1:length(ws)){
            w0 <- ws[k]
            # ------------- initial value x0: this part depends on simulation setting (design, semi_varying) 
            if(design==1)
              x0 <- c(0.5, w0^2, cos(3*w0), 0, 2*w0, -3*sin(3*w0))
            if(design==2&&semi_varying)
              x0 <- c(w0^3, sin(2*pi*w0), 1, 3*w0^2, 2*pi*cos(2*pi*w0), 0)
            if(design==2&&!semi_varying)
              x0 <- c(w0^3, sin(2*pi*w0), 3*w0^2, 2*pi*cos(2*pi*w0))
            
            Qfun <- function(x, w0){
              Z1 <- as.matrix(dat[,-c(1:3)])
              Z1 <- cbind(Z1, Z1*matrix(rep(dat$W-w0, p), ncol=p))
              e <- log(dat$V) - Z1%*%x
              u0 <- dnorm((dat$W-w0)/hn)
              v0 <- matrix(rep(e,n),ncol=n); v0 <- (v0>=t(v0))
              u <- matrix(rep(log(1-0.5), n*P), ncol=P) 
              for(j in 1:P){
                v <- u0*Z1[,j]
                u[ncind,j] <- u[ncind,j] + n*v[ncind]*(e[ncind]<=0)/colSums( matrix(rep(v,n),ncol=n)*v0 )[ncind]         
              }
              Un <- apply(u,2,mean)
              Qn <- matrix(0,P,P); for(i in 1:n) Qn <- Qn + u[i,]%*%t(u[i,])
              Qn <- Qn/n - Un%*%t(Un)
              Qn <- sum(Un * (chol2inv(chol(Qn))%*%Un))
              return(Qn)
            }
            
            fit <- optim(x0*1, Qfun, NULL, w0=w0, method='Nelder-Mead',control=list(trace=0,reltol=1e-16, maxit=1e3))
            mat[,k] <- fit$par
          }
          return(mat)
        }
        
        xhat[,,loop_hn] <- kernelMethod(dat)    #:len
        
        #      # bootstrap for variance estimation
        #      if(nboot>0) cat('start bootstrap samples...\n')
        #      for(nloop in 1:nboot){
        #        cat(nloop,' '); if(!nloop%%20) cat('\n') 
        #        
        #        dat <- dat0[sample(n,n,replace=T),]
        #        mat <- try(kernelMethod(dat, c(5,10,15) ))
        #        if(!inherits(mat, "try-error")) xhat[,,1+nloop] <- mat
        #        
        #      }
      }
      cat('\n')
    }# ============================ Yin et al (2014)'s method ends
    
  }
  
  
  if(!simuDataOnly){
    cputime <- as.numeric(proc.time()[3]-ptm)
    cputime <- cputime/60
    cat('CPUtime', round(cputime,3), 'minutes: completed!','\n')
    
    save(file=paste0('out',repl), xhat)
  }
  
  
  # return the starting folder
  setwd(dir0)
  
}

#++++++++++++++++++ util
doCombine <- FALSE  #combine the outputs from job arrays
if(doCombine){
  # make sure the following is consistent with the setting used for job arrays
  nknot <- 20 #number of knots for both methods
  myknots <- seq(1/nknot, 1-1/nknot, length=nknot)
  design <- 2
  semi_varying <- TRUE 
  nsim <- 500;
  
  #len <- length(ws <- seq(-.95, .95, by=0.1)) # 20 equal intervals
  # len <- length(ws <- seq(.1,.9,by=.1))  # 9 knots
  len <- length(ws <- myknots) 
  if(design==1)
    x0 <- t(cbind(0.5, ws^2, cos(3*ws), 0, 2*ws, -3*sin(3*ws)))
  if(design==2&&!semi_varying)
    x0 <- t(cbind(ws^3, sin(2*pi*ws), 3*ws^2, 2*pi*cos(2*pi*ws)))
  if(design==2&&semi_varying)
    x0 <- t(cbind(ws^3, sin(2*pi*ws), 1, 3*ws^2, 2*pi*cos(2*pi*ws), 0))
  
  runId <- 1
  foldNam <- paste0('./run', runId)
  ncount <- 0 
  for(repl in 1:nsim){
    fnam <- paste0(foldNam, '/out',repl)
    if(file.exists(fnam)){
      load(fnam)
      ncount <- ncount+1
      if(ncount==1){
        P <- dim(xhat)[1]; len <- dim(xhat)[2]; nbw <- dim(xhat)[3]
        #Bs <- SE <- CP <- array(NA, dim=c(P, len, nsim))
        Xs <- Bias <- MSEs <- array(NA, dim=c(P, len, nbw, nsim))
        gamma_term <- array(NA, dim=c(2, nbw, nsim))
      }
      Xs[,,,repl] <- xhat
      for(k in 1:nbw){
        Bias[,,k,repl] <- xhat[,,k] - x0
        MSEs[,,k,repl] <- (xhat[,,k] - x0)^2
        gammahat <- mean(xhat[3,,k],na.rm=T)
        gamma_term[,k,repl] <- c((gammahat-1), (gammahat-1)^2) 
      }
      #SE[,,repl] <- apply(xhat[,,-1], c(1,2), sd, na.rm=T)
      #CP[,,repl] <- (x0 >= Bs[,,repl]-1.96*SE[,,repl])*(x0 <= Bs[,,repl]+1.96*SE[,,repl])
    }
  }
  
  # for plotting
  ws0 <- seq(.1,.9,by=.1)
  len <- length(ws <- c(.1,.3,.5,.7, .9))
  inds <- c(1,3,5,7,9)
  ISE <- apply(MSEs[1:2,inds,1:3,], c(1,3,4), sum)
  findopt <- function(x) return(which(x==min(x))[1])
  optind <- matrix(NA, dim(ISE)[1],dim(ISE)[3])
  for(j in 1:nrow(optind)) optind[j,] <- apply(ISE[j,,],2,findopt)
  Xs1 <- array(NA,dim=c(dim(Xs)[2],2,dim(Xs)[4]))
  for(j in 1:nrow(optind))
    for(k in 1:dim(Xs)[4])
      Xs1[,j,k] <- Xs[j,,optind[j,k],k] 
  
  Xsm <- apply(Xs1, c(1,2), mean, na.rm=T)
  Xsl <- apply(Xs1, c(1,2), quantile,.025, na.rm=T)
  Xsu <- apply(Xs1, c(1,2), quantile,.975, na.rm=T)
  
  #  Xsm <- apply(Xs, c(1,2,3), mean)[1:3,,]
  #  Xsl <- apply(Xs, c(1,2,3), quantile,.025)[1:3,,]
  #  Xsu <- apply(Xs, c(1,2,3), quantile,.975)[1:3,,]
  p <- nrow(MSEs)/2
  MSEs <- apply(MSEs, c(1,2,3), mean, na.rm=T)[1:p,,]
  Bias <- apply(Bias, c(1,2,3), mean, na.rm=T)[1:p,,]
  gamma_term <- apply(gamma_term, c(1,2), mean, na.rm=T)
  
  #save(file=paste('tab2_s3',sep=''),Bias,MSEs,gamma_term,Xsm,Xsl,Xsu)
  # save(file=paste0('tab_s',errType,'_c',cen_pt,'_semi_',semi_varying), Bias,MSEs,gamma_term,Xsm,Xsl,Xsu)
  save(file=paste0('tab_',runId), Bias,MSEs,gamma_term,Xsm,Xsl,Xsu)
  
  #design <- 1
  #  load(paste('cout_s',design,'_y',sep=''))
  #  len <- length(ws <- seq(-.95, .95, by=0.1)) # 20 equal intervals
  #  if(design==1)
  #     x0 <- t(cbind(0.5, ws^2, cos(3*ws), 0, 2*ws, -3*sin(3*ws)))
  #  if(design==2)
  #     x0 <- t(cbind(ws^3, sin(2*pi*ws), 3*ws^2, 2*pi*cos(2*pi*ws)))
  #  Bhat <- apply(Bs, c(1,2), mean, na.rm=T)
  #  Bl <- apply(Bs, c(1,2), quantile, .025, na.rm=T)
  #  Bu <- apply(Bs, c(1,2), quantile, .975, na.rm=T)
  #  P <- dim(Bhat)[1]
  #  par(mfrow=c(2,P/2), mar=c(2,2,1,0)+.3,mgp=c(1.1,.2,0), tck=-0.01, cex.axis=.9, cex.lab=1.1, font.lab=1, cex.main=1)
  #  for(j in 1:P) {matplot(ws,cbind(Bhat[j,],Bl[j,],Bu[j,]), type='l', pch=16, lty=c(1,2,2),lwd=c(2,2,2),col=c('black','blue3','green3')) 
  #    lines(ws, x0[j,], col='red')
  #  }
}



getTable <- FALSE
if(getTable){
  
  library(R.matlab)
  getFigure <- FALSE
  
  # create the table for the manuscript
  round1 <- function(x,d=3){
    y <- as.character(round(x,d))
    for(i in 1:length(y)){
      k <- as.integer(gregexpr('.',y[i],fixed=T))
      if(k<0) {y[i] <- paste(y[i],'.',sep=''); for(j in 1:d) y[i] <- paste(y[i],'0',sep='')}
      if(k>=0){ b <- nchar(y[i])-as.integer(gregexpr('.',y[i],fixed=T))
      if(b<d) for(j in 1:(d-b))y[i] <- paste(y[i],'0',sep='')
      }
    }
    return(y)
  }
  
  ws0 <- seq(.1,.9,by=.1)
  len <- length(ws <- c(.1,.3,.5,.7, .9))
  inds <- c(1,3,5,7,9)
  rbind(ws0[inds], ws)
  hns <- c(.04, .07, .1)#, .13)
  
  
  if(getFigure){
    # postscript(file=paste('./tex/S.eps',sep=''), pointsize=16,width=8,height=7.9,horizontal=F,paper='special') 
    postscript(file=paste0('./tex/S.eps'), pointsize=16,width=8,height=7.9,horizontal=F,paper='special')
    png(paste0('./tex/S.png'), units='in', family='Times',
        height=7.9, width=8, pointsize=16, res=600) #, bg='transparent'
  }
  
  # par(mfrow=c(3,4), mar=c(2.5,1,1.5,0)+.3,mgp=c(1.1,.2,0), tck=-0.01, cex.axis=.9, cex.lab=1, font.lab=3, cex.main=1.1,font.main=3)
  par(mfrow=c(1,1), mar=c(2.5,1,1.5,0)+.3,mgp=c(1.1,.2,0), tck=-0.01, cex.axis=.9, cex.lab=1, font.lab=3, cex.main=1.1,font.main=3)
  U0 <- readMat('B1.mat')$u1
  
  for(design in 1:3){
    
    cat('\nscenario =', design,':\n')
    
    # obj <- readMat(paste0('tab1_s',design,'.mat'))
    # tab1 <- obj$tab
    # tab1 <- apply(tab1,2,round1); tab1 <- tab1[,-c(3,6)] 
    # 
    # obj2 <- readMat(paste0('tab1_s',design,'_c50.mat'))
    # tab2 <- obj2$tab
    # tab2 <- apply(tab2,2,round1); tab2 <- tab2[,-c(3,6)] 
    # tab1 <- cbind(tab1, tab2)
    # load(paste0('tab2_s',design))
    # Bias1 <-  Bias; MSEs1 <- MSEs; gamma_term1 <- gamma_term
    # Xsm1 <- Xsm; Xsl1 <- Xsl; Xsu1 <- Xsu;
    # load(paste0('tab2_s',design,'_c50'))
    runId <- 1
    load(paste0('tab_',runId))
    ns <- dim(Bias) # p, lw, nbw
    p <- ncol(obj$Fm)
    
    # table for beta
    for(i in 1:len){
      cat(  '\\multirow{4}{*}{',ws[i], '}', sep='')
      for(j in 1:length(hns)) 
      { if(j==2) cat('&$1$') else cat('&')
        cat('&$',round1(hns[j],2),'$', sep='')
        for(k in 1:p)
        {cat(paste('&',round1(Bias1[k,inds[i],j]),'&',round1(sqrt(MSEs1[k,inds[i],j])),sep='$'));cat('$') #if(k<ns[1]) cat('&') 
        }
        for(k in 1:p)
        {cat(paste('&',round1(Bias[k,inds[i],j]),'&',round1(sqrt(MSEs[k,inds[i],j])),sep='$'));cat('$') #if(k<ns[1]) cat('&') 
        }
        cat('\\\\\n')
      }
      cat('&$2$&'); cat(paste('&$',tab1[i,],'$',sep=''),'\\\\','\n',sep='')
      if(i<len) cat('&  \\\\\n') else cat('\\hline\n')
    }
    
    # table for gamma
    getTableGamma <- function(){
      for(j in 1:length(hns)) 
      { if(j==2) cat('$1$')
        cat('&$',round1(hns[j],2),'$',sep='')
        cat(paste('&',round1(gamma_term1[1,j]),'&',round1(sqrt(gamma_term1[2,j])),sep='$'));cat('$') #if(k<ns[1]) cat('&')
        cat(paste('&',round1(gamma_term[1,j]),'&',round1(sqrt(gamma_term[2,j])),sep='$'));cat('$') 
        cat('\\\\\n')
      }
      cat( '$2$&',
           paste('&$',round1(c(obj$gamma.term[1], sqrt(obj$gamma.term[2]))),'$',sep=''),
           paste('&$',round1(c(obj2$gamma.term[1], sqrt(obj2$gamma.term[2]))),'$',sep=''),  '\\\\','\n',sep='')
    }
    cat('\ntable for gamma:\n')
    getTableGamma()
    cat('\n\n')
    
    ## make plot
    # for(mycen in c(30, 50)){
    for(mycen in c(30)){
      for(j in 1:p){   
        if(mycen==30)mat <- cbind(obj$Fm[,j],obj$Fl[,j],obj$Fu[,j])
        if(mycen==50)mat <- cbind(obj2$Fm[,j],obj2$Fl[,j],obj2$Fu[,j])
        matplot(U0, mat, type='l', pch=16, lty=c(1,2,2),lwd=c(2,2,2),col=c('blue3','blue3','blue3'), 
                xlab=paste('Scenario ',design,', ', mycen,'% censoring', sep=''),#'Age',
                ylab='' #,express(paste('beta[',j, ']', sep=''))
        )
        for(k in 1:1){  #length(hns)
          lines(ws0, Xsm1[,j], col=paste('gray',k*27,sep=''), lty=3,lwd=2)
          lines(ws0, Xsl1[,j], col=paste('gray',k*27,sep=''), lty=3,lwd=2)
          lines(ws0, Xsu1[,j], col=paste('gray',k*27,sep=''), lty=3,lwd=2)
        }
        if(j==1) curve(x^3,0,1, col='red2',add=T,lwd=2,lty=2)
        if(j==2) {chippy <- function(x) sin(2*pi*x); curve(chippy,0,1, col='red2',lwd=2,lty=2,add=T) }
        #, col=c('#1a1a1a','#d7301f','#053061')
        #abline(h=0,lty=2,col='red2')
        if(j==1) title(main=substitute(paste(beta[1],a^3), list(a='(U) = U')))
        if(j==2) title(main=substitute(paste(beta[2],a1,pi,a2), list(a1='(U) = sin(2',a2='U)')))
      }
    }
    
    
  }
  if(getFigure) dev.off()
  
}

print('done')

# not run

