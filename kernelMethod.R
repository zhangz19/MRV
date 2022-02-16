
kernelMethod <- function(dat, inds, ws=seq(.1,.9,by=.1)){
  # The kernel methods in Yin et al (2014)'s paper
  # dat is as.data.frame(cbind(W, V, Delta, Z)), each variable is n-vector
  # inds is typically 1:length(ws) where ws includes knots in (0,1)
  
  len <- length(ws)
  ncind <- which(dat$Delta==1) 
  mat <- matrix(NA, P, len)
  for(k in inds){
    w0 <- ws[k]
    if(design==1)
     x0 <- c(0.5, w0^2, cos(3*w0), 0, 2*w0, -3*sin(3*w0))
    if(design==2)
     x0 <- c(w0^3, sin(2*pi*w0), 1, 3*w0^2, 2*pi*cos(2*pi*w0), 0)
    
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

# not run 