mixtwice <-
  function(thetaHat, s2, Btheta = 15, Bsigma2 = 10, df, prop = 1){
    
    ## some initial check
    
    if(length(thetaHat)==0){
      stop("the input is empty")
    }
    
    if(!length(thetaHat) == length(s2)){
      stop("the length of estimated effect sizes must be equal to the estimated squared standard errors")
    }
    
    if(sum(s2<0)>0){
      stop("estimated squared standard errors can not be negative")
    }
    
    ## give another notation for convenience
    theta0=thetaHat
    s20=s2
    
    ## let's first sample using "prop"
    ok.sample=sample(length(theta0), length(s20)*prop)
    
    theta=theta0[ok.sample]
    s2=s20[ok.sample]
    p0=length(theta0)
    p=length(theta)
    
    cc=max( abs(theta) )*1.1
    
    gridtheta=(cc/Btheta)*seq(-Btheta,Btheta,by=1)
    
    gridsigma=seq(sqrt(min(s2)), sqrt(max(s2)), by=(sqrt(max(s2))-sqrt(min(s2)))/Bsigma2)
    
    ltheta=length(gridtheta)
    lsigma=length(gridsigma)
    
    grid1=rep(gridtheta,each=length(gridsigma))
    grid2=rep(gridsigma, length(gridtheta))
    
    rbind(grid1, grid2)
    
    lik1=t(exp(-.5*( t((outer(theta,grid1,"-"))^2 )/(grid2)^2)) / (grid2*sqrt(2*pi)))
    
    ### Then, what about the likelihood of, P(S2=s2|theta_k, sigma_j)
    
    ### It is a chi-square distribution
    
    ### for simplicity, let's denote y=df*s^2/sigma^2, it is a matrix
    
    y=outer(df*s2, (1/gridsigma^2), "*")
    
    m=df/2
    lik2=y^(m -1) * exp(-0.5*y) /((2^m)*(gamma(m)))
    lik22=matrix(rep(lik2, ltheta), nrow = nrow(lik2))
    
    
    ### Then the overall likelihood is the product of lik1 and lik2
    
    lik=lik1*lik22
    
    ### Then, we write down the objective function and its gradient
    L=function(x){
      
      ## x is a combination of pointmass of theta and sigma
      ## the first entries are pointmasses of beta and rest sigma
      
      xtheta=x[1:ltheta]
      
      xsigma=x[(ltheta+1):(ltheta+lsigma)]
      
      ## I pull out as a array, the outer product of xsigma and xtheta
      
      yy=array(x[(ltheta+1):(ltheta+lsigma)] %o% x[1:ltheta])
      
      return(-sum(log(yy %*% t(lik))))
      
      
    }
    
    G=function(x){
      
      g=h=NULL
      
      xtheta=x[1:ltheta]
      
      xsigma=x[(ltheta+1):(ltheta+lsigma)]
      
      yy=array(x[(ltheta+1):(ltheta+lsigma)] %o% x[1:ltheta])
      
      ## I consider the denormator
      
      d=yy %*% t(lik)
      
      ## I then consider the numerator
      
      # first the gradient over the component of beta
      
      for (i in (1:ltheta)) {
        g[i]=-sum((xsigma %*% t(lik[,c((lsigma*(i-1)+1):(lsigma*i))]))/d)
      }
      
      # second the gradient over the component of sigma
      
      for (j in (1:lsigma)) {
        
        h[j]=-sum((xtheta %*% t(lik[,seq(j, (j+(ltheta-1)*(lsigma)), by=lsigma)]))/d)
      }
      
      return(c(g,h))
      
      
    }
    
    
    ### let's think about the constraint
    
    ## first, equality constraint
    
    heq <- function(x){
      
      h=NULL
      ## we gonna need the summation of signal and sigma to be 1
      
      h[1]=sum(x[1:ltheta])-1
      
      h[2]=sum(x[(ltheta+1):(lsigma+ltheta)])-1
      
      return(h)
    }
    
    hh1=c(rep(1, ltheta), rep(0, lsigma))
    hh2=c(rep(0, ltheta), rep(1, lsigma))
    
    heq.jac=rbind(hh1, hh2)
    
    heq.jac.fun=function(x){
      j=heq.jac
      return(j)
    }
    
    ## Then, inequality constaint
    ## let's first try the simplest case, only probability mass property
    hin <- function(x){
      
      h1=NULL
      
      for (i in 1:((ltheta)+(lsigma))) {
        
        h1[i]=x[i]
        
      }
      
      h2=NULL
      
      for (i in 1:(Btheta)) {
        
        h2[i]=x[i+1]-x[i]
        
      }
      
      for (i in (Btheta+1):(ltheta-1)) {
        
        h2[i]=x[i]-x[i+1]
      }
      
      h=c(h1, h2)
      return(h)
    }
    
    hin.jac1=diag(1, nrow = (ltheta+lsigma), ncol = (ltheta+lsigma))
    hin.jac2=matrix(0, ncol = ltheta, nrow = ltheta-1)
    for (i in 1:(Btheta)) {
      
      hin.jac2[i, i]=-1
      hin.jac2[i, i+1]=1
      
    }
    
    for (i in (Btheta+1):(ltheta-1)) {
      hin.jac2[i, i]=1
      hin.jac2[i, i+1]=-1
    }
    hin.jac3=matrix(0, nrow=ltheta-1, ncol=lsigma)
    hin.jac=rbind(hin.jac1, cbind(hin.jac2, hin.jac3))
    
    hin.jac.fun=function(x){
      j=hin.jac
      return(j)
    }
    ## let's try optimization
    
    ## the initial
    
    a1=rep(1,ltheta);a1=a1/sum(a1)
    
    a2=rep(1, lsigma); a2=a2/sum(a2)
    
    a=c(a1, a2)
    
    try1=suppressWarnings(alabama::auglag(par=a, fn=L,  gr=G, 
                                          heq=heq, hin=hin, 
                                          heq.jac = heq.jac.fun, hin.jac = hin.jac.fun,
                                          control.outer = list(trace = F)))
    
    est.theta=try1$par[1:ltheta]
    est.theta[est.theta<0]=0
    est.sigma=try1$par[(ltheta+1):(ltheta+lsigma)]
    est.sigma[est.sigma<0]=0
    
    est.matrix=outer(est.theta, est.sigma)
    
    est.array=NULL
    
    for (i in 1:ltheta) {
      
      est.array=c(est.array, est.matrix[i,])
      
    }
    
    LFDR=matrix(NA, ncol=ltheta, nrow=p0)
    theta=theta0; s2=s20
    lik1=t(exp(-.5*( t((outer(theta,grid1,"-"))^2 )/(grid2)^2)) / (grid2*sqrt(2*pi)))
    for (i in 1:p0) {
      
      ddd=lik1[i,]*est.array
      
      UUU=NULL
      for (j in 1:ltheta) {
        
        begin=(j-1)*lsigma+1
        end=j*lsigma
        
        uuu=sum(ddd[begin:end])
        
        UUU=c(UUU, uuu)
        
      }
      UUU=UUU/sum(UUU)
      LFDR[i,]=UUU
      
    }
    
    lfdr=LFDR[,(Btheta+1)]
    
    lfsr=NULL
    
    for (i in 1:p0) {
      
      xi=theta[i]
      
      if(xi>0){
        lfsr[i]=sum(LFDR[i, 1:(Btheta+1)])
      }
      
      if(xi<0){
        lfsr[i]=sum(LFDR[i, (Btheta+1):ltheta])
      }
    }
    
    return(list(grid.theta = gridtheta, 
                grid.sigma2 = gridsigma^2, 
                mix.theta = est.theta, 
                mix.sigma2 = est.sigma, 
                lfdr=lfdr, 
                lfsr=lfsr))
    
  }
