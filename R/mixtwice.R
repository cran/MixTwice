mixtwice <-
function (thetaHat, s2, 
                      Btheta = 15, Bsigma2 = 10, 
                      df, 
                      method = c("EM-pava", "AugLag"), 
                      maxit = 100,
                      prop = 1) 
{
  
  ## quality control process
  
  if(length(thetaHat) == 0){
    stop("the input is empty")
  }
  if(!length(thetaHat) == length(s2)){
    stop("the length of estimated effect sizes must be equal to the estimated squared standard errors")
  }
  if(sum(s2 < 0) > 0){
    stop("estimated squared standard errors can not be negative")
  }
  
  ## make a copy of the original data, would be used after mixing distribution estimation
  
  theta0 = thetaHat
  s20 = s2
  
  ## random sampling process
  ok.sample = sample(length(theta0), length(s20) * prop)
  theta = theta0[ok.sample]
  s2 = s20[ok.sample]
  
  p0 = length(theta0)
  p = length(theta)
  cc = max(abs(theta)) * 1.1
  gridtheta = (cc/Btheta) * seq(-Btheta, Btheta, by = 1)
  gridsigma = seq(sqrt(min(s2)), sqrt(max(s2)), by = (sqrt(max(s2)) - 
                                                        sqrt(min(s2)))/Bsigma2)
  ltheta = length(gridtheta)
  lsigma = length(gridsigma)
  grid1 = rep(gridtheta, each = length(gridsigma))
  grid2 = rep(gridsigma, length(gridtheta))
  
  likelihood.theta = matrix(NA, nrow = p, ncol = length(grid1))
  likelihood.s2 = matrix(NA, nrow = p, ncol = length(grid2))
  
  for(i in 1:length(theta)){
    
    likelihood.theta[i,] = stats::dnorm(theta[i], mean = grid1, sd = grid2)
    
    likelihood.s2[i,] = stats::dchisq((df*s2[i]/grid2^2), df = df)*df/grid2^2
    
  }
  
  liklihood.conditional = likelihood.theta * likelihood.s2
  loglik = log(liklihood.conditional)
  
  
  if(method == "EM-pava"){
    
    t1 = rep(c(1:ltheta), each = lsigma)
    t2 = rep(c(1:lsigma), ltheta)
    
    est.theta = rep(1, ltheta)
    est.theta = est.theta/sum(est.theta)
    est.sigma = rep(1, lsigma)
    est.sigma = est.sigma/sum(est.sigma)
    
    count = 0
    notdone = TRUE
    
    while (notdone) {
      
      est = as.numeric(array(est.sigma %o% est.theta))
      
      vv = t(t(loglik) + log(est))
      uu = exp(vv)
      vv = uu/rowSums(uu)
      est = colMeans(vv)
      
      ## de-compose est into est.theta and est.sigma
      
      est.theta = as.numeric(tapply(est, t1, sum))
      est.sigma = as.numeric(tapply(est, t2, sum))
      
      ## use Pava to get est.theta unimodal
      
      tmp = Iso::ufit(y = est.theta, x = 1:length(est.theta), lmode = Btheta + 1)
      est.theta = tmp$y
      
      notdone = (count < maxit )
      count = count + 1
      
    }
    
    
  }
  
  if(method == "AugLag"){
    
    ## the optimization function [likelihood, rather than the conditional sampling]
    
    L = function(x) {
      xtheta = x[1:ltheta]
      xsigma = x[(ltheta + 1):(ltheta + lsigma)]
      yy = array(x[(ltheta + 1):(ltheta + lsigma)] %o% x[1:ltheta])
      return(-sum(log(yy %*% t(liklihood.conditional))))
    }
    
    ## the gradient of the likelihood
    
    G = function(x) {
      g = h = NULL
      xtheta = x[1:ltheta]
      xsigma = x[(ltheta + 1):(ltheta + lsigma)]
      yy = array(x[(ltheta + 1):(ltheta + lsigma)] %o% x[1:ltheta])
      d = yy %*% t(liklihood.conditional)
      for (i in (1:ltheta)) {
        g[i] = -sum((xsigma %*% t(liklihood.conditional[, c((lsigma * (i - 
                                                                         1) + 1):(lsigma * i))]))/d)
      }
      for (j in (1:lsigma)) {
        h[j] = -sum((xtheta %*% t(liklihood.conditional[, seq(j, (j + (ltheta - 
                                                                         1) * (lsigma)), by = lsigma)]))/d)
      }
      return(c(g, h))
    }
    heq <- function(x) {
      h = NULL
      h[1] = sum(x[1:ltheta]) - 1
      h[2] = sum(x[(ltheta + 1):(lsigma + ltheta)]) - 1
      return(h)
    }
    hh1 = c(rep(1, ltheta), rep(0, lsigma))
    hh2 = c(rep(0, ltheta), rep(1, lsigma))
    heq.jac = rbind(hh1, hh2)
    heq.jac.fun = function(x) {
      j = heq.jac
      return(j)
    }
    hin <- function(x) {
      h1 = NULL
      for (i in 1:((ltheta) + (lsigma))) {
        h1[i] = x[i]
      }
      h2 = NULL
      for (i in 1:(Btheta)) {
        h2[i] = x[i + 1] - x[i]
      }
      for (i in (Btheta + 1):(ltheta - 1)) {
        h2[i] = x[i] - x[i + 1]
      }
      h = c(h1, h2)
      return(h)
    }
    hin.jac1 = diag(1, nrow = (ltheta + lsigma), ncol = (ltheta + 
                                                           lsigma))
    hin.jac2 = matrix(0, ncol = ltheta, nrow = ltheta - 1)
    for (i in 1:(Btheta)) {
      hin.jac2[i, i] = -1
      hin.jac2[i, i + 1] = 1
    }
    for (i in (Btheta + 1):(ltheta - 1)) {
      hin.jac2[i, i] = 1
      hin.jac2[i, i + 1] = -1
    }
    hin.jac3 = matrix(0, nrow = ltheta - 1, ncol = lsigma)
    hin.jac = rbind(hin.jac1, cbind(hin.jac2, hin.jac3))
    hin.jac.fun = function(x) {
      j = hin.jac
      return(j)
    }
    est.theta = rep(1, ltheta)
    est.theta = est.theta/sum(est.theta)
    est.sigma = rep(1, lsigma)
    est.sigma = est.sigma/sum(est.sigma)
    a = c(est.theta, est.sigma)
    try1 = suppressWarnings(alabama::auglag(par = a, fn = L, 
                                            gr = G, heq = heq, hin = hin, heq.jac = heq.jac.fun, 
                                            hin.jac = hin.jac.fun, control.outer = list(trace = F)))
    
    est.theta = try1$par[1:ltheta]
    est.theta[est.theta < 0] = 0
    est.sigma = try1$par[(ltheta + 1):(ltheta + lsigma)]
    est.sigma[est.sigma < 0] = 0
    
  }
  
  est.matrix = outer(est.theta, est.sigma)
  est.array = NULL
  for (i in 1:ltheta) {
    est.array = c(est.array, est.matrix[i, ])
  }

  
  ##get back to something before the sampling, recompute the conditional likelihood
  
  LFDR = matrix(NA, ncol = ltheta, nrow = p0)
  theta = theta0
  s2 = s20
  
  likelihood.theta = matrix(NA, nrow = p0, ncol = length(grid1))
  likelihood.s2 = matrix(NA, nrow = p0, ncol = length(grid2))
  
  for(i in 1:length(theta)){
    
    likelihood.theta[i,] = stats::dnorm(theta[i], mean = grid1, sd = grid2)
    
    likelihood.s2[i,] = stats::dchisq((df*s2[i]/grid2^2), df = df)*df/grid2^2
    
  }
  
  liklihood.conditional = likelihood.theta * likelihood.s2
  
  for (i in 1:p0) {
    
    ddd = liklihood.conditional[i, ] * est.array
    UUU = NULL
    
    for (j in 1:ltheta) {
      begin = (j - 1) * lsigma + 1
      end = j * lsigma
      uuu = sum(ddd[begin:end])
      UUU = c(UUU, uuu)
    }
    UUU = UUU/sum(UUU)
    LFDR[i, ] = UUU
  }
  
  lfdr = LFDR[, (Btheta + 1)]
  
  lfsr = NULL
  for (i in 1:p0) {
    xi = theta[i]
    if (xi > 0) {
      lfsr[i] = sum(LFDR[i, 1:(Btheta + 1)])
    }
    if (xi < 0) {
      lfsr[i] = sum(LFDR[i, (Btheta + 1):ltheta])
    }
  }
  
  ans = list(grid.theta = gridtheta, grid.sigma2 = gridsigma^2, 
             mix.theta = est.theta, mix.sigma2 = est.sigma, lfdr = lfdr, 
             lfsr = lfsr)
  return(ans)
  
}
