###----------------------------------------------###
###         Bayesian Isotonic Regression         ###
###----------------------------------------------###

# y: vector of observations
# x:vector of lacations

library("MASS")
library("truncnorm")
library("GIGrvg")

BIR = function(y,x, mc=3000, burn=500){
  n = length(y)
  D = lower.tri(diag(n), diag = T)
  w = diff(x)
  MC = mc + burn
  
  ##MCMC sample box
  eta.pos = matrix(NA, MC, n)
  tau_sqr.pos = matrix(NA, MC, n)
  nu.pos = matrix(NA, MC, n)
  lam_sqr.pos = c()
  xi.pos = c()
  sigma_sqr.pos = c()
  
  ##initial values
  eta = rep(1, n)
  tau_sqr = rep(1, n)
  nu = rep(1, n)
  lam_sqr = 1
  xi = 1
  sigma_sqr = 1
  
  ##MCMC iterations
  for(k in 1:MC){
    #sampling from eta
    e = y - D %*% eta + D[, 1] * eta[1] 
    m = t(e) %*% D[, 1]/(sum(D[, 1]^2) + 1/tau_sqr[1])
    s_sqr = sigma_sqr /(sum(D[, 1]^2) + 1/tau_sqr[1])
    eta[1] = rnorm(n=1,mean = m, sd = sqrt(s_sqr))
    for(j in 2:n){
      e = y - D %*% eta + D[, j] * eta[j] 
      m = t(e) %*% D[, j]/(sum(D[, j]^2) + 1/(lam_sqr*tau_sqr[j]*w[j-1]))
      s_sqr = sigma_sqr /(sum(D[, j]^2) + 1/(lam_sqr*tau_sqr[j]*w[j-1]))
      eta[j] = rtruncnorm(n=1, a = 0, b = Inf, mean = m, sd = sqrt(s_sqr))
    }
    eta.pos[k, ] = eta
  
  #sampling from nu
  nu[1] = rgamma(n=1, shape=3/2, rate = 1+tau_sqr[1])
  for(j in 2:n){
    nu[j] = rgamma(n=1, shape = 1, rate = 1+tau_sqr[j]) 
  }
  nu.pos[k, ] = nu
  
  #sampling from tau_sqr
  tau_sqr[1] = rgig(n=1, lambda = 1/2, chi = eta[1]^2/sigma_sqr , psi = 2 * nu[1])
  for(j in 2:n){
    tau_sqr[j] = rgig(n=1, lambda = 0, chi = eta[j]^2/(sigma_sqr*lam_sqr*w[j-1]), psi = 2 * nu[j])
  }
  tau_sqr.pos[k, ] = tau_sqr
  
  #sampling from xi
  xi = rgamma(n=1, shape = 1, rate = 1+lam_sqr)
  xi.pos[k] = xi
  
  #sampling from lambda_sqr
  lam_sqr = rgig(n=1, lambda = (-n+2)/2, chi = sum((eta^2/tau_sqr)[2:n]/w)/sigma_sqr, psi = 2*xi)
  lam_sqr.pos[k] = lam_sqr
  
  #sampling from sigma_sqr
  scale = 0.5 * (sum((y - D%*%eta)^2) + lam_sqr^(-1) * sum((eta^2/tau_sqr)[2:n]/w) + eta[1]^2/tau_sqr[1])
  sigma_sqr = 1/rgamma(n=1, shape = n, scale = 1/scale)
  sigma_sqr.pos[k] = sigma_sqr
  }
  
  ##Summary
  om = 1:burn
  eta.pos = eta.pos[-om, ]
  tau_sqr.pos = tau_sqr.pos[-om, ]
  nu.pos = nu.pos[-om, ]
  lam_sqr.pos = lam_sqr.pos[-om]
  xi.pos = xi.pos[-om]
  sigma_sqr.pos = sigma_sqr.pos[-om]
  Res = list(eta = eta.pos, tau_sqr = tau_sqr.pos, nu = nu.pos, lam_sqr = lam_sqr.pos, xi = xi.pos, sigma_sqr = sigma_sqr.pos)
  return(Res)
}