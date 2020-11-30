#### PPLC
pplc = function(yinew, yi, k) {
  var.sum = sum(apply(yinew, 2, var))
  sse = (yi-apply(yinew, 2, mean))^2
  var.sum + k/(k+1)*sum(sse)
}

#### CPO
f_ig_inv <- function(y, x, b0, b1, phi){
  mu <- 1/(b0 + b1*x)
  sqrt(2*pi*phi)*y^3/2*exp((y - mu)^2/(2*phi*y*mu^2))
}

f_gam_inv <- function(y, x, b0, b1, nu){
  mu <- 1/(b0 + b1*x)
  exp(lgamma(nu) + nu*log(mu/nu)+(1-nu)*log(y) + nu*y/mu)
}

cpo_ig <- sapply(1:length(xi), function(i){
  1/mean(f_ig_inv(yi[i], xi[i], b1.ig, b2.ig, phi.ig))
})

cvLS.ig <- mean(log(cpo_ig))

cpo_gam <- sapply(1:length(xi), function(i){
  1/mean(f_gam_inv(yi[i], xi[i], b1.inv.gamma, b2.inv.gamma, nu.inv.gamma))
})

#### Binomial RC (cloglog link, normal prior on beta)
post.cloglog.binom.norm = function(beta) {
  mle.b1 = -39.572
  mle.b2 = 22.041
  prior.sigma = 100*diag(2)
  post = dmvnorm(beta, c(mle.b1,mle.b2), prior.sigma, log = TRUE) +
    sum(dbinom(x = yi,size = mi,prob = 1-exp(-exp(beta[1] + beta[2]*xi)), log=TRUE))
  return(post)
}

#### Binomial RC (cloglog link, flat prior on beta)
post.cloglog.binom.flat = function(beta) {
  mle.b1 = -39.572
  mle.b2 = 22.041
  prior.sigma = 100*diag(2)
  post = sum(dbinom(x = yi,size = mi, 
                    prob = 1-exp(-exp(beta[1] + beta[2]*xi)), 
                    log=TRUE))
  return(post)
}

#### Binomial RC (logit link, flat prior on beta)
post.logit.binom.flat = function(beta) {
  post = sum(dbinom(x = yi, size = mi, 
                    prob = (exp(beta[1] + beta[2]*xi)) / (1 + exp(beta[1] + beta[2]*xi)), 
                    log=TRUE))
  return(post)
}

#### Binomial RC (logit link, flat prior on beta, beta prior on alpha)
post.genlink.betabinom.loglike = function(beta) {
  b1 = beta[1]; b2 = beta[2]; b3 = beta[3]
  eta = b1 + b2 * xi
  link = (exp(eta * b3)) / (1 + exp(eta))^b3
  if(any(link < 0) | any(link > 1)) {
    return(-Inf)
  }
  post = sum(yi * log(link) + (mi - yi) * log(1 - link) + 
               dbeta(beta[3],3,6,log = TRUE))
  return(post)
}

#### Binomial RC (logit link, flat prior on beta, gamma prior on alpha)
post.genlink.gammabinom.loglike = function(beta) {
  
  b1 = beta[1]; b2 = beta[2]; b3 = beta[3]
  eta = b1 + b2 * xi
  link = (exp((eta) * b3)) / (1 + exp(eta))^b3
  alpha.prior = 1/dgamma(b3, .00001, .00001, log = TRUE)
  
  if(any(link < 0) | any(link > 1)) {
    return(-Inf)
  }
  post = sum(yi * log(link) + (mi - yi) * log(1 - link) +
               alpha.prior,
             na.rm=TRUE)
  
  return(post)
}

#### Gamma Log Likelihood (inverse link)
gamma.llh.invlink = function(beta) {
  
  b1=beta[1]; b2=beta[2]; b3=beta[3]
  
  prior.mu = c(b1.gam.mle, b2.gam.mle)
  prior.Sigma = 100*diag(2)
  prior.shape = 0.0001
  prior.rate = 0.0001
  link = b1+b2*xi
  if (any((link)<0)) return(-Inf)
  
  gam.prior = 1/dgamma(exp(b3), prior.shape, prior.rate, log=TRUE)
  norm.prior = mvtnorm::dmvnorm(c(b1, b2), prior.mu, prior.Sigma, log=TRUE)
  likelihood = sum(-lgamma(exp(b3)) + exp(b3)*b3 + exp(b3)*log(link) + 
                     (exp(b3)-1)*log(yi)-yi*exp(b3)*(link))
  jacobian = b3
  density = gam.prior + norm.prior + likelihood + jacobian
  return(density)
}

#### Inverse Gaussian Random Component (inverse link)
loglike.ig = function(beta) {
  b1 = beta[1]
  b2 = beta[2]
  b3 = beta[3]
  prior.mu = c(b1.ig.mle, b2.ig.mle)
  prior.Sigma = 100 * diag(2)
  prior.shape = 0.0001
  prior.rate = 0.0001
  phi = exp(b3)
  link = b1 + b2 * xi
  if (any(link < 0))
    return(-Inf)
  gam.prior = dgamma(exp(b3), prior.shape, prior.rate, log = TRUE)
  norm.prior = mvtnorm::dmvnorm(c(b1, b2), prior.mu, prior.Sigma, log =
                                  TRUE)
  likelihood = sum(-1 / 2 * log(2 * pi * phi * yi ^ 3) - yi / (2 * phi) *
                     (link) ^ 2 +
                     1 / phi * (link) - 1 / (2 * phi * yi))
  jacobian = b3
  density = gam.prior + norm.prior + likelihood + jacobian
  return(density)
}

############# MCMC CODE ###################
  
my.mcmc = function(initial_beta,
                   log.post,
                   cov.mat,
                   niter,
                   var.tune) {
  current = initial_beta
  beta.store = matrix(NA, niter, length(current))
  accept = 0
  
  for (i in 1:niter) {
    proposed = mvtnorm::rmvnorm(1, current, cov.mat * var.tune)
    #print(proposed); #browser()
    ratio = log.post(proposed) - log.post(current)
    if (log(runif(1)) < ratio) {
      current = proposed
      accept = accept + 1
    }
    else{
      current = current
      accept = accept
    }
    beta.store[i,] = current
  }
  #print(accept/niter)
  out = list(beta.store = beta.store)
}
  
### Three Parameter MCMC (positive truncation)
mcmc.3param = function(current,
                       covar.mat,
                       var.tuning,
                       niter,
                       log.post) {
  accept = 0
  beta = c(current[1], current[2])
  alpha = current[3]
  
  param_store = matrix(nrow = niter, ncol = 3)
  
  for (i in 1:niter) {
    prop = rtmvnorm(1,
                    c(beta, alpha),
                    covar.mat * var.tuning,
                    lower = c(-Inf, -Inf, 0))
    ratio = log.post(prop) - log.post(c(beta, alpha))
    
    if (log(runif(1)) < ratio) {
      beta = prop[1:2]
      alpha = prop[3]
      accept = accept + 1
    }
    param_store[i, ] = c(beta, alpha)
  }
  
  out = list(beta.store = param_store, accept = accept)
}

binom.glm.MCMC = function(initial_beta,
                          covar.mat,
                          var.tuning,
                          niter,
                          log.post) {
  accept = 0
  beta = initial_beta
  beta.store <- matrix(NA, niter, length(beta))
  
  for (i in 1:niter) {
    prop <- mvtnorm::rmvnorm(1, beta, var.tuning * covar.mat)
    ratio <- log.post(prop) - log.post(beta)
    if (log(runif(1)) < ratio) {
      beta = prop
      accept = accept + 1
    }
    else{
      beta = beta
      
      accep = accept
    }
    beta.store[i, ] = beta
  }
  out = list(beta.store = beta.store, accept = accept)
}

# Using truncated normal
binom.genlink.glm.MCMC = function(initial_beta,
                                  covar.mat,
                                  var.tuning,
                                  niter,
                                  log.post) {
  accept = 0
  beta = c(initial_beta[1], initial_beta[2])
  alpha = initial_beta[3]
  
  param_store = matrix(nrow = niter, ncol = 3)
  
  for (i in 1:niter) {
    prop = rtmvnorm(1,
                    c(beta, alpha),
                    covar.mat * var.tuning,
                    lower = c(-Inf, -Inf, 0))
    ratio = log.post(prop) - log.post(c(beta, alpha))
    
    if (log(runif(1)) < ratio) {
      beta = prop[1:2]
      alpha = prop[3]
      accept = accept + 1
    }
    param_store[i, ] = c(beta, alpha)
  }
  
  out = list(beta.store = param_store, accept = accept)
}