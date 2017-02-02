#############################################3
#
#
# Maximum likelihood estimation
#
#
#############################################

#############################################
#' 
#' Fits distributions to data using newton raphson#' 
#' Currently implemented
#' 1) Beta distribution
#' 2) ...
#' Also included are unit tests for each distribution
#' 
#############################################

#Uses my implmentation of optimization methods in the package
#optimization
source("optimization.r")



#############################################################
#'
#'
#' 1) Beta Distribution
#'
#'
#############################################################

#'
#' Method of the moments
#' 
mm_beta <- function(mu, var) {
  alpha <-    mu  * (mu*(1-mu)/var -1)
  beta  <- (1-mu) * (mu*(1-mu)/var -1)
  return(params = list(alpha = alpha, beta = beta))
}

#'
#' MLE 
#' 
mle_beta(p){
  
  if(any(is.na(p))) error("P values vector contains NAs")
  
  pbm = mm_beta(mu = mean(p),var = var(p))
  
  if(!all(is.finite(c(pbm$alpha,pbm$beta)))) error("Method of the moments failed, aborting. Are you sure this is beta distributed?")
  
  sol=newton_raphson(fun=loglik_beta, grad=gradient_beta,hess=hessian_beta,data=p,start=c(pbm$alpha,pbm$beta),eta=0.05, epsilon=1E-12, nmax=10000)

  #Test convergence & Test gradient should be close to zero
  if(sol$conv == F | any(sol$g_hat>1E-3)){ 
    message("Beta mle failed, calculating estimates using the method of the moments")
    alpha = pbm$alpha
    beta  = pbm$beta
  }else{
    alpha = sol$param_hat[1]
    beta  = sol$param_hat[2]
  }
  
  return(list(alpha=alpha,beta=beta));
  
}

#'
#' Simulation test of the goodness of fit
#' 
gof_beta(p,alpha_hat,beta_hat){
  n= length(p)
  ks.pvalue = mean(vapply(X = 1:100,FUN.VALUE = 0.001 ,FUN = 
                            function(x){
                              testFit = ks.test(x = rbeta(n, shape1 = alpha_hat, shape2 = beta_hat),
                                                y = p,
                                                alternative = "two.sided",exact = T)
                              return(testFit$p.value)
                              }))
  return(ks.pvalue)
}


gradient_beta <-function( params, data ){
  p     = data
  n     = length(p)
  a     = params[1]
  b     = params[2]
  
  df.da = n*digamma(a+b) - n*digamma(a) + sum(log(p))   #ok
  df.db = n*digamma(a+b) - n*digamma(b) + sum(log(1-p)) #ok
  return(c(df.da,df.db))
}
hessian_beta  <-function( params, data ){
  p     = data
  n     = length(p)
  m     = length(params)
  a     = params[1]
  b     = params[2]
  
  d2f.da2 = n*trigamma(a+b) - n*trigamma(a)
  d2f.db2 = n*trigamma(a+b) - n*trigamma(b) 
  d2f.dab = n*trigamma(a+b) 
  H       = matrix(data = c(d2f.da2,d2f.dab,d2f.dab,d2f.db2),byrow = T,nrow = m,ncol = m)
  return(H)
}
loglik_beta   <-function( params, data ){
  p  = data
  n  = length(p)
  a  = params[1]
  b  = params[2]
  ll = n*log(gamma(a+b)) - n*log(gamma(a)) - n*log(gamma(b)) + (a-1)*sum(log(p)) + (b-1)*sum(log(1-p)) #ok
  return(ll)
}




#############################################################
#'
#'
#' 3) Unit tests for all distributions
#'
#'
#############################################################



unit_test <- function(){
  
  epsilon = 1E-6
  nmax    = 1000
  
  ###########################################################################################################
  #Test case 1
  a = 1
  b = 10
  p = rbeta(n = nmax,shape1 = a,shape2 = b) 
  
  sol=newton_raphson(fun=loglik_beta, grad=gradient_beta,hess=hessian_beta,data=p,start=c(a+3,b+5),eta=0.05, epsilon=1E-12, nmax=10000)
  
  print(t(sol$param_hat))
  print(c(a,b))
  
  stopifnot( sol$conv == T)                                   #convergence
  stopifnot( all(sol$g_hat<epsilon))                          #gradient should be close to zero
  stopifnot( all(abs(sol$param_hat - c(a,b)) < 0.05*c(a,b)) ) #Less than 5% error from the original solution
  stopifnot( ( loglik_beta(params=c(a,b), data=p) - sum(log(dbeta(p,a,b))) ) < epsilon )
  
  print("Unit test 1: PASSED")
  
  ###########################################################################################################
  #Test case 2
  a = 0.8
  b = 20
  p = rbeta(n = nmax,shape1 = a,shape2 = b)
  
  sol=newton_raphson(fun=loglik_beta, grad=gradient_beta,hess=hessian_beta,data=p,start=c(a+2,b+5),eta=0.05, epsilon=1E-12, nmax=10000)
  
  print(t(sol$param_hat))
  print(c(a,b))
  
  stopifnot( sol$conv == T)                                   #convergence
  stopifnot( all(sol$g_hat<epsilon))                          #gradient should be close to zero
  stopifnot( all(abs(sol$param_hat - c(a,b)) < 0.05*c(a,b)) ) #Less than 5% error from the original solution
  stopifnot( ( loglik_beta(params=c(a,b), data=p) - sum(log(dbeta(p,a,b))) ) < epsilon )
  
  print("Unit test 2: PASSED")
  
  ###########################################################################################################

  
}