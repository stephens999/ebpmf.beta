#' @title Simple SVD version of EBPMF
#' @param X a matrix of counts
#' @param nu the number of singular vectors to compute
#' @param maxiter maximum number of iterations
#' @param tol convergence tolerance
#' @details
ebpmf_svd = function(X,nu,maxiter=100,tol=1e-3){
  fit = ebpmf_init(X)
  for(i in 1:maxiter){
    fit = ebpmf_update_mu(fit)
    fit = ebpmf_update_svd(fit,nu=nu)
    fit = ebpmf_update_obj(fit)
    if(converged(fit,tol))
      break
  }
  return(fit)
}
