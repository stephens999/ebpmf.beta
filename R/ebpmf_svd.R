#' @title Simple SVD version of EBPMF
#' @param X a matrix of counts
#' @param nu the number of singular vectors to compute
#' @param maxiter maximum number of iterations
#' @param tol convergence tolerance
#' @param method string indicating method to use for update, either "als" (alternating least squares) or "naive" (which directly calls svd)
#' @details Runs a version of EBPMF where the prior families for L and F are point masses
#' The EBMF step is equivalent to performing an SVD
#' @examples
#' set.seed(1)
#' n= 10
#' p = 20
#' K = 4
#' LL = matrix(rnorm(n*K),nrow=n)
#' FF = matrix(rnorm(p*K),nrow=p)
#' mu = matrix(LL %*% t(FF),ncol=p, nrow=n)
#' X= matrix(rpois(n*p,exp(mu)), ncol=p, nrow=n)
#' fit = ebpmf_svd(X,nu=5)
#' fit$tau # should be big
#' plot(fit$obj) # should show converged
#' @export
ebpmf_svd = function(X,nu,maxiter=100,tol=1e-3,method=c("naive","als")){
  method = match.arg(method)
  fit = ebpmf_init(X)
  fit = ebpmf_init_udv(fit,nu)
  for(i in 1:maxiter){
    fit = ebpmf_update_mu(fit)
    if(method == "als")
      fit = ebpmf_update_svd_als(fit)
    else if(method=="naive")
      fit = ebpmf_update_svd_naive(fit,nu)
    fit = ebpmf_update_obj(fit)
    if(converged(fit,tol))
      break
  }
  return(fit)
}
