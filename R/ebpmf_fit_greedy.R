#' @title ebpmf_fit_greedy
#' @details adds one factor at a time using flash or svd (if ebnm_fn = NULL)
#' The first factor is fit by svd, so unconstrained.
#' @param backfit_first_factor whether to backfit the first factor after every iteration
#' I thought it might be useful if ebnm_fn is non-negative as otherwise the greedy fit
#' may stop prematurely. but it turns out the just want to keep adding factors in this case so not recommended. Could maybe be removed
#' and not recommended for now.
#' @export
ebpmf_fit_greedy = function(X,Kmax=100,ebnm_fn = NULL,tol=1e-3,maxiter.mu = 10, est.tau.dim = 0, backfit_first_factor = FALSE){

  if(is.matrix(X)){
    fit = ebpmf_init(X, est.tau.dim = est.tau.dim)
    fit = ebpmf_update_mu(fit,maxiter = maxiter.mu)
    fit = ebpmf_update_r1(fit) #uses svd method so first factor unconstrained
  } else { # X is a previous fit
    fit=X
    if(!missing(est.tau.dim)){
      warning("Initializing based on previous fit; ignoring value of est.tau.dim")
    }
    fit$obj = -Inf
  }

  for(i in 1:Kmax){
    fit = ebpmf_update_mu(fit,maxiter = maxiter.mu)
    fit = ebpmf_update_r1(fit,ebnm_fn = ebnm_fn)

    if(backfit_first_factor){
      fit = ebpmf_update_mu(fit,maxiter = maxiter.mu)
      fit = ebpmf_update_r1(fit,factor_to_update = 1)
    }

    fit = ebpmf_update_obj(fit)

    if(converged(fit,tol))
      break
  }
  return(fit)
  fit = ebpmf_init(X)
  for(i in 1:10){
    fit = ebpmf_update_mu(fit)
    fit = ebpmf_update_flash(fit) #,flash_fn) # apply flash_fn to M, V
    print(ebpmf_obj(fit))
  }
  return(fit)
}
