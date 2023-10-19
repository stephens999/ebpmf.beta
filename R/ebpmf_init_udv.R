#' @title ebpmf_init_udv
#' @param fit an existing ebpmf fit
#' @param nu number of singular vectors (columns of u,v) to use
#' @details Initializes the u element (only) of the udv list by randomly drawing any needed additional columns for u as white noise
#' (If udv is already specified, and nu is less than the number of columns of u, then the initialization removes the excess columns)
ebpmf_init_udv = function(fit,nu){
  if(is.null(fit$udv)){
    fit$udv = list(u=NULL,d=NULL,v=NULL)
  }
  if(is.null(fit$udv$u)){
    fit$udv$u = matrix(rnorm(fit$n*nu),nrow=fit$n)
  }
  ncu = ncol(fit$udv$u)
  if(nu>ncu){
    fit$udv$u = cbind(fit$udv$u,matrix(rnorm(fit$n*(nu-ncu)),nrow=fit$n))
  }
  if(nu<ncu){
    fit$udv$u = fit$udv$u[,1:nu,drop=FALSE]
  }
  return(fit)
}
