#' @param fit an existing ebpmf fit
#' @param nu number of singular vectors (columns of u,v) to use
#' @details Initializes by randomly drawing any needed additional columns for u as white noise and performing one ALS update.
#' (If udv is already specified, and nu is less than the numer of columns of u, then the initialization removes the excess columns)
ebpmf_init_udv = function(fit,nu){
  if(is.null(fit$udv)){
    fit$udv = list(u=NULL,d=NULL,v=NULL)
  }
  if(is.null(fit$udv$u)){
    fit$udv$u = matrix(rnorm(n*nu),nrow=n)
  }
  ncu = ncol(fit$udv$u)
  if(nu>ncu){
    fit$udv$u = cbind(fit$udv$u,matrix(rnorm(n*(nu-ncu)),nrow=n))
  }
  if(nu<ncu){
    fit$udv$u = fit$udv$u[,1:nu]
  }
  fit = ebpmf_update_svd_als(fit) # performs one iteration of ALS to get udv
}
