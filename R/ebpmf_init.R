ebpmf_init= function(X,M = NULL, B=NULL, est.tau.dim = 0){
  fit = list()
  fit$X = X
  fit$n = nrow(X)
  fit$p = ncol(X)
  zero_mat = matrix(0,nrow=fit$n,ncol=fit$p)

  if(is.null(M))
    fit$M = log(0.5+X) #a rough initial guess
  else
    fit$M = M

  if(is.null(B))
    fit$B = zero_mat + mean(fit$M)
  else
    fit$B = B

  fit$udv = NULL

  fit$V = zero_mat # equivalent to initializing q_mu to a point mass at M

  fit$est.tau.dim = est.tau.dim
  fit = ebpmf_update_tau(fit,sub=(X>0))

  fit$elbo = svd_elbo(fit)

  #fit = ebpmf_update_svd(fit)
  fit$obj = -Inf

  return(fit)
}
