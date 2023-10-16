ebpmf_init= function(X,M = NULL, B=NULL){
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
    fit$B = zero_mat
  else
    fit$B = B

  fit$udv = NULL

  fit$V = zero_mat # equivalent to initializing q_mu to a point mass at M

  res = fit$M-fit$B
  fit$tau = 1/mean((res)^2+fit$V)
  fit$est.tau.dim = 0 #  Because the svd assumes all tau are equal;
  # could change to column-specific, but would need to change SVD routine

  fit$elbo = sum(dnorm(res,mean=0,sd=sqrt(1/fit$tau),log=TRUE)) - sum(fit$V*fit$tau/2) #elbo for the svd (point mass on B) is E_mu log(p(mu|B,tau))

  #fit = ebpmf_update_svd(fit)
  fit$obj = -Inf

  return(fit)
}
