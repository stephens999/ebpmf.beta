ebpmf_update_svd = function(fit,nu=1){
  s = svd(fit$M,nu=nu,nv=nu)
  fit$B = s$u %*% (s$d[1:nu] * t(s$v))

  res = fit$M-fit$B

  fit$tau = 1/mean((res)^2+fit$V)
  fit$est.tau.dim = 0 #  Because the svd assumes all tau are equal;
  # could change to column-specific, but would need to change SVD routine

  fit$elbo = sum(dnorm(res,mean=0,sd=sqrt(1/fit$tau),log=TRUE)) - sum(fit$V*fit$tau/2) #elbo for the svd (point mass on B) is E_mu log(p(mu|B,tau))

  return(fit)
}
