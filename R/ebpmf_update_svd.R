ebpmf_update_svd_als = function(fit){

  if(is.null(fit$udv)){
    stop("must initialize svd using ebpmf_init_udv before updating")
  }

  # implements als update; based on code in softImpute package (file simpute.als.R)
  UtM = t(fit$udv$u) %*% fit$M
  s = svd(t(UtM))
  fit$udv$v = s$u
  fit$udv$d = s$d
  fit$udv$u = fit$udv$u %*% s$v

  MV = fit$M %*% fit$udv$v
  s = svd(MV)
  fit$udv$u = s$u
  fit$udv$d = s$d
  fit$udv$v = fit$udv$v %*% s$v

  # order in terms of decreasing singular values; not necessary but included
  # for neatness
  o = order(fit$udv$d,decreasing=TRUE)
  fit$udv$u = fit$udv$u[,o,drop=FALSE]
  fit$udv$d = fit$udv$d[o,drop=FALSE]
  fit$udv$v = fit$udv$v[,o,drop=FALSE]

  fit$B = expand_svd(fit$udv)

  res = fit$M-fit$B

  fit$tau = 1/mean((res)^2+fit$V)
  fit$est.tau.dim = 0 #  Because the svd assumes all tau are equal;
  # could change to column-specific, but would need to change SVD routine

  fit$elbo = -0.5*fit$n*fit$p*log(2*pi/fit$tau) -
    0.5*fit$tau*sum(res^2) - sum(fit$V*fit$tau/2) #elbo for the svd (point mass on B) is E_mu log(p(mu|B,tau))

  return(fit)
}

ebpmf_update_svd_naive = function(fit,nu=1){
  fit$udv = svd(fit$M,nu=nu,nv=nu)
  fit$udv$d = fit$udv$d[1:nu]

  fit$B = expand_svd(fit$udv)

  res = fit$M-fit$B

  fit$tau = 1/mean((res)^2+fit$V)
  fit$est.tau.dim = 0 #  Because the svd assumes all tau are equal;
  # could change to column-specific, but would need to change SVD routine

  fit$elbo = sum(dnorm(res,mean=0,sd=sqrt(1/fit$tau),log=TRUE)) - sum(fit$V*fit$tau/2) #elbo for the svd (point mass on B) is E_mu log(p(mu|B,tau))

  return(fit)
}
