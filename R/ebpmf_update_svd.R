ebpmf_update_svd_als = function(fit){

  if(is.null(fit$udv)){
    stop("must initialize svd using ebpmf_init_udv before updating")
  }

  # implements als update; based on code in softImpute package (file simpute.als.R)
  # To deal with column-specific tau we apply these updates to M with columns scaled by sqrt(tau) and then rescales v at the end
  UtM = t(fit$udv$u) %*% fit$M
  s = svd(t(UtM) * sqrt(fit$tau))
  fit$udv$v = s$u
  fit$udv$d = s$d
  fit$udv$u = fit$udv$u %*% s$v


  MV = fit$M %*% (fit$udv$v * sqrt(fit$tau))
  s = svd(MV)
  fit$udv$u = s$u
  fit$udv$d = s$d
  fit$udv$v = (fit$udv$v %*% s$v) / sqrt(fit$tau)


  # order in terms of decreasing singular values; not necessary but included
  # for neatness
  o = order(fit$udv$d,decreasing=TRUE)
  fit$udv$u = fit$udv$u[,o,drop=FALSE]
  fit$udv$d = fit$udv$d[o,drop=FALSE]
  fit$udv$v = fit$udv$v[,o,drop=FALSE]

  fit$B = expand_svd(fit$udv)

  fit = ebpmf_update_tau(fit)

  fit$elbo = svd_elbo(fit) #elbo for the svd (point mass on B) is E_mu log(p(mu|B,tau))

  return(fit)
}

ebpmf_update_svd_naive = function(fit,nu=1){
  fit$udv = svd(t(sqrt(fit$tau)*t(fit$M)),nu=nu,nv=nu)
  fit$udv$v = fit$udv$v/sqrt(fit$tau)

  fit$udv$d = fit$udv$d[1:nu]

  fit$B = expand_svd(fit$udv)

  fit = ebpmf_update_tau(fit)

  fit$elbo = svd_elbo(fit)

  return(fit)
}

