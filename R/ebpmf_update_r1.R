#' @details Adds a single factor to current fit using either flash (if ebnm_fn supplied) or svd otherwise.
#' This is done by applying flash/svd to M-B. Also updates tau
#' If factor to update is supplied then that factor is updates ("backfitting")
#' Otherwise a new factor is added
ebpmf_update_r1 = function(fit, factor_to_update = NULL, ebnm_fn = NULL){

  if(!is.null(factor_to_update)) #remove the factor to be updated from current fit
    fit$B = fit$B - fit$udv$d[factor_to_update] * fit$udv$u[,factor_to_update] %*% t(fit$udv$v[,factor_to_update])

  if(!is.null(ebnm_fn)){
    dat = fit$M-fit$B

    ffit = flash(dat, greedy_Kmax = 1, ebnm_fn = ebnm_fn)

    if(ffit$n_factors == 0){
      message("flash decided not to add a factor")
      return(fit)
    }

    s = list()

    #IMPORTANT: I re-estimate d to avoid overshrinking; this is because I found
    # that the current greedy approach, which does not keep a track of second
    # moments of u and v, tends to repeatedly fit the same factor
    s$u = ffit$L_pm
    s$v = ffit$F_pm
    s$d = (t(s$u) %*%  dat %*% s$v) /
      (sum(s$u^2)*sum(s$v^2))
  } else { #run svd method
    s = svd(t(sqrt(fit$tau)*t(fit$M-fit$B)),nu=1, nv=1)
    s$v = s$v[,drop=FALSE]/sqrt(fit$tau)
  }

  # update fit$udv
  if(is.null(fit$udv)){
    fit$udv$u = s$u
    fit$udv$d = s$d[1]
    fit$udv$v = s$v
  } else if(is.null(factor_to_update)){
      fit$udv$u = cbind(fit$udv$u,s$u)
      fit$udv$v = cbind(fit$udv$v,s$v)
      fit$udv$d = c(fit$udv$d, s$d[1])
  } else {
      fit$udv$u[,factor_to_update] = s$u
      fit$udv$v[,factor_to_update] = s$v
      fit$udv$d[factor_to_update] = s$d[1]
  }


  fit$B = expand_svd(fit$udv)

  fit = ebpmf_update_tau(fit)

  fit$elbo = svd_elbo(fit)

  return(fit)
}
