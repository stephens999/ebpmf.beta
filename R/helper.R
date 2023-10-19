adjust_var_shape <- function (sigma2, var_type, n, p)
{
  if (var_type == 0) {
    sigma2 = matrix(sigma2, nrow = n, ncol = p)
  }
  else if (var_type == 1) { #by row
    sigma2 = matrix(sigma2, nrow = n, ncol = p, byrow = F)
  }
  else if (var_type == 2) { #by column
    sigma2 = matrix(sigma2, nrow = n, ncol = p, byrow = T)
  }
  else {
    stop("Non-supported var type")
  }
  sigma2
}

get_sigma2 = function(fit){
  adjust_var_shape(1/fit$tau, fit$est.tau.dim, fit$n, fit$p)
}

converged = function(fit,tol=1e-3){
  o = rev(fit$obj)
  if(abs(o[1]-o[2])< tol)
    return(TRUE)
  else
    return(FALSE)
}

expand_svd = function(s){
  return(s$u %*% (s$d * t(s$v)))
}

svd_elbo = function(fit){
  res = fit$M-fit$B
  if(fit$est.tau.dim==0){ # constant tau
    return(-0.5*fit$n*fit$p*log(2*pi/fit$tau) -
             0.5*fit$tau*sum(res^2) - 0.5* fit$tau* sum(fit$V))
  } else { # column-specific tau
    return(-0.5*fit$n*sum(log(2*pi/fit$tau)) -
             0.5*sum(fit$tau*colSums(res^2)) - 0.5* sum(fit$tau* colSums(fit$V)))
  }
}
