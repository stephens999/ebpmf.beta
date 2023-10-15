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
  if(o[1]-o[2]< tol)
    return(TRUE)
  else
    return(FALSE)
}
