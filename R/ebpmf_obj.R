#' @title ebpmf_obj
#' @param fit a ebpmf fit
#' @details returns the objective function for the current fit. Note that
#' the objective function for SVD-based and flash-based fits may not be comparable.
#' @export
ebpmf_obj = function(fit){
  return(sum(fit$X * fit$M - exp(fit$M+fit$V/2) +
        0.5*(log(2*pi*fit$V)+1)) +   #the relative entropy of Normal with variance V
    +fit$elbo)
}


