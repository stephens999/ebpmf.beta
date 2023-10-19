#' @param sub indicator of which observations to use when updating tau (uses all obervations if NULL)
ebpmf_update_tau = function(fit, sub=NULL){

  if(is.null(sub)){sub=matrix(TRUE,nrow=fit$n,ncol=fit$p)}
  res = fit$M - fit$B

  if(fit$est.tau.dim == 0){ # assume all equal
      fit$tau = sum(sub)/sum(sub*(res^2 + fit$V))
  } else if(fit$est.tau.dim == 2){
      fit$tau = colSums(sub)/colSums(sub*(res^2 + fit$V))
  } else {
    stop("implemented only for est.tau.dim = 0 or 2")
  }
  return(fit)
}
