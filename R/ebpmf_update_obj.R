ebpmf_update_obj = function(fit){
  fit$obj = c(fit$obj,ebpmf_obj(fit))
  return(fit)
}
