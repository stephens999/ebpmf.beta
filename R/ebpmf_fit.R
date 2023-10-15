ebpmf_fit = function(X){
  fit = ebpmf_init(X)
  for(i in 1:10){
    fit = ebpmf_update_mu(fit)
    fit = ebpmf_update_flash(fit) #,flash_fn) # apply flash_fn to M, V
    print(ebpmf_obj(fit))
  }
  return(fit)
}
