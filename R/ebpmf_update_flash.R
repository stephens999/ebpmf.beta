ebpmf_update_flash = function(fit,useV = FALSE){

  ## maybe these first two lines should be part of the mu update?
  ## To maximize the splitting ELBO we shoudl use V, but I suspect
  ## This may lead to underfitting, so have also implemented an option
  ## not to use V
  if(useV)
    new.flash = flash_init(fit$M, V = fit$V)
  else
    new.flash = flash_init(fit$M)

  new.flash = flash_factors_init(new.flash, fit$f) # this doesn't work if fit$f has no factors; maybe should add this?

  new.flash = flash_backfit(new.flash,maxiter=1) # this is necessary for now because the above initialization is not valid until after backfit

  new.flash = flash_greedy(new.flash,Kmax = 1) #or new.flash = do.call(flash_fn,new.flash)
  new.flash = flash_nullcheck(new.flash)

  #fit$f = new.flash
  fit$B = fitted(new.flash)
  fit$elbo = new.flash$elbo # log p(M | g,tau)
  fit$tau = new.flash$flash_fit$tau
  fit$est.tau.dim = new.flash$flash_fit$est.tau.dim

  return(fit)
}



