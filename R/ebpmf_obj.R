ebpmf_obj = function(fit){
  return(sum(fit$X * fit$M - exp(fit$M+fit$V/2) +
        0.5*(log(2*pi*fit$V)+1)) +   #the relative entropy of Normal with variance V
    +fit$elbo)
  #+ sum(dnorm(fit$M-fit$B,mean=0,sd=sqrt(1/fit$tau),log=TRUE))
  #- sum(fit$V*fit$tau/2))
  #last line is only for svd; done while debugging..., should use elbo if using flash
}
