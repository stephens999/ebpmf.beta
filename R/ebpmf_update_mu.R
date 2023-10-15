ebpmf_update_mu = function(fit){
  res = vga_pois_solver_mat_newton(fit$M,
                                   fit$X,
                                   1,
                                   fit$B,
                                   get_sigma2(fit),
                                   maxiter = 1000,
                                   tol=1e-8,return_V = TRUE)
  fit$M = res$M
  fit$V = res$V

  return(fit)
}
