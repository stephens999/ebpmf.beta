#'@title a matrix version of the vga solver using Newton's method
#'@importFrom Rfast Pmin
vga_pois_solver_mat_newton = function(M,X,S,Beta,Sigma2,maxiter=1000,tol=1e-8,return_V = TRUE){

  #X = as.matrix(X)
  const0 = (Sigma2*X+Beta + 1)
  const1 = 1/Sigma2
  const2 = Sigma2/2
  const3 = Beta/Sigma2

  # make sure m < sigma2*x+beta
  M = Pmin(M,const0-1)
  # idx = (M>(const0-1))
  # if(sum(idx)>0){
  #   M[idx] =suppressWarnings(vga_pois_solver_bisection(c(X[idx]),c(S[idx]),c(Beta[idx]),c(Sigma2[idx]),maxiter = 10)$m)
  # }
  for(i in 1:maxiter){
    temp = (const0-M)
    sexp = S*exp(M+const2/temp)
    # f = X - sexp - (M-Beta)/Sigma2
    f = X - sexp - M*const1 + const3
    if(max(abs(f))<tol){
      break
    }
    # f_grad = -sexp*(1+const2/temp^2)-const1
    # direction = (X - sexp - (M-Beta)/Sigma2)/(-sexp*(1+const2/temp^2)-const1)
    M = M - f/(-sexp*(1+const2/temp^2)-const1)

  }
  #gc()
  if(return_V){
    return(list(M=M,V=Sigma2/temp))
  }else{
    return(M)
  }

}
