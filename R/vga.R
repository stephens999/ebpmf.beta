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

  #this version stops updating subsets of columns of M once they are converged to speed things up
  subset = 1:ncol(M)
  S = matrix(S, nrow=nrow(M),ncol=ncol(M))

  for(i in 1:maxiter){
    temp = (const0[,subset]-M[,subset])
    sexp = S[,subset]*exp(M[,subset]+const2[,subset]/temp)
    # f = X - sexp - (M-Beta)/Sigma2
    f = X[,subset,drop=FALSE] - sexp - M[,subset]*const1[,subset] + const3[,subset]

    # f_grad = -sexp*(1+const2/temp^2)-const1
    # direction = (X - sexp - (M-Beta)/Sigma2)/(-sexp*(1+const2/temp^2)-const1)
    M[,subset] = M[,subset] - f/(-sexp*(1+const2[,subset]/temp^2)-const1[,subset])
    converged = apply(abs(f),2,function(x){max(x)<tol})
    subset = subset[!converged]
    if(all(converged))
      break
  }
  #gc()
  if(return_V){
    return(list(M=M,V=Sigma2/(const0-M)))
  }else{
    return(M)
  }

}
