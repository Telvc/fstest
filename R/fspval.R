fspval <- function(X, kernel, rho = 0.99, p, cl_fun, rectime, toltime, lambd = 1e-5, ndraws=2000, seed = 100){
  #compute the dimension of the input dataset
  m = dim(X)[1]
  n_curve = dim(X)[2]

  eigen = c()
  for(i in seq(p)){
    eigen = c(eigen, (1-rho)*rho**(i-1))
  }
  core_mat = diag(eigen)

  # kernel ridge regression
  tol_sam = array(0, dim = c(m,n_curve,p))
  for(fea in 1:n_curve){
    datcur = X[,fea,]
    cursam = array(0, dim = c(m, p))

    for(id in 1:m){
      dcur = na.omit(datcur[id,])
      time = na.omit(rectime[id,fea,]/toltime)
      K = length(time)
      val = as.matrix(dcur)

      kercur = kernel(time,K,rho)
      kernel_f = matrix(rep(0,p*K),p,K)
      for(i in seq(p)){
        kernel_f[i,] = f(i-1,time)
      }
      mat1 = t(kernel_f)%*%core_mat
      inv1 = solve(kercur+lambd*diag(K))

      coef = t(val)%*%inv1
      cursam[id,] = coef%*%mat1
    }

    tol_sam[,fea,] = cursam
  }


  # estimate the covariance matrix and whitten the data
  whit_vec = matrix(tol_sam, m, n_curve*p)
  sam_cov = solve(sqrtm(cov(whit_vec)))
  whit_vec = whit_vec%*%sam_cov

  # compute the selective p-value

  subset = sample(1:m,400,replace = FALSE)
  cl = km_cluster(whit_vec[subset,])
  pva = test_clusters_approx_tensor(whit_vec[subset,], k1=1, k2=2, iso=TRUE, sig=NULL, SigInv=NULL, cl_fun=km_cluster, cl = cl, ndraws)

  return(pva$pval)
}
