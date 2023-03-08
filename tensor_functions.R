

svd <- function (x, nu = min(n, p), nv = min(n, p), LINPACK = TRUE)
{
  x <- as.matrix(x)
  if (any(!is.finite(x)))
    stop("infinite or missing values in 'x'")
  dx <- dim(x)
  n <- dx[1L]
  p <- dx[2L]
  if (!n || !p)
    stop("a dimension is zero")
  La.res <- La.svd(x, nu, nv)   ## your problem line
  res <- list(d = La.res$d)
  if (nu)
    res$u <- La.res$u
  if (nv) {
    if (is.complex(x))
      res$v <- Conj(t(La.res$vt))
    else res$v <- t(La.res$vt)
  }
  res
}



Tensor_Hetero_PCA_test <- function(Y, r){
  
  p = dim(Y)
  d = length(p)
  U_0 = list()
  for (i in 1:d){
    MY = k_unfold(Y, i)@data
    MY_Hetero_PCA = Hetero_PCA_test(MY %*% t(MY), r[i])
    U_0 = c(U_0, list(MY_Hetero_PCA))
  }
  return(U_0)
}

Hetero_PCA_test <- function(Y, r, tmax, vartol){
  try(if(missing("tmax")) tmax = 20)
  try(if(missing("vartol")) vartol = 1e-6)
  
  N_t = Y
  r = min(c(r, dim(N_t)))
  diag(N_t) = 0
  U_t = matrix(NA, nrow = dim(Y)[1], r)
  t = 1
  approx = -1
  
  while(t<tmax){ # Stop criterion: convergence or maximum number of iteration reached
     temp = svd(N_t) 
     U_t = temp$u[,1:r]
     tilde_N_test_t = temp$u[,1:r] %*% diag(temp$d[1:r]) %*% t(temp$v[,1:r])
     N_test_new =   N_t 
     diag(N_test_new) = diag(tilde_N_test_t)
     N_t =  N_test_new
     svector = diag(tilde_N_test_t)
  if (abs(sum(svector^2) - approx) > vartol){
    t = t+1
    approx = sum(svector^2)
  }
  else {
    break
  }
}
  return(U_t)
}

