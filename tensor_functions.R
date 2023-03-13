#### data generation ####

#' generate a tensor
#' @param n_1 dimension 1 of the tensor
#' @param n_2 dimension 2 of the tensor
#' @param d dimension of the latent space
#' @param L dimension 3 of the tensor, i.e., number of layers
#' @return (n_1, n_2, L)-shaped tensor
#' @export 
generate_tensor_dirichlet <- function(n_1, n_2, L, W,
                                      dirichlet_x, dirichlet_y){
  # n_1 = dim_[1]
  # n_2 = dim_[2]
  # d = dim_[3]
  # L = dim_[4]
  dim_ = c(n_1, n_2, L)
  A = array(NA,dim_)
  probability = array(NA,dim_)
  
  for (layer in 1: L)
  {
    temp_1 = rdirichlet(n_1, dirichlet_x)
    temp_2 = rdirichlet(n_2, dirichlet_y)
    P =  temp_1  %*% W[, , layer] %*% t(temp_2)
    probability[, , layer] = P
    
    A[, , layer] = matrix(rbinom(matrix(1,n_1,n_2),matrix(1,n_1,n_2), probability[ , , layer]),n_1,n_2)
  }
  hat.rank =c(10,10,10)
  
  Y.tensor =  as.tensor(A)
  
  return(list(Y.tensor = Y.tensor, probability = probability, hat.rank = hat.rank))
}


generate_tensor_probability <- function(n_1, n_2, L, probability){
  
  dim_ = c(n_1, n_2, L)
  A = array(NA,dim_)
  
  for (layer in 1: L)
  {
    A[, , layer] = matrix(rbinom(matrix(1,n_1,n_2),matrix(1,n_1,n_2), probability[ , , layer]),n_1,n_2)
  }
  hat.rank =c(10,10,10)
  
  Y.tensor =  as.tensor(A)
  
  return(list(Y.tensor = Y.tensor, probability = probability, hat.rank = hat.rank))
}




#### utils ####

frobenius <- function(A, B){
  return (sum((A - B)^2)^0.5)
}



my_svd <- function (x, nu = min(n, p), nv = min(n, p), LINPACK = TRUE){
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
  return(res)
}





HOSVD_test <- function(Y){
  p = dim(Y)
  d = length(p)
  U_0 = list()
  for (i in 1:d){
    MY = k_unfold(Y, i)@data
    temp = svd( MY) 
    U_t = temp$u[,1]
    U_t = matrix(U_t, nrow  = length(U_t))
    U_0 = c(U_0,list(U_t))
  }
  return(U_0)
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






#### estimation ####



estimate_flex <- function(A, dim_latent, TT = 200, eta_u = 5e-4,
                          eta_a = 5e-4,
                          eta_l = 5e-4, err_tol = 1e-3){
  
  sigmoid <- function(x){
    return (1 / (1 + exp(-x)))
  }
  
  n = dim(A)[1]
  I_n = diag(rep(1, n))
  J_n = I_n - matrix(1, nrow = n, ncol = n) / n
  one_n = matrix(1, nrow = n, ncol = 1)
  
  k = dim_latent
  R = dim(A)[3]
  U_t = matrix(rnorm(n * k, 0, 1), nrow = n, ncol = k)
  Lambda_t = array(0, c(R, k, k))
  Lambda_t1 = array(0, c(R, k, k))
  for (r in 1:R){
    Lambda_t[r, , ] = diag(rep(1, k))
  }
  alpha_t = matrix(0, nrow = n, ncol = R)
  alpha_t1 = matrix(0, nrow = n, ncol = R)
  P_t = array(0, c(n, n, R))
  for (r in 1:R){
    P_t[, , r] = sigmoid(alpha_t[, r] %*% t(one_n) + one_n %*% t(alpha_t[, r]) + U_t %*% Lambda_t[r, , ] %*% t(U_t))
  }
  P_t1 = array(0, c(n, n, R))
  
  
  # TT = 200
  # eta_u = 5e-4
  # eta_a = 5e-4
  # eta_l = 5e-4
  diff = 1e2
  
  # err_tol = 1e-1
  
  # true_error_list = rep(0, TT)
  # diff_record = rep(0, TT)
  for (t in 1:TT){
    U_sum = 0
    for (r in 1:R){
      A_sigma = A[, , r] - P_t[, , r]
      U_sum = U_sum + A_sigma %*% U_t %*% Lambda_t[r, , ]
      alpha_t1[, r] = alpha_t[, r] + 2 * eta_a * A_sigma %*% one_n
      Lambda_t1[r, , ] = Lambda_t[r, , ] + eta_l * t(U_t) %*% A_sigma %*% U_t
    }
    
    U_t1 = U_t + 2 * eta_u * U_sum
    U_t1 = J_n %*% U_t1
    
    G = t(U_t1) %*% U_t1
    tmp = my_svd(G)
    W = n^0.5 * tmp$u %*% diag(1 / (tmp$d)^0.5)
    U_t1 = U_t1 %*% W
    # t(U_t1) %*% U_t1
    
    W_inv = 1/n^0.5 * diag((tmp$d)^0.5) %*% t(tmp$u)
    # W_inv %*% W
    for (r in 1:R){
      Lambda_t1[r, , ] = t(W_inv) %*% Lambda_t1[r, , ] %*% W_inv
      P_t1[, , r] = sigmoid(alpha_t1[, r] %*% t(one_n) + one_n %*% t(alpha_t1[, r]) + U_t1 %*% Lambda_t1[r, , ] %*% t(U_t1))
    }
    
    
    diff = frobenius(P_t, P_t1)
    # diff_record[t] = diff
    if (diff < err_tol){
      break
    }
    
    #### update variables
    U_t = U_t1
    alpha_t = alpha_t1
    Lambda_t = Lambda_t1
    P_t = P_t1
    
    # true_error_list[t] = frobenius(P_t1, probability_1)
  }
  return (P_t1)
}


estimate_thpca <- function(Y.tensor, hat.rank){
  U.hat = Tensor_Hetero_PCA_test(Y.tensor, hat.rank)
  P.U1 = U.hat[[1]]%*%t(U.hat[[1]])
  P.U2 = U.hat[[2]]%*%t(U.hat[[2]])
  P.U3 = U.hat[[3]]%*%t(U.hat[[3]])
  Y.hat = ttm(ttm(ttm(Y.tensor, P.U1, 1), P.U2, 2), P.U3, 3) 
  
  P_hat  = Y.hat@data
  
  P_hat[P_hat > 1]  = 1
  P_hat[P_hat < 0]  = 0
  
  return(P_hat)
}

estimate_hosvd <- function(Y.tensor, hat.rank){
  U.hat = rTensor::hosvd(Y.tensor, hat.rank)$U
  P.U1 = U.hat[[1]]%*%t(U.hat[[1]])
  P.U2 = U.hat[[2]]%*%t(U.hat[[2]])
  P.U3 = U.hat[[3]]%*%t(U.hat[[3]])
  Y.hat.3 = ttm(ttm(ttm(Y.tensor, P.U1, 1), P.U2, 2), P.U3, 3) 
  
  P_hat_3  = Y.hat.3@data
  
  P_hat_3[P_hat_3 > 1]  = 1
  P_hat_3[P_hat_3 < 0]  = 0
  
  return(P_hat_3)
}


estimate_hooi <- function(Y.tensor, hat.rank){
  U.hat = HOOI(Y.tensor, hat.rank)
  P.U1 = U.hat[[1]]%*%t(U.hat[[1]])
  P.U2 = U.hat[[2]]%*%t(U.hat[[2]])
  P.U3 = U.hat[[3]]%*%t(U.hat[[3]])
  Y.hat.1 = ttm(ttm(ttm(Y.tensor, P.U1, 1), P.U2, 2), P.U3, 3) 
  
  P_hat_1  = Y.hat.1@data
  
  P_hat_1[P_hat_1 > 1]  = 1
  P_hat_1[P_hat_1 < 0]  = 0
  
  return(P_hat_1)
}

estimate_svd <- function(A, hat.rank){
  L = dim(A)[3]
  P_hat_4 = array(NA,dim(A))
  
  for (layer in 1: L){
    temp = my_svd(A[, , layer])
    xhat =  temp$u[,1:hat.rank[1]] %*%  diag( sqrt(temp$d[1:hat.rank[1]]) )
    yhat =  temp$v[,1:hat.rank[1]] %*%  diag( sqrt(temp$d[1:hat.rank[1]]) )
    P_hat_4[, , layer] = xhat %*% t(yhat)
  }
  
  P_hat_4[P_hat_4 > 1]  = 1
  P_hat_4[P_hat_4 < 0]  = 0
  
  return(P_hat_4)
}







uase <- function(A_tensor, rank){
  
  Y.tensor = A_tensor
  dim =dim(Y.tensor)
  Y.matrix = t(cs_unfold(Y.tensor, 1)@data)
  temp.Y = my_svd(Y.matrix)
  xhat.Y =  temp.Y$u[,1:rank] %*%  diag( sqrt(temp.Y$d[1:rank]) )
  yhat.Y =  temp.Y$v[,1:rank] %*%  diag( sqrt(temp.Y$d[1:rank]) )
  
  probability.uase = xhat.Y %*% t(yhat.Y) 
  
  probability.uase[probability.uase > 1]  = 1
  
  probability.uase[probability.uase < 0]  = 0
  
  probability.uase =  array(probability.uase, c(dim)) 
  
  return(probability.uase)
  
}



hetero_pca_estimate <- function(Y.tensor, hat.rank){
  U.hat = Tensor_Hetero_PCA_test(Y.tensor, hat.rank)
  P.U1 = U.hat[[1]] %*% t(U.hat[[1]])
  P.U2 = U.hat[[2]] %*% t(U.hat[[2]])
  P.U3 = U.hat[[3]] %*% t(U.hat[[3]])
  Y.hat.2 = ttm(ttm(ttm(Y.tensor, P.U1, 1), P.U2, 2), P.U3, 3) 
  
  probability.2 = Y.hat.2@data
  
  probability.2[probability.2 > 1] = 1
  
  probability.2[probability.2 < 0] = 0
  
  return(probability.2)
}

hosvd_pca_estimate <- function(Y.tensor){
  U.hat = HOSVD_test(Y.tensor)
  P.U1 = U.hat[[1]]%*%t(U.hat[[1]])
  P.U2 = U.hat[[2]]%*%t(U.hat[[2]])
  P.U3 = U.hat[[3]]%*%t(U.hat[[3]])
  Y.hat.3 = ttm(ttm(ttm(Y.tensor, P.U1, 1), P.U2, 2), P.U3, 3) 
  
  probability.3  = Y.hat.3@data
  
  probability.3[probability.3 > 1] = 1
  
  probability.3[probability.3 < 0] = 0
  
  return(probability.3)
}



#### multiness

multiness_adaptive <- function(A_tensor, rank){
  dimm = dim(A_tensor)
  probability.multi = array(NA, dimm)
  fit1 <- multiness_fit(A=A_tensor,
                        self_loops=TRUE,
                        refit=TRUE,
                        tuning="adaptive",
                        tuning_opts=list(layer_wise=FALSE))
  for (layer in 1: dimm[3]){
    probability.multi[, , layer] = fit1$F_hat + fit1$G_hat[[layer]]
    
  }
  probability.multi[probability.multi > 1]  = 1
  probability.multi[probability.multi < 0]  = 0
  return(probability.multi)
}

multiness <- function(A_tensor, rank){
  dimm = dim(A_tensor)
  probability.multi = array(NA, dimm)
  fit1 <- multiness_fit(A=A_tensor,
                        self_loops=TRUE,
                        refit=FALSE,
                        tuning="fixed",
                        tuning_opts=list(lambda=40,alpha=1/2),
                        optim_opts=list(max_rank=rank,verbose=FALSE))
  for (layer in 1: dimm[3]){
    probability.multi[, , layer] = fit1$F_hat + fit1$G_hat[[layer]]
    
  }
  probability.multi[probability.multi > 1]  = 1
  probability.multi[probability.multi < 0]  = 0
  return(probability.multi)
}


