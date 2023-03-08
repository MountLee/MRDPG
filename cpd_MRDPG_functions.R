#### data generation ####


#' Calculate the sine-theta distance between two subspaces
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

hosvd_pca_estimate <- function(Y.tensor)
{
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


#### cpd ####

get_D_s_t_R <- function(P_left, P_right, Z, Sigma, S_tilde_p){
  n_1 = dim(P_left)[1]
  n_2 = dim(P_left)[2]
  L = dim(P_left)[3]
  M_t = dim(Z)[1]
  
  D_K_t_left = rep(0, M_t)
  D_K_t_right = rep(0, M_t)
  
  for (p in 1:(n_1 + n_2 - 1)){
    cur_P_sp_left = rep(0, L)
    cur_P_sp_right = rep(0, L)
    if (p <= n_2){
      for (i in 1:min(n_2 + 1 - p, n_1)){
        cur_P_sp_left = cur_P_sp_left + c(P_left[i, i + p - 1, ])
        cur_P_sp_right = cur_P_sp_right + c(P_right[i, i + p - 1, ])  
      }
    }
    else{
      for (i in 1:min(n_1 + n_2 - p, n_2)){
        cur_P_sp_left = cur_P_sp_left + c(P_left[i + p - n_2, i, ])
        cur_P_sp_right = cur_P_sp_right + c(P_right[i + p - n_2, i, ])  
      }
    }
    cur_P_sp_left = cur_P_sp_left / S_tilde_p[p]
    cur_P_sp_right = cur_P_sp_right / S_tilde_p[p]
    
    D_K_t_left = D_K_t_left + dmvnorm(Z, mean = cur_P_sp_left,
                                      sigma = diag(L) / h_kernel) * (S_tilde_p[p] / (n_1 * n_2))
    D_K_t_right = D_K_t_right + dmvnorm(Z, mean = cur_P_sp_right,
                                        sigma = diag(L) / h_kernel) * (S_tilde_p[p] / (n_1 * n_2))
   
  }
  
  return (abs(D_K_t_left - D_K_t_right))
}

max_D_s_t_rescale <- function(A_list, h_kernel, verbose = TRUE){
  n_1 = dim(A_list)[1]
  n_2 = dim(A_list)[2]
  L = dim(A_list)[3]
  TT = dim(A_list)[4]
  
  C_M = 20
  M_up = 10^3
  
  Sigma = diag(L) / h_kernel
  
  D_K_t_max_rescale = rep(0, TT)
  
  
  for (t in 2:TT){
    D_K_t = rep(0, t)
    
    for (s in 1:(t - 1)){
      
      A_sum_left = apply(A_list[ , , , c(1:s)], c(1,2,3), sum)
      A_sum_right = apply(A_list[ , , , c(min(s+1, t):t)], c(1,2,3), sum)
      
      n_left = s
      n_right = t - s
      
      P_left = hosvd_pca_estimate(as.tensor(A_sum_left / n_left))
      P_right = hosvd_pca_estimate(as.tensor(A_sum_right / n_right))
      
      M_t = C_M * ((t*min(n_1, n_2)) / (2 * L^2 * log(max(c(n_1, n_2, t)))))^(L / 2)
      M_t = as.integer(M_t)
      M_t = min(M_t, M_up)
      
      # each column is a sample
      Z = matrix(runif(M_t * L, 0, 1), nrow = M_t, ncol = L)
      
      D_K_t[s] = max(get_D_s_t_R(P_left, P_right, Z, Sigma, S_tilde_p)) / (1/s^0.5 + 1/(t-s)^0.5) / log(max(c(n_1, n_2, t)))^0.5
    }
    
    D_K_t_max_rescale[t] = max(D_K_t)
    if (verbose){
      print(t) 
    }
  }
  return(max(D_K_t_max_rescale))
}

online_cpd <- function(A_list, tau_factor, h_kernel, verbose = TRUE){
  n_1 = dim(A_list)[1]
  n_2 = dim(A_list)[2]
  L = dim(A_list)[3]
  TT = dim(A_list)[4]
  
  C_M = 20
  M_up = 10^3
  
  Sigma = diag(L) / h_kernel
  
  D_K_t_max_rescale = rep(0, TT)
  
  
  for (t in 2:TT){
    
    D_K_t = rep(0, t)
    
    for (s in 1:(t - 1)){
      
      A_sum_left = apply(A_list[ , , , c(1:s)], c(1,2,3), sum)
      A_sum_right = apply(A_list[ , , , c(min(s+1, t):t)], c(1,2,3), sum)
      
      n_left = s
      n_right = t - s
      
      P_left = hosvd_pca_estimate(as.tensor(A_sum_left / n_left))
      P_right = hosvd_pca_estimate(as.tensor(A_sum_right / n_right))
      
      M_t = C_M * ((t*min(n_1, n_2)) / (2 * L^2 * log(max(c(n_1, n_2, t)))))^(L / 2)
      M_t = as.integer(M_t)
      M_t = min(M_t, M_up)
      
      # each column is a sample
      Z = matrix(runif(M_t * L, 0, 1), nrow = M_t, ncol = L)
      
      D_K_t[s] = max(get_D_s_t_R(P_left, P_right, Z, Sigma, S_tilde_p)) / (1/s^0.5 + 1/(t-s)^0.5) / log(max(c(n_1, n_2, t)))^0.5
    }
    
    
    s_max = which.max(D_K_t)
    D_K_t_max_rescale[t] = max(D_K_t) 
    
    if (D_K_t_max_rescale[t] > tau_factor){
      return (list(t = t, D_K_t_max_rescale = D_K_t_max_rescale, find = TRUE))
    }
    
    if (verbose){
      print(t)
    }
  }  
  # didn't find a t with D_t exceeding tau_factor
  return (list(t = t, D_K_t_max_rescale = D_K_t_max_rescale, find = FALSE))
}

max_D_s_t_rescale_thpca <- function(A_list, h_kernel, rank, verbose = TRUE){
  n_1 = dim(A_list)[1]
  n_2 = dim(A_list)[2]
  L = dim(A_list)[3]
  TT = dim(A_list)[4]
  
  C_M = 20
  M_up = 10^3
  
  Sigma = diag(L) / h_kernel
  
  D_K_t_max_rescale = rep(0, TT)
  
  
  for (t in 2:TT){
    D_K_t = rep(0, t)
    
    for (s in 1:(t - 1)){
      
      A_sum_left = apply(A_list[ , , , c(1:s)], c(1,2,3), sum)
      A_sum_right = apply(A_list[ , , , c(min(s+1, t):t)], c(1,2,3), sum)
      
      n_left = s
      n_right = t - s
      
      P_left = hetero_pca_estimate(as.tensor(A_sum_left / n_left), rank)
      P_right = hetero_pca_estimate(as.tensor(A_sum_right / n_right), rank)
      
      M_t = C_M * ((t*min(n_1, n_2)) / (2 * L^2 * log(max(c(n_1, n_2, t)))))^(L / 2)
      M_t = as.integer(M_t)
      M_t = min(M_t, M_up)
      
      # each column is a sample
      Z = matrix(runif(M_t * L, 0, 1), nrow = M_t, ncol = L)
      
      D_K_t[s] = max(get_D_s_t_R(P_left, P_right, Z, Sigma, S_tilde_p)) / (1/s^0.5 + 1/(t-s)^0.5) / log(max(c(n_1, n_2, t)))^0.5
    }
    
    D_K_t_max_rescale[t] = max(D_K_t)
    if (verbose){
      print(t) 
    }
  }
  return(max(D_K_t_max_rescale))
}

online_cpd_thpca <- function(A_list, tau_factor, h_kernel, rank, verbose = TRUE){
  n_1 = dim(A_list)[1]
  n_2 = dim(A_list)[2]
  L = dim(A_list)[3]
  TT = dim(A_list)[4]
  
  C_M = 20
  M_up = 10^3
  
  Sigma = diag(L) / h_kernel
  
  D_K_t_max_rescale = rep(0, TT)
  
  
  for (t in 2:TT){
    
    D_K_t = rep(0, t)
    
    for (s in 1:(t - 1)){
      
      A_sum_left = apply(A_list[ , , , c(1:s)], c(1,2,3), sum)
      A_sum_right = apply(A_list[ , , , c(min(s+1, t):t)], c(1,2,3), sum)
      
      n_left = s
      n_right = t - s
      
      P_left = hetero_pca_estimate(as.tensor(A_sum_left / n_left), rank)
      P_right =hetero_pca_estimate(as.tensor(A_sum_right / n_right), rank)
      
      M_t = C_M * ((t*min(n_1, n_2)) / (2 * L^2 * log(max(c(n_1, n_2, t)))))^(L / 2)
      M_t = as.integer(M_t)
      M_t = min(M_t, M_up)
      
      # each column is a sample
      Z = matrix(runif(M_t * L, 0, 1), nrow = M_t, ncol = L)
      
      D_K_t[s] = max(get_D_s_t_R(P_left, P_right, Z, Sigma, S_tilde_p)) / (1/s^0.5 + 1/(t-s)^0.5) / log(max(c(n_1, n_2, t)))^0.5
    }
    
    
    s_max = which.max(D_K_t)
    D_K_t_max_rescale[t] = max(D_K_t) 
    
    if (D_K_t_max_rescale[t] > tau_factor){
      return (list(t = t, D_K_t_max_rescale = D_K_t_max_rescale, find = TRUE))
    }
    
    if (verbose){
      print(t)
    }
  }  
  # didn't find a t with D_t exceeding tau_factor
  return (list(t = t, D_K_t_max_rescale = D_K_t_max_rescale, find = FALSE))
}

uase <- function(A_tensor, rank){
  
  Y.tensor = A_tensor
  dim =dim(Y.tensor)
  Y.matrix = t(cs_unfold(Y.tensor, 1)@data)
  temp.Y = svd(Y.matrix)
  xhat.Y =  temp.Y$u[,1:rank] %*%  diag( sqrt(temp.Y$d[1:rank]) )
  yhat.Y =  temp.Y$v[,1:rank] %*%  diag( sqrt(temp.Y$d[1:rank]) )
  
  probability.uase = xhat.Y %*% t(yhat.Y) 
  
  probability.uase[probability.uase > 1]  = 1
  
  probability.uase[probability.uase < 0]  = 0
  
  probability.uase =  array(probability.uase, c(dim)) 
  
  return(probability.uase)
  
}

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

max_D_s_t_rescale_uase <- function(A_list, h_kernel, rank, verbose = TRUE){
  n_1 = dim(A_list)[1]
  n_2 = dim(A_list)[2]
  L = dim(A_list)[3]
  TT = dim(A_list)[4]
  
  C_M = 20
  M_up = 10^3
  
  Sigma = diag(L) / h_kernel
  
  D_K_t_max_rescale = rep(0, TT)
  
  
  for (t in 2:TT){
    D_K_t = rep(0, t)
    
    for (s in 1:(t - 1)){
      
      A_sum_left = apply(A_list[ , , , c(1:s)], c(1,2,3), sum)
      A_sum_right = apply(A_list[ , , , c(min(s+1, t):t)], c(1,2,3), sum)
      
      n_left = s
      n_right = t - s
      
      P_left = uase(as.tensor(A_sum_left / n_left), rank)
      P_right = uase(as.tensor(A_sum_right / n_right), rank)
      
      M_t = C_M * ((t*min(n_1, n_2)) / (2 * L^2 * log(max(c(n_1, n_2, t)))))^(L / 2)
      M_t = as.integer(M_t)
      M_t = min(M_t, M_up)
      
      # each column is a sample
      Z = matrix(runif(M_t * L, 0, 1), nrow = M_t, ncol = L)
      
      D_K_t[s] = max(get_D_s_t_R(P_left, P_right, Z, Sigma, S_tilde_p)) / (1/s^0.5 + 1/(t-s)^0.5) / log(max(c(n_1, n_2, t)))^0.5
    }
    
    D_K_t_max_rescale[t] = max(D_K_t)
    if (verbose){
      print(t) 
    }
  }
  return(max(D_K_t_max_rescale))
}


online_cpd_uase <- function(A_list, tau_factor, h_kernel, rank, verbose = TRUE){
  n_1 = dim(A_list)[1]
  n_2 = dim(A_list)[2]
  L = dim(A_list)[3]
  TT = dim(A_list)[4]
  
  C_M = 20
  M_up = 10^3
  
  Sigma = diag(L) / h_kernel
  
  D_K_t_max_rescale = rep(0, TT)
  
  
  for (t in 2:TT){
    
    D_K_t = rep(0, t)
    
    for (s in 1:(t - 1)){
      
      A_sum_left = apply(A_list[ , , , c(1:s)], c(1,2,3), sum)
      A_sum_right = apply(A_list[ , , , c(min(s+1, t):t)], c(1,2,3), sum)
      
      n_left = s
      n_right = t - s
      
      P_left = uase(as.tensor(A_sum_left / n_left), rank)
      P_right = uase(as.tensor(A_sum_right / n_right), rank)
      
      M_t = C_M * ((t*min(n_1, n_2)) / (2 * L^2 * log(max(c(n_1, n_2, t)))))^(L / 2)
      M_t = as.integer(M_t)
      M_t = min(M_t, M_up)
      
      # each column is a sample
      Z = matrix(runif(M_t * L, 0, 1), nrow = M_t, ncol = L)
      
      D_K_t[s] = max(get_D_s_t_R(P_left, P_right, Z, Sigma, S_tilde_p)) / (1/s^0.5 + 1/(t-s)^0.5) / log(max(c(n_1, n_2, t)))^0.5
    }
    
    
    s_max = which.max(D_K_t)
    D_K_t_max_rescale[t] = max(D_K_t) 
    
    if (D_K_t_max_rescale[t] > tau_factor){
      return (list(t = t, D_K_t_max_rescale = D_K_t_max_rescale, find = TRUE))
    }
    
    if (verbose){
      print(t)
    }
  }  
  # didn't find a t with D_t exceeding tau_factor
  return (list(t = t, D_K_t_max_rescale = D_K_t_max_rescale, find = FALSE))
}

max_D_s_t_rescale_multi <- function(A_list, h_kernel, rank, verbose = TRUE){
  n_1 = dim(A_list)[1]
  n_2 = dim(A_list)[2]
  L = dim(A_list)[3]
  TT = dim(A_list)[4]
  
  C_M = 20
  M_up = 10^3
  
  Sigma = diag(L) / h_kernel
  
  D_K_t_max_rescale = rep(0, TT)
  
  
  for (t in 2:TT){
    D_K_t = rep(0, t)
    
    for (s in 1:(t - 1)){
      
      A_sum_left = apply(A_list[ , , , c(1:s)], c(1,2,3), sum)
      A_sum_right = apply(A_list[ , , , c(min(s+1, t):t)], c(1,2,3), sum)
      
      n_left = s
      n_right = t - s
      
      P_left = multiness_adaptive(A_sum_left / n_left, rank)
      P_right =multiness_adaptive(A_sum_right / n_right, rank)
      
      M_t = C_M * ((t*min(n_1, n_2)) / (2 * L^2 * log(max(c(n_1, n_2, t)))))^(L / 2)
      M_t = as.integer(M_t)
      M_t = min(M_t, M_up)
      
      # each column is a sample
      Z = matrix(runif(M_t * L, 0, 1), nrow = M_t, ncol = L)
      
      D_K_t[s] = max(get_D_s_t_R(P_left, P_right, Z, Sigma, S_tilde_p)) / (1/s^0.5 + 1/(t-s)^0.5) / log(max(c(n_1, n_2, t)))^0.5
    }
    
    D_K_t_max_rescale[t] = max(D_K_t)
    if (verbose){
      print(t) 
    }
  }
  return(max(D_K_t_max_rescale))
}

online_cpd_multi <- function(A_list, tau_factor, h_kernel, rank, verbose = TRUE){
  n_1 = dim(A_list)[1]
  n_2 = dim(A_list)[2]
  L = dim(A_list)[3]
  TT = dim(A_list)[4]
  
  C_M = 20
  M_up = 10^3
  
  Sigma = diag(L) / h_kernel
  
  D_K_t_max_rescale = rep(0, TT)
  
  
  for (t in 2:TT){
    
    D_K_t = rep(0, t)
    
    for (s in 1:(t - 1)){
      
      A_sum_left = apply(A_list[ , , , c(1:s)], c(1,2,3), sum)
      A_sum_right = apply(A_list[ , , , c(min(s+1, t):t)], c(1,2,3), sum)
      
      n_left = s
      n_right = t - s
      
      P_left = multiness_adaptive(A_sum_left / n_left, rank)
      P_right = multiness_adaptive(A_sum_right / n_right, rank)
      
      M_t = C_M * ((t*min(n_1, n_2)) / (2 * L^2 * log(max(c(n_1, n_2, t)))))^(L / 2)
      M_t = as.integer(M_t)
      M_t = min(M_t, M_up)
      
      # each column is a sample
      Z = matrix(runif(M_t * L, 0, 1), nrow = M_t, ncol = L)
      
      D_K_t[s] = max(get_D_s_t_R(P_left, P_right, Z, Sigma, S_tilde_p)) / (1/s^0.5 + 1/(t-s)^0.5) / log(max(c(n_1, n_2, t)))^0.5
    }
    
    
    s_max = which.max(D_K_t)
    D_K_t_max_rescale[t] = max(D_K_t) 
    
    if (D_K_t_max_rescale[t] > tau_factor){
      return (list(t = t, D_K_t_max_rescale = D_K_t_max_rescale, find = TRUE))
    }
    
    if (verbose){
      print(t)
    }
  }  
  # didn't find a t with D_t exceeding tau_factor
  return (list(t = t, D_K_t_max_rescale = D_K_t_max_rescale, find = FALSE))
}



##### functions #####

distance_Y <- function(Y1, Y2){
  dist = sum((Y1 - Y2)^2)
  return(dist)
}


get_dist_mat <- function(A_list, distance_func = distance_A){
  n_ob = dim(A_list)[4]
  
  dist_mat = matrix(0, n_ob, n_ob)
  for (i in 1:(n_ob - 1)){
    for (j in (i + 1):n_ob){
      dist_mat[i, j] = distance_func(A_list[, , , i], A_list[ , , , j])
      dist_mat[j, i] = dist_mat[i, j]
    }
  }
  
  return(dist_mat)
}





