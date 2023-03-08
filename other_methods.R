#### other cpd methods


#### multi ####


multiness <- function(A_tensor, rank){
  dimm = dim(A_tensor)
  probability.multi = array(NA, dimm)
  fit1 <- multiness_fit(A=A_tensor,
                        self_loops=TRUE,
                        refit=FALSE,
                        tuning="fixed",
                        tuning_opts=list(lambda=40, alpha=1/2),
                        optim_opts=list(max_rank=rank,verbose=FALSE))
  for (layer in 1: dimm[3]){
    probability.multi[, , layer] = fit1$F_hat + fit1$G_hat[[layer]]
    
  }
  probability.multi[probability.multi > 1]  = 1
  probability.multi[probability.multi < 0]  = 0
  return(probability.multi)
}


max_D_s_t_rescale_multi <- function(A_list, h_kernel = 0.1, rank = 10, verbose = FALSE){
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
      
      P_left = multiness(A_sum_left / n_left, rank)
      P_right = multiness(A_sum_right / n_right, rank)
      
      M_t = C_M * ((t*min(n_1, n_2)) / (2 * L^2 * log(max(c(n_1, n_2, t)))))^(L / 2)
      M_t = as.integer(M_t)
      M_t = min(M_t, M_up)
      
      # each column is a sample
      Z = matrix(runif(M_t * L, 0, 1), nrow = M_t, ncol = L)
      
      D_K_t[s] = get_D_K_t_cpp(P_left, P_right, Z, Sigma) / (1/s^0.5 + 1/(t-s)^0.5) / log(max(c(n_1, n_2, t)))^0.5
    }
    
    D_K_t_max_rescale[t] = max(D_K_t)
    if (verbose){
      print(t) 
    }
  }
  return(max(D_K_t_max_rescale))
}

online_cpd_multi <- function(A_list, tau_factor, h_kernel = 0.1, rank = 10, verbose = FALSE){
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
      
      P_left = multiness(A_sum_left / n_left, rank)
      P_right = multiness(A_sum_right / n_right, rank)
      
      M_t = C_M * ((t*min(n_1, n_2)) / (2 * L^2 * log(max(c(n_1, n_2, t)))))^(L / 2)
      M_t = as.integer(M_t)
      M_t = min(M_t, M_up)
      
      # each column is a sample
      Z = matrix(runif(M_t * L, 0, 1), nrow = M_t, ncol = L)
      
      D_K_t[s] = get_D_K_t_cpp(P_left, P_right, Z, Sigma) / (1/s^0.5 + 1/(t-s)^0.5) / log(max(c(n_1, n_2, t)))^0.5
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



#### knn ####

online_cpd_knn <- function(A_b, A_list, k_nn, alpha, n0, n1, ARL){
  n_1 = dim(A_list)[1]
  n_2 = dim(A_list)[2]
  L = dim(A_list)[3]
  TT = dim(A_list)[4]
  T_burn = dim(A_b)[4]
  
  A_all = array(0, c(n_1, n_2, L, TT + T_burn))
  A_all[, , , 1:T_burn] = A_b
  A_all[, , , (T_burn + 1):(T_burn + TT)] = A_list
  
  n_all = dim(A_all)[4]
  
  dist_mat_all = get_dist_mat(A_all, distance_Y)
  diag(dist_mat_all) = max(dist_mat_all) + 100
  
  t_hat = tryCatch(
    expr = {
      res = gstream(dist_mat_all, L = T_burn, N0 = T_burn, k_nn, statistics = "m",
                    n0, n1, ARL, alpha, skew.corr = TRUE, asymp = FALSE)
      tau_hat = res$tauhat$max.type
      if(length(tau_hat) == 0){
        tau_hat = c(TT + 1)
      }
      min(tau_hat)
    },
    error = function(cond){return (TT + 2)}
  )
  # res = gstream(dist_mat_all, L = T_burn, N0 = T_burn, k_nn, statistics = "m",
  #               n0, n1, ARL, alpha, skew.corr = TRUE, asymp = FALSE)
  # tau_hat = res$tauhat$max.type
  # if(length(tau_hat) == 0){
  #   tau_hat = c(TT + 1)
  # }
  # t_hat = min(tau_hat)
  find_knn = t_hat <= TT
  
  # didn't find a t with D_t exceeding tau_factor
  return (list(t = t_hat, find = find_knn))
}



distance_Y <- function(Y1, Y2){
  dist = sum((Y1 - Y2)^2)
  return(dist)
}


get_dist_mat <- function(A_list, distance_func = distance_Y){
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







#### uase ####


uase <- function(A_tensor, rank){
  
  dim =dim(A_tensor)
  Y.matrix = t(cs_unfold(A_tensor, 1)@data)
  temp.Y = svd(Y.matrix)
  xhat.Y =  temp.Y$u[,1:rank] %*%  diag( sqrt(temp.Y$d[1:rank]) )
  yhat.Y =  temp.Y$v[,1:rank] %*%  diag( sqrt(temp.Y$d[1:rank]) )
  
  probability.uase = xhat.Y %*% t(yhat.Y) 
  
  probability.uase[probability.uase > 1]  = 1
  
  probability.uase[probability.uase < 0]  = 0
  
  probability.uase =  array(probability.uase, c(dim)) 
  
  return(probability.uase)
  
}



max_D_s_t_rescale_uase <- function(A_list, h_kernel = 0.1, rank = 10, verbose = FALSE){
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
      
      D_K_t[s] = get_D_K_t_cpp(P_left, P_right, Z, Sigma) / (1/s^0.5 + 1/(t-s)^0.5) / log(max(c(n_1, n_2, t)))^0.5
    }
    
    D_K_t_max_rescale[t] = max(D_K_t)
    if (verbose){
      print(t) 
    }
  }
  return(max(D_K_t_max_rescale))
}


online_cpd_uase <- function(A_list, tau_factor, h_kernel = 0.1, rank = 10, verbose = FALSE){
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
      
      D_K_t[s] = get_D_K_t_cpp(P_left, P_right, Z, Sigma) / (1/s^0.5 + 1/(t-s)^0.5) / log(max(c(n_1, n_2, t)))^0.5
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






