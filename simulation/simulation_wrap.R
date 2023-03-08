#### hpca ####

hpca_simulation_sbm <- function(TT, cp_truth, n, L, B = 100, N_MC = 100){
  
  alpha = 0.01
  
  n_1 = n
  n_2 = n
  set.seed(n_1*n_2*L)
  
  params = get_sbm_params(n, L)
  probability_1=params[[1]]
  probability_2=params[[2]]
  
  
  K_max = 500
  h_kernel = (K_max*log(TT*n_1 * n_2)/(TT*n_1*n_2))^{1/L}
  
  ### burn-in
  T_burn = as.integer(cp_truth*1.5)
  
  set.seed(n_1*n_2*L)
  A_burn = get_data_burn_sbm(T_burn, n, L, probability_1)
  
  max_D_rescale_B = rep(0, B)
  tic()
  for(b in 1:B){
    set.seed(n_1*n_2*L*b)
    b_ix = sample(1:T_burn, replace = FALSE)
    A_b = A_burn[, , , b_ix]
    A_b_list = list()
    for (t in 1:T_burn){
      A_b_list[[t]] = A_b[, , , t]
    }
    
    set.seed(n_1*n_2*L*b)
    max_D_rescale_B[b] = max(max_D_s_t_rescale_cpp(A_b_list, h_kernel))
    print(paste0("b = ", b))
  }
  toc()
  
  tau_factor = quantile(max_D_rescale_B, 1 - alpha, type = 1)
  
  print(paste0("tau hpca = ", tau_factor))
  
  ### online cpd ####
  
  tic()
  
  t_hat_record = rep(NA, N_MC)
  find_record = rep(FALSE, N_MC)
  D_path_record = list()
  
  for (b in 1:N_MC){
    set.seed(n_1*n_2*L*b)
    A_list = get_data_cp_sbm(TT, cp_truth, n, L, probability_1, probability_2)
    A_list_list = list()
    for (t in 1:TT){
      A_list_list[[t]] = A_list[, , , t]
    }
    
    #### reset seed for the method
    set.seed(n_1*n_2*L*b)
    result_online_cpd = online_cpd_cpp(A_list_list, tau_factor, h_kernel)
    
    t_hat_record[b] = result_online_cpd$t
    find_record[b] = result_online_cpd$find
    D_path_record[[b]] = result_online_cpd$D_K_t_max_rescale
    
    print(b)
  }
  toc()
  
  t_hat_found = t_hat_record[find_record]
  ix_correct = t_hat_record > cp_truth
  
  
  type_I_error = 1 - mean(ix_correct)
  type_II_error = 1 - mean(find_record)
  run_length = mean(t_hat_record[ix_correct] - cp_truth - 1)
  
  return(c(type_I_error,
           type_II_error,
           run_length))
}




hpca_simulation_dirichlet <- function(TT, cp_truth, n_1, n_2, L, d, B = 100, N_MC = 100){
  
  alpha = 0.01
  
  set.seed(n_1*n_2*L*d)
  
  params = get_dirichlet_params(n_1, n_2, L, d)
  dirichlet_xy_1=params[[1]]
  dirichlet_xy_2=params[[2]]
  W_1 = params[[3]]
  W_2 = params[[4]]
  
  
  K_max = 500
  h_kernel = (K_max*log(TT*n_1 * n_2)/(TT*n_1*n_2))^{1/L}
  
  
  ### burn-in
  T_burn = as.integer(cp_truth*1.5)
  
  set.seed(n_1*n_2*L*d)
  A_burn = get_data_burn_dirichlet(T_burn, n_1, n_2, L, dirichlet_xy_1, W_1)
  
  max_D_rescale_B = rep(0, B)
  tic()
  for(b in 1:B){
    set.seed(n_1*n_2*L*d* b)
    b_ix = sample(1:T_burn, replace = FALSE)
    
    A_b = A_burn[, , , b_ix]
    A_b_list = list()
    for (t in 1:T_burn){
      A_b_list[[t]] = A_b[, , , t]
    }
    
    set.seed(n_1*n_2*L*d* b)
    max_D_rescale_B[b] = max(max_D_s_t_rescale_cpp(A_b_list, h_kernel))
    print(paste0("b = ", b))
  }
  toc()
  
  tau_factor = quantile(max_D_rescale_B, 1 - alpha, type = 1)
  
  print(paste0("tau hpca = ", tau_factor))
  
  ### online cpd ####
  
  tic()
  
  t_hat_record = rep(NA, N_MC)
  find_record = rep(FALSE, N_MC)
  D_path_record = list()
  
  for (b in 1:N_MC){
    set.seed(n_1*n_2*L*d* b)
    A_list = get_data_cp_dirichlet(TT, cp_truth, n_1, n_2, L, dirichlet_xy_1, W_1, dirichlet_xy_2, W_2)
    A_list_list = list()
    for (t in 1:TT){
      A_list_list[[t]] = A_list[, , , t]
    }
    
    #### reset seed for the method
    set.seed(n_1*n_2*L*d*b)
    result_online_cpd = online_cpd_cpp(A_list_list, tau_factor, h_kernel)
    
    t_hat_record[b] = result_online_cpd$t
    find_record[b] = result_online_cpd$find
    D_path_record[[b]] = result_online_cpd$D_K_t_max_rescale
    
    print(b)
  }
  toc()
  
  t_hat_found = t_hat_record[find_record]
  ix_correct = t_hat_record > cp_truth
  
  
  type_I_error = 1 - mean(ix_correct)
  type_II_error = 1 - mean(find_record)
  run_length = mean(t_hat_record[ix_correct] - cp_truth - 1)

  return(c(type_I_error,
            type_II_error,
            run_length))
}






#### TH-PCA ####

thpca_simulation_sbm <- function(TT, cp_truth, n, L, B = 100, N_MC = 100){
  
  alpha = 0.01
  
  n_1 = n
  n_2 = n
  set.seed(n_1*n_2*L)
  
  params = get_sbm_params(n, L)
  probability_1=params[[1]]
  probability_2=params[[2]]
  
  
  K_max = 500
  h_kernel = (K_max*log(TT*n_1 * n_2)/(TT*n_1*n_2))^{1/L}
  
  
  ### burn-in
  T_burn = as.integer(cp_truth*1.5)
  
  set.seed(n_1*n_2*L)
  A_burn = get_data_burn_sbm(T_burn, n, L, probability_1)
  
  max_D_rescale_B = rep(0, B)
  tic()
  for(b in 1:B){
    set.seed(n_1*n_2*L*b)
    b_ix = sample(1:T_burn, replace = FALSE)
    A_b = A_burn[, , , b_ix]
    A_b_list = list()
    for (t in 1:T_burn){
      A_b_list[[t]] = A_b[, , , t]
    }
    
    set.seed(n_1*n_2*L*b)
    max_D_rescale_B[b] = max(max_D_s_t_rescale_thpca_cpp(A_b_list, h_kernel))
    print(paste0("b = ", b))
  }
  toc()
  
  tau_factor = quantile(max_D_rescale_B, 1 - alpha, type = 1)
  
  print(paste0("tau hpca = ", tau_factor))
  
  ### online cpd ####
  
  tic()
  
  t_hat_record = rep(NA, N_MC)
  find_record = rep(FALSE, N_MC)
  D_path_record = list()
  
  for (b in 1:N_MC){
    set.seed(n_1*n_2*L*b)
    A_list = get_data_cp_sbm(TT, cp_truth, n, L, probability_1, probability_2)
    A_list_list = list()
    for (t in 1:TT){
      A_list_list[[t]] = A_list[, , , t]
    }
    
    #### reset seed for the method
    set.seed(n_1*n_2*L*b)
    result_online_cpd = online_cpd_thpca_cpp(A_list_list, tau_factor, h_kernel)
    
    t_hat_record[b] = result_online_cpd$t
    find_record[b] = result_online_cpd$find
    D_path_record[[b]] = result_online_cpd$D_K_t_max_rescale
    
    print(b)
  }
  toc()
  
  t_hat_found = t_hat_record[find_record]
  ix_correct = t_hat_record > cp_truth
  
  
  type_I_error = 1 - mean(ix_correct)
  type_II_error = 1 - mean(find_record)
  run_length = mean(t_hat_record[ix_correct] - cp_truth - 1)
  
  return(c(type_I_error,
           type_II_error,
           run_length))
}




thpca_simulation_dirichlet <- function(TT, cp_truth, n_1, n_2, L, d, B = 100, N_MC = 100){
  
  alpha = 0.01
  
  set.seed(n_1*n_2*L*d)
  
  params = get_dirichlet_params(n_1, n_2, L, d)
  dirichlet_xy_1=params[[1]]
  dirichlet_xy_2=params[[2]]
  W_1 = params[[3]]
  W_2 = params[[4]]
  
  
  K_max = 500
  h_kernel = (K_max*log(TT*n_1 * n_2)/(TT*n_1*n_2))^{1/L}
  
  
  ### burn-in
  T_burn = as.integer(cp_truth*1.5)
  
  set.seed(n_1*n_2*L*d)
  A_burn = get_data_burn_dirichlet(T_burn, n_1, n_2, L, dirichlet_xy_1, W_1)
  
  max_D_rescale_B = rep(0, B)
  tic()
  for(b in 1:B){
    set.seed(n_1*n_2*L*d* b)
    b_ix = sample(1:T_burn, replace = FALSE)
    
    A_b = A_burn[, , , b_ix]
    A_b_list = list()
    for (t in 1:T_burn){
      A_b_list[[t]] = A_b[, , , t]
    }
    
    set.seed(n_1*n_2*L*d* b)
    max_D_rescale_B[b] = max(max_D_s_t_rescale_thpca_cpp(A_b_list, h_kernel))
    print(paste0("b = ", b))
  }
  toc()
  
  tau_factor = quantile(max_D_rescale_B, 1 - alpha, type = 1)
  
  print(paste0("tau hpca = ", tau_factor))
  
  ### online cpd ####
  
  tic()
  
  t_hat_record = rep(NA, N_MC)
  find_record = rep(FALSE, N_MC)
  D_path_record = list()
  
  for (b in 1:N_MC){
    set.seed(n_1*n_2*L*d* b)
    A_list = get_data_cp_dirichlet(TT, cp_truth, n_1, n_2, L, dirichlet_xy_1, W_1, dirichlet_xy_2, W_2)
    A_list_list = list()
    for (t in 1:TT){
      A_list_list[[t]] = A_list[, , , t]
    }
    
    #### reset seed for the method
    set.seed(n_1*n_2*L*d*b)
    result_online_cpd = online_cpd_thpca_cpp(A_list_list, tau_factor, h_kernel)
    
    t_hat_record[b] = result_online_cpd$t
    find_record[b] = result_online_cpd$find
    D_path_record[[b]] = result_online_cpd$D_K_t_max_rescale
    
    print(b)
  }
  toc()
  
  t_hat_found = t_hat_record[find_record]
  ix_correct = t_hat_record > cp_truth
  
  
  type_I_error = 1 - mean(ix_correct)
  type_II_error = 1 - mean(find_record)
  run_length = mean(t_hat_record[ix_correct] - cp_truth - 1)
  
  return(c(type_I_error,
           type_II_error,
           run_length))
}









#### uase ####


uase_simulation_sbm_R <- function(TT, cp_truth, n, L, B = 100, N_MC = 100){
  
  alpha = 0.01
  
  n_1 = n
  n_2 = n
  set.seed(n_1*n_2*L)
  
  params = get_sbm_params(n, L)
  probability_1=params[[1]]
  probability_2=params[[2]]
  
  
  K_max = 500
  h_kernel = (K_max*log(TT*n_1 * n_2)/(TT*n_1*n_2))^{1/L}
  
  
  ### burn-in
  T_burn = as.integer(cp_truth*1.5)
  
  set.seed(n_1*n_2*L)
  A_burn = get_data_burn_sbm(T_burn, n, L, probability_1)
  
  max_D_rescale_B = rep(0, B)
  tic()
  for(b in 1:B){
    set.seed(n_1*n_2*L*b)
    b_ix = sample(1:T_burn, replace = FALSE)
    A_b = A_burn[, , , b_ix]
    
    set.seed(n_1*n_2*L*b)
    max_D_rescale_B[b] = max(max_D_s_t_rescale_uase(A_b, h_kernel))
    print(paste0("b = ", b))
  }
  toc()
  
  tau_factor = quantile(max_D_rescale_B, 1 - alpha, type = 1)
  
  print(paste0("tau uase = ", tau_factor))
  
  ### online cpd ####
  
  tic()
  
  t_hat_record = rep(NA, N_MC)
  find_record = rep(FALSE, N_MC)
  D_path_record = list()
  
  for (b in 1:N_MC){
    set.seed(n_1*n_2*L*b)
    A_list = get_data_cp_sbm(TT, cp_truth, n, L, probability_1, probability_2)
    
    #### reset seed for the method
    set.seed(n_1*n_2*L*b)
    result_online_cpd = online_cpd_uase(A_list, tau_factor, h_kernel)
    
    t_hat_record[b] = result_online_cpd$t
    find_record[b] = result_online_cpd$find
    D_path_record[[b]] = result_online_cpd$D_K_t_max_rescale
    
    print(b)
  }
  toc()
  
  t_hat_found = t_hat_record[find_record]
  ix_correct = t_hat_record > cp_truth
  
  
  type_I_error = 1 - mean(ix_correct)
  type_II_error = 1 - mean(find_record)
  run_length = mean(t_hat_record[ix_correct] - cp_truth - 1)
  
  return(c(type_I_error,
           type_II_error,
           run_length))
}





uase_simulation_dirichlet_R <- function(TT, cp_truth, n_1, n_2, L, d, B = 100, N_MC = 100){
  
  alpha = 0.01
  
  set.seed(n_1*n_2*L*d)
  
  params = get_dirichlet_params(n_1, n_2, L, d)
  dirichlet_xy_1=params[[1]]
  dirichlet_xy_2=params[[2]]
  W_1 = params[[3]]
  W_2 = params[[4]]
  
  
  K_max = 500
  h_kernel = (K_max*log(TT*n_1 * n_2)/(TT*n_1*n_2))^{1/L}
  
  
  ### burn-in
  T_burn = as.integer(cp_truth*1.5)
  
  set.seed(n_1*n_2*L*d)
  A_burn = get_data_burn_dirichlet(T_burn, n_1, n_2, L, dirichlet_xy_1, W_1)
  
  max_D_rescale_B = rep(0, B)
  tic()
  for(b in 1:B){
    set.seed(n_1*n_2*L*d*b)
    b_ix = sample(1:T_burn, replace = FALSE)
    
    A_b = A_burn[, , , b_ix]
    
    set.seed(n_1*n_2*L*d*b)
    max_D_rescale_B[b] = max(max_D_s_t_rescale_uase(A_b, h_kernel))
    print(paste0("b = ", b))
  }
  toc()
  
  tau_factor = quantile(max_D_rescale_B, 1 - alpha, type = 1)
  
  print(paste0("tau uase = ", tau_factor))
  
  ### online cpd ####
  
  tic()
  
  t_hat_record = rep(NA, N_MC)
  find_record = rep(FALSE, N_MC)
  D_path_record = list()
  
  for (b in 1:N_MC){
    set.seed(n_1*n_2*L*d*b)
    A_list = get_data_cp_dirichlet(TT, cp_truth, n_1, n_2, L, dirichlet_xy_1, W_1, dirichlet_xy_2, W_2)
    
    #### reset seed for the method
    set.seed(n_1*n_2*L*d*b)
    result_online_cpd = online_cpd_uase(A_list, tau_factor, h_kernel)
    
    t_hat_record[b] = result_online_cpd$t
    find_record[b] = result_online_cpd$find
    D_path_record[[b]] = result_online_cpd$D_K_t_max_rescale
    
    print(b)
  }
  toc()
  
  t_hat_found = t_hat_record[find_record]
  ix_correct = t_hat_record > cp_truth
  
  
  type_I_error = 1 - mean(ix_correct)
  type_II_error = 1 - mean(find_record)
  run_length = mean(t_hat_record[ix_correct] - cp_truth - 1)
  
  return(c(type_I_error,
           type_II_error,
           run_length))
}





uase_simulation_sbm <- function(TT, cp_truth, n, L, B = 100, N_MC = 100){
  
  alpha = 0.01
  
  n_1 = n
  n_2 = n
  set.seed(n_1*n_2*L)
  
  params = get_sbm_params(n, L)
  probability_1=params[[1]]
  probability_2=params[[2]]
  
  
  K_max = 500
  h_kernel = (K_max*log(TT*n_1 * n_2)/(TT*n_1*n_2))^{1/L}
  
  
  ### burn-in
  T_burn = as.integer(cp_truth*1.5)
  
  set.seed(n_1*n_2*L)
  A_burn = get_data_burn_sbm(T_burn, n, L, probability_1)
  
  max_D_rescale_B = rep(0, B)
  tic()
  for(b in 1:B){
    set.seed(n_1*n_2*L*b)
    b_ix = sample(1:T_burn, replace = FALSE)
    A_b = A_burn[, , , b_ix]
    A_b_list = list()
    for (t in 1:T_burn){
      A_b_list[[t]] = A_b[, , , t]
    }
    
    set.seed(n_1*n_2*L*b)
    max_D_rescale_B[b] = max(max_D_s_t_rescale_uase_cpp(A_b_list, h_kernel))
    print(paste0("b = ", b))
  }
  toc()
  
  tau_factor = quantile(max_D_rescale_B, 1 - alpha, type = 1)
  
  print(paste0("tau uase = ", tau_factor))
  
  ### online cpd ####
  
  tic()
  
  t_hat_record = rep(NA, N_MC)
  find_record = rep(FALSE, N_MC)
  D_path_record = list()
  
  for (b in 1:N_MC){
    set.seed(n_1*n_2*L*b)
    A_list = get_data_cp_sbm(TT, cp_truth, n, L, probability_1, probability_2)
    A_list_list = list()
    for (t in 1:TT){
      A_list_list[[t]] = A_list[, , , t]
    }
    
    #### reset seed for the method
    set.seed(n_1*n_2*L*b)
    result_online_cpd = online_cpd_uase_cpp(A_list_list, tau_factor, h_kernel)
    
    t_hat_record[b] = result_online_cpd$t
    find_record[b] = result_online_cpd$find
    D_path_record[[b]] = result_online_cpd$D_K_t_max_rescale
    
    print(b)
  }
  toc()
  
  t_hat_found = t_hat_record[find_record]
  ix_correct = t_hat_record > cp_truth
  
  
  type_I_error = 1 - mean(ix_correct)
  type_II_error = 1 - mean(find_record)
  run_length = mean(t_hat_record[ix_correct] - cp_truth - 1)
  
  return(c(type_I_error,
           type_II_error,
           run_length))
}




uase_simulation_dirichlet <- function(TT, cp_truth, n_1, n_2, L, d, B = 100, N_MC = 100){
  
  alpha = 0.01
  
  set.seed(n_1*n_2*L*d)
  
  params = get_dirichlet_params(n_1, n_2, L, d)
  dirichlet_xy_1=params[[1]]
  dirichlet_xy_2=params[[2]]
  W_1 = params[[3]]
  W_2 = params[[4]]
  
  
  K_max = 500
  h_kernel = (K_max*log(TT*n_1 * n_2)/(TT*n_1*n_2))^{1/L}
  
  
  ### burn-in
  T_burn = as.integer(cp_truth*1.5)
  
  set.seed(n_1*n_2*L*d)
  A_burn = get_data_burn_dirichlet(T_burn, n_1, n_2, L, dirichlet_xy_1, W_1)
  
  max_D_rescale_B = rep(0, B)
  tic()
  for(b in 1:B){
    set.seed(n_1*n_2*L*d*b)
    b_ix = sample(1:T_burn, replace = FALSE)
    
    A_b = A_burn[, , , b_ix]
    A_b_list = list()
    for (t in 1:T_burn){
      A_b_list[[t]] = A_b[, , , t]
    }
    
    set.seed(n_1*n_2*L*d*b)
    max_D_rescale_B[b] = max(max_D_s_t_rescale_uase_cpp(A_b_list, h_kernel))
    print(paste0("b = ", b))
  }
  toc()
  
  tau_factor = quantile(max_D_rescale_B, 1 - alpha, type = 1)
  
  print(paste0("tau uase = ", tau_factor))
  
  ### online cpd ####
  
  tic()
  
  t_hat_record = rep(NA, N_MC)
  find_record = rep(FALSE, N_MC)
  D_path_record = list()
  
  for (b in 1:N_MC){
    set.seed(n_1*n_2*L*d*b)
    A_list = get_data_cp_dirichlet(TT, cp_truth, n_1, n_2, L, dirichlet_xy_1, W_1, dirichlet_xy_2, W_2)
    A_list_list = list()
    for (t in 1:TT){
      A_list_list[[t]] = A_list[, , , t]
    }
    
    #### reset seed for the method
    set.seed(n_1*n_2*L*d*b)
    result_online_cpd = online_cpd_uase_cpp(A_list_list, tau_factor, h_kernel)
    
    t_hat_record[b] = result_online_cpd$t
    find_record[b] = result_online_cpd$find
    D_path_record[[b]] = result_online_cpd$D_K_t_max_rescale
    
    print(b)
  }
  toc()
  
  t_hat_found = t_hat_record[find_record]
  ix_correct = t_hat_record > cp_truth
  
  
  type_I_error = 1 - mean(ix_correct)
  type_II_error = 1 - mean(find_record)
  run_length = mean(t_hat_record[ix_correct] - cp_truth - 1)
  
  return(c(type_I_error,
           type_II_error,
           run_length))
}




#### multi ####

multi_simulation_sbm <- function(TT, cp_truth, n, L, B = 100, N_MC = 100){
  
  alpha = 0.01
  
  n_1 = n
  n_2 = n
  set.seed(n_1*n_2*L)
  
  params = get_sbm_params(n, L)
  probability_1=params[[1]]
  probability_2=params[[2]]
  
  
  K_max = 500
  h_kernel = (K_max*log(TT*n_1 * n_2)/(TT*n_1*n_2))^{1/L}
  
  
  ### burn-in
  T_burn = as.integer(cp_truth*1.5)
  
  set.seed(n_1*n_2*L)
  A_burn = get_data_burn_sbm(T_burn, n, L, probability_1)
  
  max_D_rescale_B = rep(0, B)
  tic()
  for(b in 1:B){
    set.seed(n_1*n_2*L*b)
    b_ix = sample(1:T_burn, replace = FALSE)
    A_b = A_burn[, , , b_ix]
    
    set.seed(n_1*n_2*L*b)
    max_D_rescale_B[b] = max(max_D_s_t_rescale_multi(A_b, h_kernel))
    print(paste0("b = ", b))
  }
  toc()
  
  tau_factor = quantile(max_D_rescale_B, 1 - alpha, type = 1)
  
  print(paste0("tau multi = ", tau_factor))
  
  ### online cpd ####
  
  tic()
  
  t_hat_record = rep(NA, N_MC)
  find_record = rep(FALSE, N_MC)
  D_path_record = list()
  
  for (b in 1:N_MC){
    set.seed(n_1*n_2*L*b)
    A_list = get_data_cp_sbm(TT, cp_truth, n, L, probability_1, probability_2)
    
    #### reset seed for the method
    set.seed(n_1*n_2*L*b)
    result_online_cpd = online_cpd_multi(A_list, tau_factor, h_kernel)
    
    t_hat_record[b] = result_online_cpd$t
    find_record[b] = result_online_cpd$find
    D_path_record[[b]] = result_online_cpd$D_K_t_max_rescale
    
    print(b)
  }
  toc()
  
  t_hat_found = t_hat_record[find_record]
  ix_correct = t_hat_record > cp_truth
  
  
  type_I_error = 1 - mean(ix_correct)
  type_II_error = 1 - mean(find_record)
  run_length = mean(t_hat_record[ix_correct] - cp_truth - 1)
  
  return(c(type_I_error,
           type_II_error,
           run_length))
}




multi_simulation_dirichlet <- function(TT, cp_truth, n_1, n_2, L, d, B = 100, N_MC = 100){
  
  alpha = 0.01
  
  set.seed(n_1*n_2*L*d)
  
  params = get_dirichlet_params(n_1, n_2, L, d)
  dirichlet_xy_1=params[[1]]
  dirichlet_xy_2=params[[2]]
  W_1 = params[[3]]
  W_2 = params[[4]]
  
  
  K_max = 500
  h_kernel = (K_max*log(TT*n_1 * n_2)/(TT*n_1*n_2))^{1/L}
  
  
  ### burn-in
  T_burn = as.integer(cp_truth*1.5)
  
  set.seed(n_1*n_2*L*d)
  A_burn = get_data_burn_dirichlet(T_burn, n_1, n_2, L, dirichlet_xy_1, W_1)
  
  max_D_rescale_B = rep(0, B)
  tic()
  for(b in 1:B){
    set.seed(n_1*n_2*L*d*b)
    b_ix = sample(1:T_burn, replace = FALSE)
    
    A_b = A_burn[, , , b_ix]
    
    set.seed(n_1*n_2*L*d*b)
    max_D_rescale_B[b] = max(max_D_s_t_rescale_multi(A_b, h_kernel))
    print(paste0("b = ", b))
  }
  toc()
  
  tau_factor = quantile(max_D_rescale_B, 1 - alpha, type = 1)
  
  print(paste0("tau multi = ", tau_factor))
  
  ### online cpd ####

  tic()
  
  t_hat_record = rep(NA, N_MC)
  find_record = rep(FALSE, N_MC)
  D_path_record = list()
  
  for (b in 1:N_MC){
    set.seed(n_1*n_2*L*d*b)
    A_list = get_data_cp_dirichlet(TT, cp_truth, n_1, n_2, L, dirichlet_xy_1, W_1, dirichlet_xy_2, W_2)
    
    #### reset seed for the method
    set.seed(n_1*n_2*L*d*b)
    result_online_cpd = online_cpd_multi(A_list, tau_factor, h_kernel)
    
    t_hat_record[b] = result_online_cpd$t
    find_record[b] = result_online_cpd$find
    D_path_record[[b]] = result_online_cpd$D_K_t_max_rescale
    
    print(b)
  }
  toc()
  
  t_hat_found = t_hat_record[find_record]
  ix_correct = t_hat_record > cp_truth
  
  
  type_I_error = 1 - mean(ix_correct)
  type_II_error = 1 - mean(find_record)
  run_length = mean(t_hat_record[ix_correct] - cp_truth - 1)
  
  return(c(type_I_error,
           type_II_error,
           run_length))
}




#### knn ####

knn_simulation_sbm <- function(TT, cp_truth, n, L, N_MC = 100){
  
  alpha = 0.01
  
  n_1 = n
  n_2 = n
  set.seed(n_1*n_2*L)
  
  params = get_sbm_params(n, L)
  probability_1=params[[1]]
  probability_2=params[[2]]

  k_nn = 3
  n0 = floor(0.3 * TT)
  n1 = floor(0.7 * TT)
  ARL = 10000
  
  ### burn-in
  T_burn = as.integer(cp_truth*1.5)
  
  set.seed(n_1*n_2*L)
  A_burn = get_data_burn_sbm(T_burn, n, L, probability_1)
  
  ### online cpd ####
  
  tic()
  
  t_hat_record = rep(NA, N_MC)
  find_record = rep(FALSE, N_MC)
  
  for (b in 1:N_MC){
    set.seed(n_1*n_2*L*b)
    b_ix = sample(1:T_burn, replace = FALSE)
    A_b = A_burn[, , , b_ix]
    
    set.seed(n_1*n_2*L*b)
    A_list = get_data_cp_sbm(TT, cp_truth, n, L, probability_1, probability_2)
    
    #### reset seed for the method
    set.seed(n_1*n_2*L*b)
    result_online_cpd = online_cpd_knn(A_b, A_list, k_nn, alpha, n0, n1, ARL)
    
    t_hat_record[b] = result_online_cpd$t
    find_record[b] = result_online_cpd$find

    print(b)
  }
  toc()
  
  t_hat_found = t_hat_record[find_record]
  ix_correct = t_hat_record > cp_truth
  
  
  type_I_error = 1 - mean(ix_correct)
  type_II_error = 1 - mean(find_record)
  run_length = mean(t_hat_record[ix_correct] - cp_truth - 1)
  
  return(c(type_I_error,
           type_II_error,
           run_length))
}




knn_simulation_dirichlet <- function(TT, cp_truth, n_1, n_2, L, d, N_MC = 100){
  
  alpha = 0.01
  
  set.seed(n_1*n_2*L*d)
  
  params = get_dirichlet_params(n_1, n_2, L, d)
  dirichlet_xy_1=params[[1]]
  dirichlet_xy_2=params[[2]]
  W_1 = params[[3]]
  W_2 = params[[4]]
  
  k_nn = 3
  n0 = floor(0.3 * TT)
  n1 = floor(0.7 * TT)
  ARL = 10000
  
  ### burn-in
  T_burn = as.integer(cp_truth*1.5)
  
  set.seed(n_1*n_2*L*d)
  A_burn = get_data_burn_dirichlet(T_burn, n_1, n_2, L, dirichlet_xy_1, W_1)
  
  ### online cpd ####
  
  tic()
  t_hat_record = rep(NA, N_MC)
  find_record = rep(FALSE, N_MC)
  D_path_record = list()
  
  for (b in 1:N_MC){
    set.seed(n_1*n_2*L*d*b)
    b_ix = sample(1:T_burn, replace = FALSE)
    A_b = A_burn[, , , b_ix]
    
    set.seed(n_1*n_2*L*d*b)
    A_list = get_data_cp_dirichlet(TT, cp_truth, n_1, n_2, L, dirichlet_xy_1, W_1, dirichlet_xy_2, W_2)
    
    #### reset seed for the method
    set.seed(n_1*n_2*L*d*b)
    result_online_cpd = online_cpd_knn(A_b, A_list, k_nn, alpha, n0, n1, ARL)
    
    t_hat_record[b] = result_online_cpd$t
    find_record[b] = result_online_cpd$find
    
    print(b)
  }
  toc()
  
  t_hat_found = t_hat_record[find_record]
  ix_correct = t_hat_record > cp_truth
  
  
  type_I_error = 1 - mean(ix_correct)
  type_II_error = 1 - mean(find_record)
  run_length = mean(t_hat_record[ix_correct] - cp_truth - 1)
  
  return(c(type_I_error,
           type_II_error,
           run_length))
}

