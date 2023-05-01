#### whole experiment ####

cpd_sbm <- function(TT, cp_truth, n_list, L_list, rand_pos = FALSE, 
                    B = 100, N_MC = 100, verbose_freq = 10){

  # type_I_error = array(NA, dim=c(length(n_list),  length(L_list)))
  # type_II_error = array(NA, dim=c(length(n_list), length(L_list)))
  # run_length = array(NA, dim=c(length(n_list), length(L_list)))
  # 
  # type_I_error_thpca = array(NA, dim=c(length(n_list),  length(L_list)))
  # type_II_error_thpca = array(NA, dim=c(length(n_list), length(L_list)))
  # run_length_thpca = array(NA, dim=c(length(n_list), length(L_list)))
  # 
  # type_I_error_thpca_r1 = array(NA, dim=c(length(n_list),  length(L_list)))
  # type_II_error_thpca_r1 = array(NA, dim=c(length(n_list), length(L_list)))
  # run_length_thpca_r1 = array(NA, dim=c(length(n_list), length(L_list)))
  # 
  # type_I_error_uase = array(NA, dim=c(length(n_list),  length(L_list)))
  # type_II_error_uase = array(NA, dim=c(length(n_list), length(L_list)))
  # run_length_uase = array(NA, dim=c(length(n_list), length(L_list)))
  
  type_I_error_multi = array(NA, dim=c(length(n_list),  length(L_list)))
  type_II_error_multi = array(NA, dim=c(length(n_list), length(L_list)))
  run_length_multi = array(NA, dim=c(length(n_list), length(L_list)))
  
  # type_I_error_knn = array(NA, dim=c(length(n_list),  length(L_list)))
  # type_II_error_knn = array(NA, dim=c(length(n_list), length(L_list)))
  # run_length_knn = array(NA, dim=c(length(n_list), length(L_list)))
  
  record = array("", dim=c(length(n_list),  length(L_list)))
  
  for (n_index in 1:length(n_list)) {
    
    n = n_list[n_index]
    
    for (L_index in 1:length(L_list)){
            
      L =  L_list[L_index]
      info = paste0("n = ", n, ", L = ", L)
      print(info)
      record[n_index, L_index] = info
      
      # print("---- hosvd ----")
      # hosvd_res = hosvd_simulation_sbm(TT, cp_truth, n, L, rank = rep(10, 3), rand_pos, B, N_MC)
      # type_I_error[n_index, L_index] = hosvd_res[1]
      # type_II_error[n_index, L_index] = hosvd_res[2]
      # run_length[n_index, L_index] = hosvd_res[3]
      # print(hosvd_res)
      # 
      # print("---- thpca ----")
      # thpca_res = thpca_simulation_sbm(TT, cp_truth, n, L, rank = rep(10, 3), rand_pos, B, N_MC)
      # type_I_error_thpca[n_index, L_index] = thpca_res[1]
      # type_II_error_thpca[n_index, L_index] = thpca_res[2]
      # run_length_thpca[n_index, L_index] = thpca_res[3]
      # print(thpca_res)
      # 
      # print("---- thpca rank 1----")
      # thpca_res_r1 = thpca_simulation_sbm(TT, cp_truth, n, L, rank = rep(1, 3), rand_pos, B, N_MC)
      # type_I_error_thpca_r1[n_index, L_index] = thpca_res_r1[1]
      # type_II_error_thpca_r1[n_index, L_index] = thpca_res_r1[2]
      # run_length_thpca_r1[n_index, L_index] = thpca_res_r1[3]
      # print(thpca_res_r1)
      # 
      # print("---- uase ----")
      # uase_res = uase_simulation_sbm(TT, cp_truth, n, L, rank = 10, rand_pos, B, N_MC)
      # type_I_error_uase[n_index, L_index] = uase_res[1]
      # type_II_error_uase[n_index, L_index] = uase_res[2]
      # run_length_uase[n_index, L_index]= uase_res[3]
      # print(uase_res)
      
      print("---- multiness ----")
      multi_res = multi_simulation_sbm(TT, cp_truth, n, L, rank = 10, rand_pos, B, N_MC)
      type_I_error_multi[n_index, L_index] = multi_res[1]
      type_II_error_multi[n_index, L_index] = multi_res[2]
      run_length_multi[n_index, L_index] = multi_res[3]
      print(multi_res)
      
      # print("---- knn ----")
      # knn_res = knn_simulation_sbm(TT, cp_truth, n, L, rand_pos, N_MC)
      # type_I_error_knn[n_index, L_index] = knn_res[1]
      # type_II_error_knn[n_index, L_index] = knn_res[2]
      # run_length_knn[n_index, L_index] = knn_res[3]
      # print(knn_res)
      
      
    }
  }
  
  print(record)
  
  # print("---- hosvd ----")
  # print(type_I_error)
  # print(type_II_error)
  # print(run_length)
  # 
  # print("---- thpca ----")
  # print(type_I_error_thpca)
  # print(type_II_error_thpca)
  # print(run_length_thpca)
  # 
  # print("---- thpca rank 1----")
  # print(type_I_error_thpca_r1)
  # print(type_II_error_thpca_r1)
  # print(run_length_thpca_r1)
  # 
  # print("---- uase ----")
  # print(type_I_error_uase)
  # print(type_II_error_uase)
  # print(run_length_uase)
  
  print("---- multiness ----")
  print(type_I_error_multi)
  print(type_II_error_multi)
  print(run_length_multi)
  
  # print("---- knn ----")
  # print(type_I_error_knn)
  # print(type_II_error_knn)
  # print(run_length_knn)

}





cpd_dirichlet <- function(TT, cp_truth, n_list, L_list, directed = TRUE,
                    B = 100, N_MC = 100, verbose_freq = 10){
  
  # type_I_error = array(NA, dim=c(length(n_list_1),  length(n_list_2), length(L_list), length(d_list)))
  # type_II_error = array(NA, dim=c(length(n_list_1), length(n_list_2), length(L_list), length(d_list)))
  # run_length = array(NA, dim=c(length(n_list_1), length(n_list_2), length(L_list), length(d_list)))
  # 
  # type_I_error_thpca =  array(NA, dim=c(length(n_list_1),  length(n_list_2), length(L_list), length(d_list)))
  # type_II_error_thpca = array(NA, dim=c(length(n_list_1),  length(n_list_2), length(L_list), length(d_list)))
  # run_length_thpca =  array(NA, dim=c(length(n_list_1),  length(n_list_2), length(L_list), length(d_list)))
  # 
  # type_I_error_thpca_r1 =  array(NA, dim=c(length(n_list_1),  length(n_list_2), length(L_list), length(d_list)))
  # type_II_error_thpca_r1 = array(NA, dim=c(length(n_list_1),  length(n_list_2), length(L_list), length(d_list)))
  # run_length_thpca_r1 =  array(NA, dim=c(length(n_list_1),  length(n_list_2), length(L_list), length(d_list)))
  # 
  # type_I_error_uase = array(NA, dim=c(length(n_list_1),  length(n_list_2), length(L_list), length(d_list)))
  # type_II_error_uase = array(NA, dim=c(length(n_list_1), length(n_list_2), length(L_list), length(d_list)))
  # run_length_uase = array(NA, dim=c(length(n_list_1), length(n_list_2), length(L_list), length(d_list)))

  type_I_error_multi = array(NA, dim=c(length(n_list_1),  length(n_list_2), length(L_list), length(d_list)))
  type_II_error_multi = array(NA, dim=c(length(n_list_1), length(n_list_2), length(L_list), length(d_list)))
  run_length_multi = array(NA, dim=c(length(n_list_1), length(n_list_2), length(L_list), length(d_list)))

  # type_I_error_knn = array(NA, dim=c(length(n_list_1),  length(n_list_2), length(L_list), length(d_list)))
  # type_II_error_knn = array(NA, dim=c(length(n_list_1), length(n_list_2), length(L_list), length(d_list)))
  # run_length_knn = array(NA, dim=c(length(n_list_1), length(n_list_2), length(L_list), length(d_list)))
  
  record = array("", dim=c(length(n_list_1),  length(n_list_2), length(L_list), length(d_list)))
  
  for (n_index_1 in 1:length(n_list_1)) {
    
    n_1 = n_list_1[n_index_1]
    
    for (n_index_2 in 1:length(n_list_2)) {  
      
      n_2 = n_list_2[n_index_2]
      
      for (L_index in 1:length(L_list)){
        
        L =  L_list[L_index]
        
        for (d_index in 1:length(d_list)){
          
          d =  d_list[d_index]
          
          info = paste0("n_1 = ", n_1, ", n_2 = ", n_2, ", L = ", L, ", d = ", d)
          print(info)
          record[n_index_1, n_index_2, L_index, d_index] = info
          
          # print("---- hosvd ----")
          # hosvd_res = hosvd_simulation_dirichlet(TT, cp_truth, n_1, n_2, L, d, rank = rep(10, 3), directed, B, N_MC)
          # type_I_error[n_index_1, n_index_2, L_index, d_index] = hosvd_res[1]
          # type_II_error[n_index_1, n_index_2, L_index, d_index] = hosvd_res[2]
          # run_length[n_index_1, n_index_2, L_index, d_index] = hosvd_res[3]
          # print(hosvd_res)
          # 
          # print("---- thpca ----")
          # thpca_res = thpca_simulation_dirichlet(TT, cp_truth, n_1, n_2, L, d, rank = rep(10, 3), directed, B, N_MC)
          # type_I_error_thpca[n_index_1, n_index_2, L_index, d_index] = thpca_res[1]
          # type_II_error_thpca[n_index_1, n_index_2, L_index, d_index] = thpca_res[2]
          # run_length_thpca[n_index_1, n_index_2, L_index, d_index] = thpca_res[3]
          # print(thpca_res)
          # 
          # print("---- thpca rank 1 ----")
          # thpca_res_r1 = thpca_simulation_dirichlet(TT, cp_truth, n_1, n_2, L, d, rank = rep(1, 3), directed, B, N_MC)
          # type_I_error_thpca_r1[n_index_1, n_index_2, L_index, d_index] = thpca_res_r1[1]
          # type_II_error_thpca_r1[n_index_1, n_index_2, L_index, d_index] = thpca_res_r1[2]
          # run_length_thpca_r1[n_index_1, n_index_2, L_index, d_index] = thpca_res_r1[3]
          # print(thpca_res_r1)
          # 
          # print("---- uase ----")
          # uase_res = uase_simulation_dirichlet(TT, cp_truth, n_1, n_2, L, d, rank = 10, directed, B, N_MC)
          # type_I_error_uase[n_index_1, n_index_2, L_index, d_index] = uase_res[1]
          # type_II_error_uase[n_index_1, n_index_2, L_index, d_index] = uase_res[2]
          # run_length_uase[n_index_1, n_index_2, L_index, d_index]= uase_res[3]
          # print(uase_res)

          print("---- multiness ----")
          multi_res = multi_simulation_dirichlet(TT, cp_truth, n_1, n_2, L, d, rank = 10, directed, B, N_MC)
          type_I_error_multi[n_index_1, n_index_2, L_index, d_index] = multi_res[1]
          type_II_error_multi[n_index_1, n_index_2, L_index, d_index] = multi_res[2]
          run_length_multi[n_index_1, n_index_2, L_index, d_index] = multi_res[3]
          print(multi_res)

          # print("---- knn ----")
          # knn_res = knn_simulation_dirichlet(TT, cp_truth, n_1, n_2, L, d, directed, N_MC)
          # type_I_error_knn[n_index_1, n_index_2, L_index, d_index] = knn_res[1]
          # type_II_error_knn[n_index_1, n_index_2, L_index, d_index] = knn_res[2]
          # run_length_knn[n_index_1, n_index_2, L_index, d_index] = knn_res[3]
          # print(knn_res)

          
        }
      }
    }
  }
  
  print(record)
  
  # print("---- hosvd ----")
  # print(type_I_error)
  # print(type_II_error)
  # print(run_length)
  # 
  # print("---- thpca ----")
  # print(type_I_error_thpca)
  # print(type_II_error_thpca)
  # print(run_length_thpca)
  # 
  # print("---- thpca rank 1----")
  # print(type_I_error_thpca_r1)
  # print(type_II_error_thpca_r1)
  # print(run_length_thpca_r1)
  # 
  # print("---- uase ----")
  # print(type_I_error_uase)
  # print(type_II_error_uase)
  # print(run_length_uase)
  
  print("---- multiness ----")
  print(type_I_error_multi)
  print(type_II_error_multi)
  print(run_length_multi)
  
  # print("---- knn ----")
  # print(type_I_error_knn)
  # print(type_II_error_knn)
  # print(run_length_knn)
  
}








#### HOSVD ####

hosvd_simulation_sbm <- function(TT, cp_truth, n, L, rank = rep(10, 3), rand_pos = FALSE, 
                                B = 100, N_MC = 100, verbose_freq = 10){
  
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
  A_burn = get_data_burn_sbm(T_burn, n, L, probability_1, rand_pos)
  
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
    if (rand_pos){
      max_D_rescale_B[b] = max(max_D_s_t_rescale_cpp(A_b_list, h_kernel, rank, directed = FALSE))
    }
    else{
      max_D_rescale_B[b] = max(max_D_s_t_rescale_fixed_cpp(A_b_list, rank))
    }
    
    if (b == 1 | b%%verbose_freq == 0){
      print(paste0("b = ", b))
    }
    
  }
  toc()
  
  tau_factor = quantile(max_D_rescale_B, 1 - alpha, type = 1)
  
  print(paste0("tau hosvd = ", tau_factor))
  
  ### online cpd ####
  
  tic()
  
  t_hat_record = rep(NA, N_MC)
  find_record = rep(FALSE, N_MC)
  D_path_record = list()
  
  for (b in 1:N_MC){
    set.seed(n_1*n_2*L*b)
    A_list = get_data_cp_sbm(TT, cp_truth, n, L, probability_1, probability_2, rand_pos)
    A_list_list = list()
    for (t in 1:TT){
      A_list_list[[t]] = A_list[, , , t]
    }
    
    #### reset seed for the method
    set.seed(n_1*n_2*L*b)
    
    if (rand_pos){
      result_online_cpd = online_cpd_cpp(A_list_list, tau_factor, h_kernel, rank, directed = FALSE)
    }
    else{
      result_online_cpd = online_cpd_fixed_cpp(A_list_list, tau_factor, rank)
    }
    
    t_hat_record[b] = result_online_cpd$t
    find_record[b] = result_online_cpd$find
    D_path_record[[b]] = result_online_cpd$D_K_t_max_rescale
    
    if (b == 1 | b%%verbose_freq == 0){
      print(paste0("b = ", b))
    }
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




hosvd_simulation_dirichlet <- function(TT, cp_truth, n_1, n_2, L, d, rank = rep(10, 3), directed = TRUE,
                                      B = 100, N_MC = 100, verbose_freq = 10){
  
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
  A_burn = get_data_burn_dirichlet(T_burn, n_1, n_2, L, dirichlet_xy_1, W_1, directed)
  
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
    max_D_rescale_B[b] = max(max_D_s_t_rescale_cpp(A_b_list, h_kernel, rank, directed))
    
    if (b == 1 | b%%verbose_freq == 0){
      print(paste0("b = ", b))
    }
    
  }
  toc()
  
  tau_factor = quantile(max_D_rescale_B, 1 - alpha, type = 1)
  
  print(paste0("tau hosvd = ", tau_factor))
  
  ### online cpd ####
  
  tic()
  
  t_hat_record = rep(NA, N_MC)
  find_record = rep(FALSE, N_MC)
  D_path_record = list()
  
  for (b in 1:N_MC){
    set.seed(n_1*n_2*L*d* b)
    A_list = get_data_cp_dirichlet(TT, cp_truth, n_1, n_2, L, dirichlet_xy_1, W_1, dirichlet_xy_2, W_2, directed)
    A_list_list = list()
    for (t in 1:TT){
      A_list_list[[t]] = A_list[, , , t]
    }
    
    #### reset seed for the method
    set.seed(n_1*n_2*L*d*b)
    result_online_cpd = online_cpd_cpp(A_list_list, tau_factor, h_kernel, rank, directed)
    
    t_hat_record[b] = result_online_cpd$t
    find_record[b] = result_online_cpd$find
    D_path_record[[b]] = result_online_cpd$D_K_t_max_rescale
    
    if (b == 1 | b%%verbose_freq == 0){
      print(paste0("b = ", b))
    }
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

thpca_simulation_sbm <- function(TT, cp_truth, n, L, rank = rep(10, 3), rand_pos = FALSE, 
                                 B = 100, N_MC = 100, verbose_freq = 10){
  
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
  A_burn = get_data_burn_sbm(T_burn, n, L, probability_1, rand_pos)
  
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
    if (rand_pos){
      max_D_rescale_B[b] = max(max_D_s_t_rescale_thpca_cpp(A_b_list, h_kernel, rank, directed = FALSE))
    }
    else{
      max_D_rescale_B[b] = max(max_D_s_t_rescale_fixed_thpca_cpp(A_b_list, rank))
    }
    
    if (b == 1 | b%%verbose_freq == 0){
      print(paste0("b = ", b))
    }
  }
  toc()
  
  tau_factor = quantile(max_D_rescale_B, 1 - alpha, type = 1)
  
  print(paste0("tau thpca = ", tau_factor))
  
  ### online cpd ####
  
  tic()
  
  t_hat_record = rep(NA, N_MC)
  find_record = rep(FALSE, N_MC)
  D_path_record = list()
  
  for (b in 1:N_MC){
    set.seed(n_1*n_2*L*b)
    A_list = get_data_cp_sbm(TT, cp_truth, n, L, probability_1, probability_2, rand_pos)
    A_list_list = list()
    for (t in 1:TT){
      A_list_list[[t]] = A_list[, , , t]
    }
    
    #### reset seed for the method
    set.seed(n_1*n_2*L*b)
    if (rand_pos){
      result_online_cpd = online_cpd_thpca_cpp(A_list_list, tau_factor, h_kernel, rank, directed = FALSE)  
    }
    else{
      result_online_cpd = online_cpd_fixed_thpca_cpp(A_list_list, tau_factor, rank)
    }
    
    t_hat_record[b] = result_online_cpd$t
    find_record[b] = result_online_cpd$find
    D_path_record[[b]] = result_online_cpd$D_K_t_max_rescale
    
    if (b == 1 | b%%verbose_freq == 0){
      print(paste0("b = ", b))
    }
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




thpca_simulation_dirichlet <- function(TT, cp_truth, n_1, n_2, L, d, rank = rep(10, 3), directed = TRUE,
                                       B = 100, N_MC = 100, verbose_freq = 10){
  
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
  A_burn = get_data_burn_dirichlet(T_burn, n_1, n_2, L, dirichlet_xy_1, W_1, directed)
  
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
    max_D_rescale_B[b] = max(max_D_s_t_rescale_thpca_cpp(A_b_list, h_kernel, rank, directed))
    if (b == 1 | b%%verbose_freq == 0){
      print(paste0("b = ", b))
    }
  }
  toc()
  
  tau_factor = quantile(max_D_rescale_B, 1 - alpha, type = 1)
  
  print(paste0("tau thpca = ", tau_factor))
  
  ### online cpd ####
  
  tic()
  
  t_hat_record = rep(NA, N_MC)
  find_record = rep(FALSE, N_MC)
  D_path_record = list()
  
  for (b in 1:N_MC){
    set.seed(n_1*n_2*L*d* b)
    A_list = get_data_cp_dirichlet(TT, cp_truth, n_1, n_2, L, dirichlet_xy_1, W_1, dirichlet_xy_2, W_2, directed)
    A_list_list = list()
    for (t in 1:TT){
      A_list_list[[t]] = A_list[, , , t]
    }
    
    #### reset seed for the method
    set.seed(n_1*n_2*L*d*b)
    result_online_cpd = online_cpd_thpca_cpp(A_list_list, tau_factor, h_kernel, rank, directed)
    
    t_hat_record[b] = result_online_cpd$t
    find_record[b] = result_online_cpd$find
    D_path_record[[b]] = result_online_cpd$D_K_t_max_rescale
    
    if (b == 1 | b%%verbose_freq == 0){
      print(paste0("b = ", b))
    }
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

uase_simulation_sbm <- function(TT, cp_truth, n, L, rank = 10, rand_pos = FALSE, 
                                B = 100, N_MC = 100, verbose_freq = 10){
  
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
  A_burn = get_data_burn_sbm(T_burn, n, L, probability_1, rand_pos)
  
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
    if (rand_pos){
      max_D_rescale_B[b] = max(max_D_s_t_rescale_uase_cpp(A_b_list, h_kernel, rank, directed = FALSE))
    }
    else{
      max_D_rescale_B[b] = max(max_D_s_t_rescale_fixed_uase_cpp(A_b_list, rank))
    }
    
    if (b == 1 | b%%verbose_freq == 0){
      print(paste0("b = ", b))
    }
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
    A_list = get_data_cp_sbm(TT, cp_truth, n, L, probability_1, probability_2, rand_pos)
    A_list_list = list()
    for (t in 1:TT){
      A_list_list[[t]] = A_list[, , , t]
    }
    
    #### reset seed for the method
    set.seed(n_1*n_2*L*b)
    if(rand_pos){
      result_online_cpd = online_cpd_uase_cpp(A_list_list, tau_factor, h_kernel, rank, directed = FALSE)
    }
    else{
      result_online_cpd = online_cpd_fixed_uase_cpp(A_list_list, tau_factor, rank)
    }
    
    t_hat_record[b] = result_online_cpd$t
    find_record[b] = result_online_cpd$find
    D_path_record[[b]] = result_online_cpd$D_K_t_max_rescale
    
    if (b == 1 | b%%verbose_freq == 0){
      print(paste0("b = ", b))
    }
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




uase_simulation_dirichlet <- function(TT, cp_truth, n_1, n_2, L, d, rank = 10, directed = TRUE,
                                      B = 100, N_MC = 100, verbose_freq = 10){
  
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
  A_burn = get_data_burn_dirichlet(T_burn, n_1, n_2, L, dirichlet_xy_1, W_1, directed)
  
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
    max_D_rescale_B[b] = max(max_D_s_t_rescale_uase_cpp(A_b_list, h_kernel, rank, directed))
    
    if (b == 1 | b%%verbose_freq == 0){
      print(paste0("b = ", b))
    }
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
    A_list = get_data_cp_dirichlet(TT, cp_truth, n_1, n_2, L, dirichlet_xy_1, W_1, dirichlet_xy_2, W_2, directed)
    A_list_list = list()
    for (t in 1:TT){
      A_list_list[[t]] = A_list[, , , t]
    }
    
    #### reset seed for the method
    set.seed(n_1*n_2*L*d*b)
    result_online_cpd = online_cpd_uase_cpp(A_list_list, tau_factor, h_kernel, rank, directed)
    
    t_hat_record[b] = result_online_cpd$t
    find_record[b] = result_online_cpd$find
    D_path_record[[b]] = result_online_cpd$D_K_t_max_rescale
    
    if (b == 1 | b%%verbose_freq == 0){
      print(paste0("b = ", b))
    }
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

multi_simulation_sbm <- function(TT, cp_truth, n, L, rank = 10, rand_pos = FALSE, 
                                 B = 100, N_MC = 100, verbose_freq = 10){
  
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
  A_burn = get_data_burn_sbm(T_burn, n, L, probability_1, rand_pos)
  
  max_D_rescale_B = rep(0, B)
  tic()
  for(b in 1:B){
    set.seed(n_1*n_2*L*b)
    b_ix = sample(1:T_burn, replace = FALSE)
    A_b = A_burn[, , , b_ix]
    
    set.seed(n_1*n_2*L*b)
    
    if (rand_pos){
      max_D_rescale_B[b] = max(max_D_s_t_rescale_multi(A_b, h_kernel, rank, directed = FALSE))
    }
    else{
      max_D_rescale_B[b] = max(max_D_s_t_rescale_fixed_multi(A_b, rank))
    }
    
    if (b == 1 | b%%verbose_freq == 0){
      print(paste0("b = ", b))
    }
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
    A_list = get_data_cp_sbm(TT, cp_truth, n, L, probability_1, probability_2, rand_pos)
    
    #### reset seed for the method
    set.seed(n_1*n_2*L*b)
    
    if (rand_pos){
      result_online_cpd = online_cpd_multi(A_list, tau_factor, h_kernel, rank, directed = FALSE)
    }
    else{
      result_online_cpd = online_cpd_fixed_multi(A_list, tau_factor, rank)
    }
    
    t_hat_record[b] = result_online_cpd$t
    find_record[b] = result_online_cpd$find
    D_path_record[[b]] = result_online_cpd$D_K_t_max_rescale
    
    if (b == 1 | b%%verbose_freq == 0){
      print(paste0("b = ", b))
    }
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




multi_simulation_dirichlet <- function(TT, cp_truth, n_1, n_2, L, d, rank = 10, directed = TRUE,
                                       B = 100, N_MC = 100, verbose_freq = 10){
  
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
  A_burn = get_data_burn_dirichlet(T_burn, n_1, n_2, L, dirichlet_xy_1, W_1, directed)
  
  max_D_rescale_B = rep(0, B)
  tic()
  for(b in 1:B){
    set.seed(n_1*n_2*L*d*b)
    b_ix = sample(1:T_burn, replace = FALSE)
    
    A_b = A_burn[, , , b_ix]
    
    set.seed(n_1*n_2*L*d*b)
    max_D_rescale_B[b] = max(max_D_s_t_rescale_multi(A_b, h_kernel, rank, directed))
    
    if (b == 1 | b%%verbose_freq == 0){
      print(paste0("b = ", b))
    }
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
    A_list = get_data_cp_dirichlet(TT, cp_truth, n_1, n_2, L, dirichlet_xy_1, W_1, dirichlet_xy_2, W_2, directed)
    
    #### reset seed for the method
    set.seed(n_1*n_2*L*d*b)
    result_online_cpd = online_cpd_multi(A_list, tau_factor, h_kernel, rank, directed)
    
    t_hat_record[b] = result_online_cpd$t
    find_record[b] = result_online_cpd$find
    D_path_record[[b]] = result_online_cpd$D_K_t_max_rescale
    
    if (b == 1 | b%%verbose_freq == 0){
      print(paste0("b = ", b))
    }
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

knn_simulation_sbm <- function(TT, cp_truth, n, L, rand_pos = FALSE, 
                               N_MC = 100, verbose_freq = 10){
  
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
  A_burn = get_data_burn_sbm(T_burn, n, L, probability_1, rand_pos)
  
  ### online cpd ####
  
  tic()
  
  t_hat_record = rep(NA, N_MC)
  find_record = rep(FALSE, N_MC)
  
  for (b in 1:N_MC){
    set.seed(n_1*n_2*L*b)
    b_ix = sample(1:T_burn, replace = FALSE)
    A_b = A_burn[, , , b_ix]
    
    set.seed(n_1*n_2*L*b)
    A_list = get_data_cp_sbm(TT, cp_truth, n, L, probability_1, probability_2, rand_pos)
    
    #### reset seed for the method
    set.seed(n_1*n_2*L*b)
    result_online_cpd = online_cpd_knn(A_b, A_list, k_nn, alpha, n0, n1, ARL)
    
    t_hat_record[b] = result_online_cpd$t
    find_record[b] = result_online_cpd$find

    if (b == 1 | b%%verbose_freq == 0){
      print(paste0("b = ", b))
    }
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




knn_simulation_dirichlet <- function(TT, cp_truth, n_1, n_2, L, d, directed = TRUE,
                                     N_MC = 100, verbose_freq = 10){
  
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
  A_burn = get_data_burn_dirichlet(T_burn, n_1, n_2, L, dirichlet_xy_1, W_1, directed)
  
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
    A_list = get_data_cp_dirichlet(TT, cp_truth, n_1, n_2, L, dirichlet_xy_1, W_1, dirichlet_xy_2, W_2, directed)
    
    #### reset seed for the method
    set.seed(n_1*n_2*L*d*b)
    result_online_cpd = online_cpd_knn(A_b, A_list, k_nn, alpha, n0, n1, ARL)
    
    t_hat_record[b] = result_online_cpd$t
    find_record[b] = result_online_cpd$find
    
    if (b == 1 | b%%verbose_freq == 0){
      print(paste0("b = ", b))
    }
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

