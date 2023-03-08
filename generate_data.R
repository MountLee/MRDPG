#' generate tensor based on dirichlet distribution
#' @param n_1 dimension 1 of the tensor
#' @param n_2 dimension 2 of the tensor
#' @param d dimension of the latent space
#' @param L dimension 3 of the tensor, i.e., number of layers
#' @return (n_1, n_2, L)-shaped tensor
#' @export 
generate_tensor_dirichlet <- function(n_1, n_2, L, W,
                                      dirichlet_x, dirichlet_y){
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
  
  # Y.tensor =  as.tensor(A)
  # hat.rank =c(10,10,10)
  # return(list(Y.tensor = Y.tensor, probability = probability, hat.rank = hat.rank))
  return(A)
}


get_dirichlet_params <- function(n_1, n_2, L, d){
  
  # set.seed(n_1*n_2*L*d)
  
  dirichlet_xy_1 = list(x = rep(1, d),
                        y = rep(10, d))
  dirichlet_xy_2 = list(x = rep(500, d),
                        y = rep(1000, d))
  
  W_1 =  array(NA,c(d, d, L))
  prob = seq(0,1,  1/(2*L))
  for (layer in 1: L)
  {
    Sigma = matrix(runif(d^2,prob[L+layer], prob[L+layer+1]), ncol=d) 
    W_1[, , layer] = Sigma
  }
  
  W_2 =  array(NA,c(d, d, L))
  prob = seq(0,1,  1/(2*L))
  for (layer in 1: L)
  {
    Sigma = matrix(runif(d^2,prob[L+layer], prob[L+layer+1]), ncol=d) 
    W_2[, , layer] = Sigma
  }
  
  return(list(dirichlet_xy_1, dirichlet_xy_2, W_1, W_2))
}



get_data_burn_dirichlet <- function(T_burn, n_1, n_2, L, dirichlet_xy_1, W_1){
  ### data for burn-in ####
  A_burn = array(0, c(n_1, n_2, L, T_burn))
  for (t in 1:T_burn){
    A_burn[, , , t] =  generate_tensor_dirichlet(n_1, n_2, L, W_1,
                                                 dirichlet_x = dirichlet_xy_1$x, 
                                                 dirichlet_y = dirichlet_xy_1$y)
  }
  
  return(A_burn)
}




get_data_cp_dirichlet <- function(TT, cp_truth, n_1, n_2, L, dirichlet_xy_1, W_1, dirichlet_xy_2, W_2){
  ### data with cp ####

  A_list = array(0, c(n_1, n_2, L, TT))
  for (t in 1:cp_truth){
    A_list[, , , t] = generate_tensor_dirichlet(n_1, n_2,  L, W_1,
                                                dirichlet_x = dirichlet_xy_1$x, 
                                                dirichlet_y = dirichlet_xy_1$y)
  }
  for (t in (cp_truth + 1):TT){
    A_list[, , , t] = generate_tensor_dirichlet(n_1, n_2, L, W_2,
                                                dirichlet_x = dirichlet_xy_2$x, 
                                                dirichlet_y = dirichlet_xy_2$y)
  }

  return(A_list)
}



get_sbm_params <- function(n, L){
  
  # set.seed(n*L)
  
  probability_1 = array(NA, c(n, n, L))
  probability_2 = array(NA, c(n, n, L))
  
  prob = seq(0,1,  1/(4*L)) 
  for (layer in 1: L)
  {
    p_1 =  runif(1, prob[2*L +layer], prob[2*L+layer+1])
    p_2 = runif(1, prob[3*L+layer], prob[3*L+layer+1])
    P =  matrix(p_1,n,n)
    P[1:floor(n/4), 1:floor(n/4)] = p_2
    P[(1+floor(n/4)):(2*floor(n/4)),(1+floor(n/4)):(2*floor(n/4)) ] = p_2
    P[(1+2*floor(n/4)):(3*floor(n/4)),(1+2*floor(n/4)):(3*floor(n/4)) ] = p_2
    P[(1+3*floor(n/4)):n,(1+3*floor(n/4)):n ] = p_2
    probability_1[, , layer] = P
    probability_2[, , L - layer + 1] = P
  }
  
  return(list(probability_1, probability_2))
}



generate_tensor_probability <- function(n_1, n_2, L, probability){
  
  dim_ = c(n_1, n_2, L)
  A = array(NA,dim_)
  
  for (layer in 1: L)
  {
    A[, , layer] = matrix(rbinom(matrix(1,n_1,n_2),matrix(1,n_1,n_2), probability[ , , layer]),n_1,n_2)
  }
  
  # Y.tensor =  as.tensor(A)
  # hat.rank =c(10,10,10)
  # return(list(Y.tensor = Y.tensor, probability = probability, hat.rank = hat.rank))
  return(A)
}


get_data_burn_sbm <- function(T_burn, n, L, probability_mat){

  A_burn = array(0, c(n, n, L, T_burn))
  for (t in 1:T_burn){
    A_burn[, , , t] = generate_tensor_probability(n, n, L, probability_mat)
  }
  
  return(A_burn)
}




get_data_cp_sbm <- function(TT, cp_truth, n, L, probability_1, probability_2){
  #### data with cp
  A_list = array(0, c(n, n, L, TT))
  for (t in 1:cp_truth){
    A_list[, , , t] = generate_tensor_probability(n, n, L, probability_1)
  }
  for (t in (cp_truth + 1):TT){
    A_list[, , , t] = generate_tensor_probability(n, n, L, probability_2)
  }
  
  return(A_list)
}






#### n_1 = n_2
# generate_sbm <- function(TT, cp_truth, n, L, N_MC = 100){
#   
#   set.seed(n*L)
#   
#   n_1 = n
#   n_2 = n
#   
#   probability_1 = array(NA, c(n_1, n_2, L))
#   probability_2 = array(NA, c(n_1, n_2, L))
#   
#   prob = seq(0,1,  1/(4*L)) 
#   for (layer in 1: L)
#   {
#     p_1 =  runif(1, prob[2*L +layer], prob[2*L+layer+1])
#     p_2 = runif(1, prob[3*L+layer], prob[3*L+layer+1])
#     P =  matrix(p_1,n,n)
#     P[1:floor(n/4), 1:floor(n/4)] = p_2
#     P[(1+floor(n/4)):(2*floor(n/4)),(1+floor(n/4)):(2*floor(n/4)) ] = p_2
#     P[(1+2*floor(n/4)):(3*floor(n/4)),(1+2*floor(n/4)):(3*floor(n/4)) ] = p_2
#     P[(1+3*floor(n/4)):n,(1+3*floor(n/4)):n ] = p_2
#     probability_1[, , layer] = P
#     probability_2[, , L - layer + 1] = P
#   }
#   
#   
#   T_burn = as.integer(cp_truth*1.5)
#   A_burn = array(0, c(n_1, n_2, L, T_burn))
#   for (t in 1:T_burn){
#     A_burn[, , , t] = generate_tensor_probability(n_1, n_2, L, probability_1)$Y.tensor@data
#   }
#   
#   #### data with cp ####
#   A_list_all = array(0, c(n_1, n_2, L, TT, N_MC))
#   
#   for (b in 1:N_MC){
#     A_list = array(0, c(n_1, n_2, L, TT))
#     for (t in 1:cp_truth){
#       A_list[, , , t] = generate_tensor_probability(n_1, n_2, L, probability_1)$Y.tensor@data
#     }
#     for (t in (cp_truth + 1):TT){
#       A_list[, , , t] = generate_tensor_probability(n_1, n_2, L, probability_2)$Y.tensor@data
#     }
#     A_list_all[, , , , b] = A_list
#   }
#   
#   return(list(A_burn, A_list_all))
# }





#### backup







# generate_dirichlet <- function(n_list_1, n_list_2, L_list, d_list){
#   
#   for (n_index_1 in 1:length(n_list_1)) {
#     
#     n_1 = n_list_1[n_index_1]
#     print(paste0("n_1 = ", n_1))
#     
#     for (n_index_2 in 1:length(n_list_2)) {  
#       
#       n_2 = n_list_2[n_index_2]
#       print(paste0("n_2 = ", n_2))
#       
#       for (L_index in 1:length(L_list)){
#         
#         L =  L_list[L_index]
#         print(paste0("L = ", L))
#         
#         for (d_index in 1:length(d_list)){
#           
#           d =  d_list[d_index]
#           print(paste0("d = ", d))
#           
#           set.seed(n_1*n_2*L*d)
#           
#           dirichlet_xy_1 = list(x = rep(1, d),
#                                 y = rep(10, d))
#           dirichlet_xy_2 = list(x = rep(500, d),
#                                 y = rep(1000, d))
#           
#           
#           
#           W_1 =  array(NA,c(d, d, L))
#           prob = seq(0,1,  1/(2*L))
#           for (layer in 1: L)
#           {
#             Sigma = matrix(runif(d^2,prob[L+layer], prob[L+layer+1]), ncol=d) 
#             W_1[, , layer] = Sigma
#           }
#           
#           W_2 =  array(NA,c(d, d, L))
#           prob = seq(0,1,  1/(2*L))
#           for (layer in 1: L)
#           {
#             Sigma = matrix(runif(d^2,prob[L+layer], prob[L+layer+1]), ncol=d) 
#             W_2[, , layer] = Sigma
#           }
#           
#           ### burn-in data ###
#           
#           T_burn = as.integer(cp_truth*1.5)
#           A_burn = array(0, c(n_1, n_2, L, T_burn))
#           for (t in 1:T_burn){
#             A_burn[, , , t] =  generate_tensor_dirichlet(n_1, n_2, L, W_1,
#                                                          dirichlet_x = dirichlet_xy_1$x, 
#                                                          dirichlet_y = dirichlet_xy_1$y)$Y.tensor@data
#           }
#           
#           
#           B = 100
#           
#           ### data with cp ####
#           N_MC = 100
#           
#           for (b in 1:N_MC){
#             A_list = array(0, c(n_1, n_2, L, TT))
#             for (t in 1:cp_truth){
#               A_list[, , , t] = generate_tensor_dirichlet(n_1, n_2,  L, W_1,
#                                                           dirichlet_x = dirichlet_xy_1$x, 
#                                                           dirichlet_y = dirichlet_xy_1$y)$Y.tensor@data
#             }
#             for (t in (cp_truth + 1):TT){
#               A_list[, , , t] = generate_tensor_dirichlet(n_1, n_2, L, W_2,
#                                                           dirichlet_x = dirichlet_xy_2$x, 
#                                                           dirichlet_y = dirichlet_xy_2$y)$Y.tensor@data
#             }
#           }
#           
#         }
#         
#       }
#     }
#   }
#   
#   return(list(A_burn, A_list))
# }



