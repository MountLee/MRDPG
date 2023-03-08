# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(mvtnorm)
library(dirmult)
library(rTensor)
library(Rcpp)
library(tictoc)

library(multiness)
library(gStream)

source("/home/u2042553/generate_data.R")
source("/home/u2042553/simulation_wrap.R")
Rcpp::sourceCpp("/home/u2042553/cpd_hpca.cpp")
Rcpp::sourceCpp("/home/u2042553/cpd_uase.cpp")
source("/home/u2042553/other_methods.R")


n_list_1 = c(25)
n_list_2 = c(50)
L_list = c(3)
d_list = c(4)


TT = 100
cp_truth = 50



type_I_error = array(NA, dim=c(length(n_list_1),  length(n_list_2), length(L_list), length(d_list)))
type_II_error = array(NA, dim=c(length(n_list_1), length(n_list_2), length(L_list), length(d_list)))
run_length = array(NA, dim=c(length(n_list_1), length(n_list_2), length(L_list), length(d_list)))

type_I_error_thpca =  array(NA, dim=c(length(n_list_1),  length(n_list_2), length(L_list), length(d_list)))
type_II_error_thpca = array(NA, dim=c(length(n_list_1),  length(n_list_2), length(L_list), length(d_list)))
run_length_thpca =  array(NA, dim=c(length(n_list_1),  length(n_list_2), length(L_list), length(d_list)))

type_I_error_uase = array(NA, dim=c(length(n_list_1),  length(n_list_2), length(L_list), length(d_list)))
type_II_error_uase = array(NA, dim=c(length(n_list_1), length(n_list_2), length(L_list), length(d_list)))
run_length_uase = array(NA, dim=c(length(n_list_1), length(n_list_2), length(L_list), length(d_list)))

#type_I_error_multi = array(NA, dim=c(length(n_list_1),  length(n_list_2), length(L_list), length(d_list)))
#type_II_error_multi = array(NA, dim=c(length(n_list_1), length(n_list_2), length(L_list), length(d_list)))
#run_length_multi = array(NA, dim=c(length(n_list_1), length(n_list_2), length(L_list), length(d_list)))

type_I_error_knn = array(NA, dim=c(length(n_list_1),  length(n_list_2), length(L_list), length(d_list)))
type_II_error_knn = array(NA, dim=c(length(n_list_1), length(n_list_2), length(L_list), length(d_list)))
run_length_knn = array(NA, dim=c(length(n_list_1), length(n_list_2), length(L_list), length(d_list)))



for (n_index_1 in 1:length(n_list_1)) {
  
  n_1 = n_list_1[n_index_1]
  print(paste0("n_1 = ", n_1))
  
  for (n_index_2 in 1:length(n_list_2)) {  
    
    n_2 = n_list_2[n_index_2]
    print(paste0("n_2 = ", n_2))
    
    for (L_index in 1:length(L_list)){
      
      L =  L_list[L_index]
      print(paste0("L = ", L))
      
      for (d_index in 1:length(d_list)){
        
        d =  d_list[d_index]
        print(paste0("d = ", d))
        
        print("---- hosvd ----")
        hpca_res = hpca_simulation_dirichlet(TT, cp_truth, n_1, n_2, L, d, B = 100, N_MC = 100)
        type_I_error[n_index_1, n_index_2, L_index, d_index] = hpca_res[1]
        type_II_error[n_index_1, n_index_2, L_index, d_index] = hpca_res[2]
        run_length[n_index_1, n_index_2, L_index, d_index] = hpca_res[3]
        print(hpca_res)
        
        print("---- thpca ----")
        thpca_res = thpca_simulation_dirichlet(TT, cp_truth, n_1, n_2, L, d, B = 100, N_MC = 100)
        type_I_error_thpca[n_index_1, n_index_2, L_index, d_index] = thpca_res[1]
        type_II_error_thpca[n_index_1, n_index_2, L_index, d_index] = thpca_res[2]
        run_length_thpca[n_index_1, n_index_2, L_index, d_index] = thpca_res[3]
        print(thpca_res)
        
        print("---- uase ----")
        uase_res = uase_simulation_dirichlet(TT, cp_truth, n_1, n_2, L, d, B = 100, N_MC = 100)
        type_I_error_uase[n_index_1, n_index_2, L_index, d_index] = uase_res[1]
        type_II_error_uase[n_index_1, n_index_2, L_index, d_index] = uase_res[2]
        run_length_uase[n_index_1, n_index_2, L_index, d_index]= uase_res[3]
        print(uase_res)
        
        #print("---- multiness ----")
        #multi_res = multi_simulation_dirichlet(TT, cp_truth, n_1, n_2, L, d, B = 100, N_MC = 100)
        #type_I_error_multi[n_index_1, n_index_2, L_index, d_index] = multi_res[1]
        #type_II_error_multi[n_index_1, n_index_2, L_index, d_index] = multi_res[2]
        #run_length_multi[n_index_1, n_index_2, L_index, d_index] = multi_res[3]
        #print(multi_res)
        
        print("---- knn ----")
        knn_res = knn_simulation_dirichlet(TT, cp_truth, n_1, n_2, L, d, N_MC = 100)
        type_I_error_knn[n_index_1, n_index_2, L_index, d_index] = knn_res[1]
        type_II_error_knn[n_index_1, n_index_2, L_index, d_index] = knn_res[2]
        run_length_knn[n_index_1, n_index_2, L_index, d_index] = knn_res[3]
        print(knn_res)
        # 
        
      }
    }
  }
}


type_I_error
type_I_error_thpca
type_I_error_uase
#type_I_error_multi
type_I_error_knn
# 
type_II_error
type_II_error_thpca
type_II_error_uase
#type_II_error_multi
type_II_error_knn
# 
run_length
run_length_thpca
run_length_uase
#run_length_multi
run_length_knn
