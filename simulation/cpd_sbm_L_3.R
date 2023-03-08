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


n_list = c(50)
L_list = c(3)


TT = 100
cp_truth = 50


type_I_error = array(NA, dim=c(length(n_list),  length(L_list)))
type_II_error = array(NA, dim=c(length(n_list), length(L_list)))
run_length = array(NA, dim=c(length(n_list), length(L_list)))

type_I_error_thpca = array(NA, dim=c(length(n_list),  length(L_list)))
type_II_error_thpca = array(NA, dim=c(length(n_list), length(L_list)))
run_length_thpca = array(NA, dim=c(length(n_list), length(L_list)))

type_I_error_uase = array(NA, dim=c(length(n_list),  length(L_list)))
type_II_error_uase = array(NA, dim=c(length(n_list), length(L_list)))
run_length_uase = array(NA, dim=c(length(n_list), length(L_list)))

type_I_error_multi = array(NA, dim=c(length(n_list),  length(L_list)))
type_II_error_multi = array(NA, dim=c(length(n_list), length(L_list)))
run_length_multi = array(NA, dim=c(length(n_list), length(L_list)))

type_I_error_knn = array(NA, dim=c(length(n_list),  length(L_list)))
type_II_error_knn = array(NA, dim=c(length(n_list), length(L_list)))
run_length_knn = array(NA, dim=c(length(n_list), length(L_list)))

for (n_index in 1:length(n_list)) {
  
  n = n_list[n_index]
  print(paste0("n = ", n))
  
  for (L_index in 1:length(L_list)){
    
    L =  L_list[L_index]
    print(paste0("L = ", L))
    
    print("---- hpca ----")
    hpca_res = hpca_simulation_sbm(TT, cp_truth, n, L, B = 100, N_MC = 100)
    type_I_error[n_index, L_index] = hpca_res[1]
    type_II_error[n_index, L_index] = hpca_res[2]
    run_length[n_index, L_index] = hpca_res[3]
    print(hpca_res)
    
    print("---- thpca ----")
    thpca_res = thpca_simulation_sbm(TT, cp_truth, n, L, B = 100, N_MC = 100)
    type_I_error_thpca[n_index, L_index] = thpca_res[1]
    type_II_error_thpca[n_index, L_index] = thpca_res[2]
    run_length_thpca[n_index, L_index] = thpca_res[3]
    print(thpca_res)
    
    print("---- uase ----")
    uase_res = uase_simulation_sbm(TT, cp_truth, n, L, B = 100, N_MC = 100)
    type_I_error_uase[n_index, L_index] = uase_res[1]
    type_II_error_uase[n_index, L_index] = uase_res[2]
    run_length_uase[n_index, L_index]= uase_res[3]
    print(uase_res)
    
    print("---- multiness ----")
    multi_res = multi_simulation_sbm(TT, cp_truth, n, L, B = 100, N_MC = 100)
    type_I_error_multi[n_index, L_index] = multi_res[1]
    type_II_error_multi[n_index, L_index] = multi_res[2]
    run_length_multi[n_index, L_index] = multi_res[3]
    print(multi_res)
    
    print("---- knn ----")
    knn_res = knn_simulation_sbm(TT, cp_truth, n, L, N_MC = 100)
    type_I_error_knn[n_index, L_index] = knn_res[1]
    type_II_error_knn[n_index, L_index] = knn_res[2]
    run_length_knn[n_index, L_index] = knn_res[3]
    print(knn_res)
    
    
  }
}

type_I_error
type_I_error_thpca
type_I_error_uase
type_I_error_multi
type_I_error_knn
# 
type_II_error
type_II_error_thpca
type_II_error_uase
type_II_error_multi
type_II_error_knn
# 
run_length
run_length_thpca
run_length_uase
run_length_multi
run_length_knn
