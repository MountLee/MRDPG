# args = (commandArgs(TRUE))
# 
# for (i in 1:length(args)){
#   eval(parse(text = args[[i]]))
# }


#### set the method to run ####
# method_list = c("hosvd", "thpca", "thpca_r1", "uase", "multi", "knn", "twist")
# method_list = c("hosvd")
method_list = c("twist")
# method_list = c("hosvd", "thpca", "uase", "multi", "twist")

#### server ####
path = "~/R_code/MRDPG-main/"
source(paste0(path, "generate_data.R"))
source(paste0(path, "simulation_wrap.R"))
Rcpp::sourceCpp(paste0(path, "cpd_hpca.cpp"))
Rcpp::sourceCpp(paste0(path, "cpd_uase.cpp"))
source(paste0(path, "other_methods.R"))

#### local ####
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# source("generate_data.R")
# source("simulation_wrap.R")
# Rcpp::sourceCpp("cpd_hpca.cpp")
# Rcpp::sourceCpp("cpd_uase.cpp")
# source("other_methods.R")


#### main ####

library(mvtnorm)
library(dirmult)
library(rTensor)
library(Rcpp)
library(tictoc)

library(multiness)
library(gStream)




rand_pos = FALSE
rand_method = FALSE


TT = 100
cp_truth = 50

cpd_sbm(TT, cp_truth, n_list, L_list, rand_pos, 
        B = 10, N_MC = 10, verbose_freq = 10, alpha = 0.01,
        rand_method = rand_method, c_tau = c_tau,
        method_list = method_list)


n_list = c(50)
L_list = c(2)
c_tau = 0.15
n = n_list[1]
L = L_list[1]
rank = rep(10, 3)
tau_factor = c_tau * (rank[1]^2 * rank[3] + n * rank[1] + L * rank[3])^0.5
tau_factor

twist_res = twist_simulation_sbm(TT, cp_truth, n, L, rank = rep(10, 3), rand_pos, 
                                 B = 10, N_MC = 10, verbose_freq = 10, alpha = 0.01, 
                                 rand_method=rand_method, c_tau=c_tau)

twist_res[[4]]




a=1
b=c(1,2)
x=list(a, b)
x[[1]]


