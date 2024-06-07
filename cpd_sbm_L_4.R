n_list = c(50)
L_list = c(4)
rand_pos = FALSE


TT = 100
cp_truth = 50


#### set the method to run ####
# method_list = c("hosvd", "thpca", "thpca_r1", "uase", "multi", "knn", "twist")
method_list = c("hosvd")
# method_list = c("twist")

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



method_list = c("thpca", "uase", "twist")
cpd_sbm(TT, cp_truth, n_list, L_list, rand_pos, 
        B = 100, N_MC = 100, verbose_freq = 10, 
        method_list = method_list)
