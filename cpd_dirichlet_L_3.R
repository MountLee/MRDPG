n_list_1 = c(50)
n_list_2 = c(50)
L_list = c(3)
d_list = c(4)


TT = 100
cp_truth = 50


#### server ####
# path = "/home/u2042553/"
# source(paste0(path, "generate_data.R"))
# source(paste0(path, "simulation_wrap.R"))
# Rcpp::sourceCpp(paste0(path, "cpd_hpca.cpp"))
# Rcpp::sourceCpp(paste0(path, "cpd_uase.cpp"))
# source(paste0(path, "other_methods.R"))

#### local ####
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("generate_data.R")
source("simulation_wrap.R")
Rcpp::sourceCpp("cpd_hpca.cpp")
Rcpp::sourceCpp("cpd_uase.cpp")
source("other_methods.R")


#### main ####

library(mvtnorm)
library(dirmult)
library(rTensor)
library(Rcpp)
library(tictoc)

library(multiness)
library(gStream)



cpd_dirichlet(TT, cp_truth, n_list, L_list, directed = FALSE,
              B = 100, N_MC = 100, verbose_freq = 10)