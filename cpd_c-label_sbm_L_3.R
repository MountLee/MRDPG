n_list = c(50)
L_list = c(3)

TT = 100
cp_truth = 50

#### set the type of the change point ####
## only community labels change at the change point
rand_pos = "given"
prob_change = FALSE
n_c = c(4, 4)
flip_layer = FALSE


# ## only the number of communities changes at the change point
# rand_pos = FALSE
# prob_change = TRUE
# n_c=c(4, 8)
# flip_layer=FALSE


#### set the method to run ####
method_list = c("hosvd", "thpca", "thpca_r1", "uase", "multi", "knn", "twist")
# method_list = c("hosvd")
# method_list = c("twist")

#### server ####
# path = "/home/u2042553/"
# source(paste0(path, "generate_data.R"))
# source(paste0(path, "simulation_wrap.R"))
# Rcpp::sourceCpp(paste0(path, "cpd_hpca.cpp"))
# Rcpp::sourceCpp(paste0(path, "cpd_uase.cpp"))
# source(paste0(path, "other_methods.R"))

#### local ####
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
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



cpd_sbm(TT, cp_truth, n_list, L_list, rand_pos, 
        B = 100, N_MC = 100, verbose_freq = 10, 
        prob_change = prob_change, n_c = n_c, flip_layer = flip_layer,
        method_list = method_list)
