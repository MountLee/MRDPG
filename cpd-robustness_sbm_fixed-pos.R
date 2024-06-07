# n_list = c(25)
# L_list = c(3)
args = (commandArgs(TRUE))

for (i in 1:length(args)){
  eval(parse(text = args[[i]]))
}


TT = 100
cp_truth = 50

#### set the type of the change point and cpd method ####
# ## case I: random position but fixed method
# rand_pos = TRUE
# prob_change = TRUE
# n_c = c(4, 4)
# flip_layer = TRUE
# rand_method = FALSE

## case II: fixed position but random method
rand_pos = FALSE
prob_change = TRUE
n_c = c(4, 4)
flip_layer = TRUE
rand_method = TRUE


#### set the method to run ####
# method_list = c("hosvd", "thpca", "thpca_r1", "uase", "multi", "knn")
# method_list = c("hosvd", "thpca")
method_list = c("thpca_r1")

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



cpd_sbm(TT, cp_truth, n_list, L_list, rand_pos, 
        B = 100, N_MC = 100, verbose_freq = 10, 
        prob_change = prob_change, n_c = n_c, flip_layer = flip_layer, rand_method = rand_method,
        method_list = method_list)
