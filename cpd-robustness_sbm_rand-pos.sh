R CMD BATCH "--args n_list=c(25) L_list=c(3)" cpd-robustness_sbm_rand-pos.R robust-rand-pos_n25.Rout &
R CMD BATCH "--args n_list=c(75) L_list=c(3)" cpd-robustness_sbm_rand-pos.R robust-rand-pos_n75.Rout &
R CMD BATCH "--args n_list=c(50) L_list=c(2)" cpd-robustness_sbm_rand-pos.R robust-rand-pos_L2.Rout &
R CMD BATCH "--args n_list=c(50) L_list=c(3)" cpd-robustness_sbm_rand-pos.R robust-rand-pos_L3.Rout &
R CMD BATCH "--args n_list=c(50) L_list=c(4)" cpd-robustness_sbm_rand-pos.R robust-rand-pos_L4.Rout &
