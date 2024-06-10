# MRDPG
This repository contains code and data used in the paper *Multilayer Random Dot Product Graphs: Dynamic, Estimation and Online Change Point Detection* by Fan Wang, Wanshan Li, Oscar Hernan Madrid Padilla, Yi Yu, and Alessandro Rinaldo.

All experiments and analysis are conducted in R (version >= 4.3).

### Experiments
#### Estimation
In folder "estimation".
+ Scripts: Simulation_Scenario_1.R, Simulation_Scenario_2.R
+ Utilities: tensor_functions.R

#### Change point detection: simulation
Standard online-CPD setting in the paper
+ SBM: cpd_sbm.sh
+ Dirichlet: cpd_dirichlet_direct.sh, cpd_dirichlet_undirect.sh

SBM with various types of change points:
+ Change in the number of communities: cpd_sbm_c-k.sh
+ Change in node labels of communities: cpd_sbm_c-labels.sh

Other settings:
+ Robustness check: cpd-robustness_sbm_fixed-pos.sh, cpd-robustness_sbm_rand-pos.sh
+ Fixed $C_{\tau}$: cpd_sbm_c-tau_fixed.sh

#### Change point detection: read data
In folder "data".

### Dependencies
Some R-packages are required to reproduce all results. To install them, one can run the following lines:

```
list.of.packages <- c("dirmult", "mvtnorm", "Rcpp", "devtools", "gStream", "multiness",
                      "changepoints", "MASS", "rTensor", "tictoc")
new_packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
for (package in new_packages){
  install.packages(package, repos='http://cran.us.r-project.org', quiet = TRUE)
}


# for the method HOOI
devtools::install_github("Rungang/STATSVD")
install.packages("https://cran.r-project.org/src/contrib/Archive/ssvd/ssvd_1.0.tar.gz", 
                 repos = NULL, type="source")

```

To run the experiments related to the main method, TH-PCA, in the paper only, it suffices to run:
```
list.of.packages <- c("dirmult", "mvtnorm", "Rcpp", "tictoc")
new_packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
for (package in new_packages){
  install.packages(package, repos='http://cran.us.r-project.org', quiet = TRUE)
}
```
