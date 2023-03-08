# MRDPG
This repository contains code and data used in the paper *Multilayer Random Dot Product Graphs: Dynamic, Nonparametric Estimation and Online Change Point Detection* by Fan Wang, Wanshan Li, Oscar Hernan Madrid Padilla, Yi Yu, and Alessandro Rinaldo.

All experiments and analysis are conducted in R (version >= 4.1).

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
