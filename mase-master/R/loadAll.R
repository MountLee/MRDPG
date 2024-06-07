library(igraph)
library(lattice) #github jesusdaniel/graphclass (only needed for the plots)
library(ggplot2)
library(Matrix)
path = "E:/R_code/cpd_MRDPG/MRDPG-main/mase-master/R/"
source(paste0(path, "mase.R"))
source(paste0(path, "parametric-bootstrap.R"))
source(paste0(path, "sample_graphs.R"))


# read getElbows function from github
library(RCurl)
# script <- getURL("https://raw.githubusercontent.com/youngser/gmmase/master/R/getElbows.R", ssl.verifypeer = FALSE)
# eval(parse(text = script))
