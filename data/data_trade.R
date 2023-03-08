

library("rTensor")
library("STATSVD")
library("Rlab")
library("ggplot2")
library("egg")

library(multiness)
library(mvtnorm)
library(dirmult)
library(STATSVD)

library(devtools)
library(gStream)
library(changepoints)
source("/Users/chris/MRDPG/tensor_functions.R")
source("/Users/chris/MRDPG/cpd_MRDPG_functions.R")

library(ggplot2)
library(tictoc)



require(data.table)
data_whole = fread("/Users/chris/MRDPG/Trade_DetailedTradeMatrix_E_All_Data.csv", header = T)

data_whole = as.matrix(data_whole)
unique(data_whole[, 8])

#2: export; 4: import
data_whole_ex = data_whole[which(data_whole[, 8]== unique(data_whole[, 8])[4]), ]

country_report = rle(sort(data_whole_ex[,2]))
country_report_100 = country_report$value[order(country_report$lengths, decreasing = T)[1:75]]

country_partner= rle(sort(data_whole_ex[,4]))
country_partner_100 = country_partner$value[order(country_partner$lengths, decreasing = T)[1:75]]

product_sort = rle(sort(data_whole_ex[,6]))
product_sort_4 = product_sort$value[order(product_sort$lengths, decreasing = T)[1:4]]


data_names = c()


for (i in 1:35 ) {
  data_names = c(data_names, paste("data_year", as.character(i+1985), sep="."))
}
                 



dim = c(length(country_report_100), length(country_partner_100), length(product_sort_4), length(data_names) )



data_year_product_tensor = array(rep(NA, length(country_report_100)*length(country_partner_100)*length(product_sort_4))*length(data_names), dim)



for (u in 1:length(data_names)){

year = colnames(data_whole_ex)[ 10+2*(u-1)]

data_year = data_whole[, c(2, 4, 6, 8, 10+2*(u-1))]

data_year[is.na(data_year[, year]) , year] = 0

for(l in 1:length(product_sort_4)){
  
  data_year_product = data_year[which( data_year[, "Item"] == product_sort_4[l] ), ]

  data_year_matrix =  matrix(rep(NA, length(country_report_100)*length(country_partner_100)), ncol = length(country_partner_100))

  for(i in 1:length(country_report_100)){
    for(j in 1:length(country_partner_100)){
      x = which(data_year_product[, "Reporter Countries" ]== country_report_100[i] & 
                  data_year_product[, "Partner Countries"] ==  country_partner_100[j] )
      if ( length(x)!=0 ){
        data_year_matrix[i,j] = 1
      }
      else{
        data_year_matrix[i,j] = 0
      }
    }
  }
  
  
  data_year_product_tensor[, , l, u] = data_year_matrix 
}

}


B = 100

K_max = 500

n_1 = dim(data_year_product_tensor)[1]

n_2 = dim(data_year_product_tensor)[2]

L = dim(data_year_product_tensor)[3]

time = dim(data_year_product_tensor)[4]

T_burn = 10

TT = time - T_burn

h_kernel = (K_max*log(TT*n_1* n_2)/(TT*n_1*n_2))^{1/L}
h_kernel

hat.rank = c(10, 10, 10)
max_D_rescale_B = rep(0, B)
max_D_rescale_B_thpca = rep(0, B)
max_D_rescale_B_uase = rep(0, B)
max_D_rescale_B_multi = rep(0, B)
alpha = 0.05

tic()



S_tilde_p = c(rep(0, n_1 + n_2 - 1))

for (p in 1:(n_2)){
  S_tilde_p[p] = min(n_1, n_2 + 1 - p)
}
for (p in (n_2 + 1):(n_2 + n_1 - 1)){
  S_tilde_p[p] = min(n_2, n_2 + n_1 - p)
}

set.seed(0) 

for(b in 1:B){
  b_ix = sample(1:T_burn, replace = FALSE)
  A_b =  data_year_product_tensor[, , , b_ix]
  max_D_rescale_B[b] = max_D_s_t_rescale(A_b, h_kernel, verbose = FALSE)
  max_D_rescale_B_thpca[b] = max_D_s_t_rescale_thpca(A_b, h_kernel, hat.rank, verbose = FALSE)
  max_D_rescale_B_uase[b] = max_D_s_t_rescale_uase(A_b, h_kernel, hat.rank[1], verbose = FALSE)
  max_D_rescale_B_multi[b] = max_D_s_t_rescale_multi(A_b, h_kernel, hat.rank[1], verbose = FALSE)
  print(paste0("b = ", b))
}
toc()


tau_factor = quantile(max_D_rescale_B, 1 - alpha, type = 1)
tau_factor_thpca = quantile(max_D_rescale_B_thpca, 1 - alpha, type = 1)
tau_factor_uase = quantile(max_D_rescale_B_uase, 1 - alpha, type = 1)
tau_factor_multi = quantile(max_D_rescale_B_multi, 1 - alpha, type = 1)


A_list  =  data_year_product_tensor[, , , (T_burn+1) : time]

result_online_cpd = online_cpd(A_list, tau_factor, h_kernel,  verbose = FALSE)
result_online_cpd_thpca = online_cpd_thpca(A_list, tau_factor, h_kernel, hat.rank, verbose = FALSE)
result_online_cpd_uase = online_cpd_uase(A_list, tau_factor, h_kernel, hat.rank[1], verbose = FALSE)
result_online_cpd_multi = online_cpd_uase(A_list, tau_factor, h_kernel, hat.rank[1], verbose = FALSE)


k_nn = 3
n0 = floor(0.3 * T_burn)
n1 = floor(0.7 * T_burn)
ARL = 2000

A_all = data_year_product_tensor
n_all = dim(A_all)[4]

dist_mat_all = get_dist_mat(A_all, distance_Y)
diag(dist_mat_all) = max(dist_mat_all) + 100

t_hat = tryCatch(
  expr = {
    res = gstream(dist_mat_all, L = T_burn, N0 = T_burn, k_nn, statistics = "m",
                  n0, n1, ARL, alpha, skew.corr = TRUE, asymp = FALSE)
    tauhat = res$tauhat$max.type
    tauhat
  },
  error = function(cond){return (Inf)}
)
t_hat
length(t_hat)
min(t_hat)

data_names[T_burn + result_online_cpd$t]
data_names[T_burn + result_online_cpd_thpca$t]
data_names[T_burn + result_online_cpd_uase$t]
data_names[T_burn + result_online_cpd_multi$t]

