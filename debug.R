
n = 50
L = 2

n = 25
L = 3
rank = rep(10, 3)
c_tau = 0.16
tau_factor = c_tau * (rank[1]^2 * rank[3] + n * rank[1] + L * rank[3])^0.5
tau_factor

tau_factor = c_tau * (rank[1]^2 * rank[3] + n * rank[1] + L * rank[3])^0.5
tau_factor









save("df_1", "df_2", "plot_1", "plot_2", 
     file = "sbm_scenario-1_df-plot.RData")



a = get_blockwise_const_mat(11, 3, 0.1, 0.5)
a



a = c("a", "b", "c")

"a" %in% a

if ("a" %in% a){
  print("nihao")
}


a = NULL

is.null(a)

n = 3
L = 2
dim_ = c(n, n, L)
A = array(NA, dim_)
for (i in 1:L){
  A[,,i] = matrix(c(1:(n*n)), ncol=n)
}


A
ix = sample(1:n, n, replace = FALSE)
A[ix, ix, ]



