####setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# install.packages("rARPACK")

library(mvtnorm)
library(STATSVD)
library(changepoints)
library(egg)
library(ggplot2)
library(multiness)
library(rTensor)
# library(rMultiNet)
source("tensor_functions.R")

#### fix n, vary L ####

n_grid  =  c(100)

L_grid = c(10, 20, 30, 40)

NMC = 100

hat.rank =c(10,10,10)

result_hooi = array(NA, c( length(n_grid), length(L_grid), NMC))
result_hetero = array(NA, c( length(n_grid), length(L_grid), NMC))
result_hosvd = array(NA, c( length(n_grid), length(L_grid), NMC))
result_svd = array(NA, c( length(n_grid), length(L_grid), NMC))
result_uase = array(NA, c( length(n_grid), length(L_grid), NMC))
result_MultiNeSS_ada = array(NA, c( length(n_grid), length(L_grid), NMC))
result_flex = array(NA, c( length(n_grid), length(L_grid), NMC))
result_mase = array(NA, c( length(n_grid), length(L_grid), NMC))
result_twist = array(NA, c( length(n_grid), length(L_grid), NMC))

for(ind_n in 1:length(n_grid))
{
  n = n_grid[ind_n]
  print(paste0("n = ", n))
  
  for(ind_L in 1: length(L_grid) )
  {
    L = L_grid[ind_L]
    print(paste0("L = ", L))
        
    dim = c(n, n, L)
    prob = seq(0,1,  1/(4*L))  
    
    set.seed(n*L)  
    
    for(iter  in 1:NMC)
    {
      
      A = array(NA,dim)
      P_true = array(NA,dim)
      A_list = list()
      for (layer in 1: L)
      {
        p_1 =  runif(1, prob[2*L +layer], prob[2*L+layer+1])
        p_2 = runif(1, prob[3*L+layer], prob[3*L+layer+1])
        P =  matrix(p_1,n,n)
        P[1:floor(n/4), 1:floor(n/4)] = p_2
        P[(1+floor(n/4)):(2*floor(n/4)),(1+floor(n/4)):(2*floor(n/4)) ] = p_2
        P[(1+2*floor(n/4)):(3*floor(n/4)),(1+2*floor(n/4)):(3*floor(n/4)) ] = p_2
        P[(1+3*floor(n/4)):n,(1+3*floor(n/4)):n ] = p_2
        P_true[, , layer] = P
      }
      
      for (layer in 1: L){
        Al = matrix(rbinom(matrix(1,n,n), matrix(1,n,n), P_true[ , , layer]),n,n)
        Al[upper.tri(Al)] = t(Al)[upper.tri(Al)]
        A[, , layer] = Al
        
        A_list[[layer]] = Al
      }
      
      Y.tensor =  as.tensor(A)

      
      P_hat_1 = estimate_hooi(Y.tensor, hat.rank)
      result_hooi[ind_n, ind_L,  iter] = frobenius(P_hat_1, P_true)
      
      P_hat_2 = estimate_thpca(Y.tensor, hat.rank)
      result_hetero[ind_n, ind_L,  iter] = frobenius(P_hat_2, P_true)
      
      P_hat_3 = estimate_hosvd(Y.tensor, hat.rank)
      result_hosvd[ind_n, ind_L,  iter] = frobenius(P_hat_3, P_true)
      
      P_hat_4 = estimate_svd(A, hat.rank)
      result_svd[ind_n, ind_L,  iter] = frobenius(P_hat_4, P_true)
      
      P_hat_5 = uase(Y.tensor, hat.rank[1] ) 
      result_uase[ind_n, ind_L,  iter] = frobenius(P_hat_5, P_true)
      
      P_hat_6 = multiness_adaptive(A, hat.rank[1])
      result_MultiNeSS_ada[ind_n, ind_L,  iter] = frobenius(P_hat_6, P_true)
      
      P_hat_7 = estimate_flex(A, hat.rank[1], TT = 500, eta_u = 1e-4, eta_a = 1e-4, eta_l = 1e-4)
      result_flex[ind_n, ind_L,  iter] = frobenius(P_hat_7, P_true)
      
      P_hat_8 = estimate_mase(A, A_list, hat.rank[1])
      result_mase[ind_n, ind_L,  iter] = frobenius(P_hat_8, P_true)
      
      P_hat_9 = estimate_twist(Y.tensor, hat.rank, L)
      result_twist[ind_n, ind_L,  iter] = frobenius(P_hat_9, P_true)
      
      }
         
    }
  
}




value= c( apply(result_hetero, c(1,2), mean),
          apply(result_hosvd, c(1,2), mean),
          apply(result_hooi, c(1,2), mean),
          apply(result_twist, c(1,2), mean),
          apply(result_svd, c(1,2), mean),
          apply(result_uase, c(1,2), mean), 
          apply(result_mase, c(1,2), mean),
          apply(result_flex, c(1,2), mean),
          apply(result_MultiNeSS_ada, c(1,2), mean)
)
# value


sd = c( apply(result_hetero, c(1,2), sd),
        apply(result_hosvd, c(1,2), sd),
        apply(result_hooi, c(1,2), sd),
        apply(result_twist, c(1,2), sd),
        apply(result_svd, c(1,2), sd),
        apply(result_uase, c(1,2), sd),
        apply(result_mase, c(1,2), sd),
        apply(result_flex, c(1,2), sd),
        apply(result_MultiNeSS_ada, c(1,2), sd)
        )
# sd

layer = rep(c(10,20,30,40), times = 9)

Method = rep(c('TH-PCA', 'HOSVD', 'HOOI','TWIST', 'SASE', 'UASE','MASE', 'MLE','COA' ),each = 4) 
df_1 = data.frame( layer = layer, Method = Method, value = value, sd = sd)

df_1$Method <- factor(df_1$Method, levels = c('TH-PCA', 'HOSVD', 'HOOI','TWIST', 'SASE', 'UASE','MASE', 'MLE','COA'))

# Define the colors, making sure to name them in the same order as the factor levels
color_list <- setNames(c("#66c2a5", "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f", "#e5c494", "#737373", "#484B8A", "#C35817"),
                       levels(df_1$Method))

plot_1 = ggplot(data = df_1, mapping = aes(x = layer, y = value, colour = Method)) +
    geom_line() + 
    geom_point()+ 
   geom_errorbar(aes(ymin=value-sd, ymax= value+sd), width= 2)+
  labs( y ="Estimation error", x = "Number of layers")+
  scale_color_manual(values = color_list)+
  theme_classic()
                                                                                                                  
# plot_1

#### fix L, vary n ####



n_grid  =  c(50, 100, 150, 200)

L_grid = c(20)

NMC = 100

hat.rank =c(10,10,10)

result_hooi = array(NA, c( length(n_grid), length(L_grid), NMC))
result_hetero = array(NA, c( length(n_grid), length(L_grid), NMC))
result_hosvd = array(NA, c( length(n_grid), length(L_grid), NMC))
result_svd =   array(NA, c( length(n_grid), length(L_grid), NMC))
result_uase =   array(NA, c( length(n_grid), length(L_grid), NMC))
result_MultiNeSS_ada =   array(NA, c( length(n_grid), length(L_grid), NMC))
result_flex =   array(NA, c( length(n_grid), length(L_grid), NMC))
result_mase =   array(NA, c( length(n_grid), length(L_grid), NMC))
result_twist =   array(NA, c( length(n_grid), length(L_grid), NMC))

for(ind_n in 1:length(n_grid))
{
  n = n_grid[ind_n]
  print(paste0("n = ", n))
  
  for(ind_L in 1: length(L_grid) )
  {
    L = L_grid[ind_L]
    print(paste0("L = ", L))
    
    dim = c(n, n, L)
    
    prob = seq(0,1,  1/(4*L))  
    
    
    for(iter  in 1:NMC){
      A = array(NA,dim)
      P_true = array(NA,dim)
      A_list = list()
      
      set.seed(n*L*iter)
      for (layer in 1: L){
        p_1 =  runif(1, prob[2*L +layer], prob[2*L+layer+1])
        p_2 = runif(1, prob[3*L+layer], prob[3*L+layer+1])
        P =  matrix(p_1,n,n)
        P[1:floor(n/4), 1:floor(n/4)] = p_2
        P[(1+floor(n/4)):(2*floor(n/4)),(1+floor(n/4)):(2*floor(n/4)) ] = p_2
        P[(1+2*floor(n/4)):(3*floor(n/4)),(1+2*floor(n/4)):(3*floor(n/4)) ] = p_2
        P[(1+3*floor(n/4)):n,(1+3*floor(n/4)):n ] = p_2
        P_true[, , layer] = P
      }
      
      for (layer in 1: L){
        Al = matrix(rbinom(matrix(1,n,n), matrix(1,n,n), P_true[ , , layer]),n,n)
        Al[upper.tri(Al)] = t(Al)[upper.tri(Al)]
        A[, , layer] = Al
        
        A_list[[layer]] = Al
      }
      
      Y.tensor =  as.tensor(A)
      
      
      P_hat_1 = estimate_hooi(Y.tensor, hat.rank)
      result_hooi[ind_n, ind_L,  iter] = frobenius(P_hat_1, P_true)
      
      P_hat_2 = estimate_thpca(Y.tensor, hat.rank)
      result_hetero[ind_n, ind_L,  iter] = frobenius(P_hat_2, P_true)
      
      P_hat_3 = estimate_hosvd(Y.tensor, hat.rank)
      result_hosvd[ind_n, ind_L,  iter] = frobenius(P_hat_3, P_true)
      
      P_hat_4 = estimate_svd(A, hat.rank)
      result_svd[ind_n, ind_L,  iter] = frobenius(P_hat_4, P_true)
      
      P_hat_5 = uase(Y.tensor, hat.rank[1] ) 
      result_uase[ind_n, ind_L,  iter] = frobenius(P_hat_5, P_true)
      
      P_hat_6 = multiness_adaptive(A, hat.rank[1])
      result_MultiNeSS_ada[ind_n, ind_L,  iter] = frobenius(P_hat_6, P_true)
      
      P_hat_7 = estimate_flex(A, hat.rank[1], TT = 500, eta_u = 1e-4, eta_a = 1e-4, eta_l = 1e-4)
      result_flex[ind_n, ind_L,  iter] = frobenius(P_hat_7, P_true)
      
      P_hat_8 = estimate_mase(A, A_list, hat.rank[1])
      result_mase[ind_n, ind_L,  iter] = frobenius(P_hat_8, P_true)
      
      P_hat_9 = estimate_twist(Y.tensor, hat.rank, L)
      result_twist[ind_n, ind_L,  iter] = frobenius(P_hat_9, P_true)
      
    }
    
  }
  
}


value= c( apply(result_hetero, c(1,2), mean),
          apply(result_hosvd, c(1,2), mean),
          apply(result_hooi, c(1,2), mean),
          apply(result_twist, c(1,2), mean),
          apply(result_svd, c(1,2), mean),
          apply(result_uase, c(1,2), mean), 
          apply(result_mase, c(1,2), mean),
          apply(result_flex, c(1,2), mean),
          apply(result_MultiNeSS_ada, c(1,2), mean)
)
# value


sd = c( apply(result_hetero, c(1,2), sd),
        apply(result_hosvd, c(1,2), sd),
        apply(result_hooi, c(1,2), sd),
        apply(result_twist, c(1,2), sd),
        apply(result_svd, c(1,2), sd),
        apply(result_uase, c(1,2), sd),
        apply(result_mase, c(1,2), sd),
        apply(result_flex, c(1,2), sd),
        apply(result_MultiNeSS_ada, c(1,2), sd)
)
# sd

nodes = rep(c(50, 100, 150, 200), times = 9)

Method = rep(c('TH-PCA', 'HOSVD', 'HOOI','TWIST', 'SASE', 'UASE','MASE', 'MLE','COA' ),each = 4) 

df_2 = data.frame( nodes = nodes, Method = Method, value = value, sd = sd)

df_2$Method <- factor(df_1$Method, levels = c('TH-PCA', 'HOSVD', 'HOOI','TWIST', 'SASE', 'UASE','MASE', 'MLE','COA'))

# Define the colors, making sure to name them in the same order as the factor levels
color_list <- setNames(c("#66c2a5", "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f", "#e5c494", "#737373", "#484B8A", "#C35817"),
                       levels(df_1$Method))



plot_2 = ggplot(data = df_2, mapping = aes(x = nodes, y = value, colour = Method)) +
  geom_line() + 
  geom_point()+ 
  geom_errorbar(aes(ymin=value-sd, ymax= value+sd), width= 10)+
  labs( y ="Estimation error", x = "Number of nodes")+
  scale_color_manual(values = color_list)+
  theme_classic()




#### final plot ####

ggarrange(plot_1, plot_2, ncol = 2)

dev.copy2pdf(file=paste("/Users/chris/Desktop/MRDPG-main/estimation/scenario_1_sbm.pdf",sep = ''),width = 8,height = 2.75)
