setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(mvtnorm)
library(dirmult)
library(STATSVD)
library(changepoints)
library(egg)
library(ggplot2)
library(multiness)
source("tensor_functions.R")



#### fix n_x, n_y, d=4, vary L, hat.rank = c(10,10,10) ####

n_grid_X  = c(100)
n_grid_Y  = c(100)
L_grid = c(10, 20, 30, 40)
d_grid = c(4)

NMC = 100

hat.rank =c(10,10,10)

result_hooi = array(NA, c( length(n_grid_X), length(n_grid_Y), length(L_grid), length(d_grid), NMC))
result_hetero = array(NA, c( length(n_grid_X),length(n_grid_Y),  length(L_grid), length(d_grid),  NMC))
result_hosvd = array(NA, c( length(n_grid_X), length(n_grid_Y), length(L_grid), length(d_grid), NMC))
result_svd =   array(NA, c( length(n_grid_X), length(n_grid_Y),  length(L_grid), length(d_grid),  NMC))
result_uase =   array(NA, c( length(n_grid_X), length(n_grid_Y), length(L_grid),length(d_grid), NMC))
result_MultiNeSS_ada =   array(NA, c( length(n_grid_X), length(n_grid_Y), length(L_grid),length(d_grid), NMC))
result_flex =   array(NA, c( length(n_grid_X), length(n_grid_Y), length(L_grid),length(d_grid), NMC))

for(ind_n_1 in 1:length(n_grid_X)){
  n_1  = n_grid_X[ind_n_1]
  print(paste0("n1 = ", n_1))
  
  for(ind_n_2 in 1:length(n_grid_Y)){
    n_2  = n_grid_Y[ind_n_2]
    print(paste0("n2 = ", n_2))
    
    for(ind_L in 1: length(L_grid) ){
      L = L_grid[ind_L]
      print(paste0("L = ", L))
      
      for(ind_d in 1: length(d_grid) ){
        d = d_grid[ind_d]
        print(paste0("d = ", d))
        
        dim = c(n_1, n_2, L)
        Sigma_tensor =  array(NA,c(d, d, L))
        prob = seq(0,1,  1/(2*L))
        
        set.seed(n_1*n_2*L*d)
        
        for (layer in 1: L){
          Sigma = matrix(runif(d^2,prob[L+layer], prob[L+layer+1]), ncol=d) 
          Sigma_tensor[, , layer] = Sigma
        }
        
        for(iter  in 1:NMC){ 
          A = array(NA,dim)
          P_true = array(NA,dim)
          A_list = list()
          
          for (layer in 1: L){
            temp_1 = rdirichlet(n_1, rep(1, d))
            temp_2 = rdirichlet(n_2, rep(10, d))
            P =  temp_1  %*% Sigma_tensor[, , layer] %*% t(temp_2)
            P_true[, , layer] = P
            
            A[, , layer] = matrix(rbinom(matrix(1,n_1,n_2),matrix(1,n_1,n_2), P_true[ , , layer]),n_1,n_2)
            A_list[[layer]] = matrix(rbinom(matrix(1,n_1,n_2),matrix(1,n_1,n_2), P_true[ , , layer]),n_1,n_2)
          }
          
          Y.tensor =  as.tensor(A)
          
          P_hat_1 = estimate_hooi(Y.tensor, hat.rank)
          result_hooi[ind_n_1, ind_n_2, ind_L, ind_d, iter] = frobenius(P_hat_1, P_true)
          
          P_hat_2 = estimate_thpca(Y.tensor, hat.rank)
          result_hetero[ind_n_1, ind_n_2, ind_L, ind_d, iter] = frobenius(P_hat_2, P_true)
          
          P_hat_3 = estimate_hosvd(Y.tensor, hat.rank)
          result_hosvd[ind_n_1, ind_n_2, ind_L, ind_d, iter] = frobenius(P_hat_3, P_true)
          
          P_hat_4 = estimate_svd(A, hat.rank)
          result_svd[ind_n_1, ind_n_2, ind_L, ind_d, iter] = frobenius(P_hat_4, P_true)
          
          P_hat_5 = uase(Y.tensor, hat.rank[1] ) 
          result_uase[ind_n_1, ind_n_2, ind_L, ind_d, iter] = frobenius(P_hat_5, P_true)
          
          P_hat_6 = multiness_adaptive(A, hat.rank[1])
          result_MultiNeSS_ada[ind_n_1, ind_n_2, ind_L, ind_d, iter] = frobenius(P_hat_6, P_true)
          
          P_hat_7 = estimate_flex(A, hat.rank[1])
          result_flex[ind_n_1, ind_n_2, ind_L, ind_d, iter] = frobenius(P_hat_7, P_true)
          
        }
        
      }
    }
  }
}


value= c( apply(result_hooi, c(1,2,3,4), mean)[,,,],
          apply(result_hetero, c(1,2,3,4), mean)[,,,],
          apply(result_hosvd, c(1,2,3,4), mean)[,,,],
          apply(result_svd, c(1,2,3,4), mean)[,,,],
          apply(result_uase, c(1,2,3,4), mean)[,,,],
          apply(result_MultiNeSS_ada, c(1,2,3, 4), mean)[,,,],
          apply(result_flex, c(1,2,3, 4), mean)[,,,]
)

sd = c( apply(result_hooi, c(1,2,3,4), sd)[,,,],
        apply(result_hetero, c(1,2,3,4), sd)[,,,],
        apply(result_hosvd, c(1,2,3,4), sd)[,,,],
        apply(result_svd, c(1,2,3,4), sd)[,,,],
        apply(result_uase, c(1,2,3,4), sd)[,,,],
        apply(result_MultiNeSS_ada, c(1,2,3, 4), sd)[,,,],
        apply(result_flex, c(1,2,3, 4), sd)[,,,]
)

layer = rep(c(10,20,30,40), times = 7)

Method = rep(c('HOOI','TH-PCA','HOSVD','SASE','UASE','COA','MLE'),each = 4) 

dff_1 = data.frame( layer = layer, Method = Method, value = value, sd = sd)


plott_1= ggplot(data = dff_1, mapping = aes(x = layer, y = value, colour = Method)) +
  geom_line() + 
  geom_point()+ 
  geom_errorbar(aes(ymin=value-sd, ymax= value+sd), width= 2)+
  labs( y ="Estimation error", x = "Number of layers")+
  scale_color_manual(values = c("#66c2a5", "#8da0cb", "#e78ac3","#a6d854", "#ffd92f","#e5c494","#737373"))+
  theme_classic()









#### fix n_x, n_y, d, vary L, hat.rank = c(10,10,10) ####
n_grid_X  = c(100)
n_grid_Y  = c(100)
L_grid = c(10, 20, 30, 40)
d_grid = c(4)

# NMC = 100

hat.rank =c(10,10,10)

result_hooi = array(NA, c( length(n_grid_X), length(n_grid_Y), length(L_grid), length(d_grid), NMC))
result_hetero = array(NA, c( length(n_grid_X),length(n_grid_Y),  length(L_grid), length(d_grid),  NMC))
result_hosvd = array(NA, c( length(n_grid_X), length(n_grid_Y), length(L_grid), length(d_grid), NMC))
result_svd =   array(NA, c( length(n_grid_X), length(n_grid_Y),  length(L_grid), length(d_grid),  NMC))
result_uase =   array(NA, c( length(n_grid_X), length(n_grid_Y), length(L_grid),length(d_grid), NMC))
result_MultiNeSS_ada =   array(NA, c( length(n_grid_X), length(n_grid_Y), length(L_grid),length(d_grid), NMC))
result_flex =   array(NA, c( length(n_grid_X), length(n_grid_Y), length(L_grid),length(d_grid), NMC))

for(ind_n_1 in 1:length(n_grid_X)){
  n_1  = n_grid_X[ind_n_1]
  print(paste0("n1 = ", n_1))
  
  for(ind_n_2 in 1:length(n_grid_Y)){
    n_2  = n_grid_Y[ind_n_2]
    print(paste0("n2 = ", n_2))
    
    for(ind_L in 1: length(L_grid) ){
      L = L_grid[ind_L]
      print(paste0("L = ", L))
      
      for(ind_d in 1: length(d_grid) ){
        d = d_grid[ind_d]
        print(paste0("d = ", d))
        
        dim = c(n_1, n_2, L)
        Sigma_tensor =  array(NA,c(d, d, L))
        prob = seq(0,1,  1/(2*L))
        
        set.seed(n_1*n_2*L*d)
        
        for (layer in 1: L){
          Sigma = matrix(runif(d^2,0, 1), ncol=d) 
          Sigma_tensor[, , layer] = Sigma
        }
        
        for(iter  in 1:NMC){ 
          A = array(NA,dim)
          P_true = array(NA,dim)
          
          for (layer in 1: L){
            temp_1 = rdirichlet(n_1, rep(1, d))
            temp_2 = rdirichlet(n_2, rep(10, d))
            P =  temp_1  %*% Sigma_tensor[, , layer] %*% t(temp_2)
            P_true[, , layer] = P
            
            A[, , layer] = matrix(rbinom(matrix(1,n_1,n_2),matrix(1,n_1,n_2), P_true[ , , layer]),n_1,n_2)
          }
          
          Y.tensor =  as.tensor(A)
          
          P_hat_1 = estimate_hooi(Y.tensor, hat.rank)
          result_hooi[ind_n_1, ind_n_2, ind_L, ind_d, iter] = frobenius(P_hat_1, P_true)
          
          P_hat_2 = estimate_thpca(Y.tensor, hat.rank)
          result_hetero[ind_n_1, ind_n_2, ind_L, ind_d, iter] = frobenius(P_hat_2, P_true)
          
          P_hat_3 = estimate_hosvd(Y.tensor, hat.rank)
          result_hosvd[ind_n_1, ind_n_2, ind_L, ind_d, iter] = frobenius(P_hat_3, P_true)
          
          P_hat_4 = estimate_svd(A, hat.rank)
          result_svd[ind_n_1, ind_n_2, ind_L, ind_d, iter] = frobenius(P_hat_4, P_true)
          
          P_hat_5 = uase(Y.tensor, hat.rank[1] ) 
          result_uase[ind_n_1, ind_n_2, ind_L, ind_d, iter] = frobenius(P_hat_5, P_true)
          
          P_hat_7 = multiness_adaptive(A, hat.rank[1])
          result_MultiNeSS_ada[ind_n_1, ind_n_2, ind_L, ind_d, iter] = frobenius(P_hat_7, P_true)
          
          P_hat_8 = estimate_flex(A, hat.rank[1])
          result_flex[ind_n_1, ind_n_2, ind_L, ind_d, iter] = frobenius(P_hat_8, P_true)
          
        }
        
      }
    }
  }
}


value= c( apply(result_hooi, c(1,2,3,4), mean)[,,,],
          apply(result_hetero, c(1,2,3,4), mean)[,,,],
          apply(result_hosvd, c(1,2,3,4), mean)[,,,],
          apply(result_svd, c(1,2,3,4), mean)[,,,],
          apply(result_uase, c(1,2,3,4), mean)[,,,],
          apply(result_MultiNeSS_ada, c(1,2,3, 4), mean)[,,,],
          apply(result_flex, c(1,2,3, 4), mean)[,,,]
)

sd = c( apply(result_hooi, c(1,2,3,4), sd)[,,,],
        apply(result_hetero, c(1,2,3,4), sd)[,,,],
        apply(result_hosvd, c(1,2,3,4), sd)[,,,],
        apply(result_svd, c(1,2,3,4), sd)[,,,],
        apply(result_uase, c(1,2,3,4), sd)[,,,],
        apply(result_MultiNeSS_ada, c(1,2,3, 4), sd)[,,,],
        apply(result_flex, c(1,2,3, 4), sd)[,,,]
)

layer = rep(c(10,20,30,40), times = 7)

Method = rep(c('HOOI','TH-PCA','HOSVD','SASE','UASE', 'COA','MLE'), each = 4) 

dff_1_1 = data.frame( layer = layer, Method = Method, value = value, sd = sd)

plott_1_1 = ggplot(data = dff_1_1, mapping = aes(x = layer, y = value, colour = Method)) +
  geom_line() + 
  geom_point()+ 
  geom_errorbar(aes(ymin=value-sd, ymax= value+sd), width= 2)+
  labs( y ="Estimation error", x = "Number of layers")+
  scale_color_manual(values = c("#66c2a5", "#8da0cb", "#e78ac3","#a6d854", "#ffd92f","#e5c494","#737373"))+
  theme_classic()



#### fix n_x, n_y, L, vary d, hat.rank = c(10,10,10) ####


n_grid_X  = c(100)
n_grid_Y  = c(100)
L_grid = c(20)
d_grid = c(2, 4, 6, 8)

# NMC = 100

hat.rank =c(10,10,10)

result_hooi = array(NA, c( length(n_grid_X), length(n_grid_Y), length(L_grid), length(d_grid), NMC))
result_hetero = array(NA, c( length(n_grid_X),length(n_grid_Y),  length(L_grid), length(d_grid),  NMC))
result_hosvd = array(NA, c( length(n_grid_X), length(n_grid_Y), length(L_grid), length(d_grid), NMC))
result_svd =   array(NA, c( length(n_grid_X), length(n_grid_Y),  length(L_grid), length(d_grid),  NMC))
result_uase =   array(NA, c( length(n_grid_X), length(n_grid_Y), length(L_grid),length(d_grid), NMC))
result_MultiNeSS_ada =   array(NA, c( length(n_grid_X), length(n_grid_Y), length(L_grid),length(d_grid), NMC))
result_flex =   array(NA, c( length(n_grid_X), length(n_grid_Y), length(L_grid),length(d_grid), NMC))

for(ind_n_1 in 1:length(n_grid_X)){
  n_1  = n_grid_X[ind_n_1]
  print(paste0("n1 = ", n_1))
  
  for(ind_n_2 in 1:length(n_grid_Y)){
    n_2  = n_grid_Y[ind_n_2]
    print(paste0("n2 = ", n_2))
    
    for(ind_L in 1: length(L_grid) ){
      L = L_grid[ind_L]
      print(paste0("L = ", L))
      
      for(ind_d in 1: length(d_grid) ){
        d = d_grid[ind_d]
        print(paste0("d = ", d))
        
        dim = c(n_1, n_2, L)
        Sigma_tensor =  array(NA,c(d, d, L))
        prob = seq(0,1,  1/(2*L))
        
        set.seed(n_1*n_2*L*d)
        
        for (layer in 1: L){
          Sigma = matrix(runif(d^2,prob[L+layer], prob[L+layer+1]), ncol=d) 
          Sigma_tensor[, , layer] = Sigma
        }
        
        for(iter  in 1:NMC){ 
          A = array(NA,dim)
          P_true = array(NA,dim)
          
          for (layer in 1: L){
            temp_1 = rdirichlet(n_1, rep(1, d))
            temp_2 = rdirichlet(n_2, rep(10, d))
            P =  temp_1  %*% Sigma_tensor[, , layer] %*% t(temp_2)
            P_true[, , layer] = P
            
            A[, , layer] = matrix(rbinom(matrix(1,n_1,n_2),matrix(1,n_1,n_2), P_true[ , , layer]),n_1,n_2)
          }
          
          Y.tensor =  as.tensor(A)
          
          P_hat_1 = estimate_hooi(Y.tensor, hat.rank)
          result_hooi[ind_n_1, ind_n_2, ind_L, ind_d, iter] = frobenius(P_hat_1, P_true)
          
          P_hat_2 = estimate_thpca(Y.tensor, hat.rank)
          result_hetero[ind_n_1, ind_n_2, ind_L, ind_d, iter] = frobenius(P_hat_2, P_true)
          
          P_hat_3 = estimate_hosvd(Y.tensor, hat.rank)
          result_hosvd[ind_n_1, ind_n_2, ind_L, ind_d, iter] = frobenius(P_hat_3, P_true)
          
          P_hat_4 = estimate_svd(A, hat.rank)
          result_svd[ind_n_1, ind_n_2, ind_L, ind_d, iter] = frobenius(P_hat_4, P_true)
          
          P_hat_5 = uase(Y.tensor, hat.rank[1] ) 
          result_uase[ind_n_1, ind_n_2, ind_L, ind_d, iter] = frobenius(P_hat_5, P_true)
          
          P_hat_7 = multiness_adaptive(A, hat.rank[1])
          result_MultiNeSS_ada[ind_n_1, ind_n_2, ind_L, ind_d, iter] = frobenius(P_hat_7, P_true)
          
          P_hat_8 = estimate_flex(A, hat.rank[1])
          result_flex[ind_n_1, ind_n_2, ind_L, ind_d, iter] = frobenius(P_hat_8, P_true)
          
        }
        
      }
    }
  }
}


value= c( apply(result_hooi, c(1,2,3,4), mean)[,,,],
          apply(result_hetero, c(1,2,3,4), mean)[,,,],
          apply(result_hosvd, c(1,2,3,4), mean)[,,,],
          apply(result_svd, c(1,2,3,4), mean)[,,,],
          apply(result_uase, c(1,2,3,4), mean)[,,,],
          apply(result_MultiNeSS_ada, c(1,2,3, 4), mean)[,,,],
          apply(result_flex, c(1,2,3, 4), mean)[,,,]
)

sd = c( apply(result_hooi, c(1,2,3,4), sd)[,,,],
        apply(result_hetero, c(1,2,3,4), sd)[,,,],
        apply(result_hosvd, c(1,2,3,4), sd)[,,,],
        apply(result_svd, c(1,2,3,4), sd)[,,,],
        apply(result_uase, c(1,2,3,4), sd)[,,,],
        apply(result_MultiNeSS_ada, c(1,2,3, 4), sd)[,,,],
        apply(result_flex, c(1,2,3, 4), sd)[,,,]
)


dimension = rep( c(2, 4, 6, 8), times = 7)

Method = rep(c('HOOI','TH-PCA','HOSVD','SASE','UASE', 'COA','MLE'),each = 4)

dff_2 <- data.frame( dimension = dimension, Method = Method, value = value, sd = sd)

plott_2 = ggplot(data = dff_2, mapping = aes(x = dimension, y = value, colour = Method)) + geom_line() + geom_point() +
  geom_errorbar(aes(ymin=value-sd, ymax= value+sd), width= 0.4)+
  labs( y ="Estimation error", x = "Dimension of latent positions")+
  scale_color_manual(values = c("#66c2a5", "#8da0cb", "#e78ac3","#a6d854", "#ffd92f","#e5c494","#737373"))+
  theme_classic()



#### fix n_x, n_y, L, vary d, hat.rank = c(d,d,d) ####


n_grid_X  = c(100)

n_grid_Y  = c(100)

L_grid = c(20)

d_grid = c(2, 4, 6, 8)

# NMC = 100



result_hooi = array(NA, c( length(n_grid_X), length(n_grid_Y), length(L_grid), length(d_grid), NMC))
result_hetero = array(NA, c( length(n_grid_X),length(n_grid_Y),  length(L_grid), length(d_grid),  NMC))
result_hosvd = array(NA, c( length(n_grid_X), length(n_grid_Y), length(L_grid), length(d_grid), NMC))
result_svd =   array(NA, c( length(n_grid_X), length(n_grid_Y),  length(L_grid), length(d_grid),  NMC))
result_uase =   array(NA, c( length(n_grid_X), length(n_grid_Y), length(L_grid),length(d_grid), NMC))
result_MultiNeSS_ada =   array(NA, c( length(n_grid_X), length(n_grid_Y), length(L_grid),length(d_grid), NMC))
result_flex =   array(NA, c( length(n_grid_X), length(n_grid_Y), length(L_grid),length(d_grid), NMC))

for(ind_n_1 in 1:length(n_grid_X)){
  n_1  = n_grid_X[ind_n_1]
  print(paste0("n1 = ", n_1))
  
  for(ind_n_2 in 1:length(n_grid_Y)){
    n_2  = n_grid_Y[ind_n_2]
    print(paste0("n2 = ", n_2))
    
    for(ind_L in 1: length(L_grid) ){
      L = L_grid[ind_L]
      print(paste0("L = ", L))
      
      for(ind_d in 1: length(d_grid) ){
        d = d_grid[ind_d]
        print(paste0("d = ", d))
        
        hat.rank =c(d, d, d)
        
        dim = c(n_1, n_2, L)
        Sigma_tensor =  array(NA,c(d, d, L))
        prob = seq(0,1,  1/(2*L))
        
        set.seed(n_1*n_2*L*d)
        
        for (layer in 1: L){
          Sigma = matrix(runif(d^2,prob[L+layer], prob[L+layer+1]), ncol=d) 
          Sigma_tensor[, , layer] = Sigma
        }
        
        for(iter  in 1:NMC){ 
          A = array(NA,dim)
          P_true = array(NA,dim)
          
          for (layer in 1: L){
            temp_1 = rdirichlet(n_1, rep(1, d))
            temp_2 = rdirichlet(n_2, rep(10, d))
            P =  temp_1  %*% Sigma_tensor[, , layer] %*% t(temp_2)
            P_true[, , layer] = P
            
            A[, , layer] = matrix(rbinom(matrix(1,n_1,n_2),matrix(1,n_1,n_2), P_true[ , , layer]),n_1,n_2)
          }
          
          Y.tensor =  as.tensor(A)
          
          P_hat_1 = estimate_hooi(Y.tensor, hat.rank)
          result_hooi[ind_n_1, ind_n_2, ind_L, ind_d, iter] = frobenius(P_hat_1, P_true)
          
          P_hat_2 = estimate_thpca(Y.tensor, hat.rank)
          result_hetero[ind_n_1, ind_n_2, ind_L, ind_d, iter] = frobenius(P_hat_2, P_true)
          
          P_hat_3 = estimate_hosvd(Y.tensor, hat.rank)
          result_hosvd[ind_n_1, ind_n_2, ind_L, ind_d, iter] = frobenius(P_hat_3, P_true)
          
          P_hat_4 = estimate_svd(A, hat.rank)
          result_svd[ind_n_1, ind_n_2, ind_L, ind_d, iter] = frobenius(P_hat_4, P_true)
          
          P_hat_5 = uase(Y.tensor, hat.rank[1] ) 
          result_uase[ind_n_1, ind_n_2, ind_L, ind_d, iter] = frobenius(P_hat_5, P_true)
          
          P_hat_7 = multiness_adaptive(A, hat.rank[1])
          result_MultiNeSS_ada[ind_n_1, ind_n_2, ind_L, ind_d, iter] = frobenius(P_hat_7, P_true)
          
          P_hat_8 = estimate_flex(A, hat.rank[1])
          result_flex[ind_n_1, ind_n_2, ind_L, ind_d, iter] = frobenius(P_hat_8, P_true)
          
        }
        
      }
    }
  }
}


value= c( apply(result_hooi, c(1,2,3,4), mean)[,,,],
          apply(result_hetero, c(1,2,3,4), mean)[,,,],
          apply(result_hosvd, c(1,2,3,4), mean)[,,,],
          apply(result_svd, c(1,2,3,4), mean)[,,,],
          apply(result_uase, c(1,2,3,4), mean)[,,,],
          apply(result_MultiNeSS_ada, c(1,2,3, 4), mean)[,,,],
          apply(result_flex, c(1,2,3, 4), mean)[,,,]
)

sd = c( apply(result_hooi, c(1,2,3,4), sd)[,,,],
        apply(result_hetero, c(1,2,3,4), sd)[,,,],
        apply(result_hosvd, c(1,2,3,4), sd)[,,,],
        apply(result_svd, c(1,2,3,4), sd)[,,,],
        apply(result_uase, c(1,2,3,4), sd)[,,,],
        apply(result_MultiNeSS_ada, c(1,2,3, 4), sd)[,,,],
        apply(result_flex, c(1,2,3, 4), sd)[,,,]
)

dimension = rep( c(2, 4, 6, 8), times = 7)

Method = rep(c('HOOI','TH-PCA','HOSVD','SASE','UASE','COA','MLE'),each = 4)

dff_2_1 <- data.frame( dimension = dimension, Method = Method, value = value, sd = sd)

plott_2_1 = ggplot(data = dff_2_1, mapping = aes(x = dimension, y = value, colour = Method)) + geom_line() + geom_point() +
  geom_errorbar(aes(ymin=value-sd, ymax= value+sd), width= 0.4)+
  labs( y ="Estimation error", x = "Dimension of latent positions")+
  scale_color_manual(values = c("#66c2a5", "#8da0cb", "#e78ac3","#a6d854", "#ffd92f","#e5c494","#737373"))+
  theme_classic()





#### fix L, d=4, vary n_x=n_y, hat.rank = c(10, 10, 10) ####


n_grid_X  = c(50, 100, 150, 200)
n_grid_Y  = c(0)

L_grid = c(20)

d_grid = c(4)

# NMC = 100

hat.rank =c(10,10,10)

result_hooi = array(NA, c( length(n_grid_X), length(n_grid_Y), length(L_grid), length(d_grid), NMC))
result_hetero = array(NA, c( length(n_grid_X),length(n_grid_Y),  length(L_grid), length(d_grid),  NMC))
result_hosvd = array(NA, c( length(n_grid_X), length(n_grid_Y), length(L_grid), length(d_grid), NMC))
result_svd =   array(NA, c( length(n_grid_X), length(n_grid_Y),  length(L_grid), length(d_grid),  NMC))
result_uase =   array(NA, c( length(n_grid_X), length(n_grid_Y), length(L_grid),length(d_grid), NMC))
result_MultiNeSS_ada =   array(NA, c( length(n_grid_X), length(n_grid_Y), length(L_grid),length(d_grid), NMC))
result_flex =   array(NA, c( length(n_grid_X), length(n_grid_Y), length(L_grid),length(d_grid), NMC))

for(ind_n_1 in 1:length(n_grid_X)){
  n_1  = n_grid_X[ind_n_1]
  print(paste0("n1 = ", n_1))
  
  for(ind_n_2 in 1:length(n_grid_Y)){
    n_2  = n_1
    print(paste0("n2 = ", n_2))
    
    for(ind_L in 1: length(L_grid) ){
      L = L_grid[ind_L]
      print(paste0("L = ", L))
      
      for(ind_d in 1: length(d_grid) ){
        d = d_grid[ind_d]
        print(paste0("d = ", d))
        
        dim = c(n_1, n_2, L)
        Sigma_tensor =  array(NA,c(d, d, L))
        prob = seq(0,1,  1/(2*L))
        
        set.seed(n_1*n_2*L*d)
        
        for (layer in 1: L){
          Sigma = matrix(runif(d^2,prob[L+layer], prob[L+layer+1]), ncol=d) 
          Sigma_tensor[, , layer] = Sigma
        }
        
        for(iter  in 1:NMC){ 
          A = array(NA,dim)
          P_true = array(NA,dim)
          
          for (layer in 1: L){
            temp_1 = rdirichlet(n_1, rep(1, d))
            temp_2 = rdirichlet(n_2, rep(10, d))
            P =  temp_1  %*% Sigma_tensor[, , layer] %*% t(temp_2)
            P_true[, , layer] = P
            
            A[, , layer] = matrix(rbinom(matrix(1,n_1,n_2),matrix(1,n_1,n_2), P_true[ , , layer]),n_1,n_2)
          }
          
          Y.tensor =  as.tensor(A)
          
          P_hat_1 = estimate_hooi(Y.tensor, hat.rank)
          result_hooi[ind_n_1, ind_n_2, ind_L, ind_d, iter] = frobenius(P_hat_1, P_true)
          
          P_hat_2 = estimate_thpca(Y.tensor, hat.rank)
          result_hetero[ind_n_1, ind_n_2, ind_L, ind_d, iter] = frobenius(P_hat_2, P_true)
          
          P_hat_3 = estimate_hosvd(Y.tensor, hat.rank)
          result_hosvd[ind_n_1, ind_n_2, ind_L, ind_d, iter] = frobenius(P_hat_3, P_true)
          
          P_hat_4 = estimate_svd(A, hat.rank)
          result_svd[ind_n_1, ind_n_2, ind_L, ind_d, iter] = frobenius(P_hat_4, P_true)
          
          P_hat_5 = uase(Y.tensor, hat.rank[1] ) 
          result_uase[ind_n_1, ind_n_2, ind_L, ind_d, iter] = frobenius(P_hat_5, P_true)
          
          P_hat_7 = multiness_adaptive(A, hat.rank[1])
          result_MultiNeSS_ada[ind_n_1, ind_n_2, ind_L, ind_d, iter] = frobenius(P_hat_7, P_true)
          
          P_hat_8 = estimate_flex(A, hat.rank[1], TT = 500, eta_u = 1e-4, eta_a = 1e-4, eta_l = 1e-4)
          result_flex[ind_n_1, ind_n_2, ind_L, ind_d, iter] = frobenius(P_hat_8, P_true)
          
        }
        
      }
    }
  }
}


value= c( apply(result_hooi, c(1,2,3,4), mean)[,,,],
          apply(result_hetero, c(1,2,3,4), mean)[,,,],
          apply(result_hosvd, c(1,2,3,4), mean)[,,,],
          apply(result_svd, c(1,2,3,4), mean)[,,,],
          apply(result_uase, c(1,2,3,4), mean)[,,,],
          apply(result_MultiNeSS_ada, c(1,2,3, 4), mean)[,,,],
          apply(result_flex, c(1,2,3, 4), mean)[,,,]
)

sd = c( apply(result_hooi, c(1,2,3,4), sd)[,,,],
        apply(result_hetero, c(1,2,3,4), sd)[,,,],
        apply(result_hosvd, c(1,2,3,4), sd)[,,,],
        apply(result_svd, c(1,2,3,4), sd)[,,,],
        apply(result_uase, c(1,2,3,4), sd)[,,,],
        apply(result_MultiNeSS_ada, c(1,2,3, 4), sd)[,,,],
        apply(result_flex, c(1,2,3, 4), sd)[,,,]
)

node = rep( c(50, 100, 150, 200), times = 7)

Method = rep(c('HOOI','TH-PCA','HOSVD','SASE','UASE', 'COA', 'MLE'), each = 4)

dff_3 <- data.frame( node = node, Method = Method, value = value, sd = sd)

plott_3 = ggplot(data = dff_3, mapping = aes(x = node, y = value, colour = Method)) + geom_line() + geom_point() +
  geom_errorbar(aes(ymin=value-sd, ymax= value+sd), width= 10)+
  labs( y ="Estimation error",  x =expression(paste( "Number of nodes: n with n=", n[1],"=",n[2])))+
  scale_color_manual(values = c("#66c2a5", "#8da0cb", "#e78ac3","#a6d854", "#ffd92f","#e5c494","#737373"))+
  theme_classic()



#### fix L, d=4, n_x, vary n_y, hat.rank = c(10, 10, 10) ####

n_grid_X  = c(100)

n_grid_Y  = c(50, 100, 150, 200)

L_grid = c(20)

d_grid = c(4)

# NMC = 100

hat.rank =c(10,10,10)

result_hooi = array(NA, c( length(n_grid_X), length(n_grid_Y), length(L_grid), length(d_grid), NMC))
result_hetero = array(NA, c( length(n_grid_X),length(n_grid_Y),  length(L_grid), length(d_grid),  NMC))
result_hosvd = array(NA, c( length(n_grid_X), length(n_grid_Y), length(L_grid), length(d_grid), NMC))
result_svd =   array(NA, c( length(n_grid_X), length(n_grid_Y),  length(L_grid), length(d_grid),  NMC))
result_uase =   array(NA, c( length(n_grid_X), length(n_grid_Y), length(L_grid),length(d_grid), NMC))


for(ind_n_1 in 1:length(n_grid_X)){
  n_1  = n_grid_X[ind_n_1]
  print(paste0("n1 = ", n_1))
  
  for(ind_n_2 in 1:length(n_grid_Y)){
    n_2  = n_grid_Y[ind_n_2]
    print(paste0("n2 = ", n_2))
    
    for(ind_L in 1: length(L_grid) ){
      L = L_grid[ind_L]
      print(paste0("L = ", L))
      
      for(ind_d in 1: length(d_grid) ){
        d = d_grid[ind_d]
        print(paste0("d = ", d))
        
        dim = c(n_1, n_2, L)
        Sigma_tensor =  array(NA,c(d, d, L))
        prob = seq(0,1,  1/(2*L))
        
        set.seed(n_1*n_2*L*d)
        
        for (layer in 1: L){
          Sigma = matrix(runif(d^2,prob[L+layer], prob[L+layer+1]), ncol=d) 
          Sigma_tensor[, , layer] = Sigma
        }
        
        for(iter  in 1:NMC){ 
          A = array(NA,dim)
          
          P_true = array(NA,dim)
          
          for (layer in 1: L){
            temp_1 = rdirichlet(n_1, rep(1, d))
            temp_2 = rdirichlet(n_2, rep(10, d))
            P =  temp_1  %*% Sigma_tensor[, , layer] %*% t(temp_2)
            P_true[, , layer] = P
            
            A[, , layer] = matrix(rbinom(matrix(1,n_1,n_2),matrix(1,n_1,n_2), P_true[ , , layer]),n_1,n_2)
          }
          
          Y.tensor =  as.tensor(A)
          
          P_hat_1 = estimate_hooi(Y.tensor, hat.rank)
          result_hooi[ind_n_1, ind_n_2, ind_L, ind_d, iter] = frobenius(P_hat_1, P_true)
          
          P_hat_2 = estimate_thpca(Y.tensor, hat.rank)
          result_hetero[ind_n_1, ind_n_2, ind_L, ind_d, iter] = frobenius(P_hat_2, P_true)
          
          P_hat_3 = estimate_hosvd(Y.tensor, hat.rank)
          result_hosvd[ind_n_1, ind_n_2, ind_L, ind_d, iter] = frobenius(P_hat_3, P_true)
          
          P_hat_4 = estimate_svd(A, hat.rank)
          result_svd[ind_n_1, ind_n_2, ind_L, ind_d, iter] = frobenius(P_hat_4, P_true)
          
          P_hat_5 = uase(Y.tensor, hat.rank[1] ) 
          result_uase[ind_n_1, ind_n_2, ind_L, ind_d, iter] = frobenius(P_hat_5, P_true)
          
        }
        
      }
    }
  }
}



value= c( apply(result_hooi, c(1,2,3,4), mean)[,,,],
          apply(result_hetero, c(1,2,3,4), mean)[,,,],
          apply(result_hosvd, c(1,2,3,4), mean)[,,,],
          apply(result_svd, c(1,2,3,4), mean)[,,,],
          apply(result_uase, c(1,2,3,4), mean)[,,,])

sd = c( apply(result_hooi, c(1,2,3,4), sd)[,,,],
        apply(result_hetero, c(1,2,3,4), sd)[,,,],
        apply(result_hosvd, c(1,2,3,4), sd)[,,,],
        apply(result_svd, c(1,2,3,4), sd)[,,,],
        apply(result_uase, c(1,2,3,4), sd)[,,,])

node = rep( c(50, 100, 150, 200), times = 5)

Method = rep(c('HOOI','TH-PCA','HOSVD','SASE','UASE'),each = 4)

dff_4 <- data.frame( node = node, Method = Method, value = value, sd = sd)

plott_4 = ggplot(data = dff_4, mapping = aes(x = node, y = value, colour = Method)) + geom_line() + geom_point() +
  geom_errorbar(aes(ymin=value-sd, ymax= value+sd), width= 10)+
  labs( y ="Estimation error", x =expression(paste( "Number of tail nodes: ", n[2])))+
  scale_color_manual(breaks = c('HOOI','HOSVD','SASE','TH-PCA','UASE'),
                     values = c("#8da0cb", "#e78ac3","#ffd92f","#e5c494","#737373"))+
  theme_classic()




#### final plot ####

# ggarrange(plott_1, plott_1_1, plott_2, plott_2_1, plott_3, plott_4, ncol = 2)



ggarrange(plott_1, plott_1_1, ncol = 2)
dev.copy2pdf(file=paste("scenario_2_dirichlet_L.pdf",sep = ''),width = 8,height = 2.5)

ggarrange(plott_2, plott_2_1, ncol = 2)
dev.copy2pdf(file=paste("scenario_2_dirichlet_d.pdf",sep = ''),width = 8,height = 2.5)

ggarrange(plott_3, plott_4, ncol = 2)
dev.copy2pdf(file=paste("scenario_2_dirichlet_n.pdf",sep = ''),width = 8,height = 2.5)

