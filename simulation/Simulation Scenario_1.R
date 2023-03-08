library(mvtnorm)
library(STATSVD)
library(changepoints)
library(egg)
library(ggplot2)
library(multiness)
source("/Users/chris/MRDPG/tensor_functions.R")
source("/Users/chris/MRDPG/cpd_MRDPG_functions.R")

n_grid  =  c(100)

L_grid = c(10, 20, 30, 40)

NMC = 100

hat.rank =c(10,10,10)

result_hooi = array(NA, c( length(n_grid), length(L_grid), NMC))
result_hetero = array(NA, c( length(n_grid), length(L_grid), NMC))
result_hosvd = array(NA, c( length(n_grid), length(L_grid), NMC))
result_svd =   array(NA, c( length(n_grid), length(L_grid), NMC))
result_uase =   array(NA, c( length(n_grid), length(L_grid), NMC))
result_MultiNeSS_ada =   array(NA, c( length(n_grid), length(L_grid), NMC))
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
      probability = array(NA,dim)
      probability.4 = array(NA,dim)
      probability.7 = array(NA,dim)
      
      for (layer in 1: L)
      {
        p_1 =  runif(1, prob[2*L +layer], prob[2*L+layer+1])
        p_2 = runif(1, prob[3*L+layer], prob[3*L+layer+1])
        P =  matrix(p_1,n,n)
        P[1:floor(n/4), 1:floor(n/4)] = p_2
        P[(1+floor(n/4)):(2*floor(n/4)),(1+floor(n/4)):(2*floor(n/4)) ] = p_2
        P[(1+2*floor(n/4)):(3*floor(n/4)),(1+2*floor(n/4)):(3*floor(n/4)) ] = p_2
        P[(1+3*floor(n/4)):n,(1+3*floor(n/4)):n ] = p_2
        probability[, , layer] = P
      }
      
      for (layer in 1: L)
      {
        
        A[, , layer] = matrix(rbinom(matrix(1,n,n),matrix(1,n,n),probability[ , , layer]),n,n)
        
        temp = svd(A[, , layer])
        
        xhat =  temp$u[,1:hat.rank[1]] %*%  diag( sqrt(temp$d[1:hat.rank[1]]) )
        
        yhat =  temp$v[,1:hat.rank[1]] %*%  diag( sqrt(temp$d[1:hat.rank[1]]) )
        
        probability.4[, , layer] = xhat %*% t(yhat)
        
      }
      
      
      Y.tensor =  as.tensor(A)
      
      probability.5 =   uase(Y.tensor, hat.rank[1] ) 
      
      result_uase[ind_n, ind_L,  iter] = fnorm(as.tensor(probability.5) - as.tensor(probability))
      
      U.hat = HOOI(Y.tensor, hat.rank)
      
      P.U1 = U.hat[[1]]%*%t(U.hat[[1]])
      P.U2 = U.hat[[2]]%*%t(U.hat[[2]])
      P.U3 = U.hat[[3]]%*%t(U.hat[[3]])
      Y.hat.1 = ttm(ttm(ttm(Y.tensor, P.U1, 1), P.U2, 2), P.U3, 3) 
      
      probability.1  = Y.hat.1@data
      
      probability.1[probability.1 > 1]  = 1
      
      probability.1[probability.1 < 0]  = 0
      
      result_hooi[ind_n, ind_L,  iter] = fnorm(as.tensor(probability.1) - as.tensor(probability))
      
      
      U.hat = Tensor_Hetero_PCA_test(Y.tensor, hat.rank)
      P.U1 = U.hat[[1]]%*%t(U.hat[[1]])
      P.U2 = U.hat[[2]]%*%t(U.hat[[2]])
      P.U3 = U.hat[[3]]%*%t(U.hat[[3]])
      Y.hat.2 = ttm(ttm(ttm(Y.tensor, P.U1, 1), P.U2, 2), P.U3, 3) 
      
      probability.2  = Y.hat.2@data
      
      probability.2[probability.2 > 1]  = 1
      
      probability.2[probability.2 < 0]  = 0
      
      result_hetero[ind_n, ind_L,  iter] = fnorm(as.tensor(probability.2) - as.tensor(probability))
      
      
      U.hat = hosvd(Y.tensor, hat.rank)$U
      P.U1 = U.hat[[1]]%*%t(U.hat[[1]])
      P.U2 = U.hat[[2]]%*%t(U.hat[[2]])
      P.U3 = U.hat[[3]]%*%t(U.hat[[3]])
      Y.hat.3 = ttm(ttm(ttm(Y.tensor, P.U1, 1), P.U2, 2), P.U3, 3) 
      
      probability.3  = Y.hat.3@data
      
      probability.3[probability.3 > 1]  = 1
      
      probability.3[probability.3 < 0]  = 0
      
      result_hosvd[ind_n, ind_L,  iter] = fnorm(as.tensor(probability.3) - as.tensor(probability))
      
      probability.4[probability.4 > 1]  = 1
      
      probability.4[probability.4 < 0]  = 0
      
      result_svd[ind_n, ind_L,  iter] = fnorm(as.tensor(probability.4) - as.tensor(probability))
      
      probability.7 = multiness_adaptive(A, hat.rank[1])
      
      result_MultiNeSS_ada[ind_n, ind_L,  iter] = fnorm(as.tensor(probability.7) - as.tensor(probability))
      
       }
         
    }
  
}


apply(result_hooi, c(1,2), mean)
apply(result_hetero, c(1,2), mean)
apply(result_hosvd, c(1,2), mean)
apply(result_svd, c(1,2), mean)
apply(result_uase, c(1,2), mean)
apply(result_MultiNeSS_ada, c(1,2), mean)

apply(result_hooi, c(1,2), sd)
apply(result_hetero, c(1,2), sd)
apply(result_hosvd, c(1,2), sd)
apply(result_svd, c(1,2), sd)
apply(result_uase, c(1,2), sd)
apply(result_MultiNeSS_ada, c(1,2), sd)

value= c( apply(result_hooi, c(1,2), mean),
          apply(result_hetero, c(1,2), mean),
          apply(result_hosvd, c(1,2), mean),
          apply(result_svd, c(1,2), mean),
          apply(result_uase, c(1,2), mean), 
          apply(result_MultiNeSS_ada, c(1,2), mean))


sd = c( apply(result_hooi, c(1,2), sd),
        apply(result_hetero, c(1,2), sd),
        apply(result_hosvd, c(1,2), sd),
        apply(result_svd, c(1,2), sd),
        apply(result_uase, c(1,2), sd),
        apply(result_MultiNeSS_ada, c(1,2), sd))

layer = rep(c(10,20,30,40), times = 6)

Method = rep(c('HOOI','TH-PCA','HOSVD','SASE','UASE', 'COA'),each = 4) 

df_1 = data.frame( layer = layer, Method = Method, value = value, sd = sd)

plot_1 = ggplot(data = df_1, mapping = aes(x = layer, y = value, colour = Method)) +
    geom_line() + 
    geom_point()+ 
   geom_errorbar(aes(ymin=value-sd, ymax= value+sd), width= 2)+
  labs( y ="Estimation error", x = "Number of layers")+
  scale_color_manual(values = c("#66c2a5", "#8da0cb", "#e78ac3","#a6d854", "#ffd92f","#e5c494"))+
  theme_classic()
                                                                                                                  


n_grid  =  c(50, 100, 150, 200)

L_grid = c(20)

NMC = 100

hat.rank =c(10,10,10)

result_hooi = array(NA, c( length(n_grid), length(L_grid), NMC))
result_hetero = array(NA, c( length(n_grid), length(L_grid), NMC))
result_hosvd = array(NA, c( length(n_grid), length(L_grid), NMC))
result_svd =   array(NA, c( length(n_grid), length(L_grid), NMC))
result_uase =   array(NA, c( length(n_grid), length(L_grid), NMC))
result_uase =   array(NA, c( length(n_grid), length(L_grid), NMC))
result_MultiNeSS_ada =   array(NA, c( length(n_grid), length(L_grid), NMC))
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
    
    
    for(iter  in 1:NMC)
    {
    
      
      A = array(NA,dim)
      probability = array(NA,dim)
      probability.4 = array(NA,dim)
      probability.7 = array(NA,dim)
      
      set.seed(n*L*iter)
      for (layer in 1: L)
      {
        
        p_1 =  runif(1, prob[2*L +layer], prob[2*L+layer+1])
        p_2 = runif(1, prob[3*L+layer], prob[3*L+layer+1])
        P =  matrix(p_1,n,n)
        P[1:floor(n/4), 1:floor(n/4)] = p_2
        P[(1+floor(n/4)):(2*floor(n/4)),(1+floor(n/4)):(2*floor(n/4)) ] = p_2
        P[(1+2*floor(n/4)):(3*floor(n/4)),(1+2*floor(n/4)):(3*floor(n/4)) ] = p_2
        P[(1+3*floor(n/4)):n,(1+3*floor(n/4)):n ] = p_2
        probability[, , layer] = P
      }
      for (layer in 1: L)
      {
        set.seed(n*L*layer*iter*0.05) 
        
        A[, , layer] = matrix(rbinom(matrix(1,n,n),matrix(1,n,n),probability[ , , layer]),n,n)
        
        temp = svd(A[, , layer])
        
        xhat =  temp$u[,1:hat.rank[1]] %*%  diag( sqrt(temp$d[1:hat.rank[1]]) )
        
        yhat =  temp$v[,1:hat.rank[1]] %*%  diag( sqrt(temp$d[1:hat.rank[1]]) )
        
        probability.4[, , layer] = xhat %*% t(yhat)
        
      }
      
      
      Y.tensor =  as.tensor(A)
      
      probability.5 =   uase(Y.tensor, hat.rank[1] ) 
      
      result_uase[ind_n, ind_L,  iter] = fnorm(as.tensor(probability.5) - as.tensor(probability))
      
      U.hat = HOOI(Y.tensor, hat.rank)
      
      P.U1 = U.hat[[1]]%*%t(U.hat[[1]])
      P.U2 = U.hat[[2]]%*%t(U.hat[[2]])
      P.U3 = U.hat[[3]]%*%t(U.hat[[3]])
      Y.hat.1 = ttm(ttm(ttm(Y.tensor, P.U1, 1), P.U2, 2), P.U3, 3) 
      
      probability.1  = Y.hat.1@data
      
      probability.1[probability.1 > 1]  = 1
      
      probability.1[probability.1 < 0]  = 0
      
      result_hooi[ind_n, ind_L,  iter] = fnorm(as.tensor(probability.1) - as.tensor(probability))
      
      
      U.hat = Tensor_Hetero_PCA_test(Y.tensor, hat.rank)
      P.U1 = U.hat[[1]]%*%t(U.hat[[1]])
      P.U2 = U.hat[[2]]%*%t(U.hat[[2]])
      P.U3 = U.hat[[3]]%*%t(U.hat[[3]])
      Y.hat.2 = ttm(ttm(ttm(Y.tensor, P.U1, 1), P.U2, 2), P.U3, 3) 
      
      probability.2  = Y.hat.2@data
      
      probability.2[probability.2 > 1]  = 1
      
      probability.2[probability.2 < 0]  = 0
      
      result_hetero[ind_n, ind_L,  iter] = fnorm(as.tensor(probability.2) - as.tensor(probability))
      
      
      U.hat = hosvd(Y.tensor, hat.rank)$U
      P.U1 = U.hat[[1]]%*%t(U.hat[[1]])
      P.U2 = U.hat[[2]]%*%t(U.hat[[2]])
      P.U3 = U.hat[[3]]%*%t(U.hat[[3]])
      Y.hat.3 = ttm(ttm(ttm(Y.tensor, P.U1, 1), P.U2, 2), P.U3, 3) 
      
      probability.3  = Y.hat.3@data
      
      probability.3[probability.3 > 1]  = 1
      
      probability.3[probability.3 < 0]  = 0
      
      result_hosvd[ind_n, ind_L,  iter] = fnorm(as.tensor(probability.3) - as.tensor(probability))
      
      probability.4[probability.4 > 1]  = 1
      
      probability.4[probability.4 < 0]  = 0
      
      result_svd[ind_n, ind_L,  iter] = fnorm(as.tensor(probability.4) - as.tensor(probability))
      
      probability.7 = multiness_adaptive(A, hat.rank[1])
      
      result_MultiNeSS_ada[ind_n, ind_L,  iter] = fnorm(as.tensor(probability.7) - as.tensor(probability))
      
      
    }
    
  }
  
}


apply(result_hooi, c(1,2), mean)
apply(result_hetero, c(1,2), mean)
apply(result_hosvd, c(1,2), mean)
apply(result_svd, c(1,2), mean)
apply(result_uase, c(1,2), mean)
apply(result_MultiNeSS_ada, c(1,2), mean)

apply(result_hooi, c(1,2), sd)
apply(result_hetero, c(1,2), sd)
apply(result_hosvd, c(1,2), sd)
apply(result_svd, c(1,2), sd)
apply(result_uase, c(1,2), sd)
apply(result_MultiNeSS_ada, c(1,2), sd)

value= c( apply(result_hooi, c(1,2), mean),
          apply(result_hetero, c(1,2), mean),
          apply(result_hosvd, c(1,2), mean),
          apply(result_svd, c(1,2), mean),
          apply(result_uase, c(1,2), mean), 
          apply(result_MultiNeSS_ada, c(1,2), mean))


sd = c( apply(result_hooi, c(1,2), sd),
        apply(result_hetero, c(1,2), sd),
        apply(result_hosvd, c(1,2), sd),
        apply(result_svd, c(1,2), sd),
        apply(result_uase, c(1,2), sd),
        apply(result_MultiNeSS_ada, c(1,2), sd))

nodes = rep(c(50,100,150,200), times = 6)

Method = rep(c('HOOI','TH-PCA','HOSVD','SASE','UASE',  'COA'),each = 4) 

df_2 = data.frame( nodes = nodes, Method = Method, value = value, sd = sd)

plot_2 = ggplot(data = df_2, mapping = aes(x = nodes, y = value, colour = Method)) +
  geom_line() + 
  geom_point()+ 
  geom_errorbar(aes(ymin=value-sd, ymax= value+sd), width= 10)+
  labs( y ="Estimation error", x = "Number of nodes")+
  scale_color_manual(values = c("#66c2a5", "#8da0cb", "#e78ac3","#a6d854", "#ffd92f","#e5c494"))+
  theme_classic()



ggarrange(plot_1, plot_2, ncol = 2)

