library(mvtnorm)
library(dirmult)
library(STATSVD)
library(changepoints)
library(egg)
library(ggplot2)
source("/Users/chris/MRDPG/tensor_functions.R")
source("/Users/chris/MRDPG/cpd_MRDPG_functions.R")
library(multiness)

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


for(ind_n_1 in 1:length(n_grid_X))
{
  n_1  = n_grid_X[ind_n_1]
  
  print(paste0("n1 = ", n_1))
  
  for(ind_n_2 in 1:length(n_grid_Y))
  
  {
    n_2  = n_grid_Y[ind_n_2]
    
  print(paste0("n2 = ", n_2))
  
  for(ind_L in 1: length(L_grid) )
  {
    L = L_grid[ind_L]
    
    print(paste0("L = ", L))
    
    for(ind_d in 1: length(d_grid) )
    {
      d = d_grid[ind_d]
      
      print(paste0("d = ", d))
      
      dim = c(n_1, n_2, L)
      
      probability = array(NA,dim)
      
      Sigma_tensor =  array(NA,c(d, d, L))
      
      prob = seq(0,1,  1/(2*L))
      
      set.seed(n_1*n_2*L*d)
      
      for (layer in 1: L)
      {
        
        Sigma = matrix(runif(d^2,prob[L+layer], prob[L+layer+1]), ncol=d) 
        Sigma_tensor[, , layer] = Sigma
      }
      
      for(iter  in 1:NMC)
      { 
        A = array(NA,dim)
        
        probability.4 = array(NA,dim)
      
        
        for (layer in 1: L)
        {

          temp_1 = rdirichlet(n_1, rep(1, d))
          temp_2 = rdirichlet(n_2, rep(10, d))
          P =  temp_1  %*% Sigma_tensor[, , layer] %*% t(temp_2)
          probability[, , layer] = P
          
          A[, , layer] = matrix(rbinom(matrix(1,n_1,n_2),matrix(1,n_1,n_2), probability[ , , layer]),n_1,n_2)
          
          temp = svd(A[, , layer])
          
          xhat =  temp$u[,1:hat.rank[1]] %*%  diag( sqrt(temp$d[1:hat.rank[1]]) )
          
          yhat =  temp$v[,1:hat.rank[1]] %*%  diag( sqrt(temp$d[1:hat.rank[1]]) )
          
          probability.4[, , layer] = xhat %*% t(yhat)
          
        }
      
        Y.tensor =  as.tensor(A)

        probability.5 =   uase(Y.tensor, hat.rank[1] ) 
        
        result_uase[ind_n_1, ind_n_2, ind_L, ind_d, iter] = fnorm(as.tensor(probability.5) - as.tensor(probability))
        
        U.hat = HOOI(Y.tensor, hat.rank)
        
        P.U1 = U.hat[[1]]%*%t(U.hat[[1]])
        P.U2 = U.hat[[2]]%*%t(U.hat[[2]])
        P.U3 = U.hat[[3]]%*%t(U.hat[[3]])
        Y.hat.1 = ttm(ttm(ttm(Y.tensor, P.U1, 1), P.U2, 2), P.U3, 3) 
        
        probability.1  = Y.hat.1@data
        
        probability.1[probability.1 > 1]  = 1
        
        probability.1[probability.1 < 0]  = 0
        
        result_hooi[ind_n_1, ind_n_2, ind_L, ind_d, iter] = fnorm(as.tensor(probability.1) - as.tensor(probability))
        
        
        U.hat = Tensor_Hetero_PCA_test(Y.tensor, hat.rank)
        P.U1 = U.hat[[1]]%*%t(U.hat[[1]])
        P.U2 = U.hat[[2]]%*%t(U.hat[[2]])
        P.U3 = U.hat[[3]]%*%t(U.hat[[3]])
        Y.hat.2 = ttm(ttm(ttm(Y.tensor, P.U1, 1), P.U2, 2), P.U3, 3) 
        
        probability.2  = Y.hat.2@data
        
        probability.2[probability.2 > 1]  = 1
        
        probability.2[probability.2 < 0]  = 0
        
        result_hetero[ind_n_1, ind_n_2, ind_L, ind_d, iter] = fnorm(as.tensor(probability.2) - as.tensor(probability))
        
        
        U.hat = hosvd(Y.tensor, hat.rank)$U
        P.U1 = U.hat[[1]]%*%t(U.hat[[1]])
        P.U2 = U.hat[[2]]%*%t(U.hat[[2]])
        P.U3 = U.hat[[3]]%*%t(U.hat[[3]])
        Y.hat.3 = ttm(ttm(ttm(Y.tensor, P.U1, 1), P.U2, 2), P.U3, 3) 
        
        probability.3  = Y.hat.3@data
        
        probability.3[probability.3 > 1]  = 1
        
        probability.3[probability.3 < 0]  = 0
        
        result_hosvd[ind_n_1, ind_n_2, ind_L, ind_d, iter] = fnorm(as.tensor(probability.3) - as.tensor(probability))
        
        probability.4[probability.4 > 1]  = 1
        
        probability.4[probability.4 < 0]  = 0
        
        result_svd[ind_n_1, ind_n_2, ind_L, ind_d, iter] = fnorm(as.tensor(probability.4) - as.tensor(probability))
     
        probability.7 = multiness_adaptive(A, hat.rank[1])
        result_MultiNeSS_ada[ind_n_1, ind_n_2, ind_L, ind_d, iter] = fnorm(as.tensor(probability.7) - as.tensor(probability))
        
      }
      
    }
  }
}
}



apply(result_hooi, c(1,2,3, 4), mean)
apply(result_hetero, c(1,2,3, 4), mean)
apply(result_hosvd, c(1,2,3, 4), mean)
apply(result_svd, c(1,2,3, 4), mean)
apply(result_uase, c(1,2,3, 4), mean)
apply(result_MultiNeSS_ada, c(1,2,3, 4), mean)

apply(result_hooi, c(1,2,3, 4), sd)
apply(result_hetero, c(1,2,3, 4), sd)
apply(result_hosvd, c(1,2,3, 4), sd)
apply(result_svd, c(1,2,3, 4), sd)
apply(result_uase, c(1,2,3, 4), sd)
apply(result_MultiNeSS_ada, c(1,2,3, 4), sd)

value= c( apply(result_hooi, c(1,2,3,4), mean)[,,,],
          apply(result_hetero, c(1,2,3,4), mean)[,,,],
          apply(result_hosvd, c(1,2,3,4), mean)[,,,],
          apply(result_svd, c(1,2,3,4), mean)[,,,],
          apply(result_uase, c(1,2,3,4), mean)[,,,],
          apply(result_MultiNeSS_ada, c(1,2,3, 4), mean)[,,,]
          )

sd = c( apply(result_hooi, c(1,2,3,4), sd)[,,,],
        apply(result_hetero, c(1,2,3,4), sd)[,,,],
        apply(result_hosvd, c(1,2,3,4), sd)[,,,],
        apply(result_svd, c(1,2,3,4), sd)[,,,],
        apply(result_uase, c(1,2,3,4), sd)[,,,],
        apply(result_MultiNeSS_ada, c(1,2,3, 4), sd)[,,,]
        )

layer = rep(c(10,20,30,40), times = 6)

Method = rep(c('HOOI','TH-PCA','HOSVD','SASE','UASE','COA'),each = 4) 

dff_1 = data.frame( layer = layer, Method = Method, value = value, sd = sd)


plott_1= ggplot(data = dff_1, mapping = aes(x = layer, y = value, colour = Method)) +
  geom_line() + 
  geom_point()+ 
  geom_errorbar(aes(ymin=value-sd, ymax= value+sd), width= 2)+
  labs( y ="Estimation error", x = "Number of layers")+
  scale_color_manual(values = c("#66c2a5", "#8da0cb", "#e78ac3","#a6d854", "#ffd92f","#e5c494"))+
  theme_classic()




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


for(ind_n_1 in 1:length(n_grid_X))
{
  n_1  = n_grid_X[ind_n_1]
  
  print(paste0("n1 = ", n_1))
  
  for(ind_n_2 in 1:length(n_grid_Y))
    
  {
    n_2  = n_grid_Y[ind_n_2]
    
    print(paste0("n2 = ", n_2))
    
    for(ind_L in 1: length(L_grid) )
    {
      L = L_grid[ind_L]
      
      print(paste0("L = ", L))
      
      for(ind_d in 1: length(d_grid) )
      {
        d = d_grid[ind_d]
        
        print(paste0("d = ", d))
        
        dim = c(n_1, n_2, L)
        
        probability = array(NA,dim)
        
        Sigma_tensor =  array(NA,c(d, d, L))
        
        prob = seq(0,1,  1/(2*L))
        
        set.seed(n_1*n_2*L*d)
        
        for (layer in 1: L)
        {
          
          Sigma = matrix(runif(d^2,0, 1), ncol=d) 
          Sigma_tensor[, , layer] = Sigma
        }
        
        for(iter  in 1:NMC)
        { 
          A = array(NA,dim)
          
          probability.4 = array(NA,dim)
          
          for (layer in 1: L)
          {
            
            temp_1 = rdirichlet(n_1, rep(1, d))
            temp_2 = rdirichlet(n_2, rep(10, d))
            P =  temp_1  %*% Sigma_tensor[, , layer] %*% t(temp_2)
            probability[, , layer] = P
            
            A[, , layer] = matrix(rbinom(matrix(1,n_1,n_2),matrix(1,n_1,n_2), probability[ , , layer]),n_1,n_2)
            
            temp = svd(A[, , layer])
            
            xhat =  temp$u[,1:hat.rank[1]] %*%  diag( sqrt(temp$d[1:hat.rank[1]]) )
            
            yhat =  temp$v[,1:hat.rank[1]] %*%  diag( sqrt(temp$d[1:hat.rank[1]]) )
            
            probability.4[, , layer] = xhat %*% t(yhat)
            
          }
          
          Y.tensor =  as.tensor(A)
          
          probability.5 =   uase(Y.tensor, hat.rank[1] ) 
          
          result_uase[ind_n_1, ind_n_2, ind_L, ind_d, iter] = fnorm(as.tensor(probability.5) - as.tensor(probability))
          
          U.hat = HOOI(Y.tensor, hat.rank)
          
          P.U1 = U.hat[[1]]%*%t(U.hat[[1]])
          P.U2 = U.hat[[2]]%*%t(U.hat[[2]])
          P.U3 = U.hat[[3]]%*%t(U.hat[[3]])
          Y.hat.1 = ttm(ttm(ttm(Y.tensor, P.U1, 1), P.U2, 2), P.U3, 3) 
          
          probability.1  = Y.hat.1@data
          
          probability.1[probability.1 > 1]  = 1
          
          probability.1[probability.1 < 0]  = 0
          
          result_hooi[ind_n_1, ind_n_2, ind_L, ind_d, iter] = fnorm(as.tensor(probability.1) - as.tensor(probability))
          
          
          U.hat = Tensor_Hetero_PCA_test(Y.tensor, hat.rank)
          P.U1 = U.hat[[1]]%*%t(U.hat[[1]])
          P.U2 = U.hat[[2]]%*%t(U.hat[[2]])
          P.U3 = U.hat[[3]]%*%t(U.hat[[3]])
          Y.hat.2 = ttm(ttm(ttm(Y.tensor, P.U1, 1), P.U2, 2), P.U3, 3) 
          
          probability.2  = Y.hat.2@data
          
          probability.2[probability.2 > 1]  = 1
          
          probability.2[probability.2 < 0]  = 0
          
          result_hetero[ind_n_1, ind_n_2, ind_L, ind_d, iter] = fnorm(as.tensor(probability.2) - as.tensor(probability))
          
          
          U.hat = hosvd(Y.tensor, hat.rank)$U
          P.U1 = U.hat[[1]]%*%t(U.hat[[1]])
          P.U2 = U.hat[[2]]%*%t(U.hat[[2]])
          P.U3 = U.hat[[3]]%*%t(U.hat[[3]])
          Y.hat.3 = ttm(ttm(ttm(Y.tensor, P.U1, 1), P.U2, 2), P.U3, 3) 
          
          probability.3  = Y.hat.3@data
          
          probability.3[probability.3 > 1]  = 1
          
          probability.3[probability.3 < 0]  = 0
          
          result_hosvd[ind_n_1, ind_n_2, ind_L, ind_d, iter] = fnorm(as.tensor(probability.3) - as.tensor(probability))
          
          probability.4[probability.4 > 1]  = 1
          
          probability.4[probability.4 < 0]  = 0
          
          result_svd[ind_n_1, ind_n_2, ind_L, ind_d, iter] = fnorm(as.tensor(probability.4) - as.tensor(probability))
          
          probability.7 = multiness_adaptive(A, hat.rank[1])
          
          result_MultiNeSS_ada[ind_n_1, ind_n_2, ind_L, ind_d, iter] = fnorm(as.tensor(probability.7) - as.tensor(probability))
          
          
        }
        
      }
    }
  }
}



apply(result_hooi, c(1,2,3, 4), mean)
apply(result_hetero, c(1,2,3, 4), mean)
apply(result_hosvd, c(1,2,3, 4), mean)
apply(result_svd, c(1,2,3, 4), mean)
apply(result_uase, c(1,2,3, 4), mean)
apply(result_MultiNeSS_ada, c(1,2,3, 4), mean)

apply(result_hooi, c(1,2,3, 4), sd)
apply(result_hetero, c(1,2,3, 4), sd)
apply(result_hosvd, c(1,2,3, 4), sd)
apply(result_svd, c(1,2,3, 4), sd)
apply(result_uase, c(1,2,3, 4), sd)
apply(result_MultiNeSS_ada, c(1,2,3, 4), sd)

value= c( apply(result_hooi, c(1,2,3,4), mean)[,,,],
          apply(result_hetero, c(1,2,3,4), mean)[,,,],
          apply(result_hosvd, c(1,2,3,4), mean)[,,,],
          apply(result_svd, c(1,2,3,4), mean)[,,,],
          apply(result_uase, c(1,2,3,4), mean)[,,,],
          apply(result_MultiNeSS_ada, c(1,2,3, 4), mean)[,,,])

sd = c( apply(result_hooi, c(1,2,3,4), sd)[,,,],
        apply(result_hetero, c(1,2,3,4), sd)[,,,],
        apply(result_hosvd, c(1,2,3,4), sd)[,,,],
        apply(result_svd, c(1,2,3,4), sd)[,,,],
        apply(result_uase, c(1,2,3,4), sd)[,,,],
        apply(result_MultiNeSS_ada, c(1,2,3, 4), sd)[,,,])

layer = rep(c(10,20,30,40), times = 6)

Method = rep(c('HOOI','TH-PCA','HOSVD','SASE','UASE', 'COA'), each = 4) 

dff_1_1 = data.frame( layer = layer, Method = Method, value = value, sd = sd)

plott_1_1 = ggplot(data = dff_1_1, mapping = aes(x = layer, y = value, colour = Method)) +
  geom_line() + 
  geom_point()+ 
  geom_errorbar(aes(ymin=value-sd, ymax= value+sd), width= 2)+
  labs( y ="Estimation error", x = "Number of layers")+
  scale_color_manual(values = c("#66c2a5", "#8da0cb", "#e78ac3","#a6d854", "#ffd92f","#e5c494"))+
  theme_classic()



n_grid_X  = c(100)

n_grid_Y  = c(100)

L_grid = c(20)

d_grid = c(2, 4, 6, 8)

NMC = 100

hat.rank =c(10,10,10)

result_hooi = array(NA, c( length(n_grid_X), length(n_grid_Y), length(L_grid), length(d_grid), NMC))
result_hetero = array(NA, c( length(n_grid_X),length(n_grid_Y),  length(L_grid), length(d_grid),  NMC))
result_hosvd = array(NA, c( length(n_grid_X), length(n_grid_Y), length(L_grid), length(d_grid), NMC))
result_svd =   array(NA, c( length(n_grid_X), length(n_grid_Y),  length(L_grid), length(d_grid),  NMC))
result_uase =   array(NA, c( length(n_grid_X), length(n_grid_Y), length(L_grid),length(d_grid), NMC))
result_MultiNeSS_ada =   array(NA, c( length(n_grid_X), length(n_grid_Y), length(L_grid),length(d_grid), NMC))


for(ind_n_1 in 1:length(n_grid_X))
{
  n_1  = n_grid_X[ind_n_1]
  
  print(paste0("n1 = ", n_1))
  
  for(ind_n_2 in 1:length(n_grid_Y))
    
  {
    n_2  = n_grid_Y[ind_n_2]
    
    print(paste0("n2 = ", n_2))
    
    for(ind_L in 1: length(L_grid) )
    {
      L = L_grid[ind_L]
      
      print(paste0("L = ", L))
      
      for(ind_d in 1: length(d_grid) )
      {
        d = d_grid[ind_d]
        
        print(paste0("d = ", d))
        
        dim = c(n_1, n_2, L)
        
        probability = array(NA,dim)
        
        Sigma_tensor =  array(NA,c(d, d, L))
        
        prob = seq(0,1,  1/(2*L))
        
        set.seed(n_1*n_2*L*d)
        
        for (layer in 1: L)
        {
          
          Sigma = matrix(runif(d^2,prob[L+layer], prob[L+layer+1]), ncol=d) 
          Sigma_tensor[, , layer] = Sigma
        }
        
        for(iter  in 1:NMC)
        { 
          A = array(NA,dim)
          
          probability.4 = array(NA,dim)
          
          for (layer in 1: L)
          {
            
            temp_1 = rdirichlet(n_1, rep(1, d))
            temp_2 = rdirichlet(n_2, rep(10, d))
            P =  temp_1  %*% Sigma_tensor[, , layer] %*% t(temp_2)
            probability[, , layer] = P
            
            A[, , layer] = matrix(rbinom(matrix(1,n_1,n_2),matrix(1,n_1,n_2), probability[ , , layer]),n_1,n_2)
            
            temp = svd(A[, , layer])
            
            xhat =  temp$u[,1:hat.rank[1]] %*%  diag( sqrt(temp$d[1:hat.rank[1]]) )
            
            yhat =  temp$v[,1:hat.rank[1]] %*%  diag( sqrt(temp$d[1:hat.rank[1]]) )
            
            probability.4[, , layer] = xhat %*% t(yhat)
            
          }
          Y.tensor =  as.tensor(A)
          
          probability.5 =   uase(Y.tensor, hat.rank[1] ) 
          
          result_uase[ind_n_1, ind_n_2, ind_L, ind_d, iter] = fnorm(as.tensor(probability.5) - as.tensor(probability))
          
          U.hat = HOOI(Y.tensor, hat.rank)
          
          P.U1 = U.hat[[1]]%*%t(U.hat[[1]])
          P.U2 = U.hat[[2]]%*%t(U.hat[[2]])
          P.U3 = U.hat[[3]]%*%t(U.hat[[3]])
          Y.hat.1 = ttm(ttm(ttm(Y.tensor, P.U1, 1), P.U2, 2), P.U3, 3) 
          
          probability.1  = Y.hat.1@data
          
          probability.1[probability.1 > 1]  = 1
          
          probability.1[probability.1 < 0]  = 0
          
          result_hooi[ind_n_1, ind_n_2, ind_L, ind_d, iter] = fnorm(as.tensor(probability.1) - as.tensor(probability))
          
          
          U.hat = Tensor_Hetero_PCA_test(Y.tensor, hat.rank)
          P.U1 = U.hat[[1]]%*%t(U.hat[[1]])
          P.U2 = U.hat[[2]]%*%t(U.hat[[2]])
          P.U3 = U.hat[[3]]%*%t(U.hat[[3]])
          Y.hat.2 = ttm(ttm(ttm(Y.tensor, P.U1, 1), P.U2, 2), P.U3, 3) 
          
          probability.2  = Y.hat.2@data
          
          probability.2[probability.2 > 1]  = 1
          
          probability.2[probability.2 < 0]  = 0
          
          result_hetero[ind_n_1, ind_n_2, ind_L, ind_d, iter] = fnorm(as.tensor(probability.2) - as.tensor(probability))
          
          
          U.hat = hosvd(Y.tensor, hat.rank)$U
          P.U1 = U.hat[[1]]%*%t(U.hat[[1]])
          P.U2 = U.hat[[2]]%*%t(U.hat[[2]])
          P.U3 = U.hat[[3]]%*%t(U.hat[[3]])
          Y.hat.3 = ttm(ttm(ttm(Y.tensor, P.U1, 1), P.U2, 2), P.U3, 3) 
          
          probability.3  = Y.hat.3@data
          
          probability.3[probability.3 > 1]  = 1
          
          probability.3[probability.3 < 0]  = 0
          
          result_hosvd[ind_n_1, ind_n_2, ind_L, ind_d, iter] = fnorm(as.tensor(probability.3) - as.tensor(probability))
          
          probability.4[probability.4 > 1]  = 1
          
          probability.4[probability.4 < 0]  = 0
          
          result_svd[ind_n_1, ind_n_2, ind_L, ind_d, iter] = fnorm(as.tensor(probability.4) - as.tensor(probability))
          
          probability.7 = multiness_adaptive(A, hat.rank[1])
          
          result_MultiNeSS_ada[ind_n_1, ind_n_2, ind_L, ind_d, iter] = fnorm(as.tensor(probability.7) - as.tensor(probability))
          
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
          apply(result_MultiNeSS_ada, c(1,2,3, 4), mean)[,,,]
          )

sd = c( apply(result_hooi, c(1,2,3,4), sd)[,,,],
        apply(result_hetero, c(1,2,3,4), sd)[,,,],
        apply(result_hosvd, c(1,2,3,4), sd)[,,,],
        apply(result_svd, c(1,2,3,4), sd)[,,,],
        apply(result_uase, c(1,2,3,4), sd)[,,,],
        apply(result_MultiNeSS_ada, c(1,2,3, 4), sd)[,,,]
        )


dimension = rep( c(2, 4, 6, 8), times = 6)

Method = rep(c('HOOI','TH-PCA','HOSVD','SASE','UASE', 'COA'),each = 4)

dff_2 <- data.frame( dimension = dimension, Method = Method, value = value, sd = sd)

plott_2 = ggplot(data = dff_2, mapping = aes(x = dimension, y = value, colour = Method)) + geom_line() + geom_point() +
  geom_errorbar(aes(ymin=value-sd, ymax= value+sd), width= 0.4)+
  labs( y ="Estimation error", x = "Dimension of latent positions")+
  scale_color_manual(values = c("#66c2a5", "#8da0cb", "#e78ac3","#a6d854", "#ffd92f","#e5c494"))+
  theme_classic()




n_grid_X  = c(100)

n_grid_Y  = c(100)

L_grid = c(20)

d_grid = c(2, 4, 6, 8)

NMC = 100



result_hooi = array(NA, c( length(n_grid_X), length(n_grid_Y), length(L_grid), length(d_grid), NMC))
result_hetero = array(NA, c( length(n_grid_X),length(n_grid_Y),  length(L_grid), length(d_grid),  NMC))
result_hosvd = array(NA, c( length(n_grid_X), length(n_grid_Y), length(L_grid), length(d_grid), NMC))
result_svd =   array(NA, c( length(n_grid_X), length(n_grid_Y),  length(L_grid), length(d_grid),  NMC))
result_uase =   array(NA, c( length(n_grid_X), length(n_grid_Y), length(L_grid),length(d_grid), NMC))
result_MultiNeSS_ada =   array(NA, c( length(n_grid_X), length(n_grid_Y), length(L_grid),length(d_grid), NMC))

for(ind_n_1 in 1:length(n_grid_X))
{
  n_1  = n_grid_X[ind_n_1]
  
  print(paste0("n1 = ", n_1))
  
  for(ind_n_2 in 1:length(n_grid_Y))
    
  {
    n_2  = n_grid_Y[ind_n_2]
    
    print(paste0("n2 = ", n_2))
    
    for(ind_L in 1: length(L_grid) )
    {
      L = L_grid[ind_L]
      
      print(paste0("L = ", L))
      
      for(ind_d in 1: length(d_grid) )
      {
        d = d_grid[ind_d]
        
        hat.rank =c(d, d, d)
        
        print(paste0("d = ", d))
        
        dim = c(n_1, n_2, L)
        
        probability = array(NA,dim)
        
        Sigma_tensor =  array(NA,c(d, d, L))
        
        prob = seq(0,1,  1/(2*L))
        
        set.seed(n_1*n_2*L*d)
        
        for (layer in 1: L)
        {
          
          Sigma = matrix(runif(d^2,prob[L+layer], prob[L+layer+1]), ncol=d) 
          Sigma_tensor[, , layer] = Sigma
        }
        
        for(iter  in 1:NMC)
        { 
          A = array(NA,dim)
          
          probability.4 = array(NA,dim)
          
          for (layer in 1: L)
          {
            
            temp_1 = rdirichlet(n_1, rep(1, d))
            temp_2 = rdirichlet(n_2, rep(10, d))
            P =  temp_1  %*% Sigma_tensor[, , layer] %*% t(temp_2)
            probability[, , layer] = P
            
            A[, , layer] = matrix(rbinom(matrix(1,n_1,n_2),matrix(1,n_1,n_2), probability[ , , layer]),n_1,n_2)
            
            temp = svd(A[, , layer])
            
            xhat =  temp$u[,1:hat.rank[1]] %*%  diag( sqrt(temp$d[1:hat.rank[1]]) )
            
            yhat =  temp$v[,1:hat.rank[1]] %*%  diag( sqrt(temp$d[1:hat.rank[1]]) )
            
            probability.4[, , layer] = xhat %*% t(yhat)
            
          }
          
          Y.tensor =  as.tensor(A)
          
          probability.5 =   uase(Y.tensor, hat.rank[1] ) 
          
          result_uase[ind_n_1, ind_n_2, ind_L, ind_d, iter] = fnorm(as.tensor(probability.5) - as.tensor(probability))
          
          U.hat = HOOI(Y.tensor, hat.rank)
          
          P.U1 = U.hat[[1]]%*%t(U.hat[[1]])
          P.U2 = U.hat[[2]]%*%t(U.hat[[2]])
          P.U3 = U.hat[[3]]%*%t(U.hat[[3]])
          Y.hat.1 = ttm(ttm(ttm(Y.tensor, P.U1, 1), P.U2, 2), P.U3, 3) 
          
          probability.1  = Y.hat.1@data
          
          probability.1[probability.1 > 1]  = 1
          
          probability.1[probability.1 < 0]  = 0
          
          result_hooi[ind_n_1, ind_n_2, ind_L, ind_d, iter] = fnorm(as.tensor(probability.1) - as.tensor(probability))
          
          
          U.hat = Tensor_Hetero_PCA_test(Y.tensor, hat.rank)
          P.U1 = U.hat[[1]]%*%t(U.hat[[1]])
          P.U2 = U.hat[[2]]%*%t(U.hat[[2]])
          P.U3 = U.hat[[3]]%*%t(U.hat[[3]])
          Y.hat.2 = ttm(ttm(ttm(Y.tensor, P.U1, 1), P.U2, 2), P.U3, 3) 
          
          probability.2  = Y.hat.2@data
          
          probability.2[probability.2 > 1]  = 1
          
          probability.2[probability.2 < 0]  = 0
          
          result_hetero[ind_n_1, ind_n_2, ind_L, ind_d, iter] = fnorm(as.tensor(probability.2) - as.tensor(probability))
          
          
          U.hat = hosvd(Y.tensor, hat.rank)$U
          P.U1 = U.hat[[1]]%*%t(U.hat[[1]])
          P.U2 = U.hat[[2]]%*%t(U.hat[[2]])
          P.U3 = U.hat[[3]]%*%t(U.hat[[3]])
          Y.hat.3 = ttm(ttm(ttm(Y.tensor, P.U1, 1), P.U2, 2), P.U3, 3) 
          
          probability.3  = Y.hat.3@data
          
          probability.3[probability.3 > 1]  = 1
          
          probability.3[probability.3 < 0]  = 0
          
          result_hosvd[ind_n_1, ind_n_2, ind_L, ind_d, iter] = fnorm(as.tensor(probability.3) - as.tensor(probability))
          
          probability.4[probability.4 > 1]  = 1
          
          probability.4[probability.4 < 0]  = 0
          
          result_svd[ind_n_1, ind_n_2, ind_L, ind_d, iter] = fnorm(as.tensor(probability.4) - as.tensor(probability))
          
          probability.7 = multiness_adaptive(A, hat.rank[1])
          
          result_MultiNeSS_ada[ind_n_1, ind_n_2, ind_L, ind_d, iter] = fnorm(as.tensor(probability.7) - as.tensor(probability))
          
          
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
          apply(result_MultiNeSS_ada, c(1,2,3, 4), mean)[,,,]
          )

sd = c( apply(result_hooi, c(1,2,3,4), sd)[,,,],
        apply(result_hetero, c(1,2,3,4), sd)[,,,],
        apply(result_hosvd, c(1,2,3,4), sd)[,,,],
        apply(result_svd, c(1,2,3,4), sd)[,,,],
        apply(result_uase, c(1,2,3,4), sd)[,,,],
        apply(result_MultiNeSS_ada, c(1,2,3, 4), sd)[,,,])

dimension = rep( c(2, 4, 6, 8), times = 6)

Method = rep(c('HOOI','TH-PCA','HOSVD','SASE','UASE','COA'),each = 4)

dff_2_1 <- data.frame( dimension = dimension, Method = Method, value = value, sd = sd)

plott_2_1 = ggplot(data = dff_2_1, mapping = aes(x = dimension, y = value, colour = Method)) + geom_line() + geom_point() +
  geom_errorbar(aes(ymin=value-sd, ymax= value+sd), width= 0.4)+
  labs( y ="Estimation error", x = "Number of latent dimension")+
  scale_color_manual(values = c("#66c2a5", "#8da0cb", "#e78ac3","#a6d854", "#ffd92f","#e5c494"))+
  theme_classic()




n_grid_X  = c(100)

n_grid_Y  = c(100)

L_grid = c(20)

d_grid = c(2, 4, 6, 8)

NMC = 100

hat.rank =c(2,2,2)

result_hooi = array(NA, c( length(n_grid_X), length(n_grid_Y), length(L_grid), length(d_grid), NMC))
result_hetero = array(NA, c( length(n_grid_X),length(n_grid_Y),  length(L_grid), length(d_grid),  NMC))
result_hosvd = array(NA, c( length(n_grid_X), length(n_grid_Y), length(L_grid), length(d_grid), NMC))
result_svd =   array(NA, c( length(n_grid_X), length(n_grid_Y),  length(L_grid), length(d_grid),  NMC))
result_uase =   array(NA, c( length(n_grid_X), length(n_grid_Y), length(L_grid),length(d_grid), NMC))
result_MultiNeSS_ada =   array(NA, c( length(n_grid_X), length(n_grid_Y), length(L_grid),length(d_grid), NMC))

for(ind_n_1 in 1:length(n_grid_X))
{
  n_1  = n_grid_X[ind_n_1]
  
  print(paste0("n1 = ", n_1))
  
  for(ind_n_2 in 1:length(n_grid_Y))
    
  {
    n_2  = n_grid_Y[ind_n_2]
    
    print(paste0("n2 = ", n_2))
    
    for(ind_L in 1: length(L_grid) )
    {
      L = L_grid[ind_L]
      
      print(paste0("L = ", L))
      
      for(ind_d in 1: length(d_grid) )
      {
        d = d_grid[ind_d]
        
        print(paste0("d = ", d))
        
        dim = c(n_1, n_2, L)
        
        probability = array(NA,dim)
        
        Sigma_tensor =  array(NA,c(d, d, L))
        
        prob = seq(0,1,  1/(2*L))
        
        set.seed(n_1*n_2*L*d)
        
        for (layer in 1: L)
        {
          
          Sigma = matrix(runif(d^2,prob[L+layer], prob[L+layer+1]), ncol=d) 
          Sigma_tensor[, , layer] = Sigma
        }
        
        for(iter  in 1:NMC)
        { 
          A = array(NA,dim)
          
          probability.4 = array(NA,dim)
          
          for (layer in 1: L)
          {
            
            temp_1 = rdirichlet(n_1, rep(1, d))
            temp_2 = rdirichlet(n_2, rep(10, d))
            P =  temp_1  %*% Sigma_tensor[, , layer] %*% t(temp_2)
            probability[, , layer] = P
            
            A[, , layer] = matrix(rbinom(matrix(1,n_1,n_2),matrix(1,n_1,n_2), probability[ , , layer]),n_1,n_2)
            
            temp = svd(A[, , layer])
            
            xhat =  temp$u[,1:hat.rank[1]] %*%  diag( sqrt(temp$d[1:hat.rank[1]]) )
            
            yhat =  temp$v[,1:hat.rank[1]] %*%  diag( sqrt(temp$d[1:hat.rank[1]]) )
            
            probability.4[, , layer] = xhat %*% t(yhat)
            
          }
          
          Y.tensor =  as.tensor(A)
          
          probability.5 =   uase(Y.tensor, hat.rank[1] ) 
          
          result_uase[ind_n_1, ind_n_2, ind_L, ind_d, iter] = fnorm(as.tensor(probability.5) - as.tensor(probability))
          
          U.hat = HOOI(Y.tensor, hat.rank)
          
          P.U1 = U.hat[[1]]%*%t(U.hat[[1]])
          P.U2 = U.hat[[2]]%*%t(U.hat[[2]])
          P.U3 = U.hat[[3]]%*%t(U.hat[[3]])
          Y.hat.1 = ttm(ttm(ttm(Y.tensor, P.U1, 1), P.U2, 2), P.U3, 3) 
          
          probability.1  = Y.hat.1@data
          
          probability.1[probability.1 > 1]  = 1
          
          probability.1[probability.1 < 0]  = 0
          
          result_hooi[ind_n_1, ind_n_2, ind_L, ind_d, iter] = fnorm(as.tensor(probability.1) - as.tensor(probability))
          
          
          U.hat = Tensor_Hetero_PCA_test(Y.tensor, hat.rank)
          P.U1 = U.hat[[1]]%*%t(U.hat[[1]])
          P.U2 = U.hat[[2]]%*%t(U.hat[[2]])
          P.U3 = U.hat[[3]]%*%t(U.hat[[3]])
          Y.hat.2 = ttm(ttm(ttm(Y.tensor, P.U1, 1), P.U2, 2), P.U3, 3) 
          
          probability.2  = Y.hat.2@data
          
          probability.2[probability.2 > 1]  = 1
          
          probability.2[probability.2 < 0]  = 0
          
          result_hetero[ind_n_1, ind_n_2, ind_L, ind_d, iter] = fnorm(as.tensor(probability.2) - as.tensor(probability))
          
          
          U.hat = hosvd(Y.tensor, hat.rank)$U
          P.U1 = U.hat[[1]]%*%t(U.hat[[1]])
          P.U2 = U.hat[[2]]%*%t(U.hat[[2]])
          P.U3 = U.hat[[3]]%*%t(U.hat[[3]])
          Y.hat.3 = ttm(ttm(ttm(Y.tensor, P.U1, 1), P.U2, 2), P.U3, 3) 
          
          probability.3  = Y.hat.3@data
          
          probability.3[probability.3 > 1]  = 1
          
          probability.3[probability.3 < 0]  = 0
          
          result_hosvd[ind_n_1, ind_n_2, ind_L, ind_d, iter] = fnorm(as.tensor(probability.3) - as.tensor(probability))
          
          probability.4[probability.4 > 1]  = 1
          
          probability.4[probability.4 < 0]  = 0
          
          result_svd[ind_n_1, ind_n_2, ind_L, ind_d, iter] = fnorm(as.tensor(probability.4) - as.tensor(probability))
          
          probability.7 = multiness_adaptive(A, hat.rank[1])
          
          result_MultiNeSS_ada[ind_n_1, ind_n_2, ind_L, ind_d, iter] = fnorm(as.tensor(probability.7) - as.tensor(probability))
          
          
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
          apply(result_MultiNeSS_ada, c(1,2,3, 4), mean)[,,,]
)

sd = c( apply(result_hooi, c(1,2,3,4), sd)[,,,],
        apply(result_hetero, c(1,2,3,4), sd)[,,,],
        apply(result_hosvd, c(1,2,3,4), sd)[,,,],
        apply(result_svd, c(1,2,3,4), sd)[,,,],
        apply(result_uase, c(1,2,3,4), sd)[,,,],
        apply(result_MultiNeSS_ada, c(1,2,3, 4), sd)[,,,])

dimension = rep( c(2, 4, 6, 8), times = 6)

Method = rep(c('HOOI','TH-PCA','HOSVD','SASE','UASE','COA'),each = 4)

dff_2_2 <- data.frame( dimension = dimension, Method = Method, value = value, sd = sd)

plott_2_2 = ggplot(data = dff_2_2, mapping = aes(x = dimension, y = value, colour = Method)) + geom_line() + geom_point() +
  geom_errorbar(aes(ymin=value-sd, ymax= value+sd), width= 0.4)+
  labs( y ="Estimation error", x = "Dimension of latent positions")+
  scale_color_manual(values = c("#66c2a5", "#8da0cb", "#e78ac3","#a6d854", "#ffd92f","#e5c494"))+
  theme_classic()



n_grid_X  = c(50, 100, 150, 200)



L_grid = c(20)

d_grid = c(4)

NMC = 100

hat.rank =c(10,10,10)

result_hooi = array(NA, c( length(n_grid_X), length(n_grid_Y), length(L_grid), length(d_grid), NMC))
result_hetero = array(NA, c( length(n_grid_X),length(n_grid_Y),  length(L_grid), length(d_grid),  NMC))
result_hosvd = array(NA, c( length(n_grid_X), length(n_grid_Y), length(L_grid), length(d_grid), NMC))
result_svd =   array(NA, c( length(n_grid_X), length(n_grid_Y),  length(L_grid), length(d_grid),  NMC))
result_uase =   array(NA, c( length(n_grid_X), length(n_grid_Y), length(L_grid),length(d_grid), NMC))
result_MultiNeSS_ada =   array(NA, c( length(n_grid_X), length(n_grid_Y), length(L_grid),length(d_grid), NMC))

for(ind_n_1 in 1:length(n_grid_X))
{
  n_1  = n_grid_X[ind_n_1]
  
  print(paste0("n1 = ", n_1))
  
  for(ind_n_2 in 1:length(n_grid_Y))
    
  {
    n_2  = n_1
    
    print(paste0("n2 = ", n_2))
    
    for(ind_L in 1: length(L_grid) )
    {
      L = L_grid[ind_L]
      
      print(paste0("L = ", L))
      
      for(ind_d in 1: length(d_grid) )
      {
        d = d_grid[ind_d]
        
        print(paste0("d = ", d))
        
        dim = c(n_1, n_2, L)
        
        probability = array(NA,dim)
        
        Sigma_tensor =  array(NA,c(d, d, L))
        
        prob = seq(0,1,  1/(2*L))
        
        set.seed(n_1*n_2*L*d)
        
        for (layer in 1: L)
        {
          
          Sigma = matrix(runif(d^2,prob[L+layer], prob[L+layer+1]), ncol=d) 
          Sigma_tensor[, , layer] = Sigma
        }
        
        for(iter  in 1:NMC)
        { 
          A = array(NA,dim)
          
          probability.4 = array(NA,dim)
          
          for (layer in 1: L)
          {
            
            temp_1 = rdirichlet(n_1, rep(1, d))
            temp_2 = rdirichlet(n_2, rep(10, d))
            P =  temp_1  %*% Sigma_tensor[, , layer] %*% t(temp_2)
            probability[, , layer] = P
            
            A[, , layer] = matrix(rbinom(matrix(1,n_1,n_2),matrix(1,n_1,n_2), probability[ , , layer]),n_1,n_2)
            
            temp = svd(A[, , layer])
            
            xhat =  temp$u[,1:hat.rank[1]] %*%  diag( sqrt(temp$d[1:hat.rank[1]]) )
            
            yhat =  temp$v[,1:hat.rank[1]] %*%  diag( sqrt(temp$d[1:hat.rank[1]]) )
            
            probability.4[, , layer] = xhat %*% t(yhat)
            
          }
          
          Y.tensor =  as.tensor(A)
          
          probability.5 =   uase(Y.tensor, hat.rank[1] ) 
          
          result_uase[ind_n_1, ind_n_2, ind_L, ind_d, iter] = fnorm(as.tensor(probability.5) - as.tensor(probability))
          
          U.hat = HOOI(Y.tensor, hat.rank)
          
          P.U1 = U.hat[[1]]%*%t(U.hat[[1]])
          P.U2 = U.hat[[2]]%*%t(U.hat[[2]])
          P.U3 = U.hat[[3]]%*%t(U.hat[[3]])
          Y.hat.1 = ttm(ttm(ttm(Y.tensor, P.U1, 1), P.U2, 2), P.U3, 3) 
          
          probability.1  = Y.hat.1@data
          
          probability.1[probability.1 > 1]  = 1
          
          probability.1[probability.1 < 0]  = 0
          
          result_hooi[ind_n_1, ind_n_2, ind_L, ind_d, iter] = fnorm(as.tensor(probability.1) - as.tensor(probability))
          
          
          U.hat = Tensor_Hetero_PCA_test(Y.tensor, hat.rank)
          P.U1 = U.hat[[1]]%*%t(U.hat[[1]])
          P.U2 = U.hat[[2]]%*%t(U.hat[[2]])
          P.U3 = U.hat[[3]]%*%t(U.hat[[3]])
          Y.hat.2 = ttm(ttm(ttm(Y.tensor, P.U1, 1), P.U2, 2), P.U3, 3) 
          
          probability.2  = Y.hat.2@data
          
          probability.2[probability.2 > 1]  = 1
          
          probability.2[probability.2 < 0]  = 0
          
          result_hetero[ind_n_1, ind_n_2, ind_L, ind_d, iter] = fnorm(as.tensor(probability.2) - as.tensor(probability))
          
          
          U.hat = hosvd(Y.tensor, hat.rank)$U
          P.U1 = U.hat[[1]]%*%t(U.hat[[1]])
          P.U2 = U.hat[[2]]%*%t(U.hat[[2]])
          P.U3 = U.hat[[3]]%*%t(U.hat[[3]])
          Y.hat.3 = ttm(ttm(ttm(Y.tensor, P.U1, 1), P.U2, 2), P.U3, 3) 
          
          probability.3  = Y.hat.3@data
          
          probability.3[probability.3 > 1]  = 1
          
          probability.3[probability.3 < 0]  = 0
          
          result_hosvd[ind_n_1, ind_n_2, ind_L, ind_d, iter] = fnorm(as.tensor(probability.3) - as.tensor(probability))
          
          probability.4[probability.4 > 1]  = 1
          
          probability.4[probability.4 < 0]  = 0
          
          result_svd[ind_n_1, ind_n_2, ind_L, ind_d, iter] = fnorm(as.tensor(probability.4) - as.tensor(probability))
          
          probability.7 = multiness_adaptive(A, hat.rank[1])
          
          result_MultiNeSS_ada[ind_n_1, ind_n_2, ind_L, ind_d, iter] = fnorm(as.tensor(probability.7) - as.tensor(probability))
          
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
          apply(result_MultiNeSS_ada, c(1,2,3, 4), mean)[,,,])

sd = c( apply(result_hooi, c(1,2,3,4), sd)[,,,],
        apply(result_hetero, c(1,2,3,4), sd)[,,,],
        apply(result_hosvd, c(1,2,3,4), sd)[,,,],
        apply(result_svd, c(1,2,3,4), sd)[,,,],
        apply(result_uase, c(1,2,3,4), sd)[,,,],
        apply(result_MultiNeSS_ada, c(1,2,3, 4), sd)[,,,])

node = rep( c(50, 100, 150, 200), times = 6)

Method = rep(c('HOOI','TH-PCA','HOSVD','SASE','UASE', 'COA'), each = 4)

dff_3 <- data.frame( node = node, Method = Method, value = value, sd = sd)

plott_3 = ggplot(data = dff_3, mapping = aes(x = node, y = value, colour = Method)) + geom_line() + geom_point() +
  geom_errorbar(aes(ymin=value-sd, ymax= value+sd), width= 10)+
  labs( y ="Estimation error",  x =expression(paste( "Number of nodes: n with n=", n[1],"=",n[2])))+
  scale_color_manual(values = c("#66c2a5", "#8da0cb", "#e78ac3","#a6d854", "#ffd92f","#e5c494"))+
  theme_classic()





n_grid_X  = c(100)

n_grid_Y  = c(50, 100, 150, 200)

L_grid = c(20)

d_grid = c(4)

NMC = 100

hat.rank =c(10,10,10)

result_hooi = array(NA, c( length(n_grid_X), length(n_grid_Y), length(L_grid), length(d_grid), NMC))
result_hetero = array(NA, c( length(n_grid_X),length(n_grid_Y),  length(L_grid), length(d_grid),  NMC))
result_hosvd = array(NA, c( length(n_grid_X), length(n_grid_Y), length(L_grid), length(d_grid), NMC))
result_svd =   array(NA, c( length(n_grid_X), length(n_grid_Y),  length(L_grid), length(d_grid),  NMC))
result_uase =   array(NA, c( length(n_grid_X), length(n_grid_Y), length(L_grid),length(d_grid), NMC))


for(ind_n_1 in 1:length(n_grid_X))
{
  n_1  = n_grid_X[ind_n_1]
  
  print(paste0("n1 = ", n_1))
  
  for(ind_n_2 in 1:length(n_grid_Y))
    
  {
    n_2  = n_grid_Y[ind_n_2]
    
    print(paste0("n2 = ", n_2))
    
    for(ind_L in 1: length(L_grid) )
    {
      L = L_grid[ind_L]
      
      print(paste0("L = ", L))
      
      for(ind_d in 1: length(d_grid) )
      {
        d = d_grid[ind_d]
        
        print(paste0("d = ", d))
        
        dim = c(n_1, n_2, L)
        
        probability = array(NA,dim)
        
        Sigma_tensor =  array(NA,c(d, d, L))
        
        prob = seq(0,1,  1/(2*L))
        
        set.seed(n_1*n_2*L*d)
        
        for (layer in 1: L)
        {
          
          Sigma = matrix(runif(d^2,prob[L+layer], prob[L+layer+1]), ncol=d) 
          Sigma_tensor[, , layer] = Sigma
        }
        
        for(iter  in 1:NMC)
        { 
          A = array(NA,dim)
          
          probability.4 = array(NA,dim)
          
          for (layer in 1: L)
          {
            
            temp_1 = rdirichlet(n_1, rep(1, d))
            temp_2 = rdirichlet(n_2, rep(10, d))
            P =  temp_1  %*% Sigma_tensor[, , layer] %*% t(temp_2)
            probability[, , layer] = P
            
            A[, , layer] = matrix(rbinom(matrix(1,n_1,n_2),matrix(1,n_1,n_2), probability[ , , layer]),n_1,n_2)
            
            temp = svd(A[, , layer])
            
            xhat =  temp$u[,1:hat.rank[1]] %*%  diag( sqrt(temp$d[1:hat.rank[1]]) )
            
            yhat =  temp$v[,1:hat.rank[1]] %*%  diag( sqrt(temp$d[1:hat.rank[1]]) )
            
            probability.4[, , layer] = xhat %*% t(yhat)
            
          }
          
          Y.tensor =  as.tensor(A)
          
          probability.5 =   uase(Y.tensor, hat.rank[1] ) 
          
          result_uase[ind_n_1, ind_n_2, ind_L, ind_d, iter] = fnorm(as.tensor(probability.5) - as.tensor(probability))
          
          U.hat = HOOI(Y.tensor, hat.rank)
          
          P.U1 = U.hat[[1]]%*%t(U.hat[[1]])
          P.U2 = U.hat[[2]]%*%t(U.hat[[2]])
          P.U3 = U.hat[[3]]%*%t(U.hat[[3]])
          Y.hat.1 = ttm(ttm(ttm(Y.tensor, P.U1, 1), P.U2, 2), P.U3, 3) 
          
          probability.1  = Y.hat.1@data
          
          probability.1[probability.1 > 1]  = 1
          
          probability.1[probability.1 < 0]  = 0
          
          result_hooi[ind_n_1, ind_n_2, ind_L, ind_d, iter] = fnorm(as.tensor(probability.1) - as.tensor(probability))
          
          
          U.hat = Tensor_Hetero_PCA_test(Y.tensor, hat.rank)
          P.U1 = U.hat[[1]]%*%t(U.hat[[1]])
          P.U2 = U.hat[[2]]%*%t(U.hat[[2]])
          P.U3 = U.hat[[3]]%*%t(U.hat[[3]])
          Y.hat.2 = ttm(ttm(ttm(Y.tensor, P.U1, 1), P.U2, 2), P.U3, 3) 
          
          probability.2  = Y.hat.2@data
          
          probability.2[probability.2 > 1]  = 1
          
          probability.2[probability.2 < 0]  = 0
          
          result_hetero[ind_n_1, ind_n_2, ind_L, ind_d, iter] = fnorm(as.tensor(probability.2) - as.tensor(probability))
          
          
          U.hat = hosvd(Y.tensor, hat.rank)$U
          P.U1 = U.hat[[1]]%*%t(U.hat[[1]])
          P.U2 = U.hat[[2]]%*%t(U.hat[[2]])
          P.U3 = U.hat[[3]]%*%t(U.hat[[3]])
          Y.hat.3 = ttm(ttm(ttm(Y.tensor, P.U1, 1), P.U2, 2), P.U3, 3) 
          
          probability.3  = Y.hat.3@data
          
          probability.3[probability.3 > 1]  = 1
          
          probability.3[probability.3 < 0]  = 0
          
          result_hosvd[ind_n_1, ind_n_2, ind_L, ind_d, iter] = fnorm(as.tensor(probability.3) - as.tensor(probability))
          
          probability.4[probability.4 > 1]  = 1
          
          probability.4[probability.4 < 0]  = 0
          
          result_svd[ind_n_1, ind_n_2, ind_L, ind_d, iter] = fnorm(as.tensor(probability.4) - as.tensor(probability))

          
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
  scale_color_manual(values = c("#8da0cb", "#e78ac3","#a6d854", "#ffd92f","#e5c494"))+
  theme_classic()




ggarrange(plott_1, plott_1_1, plott_2, plott_2_1, plott_3, plott_4, ncol = 2)






