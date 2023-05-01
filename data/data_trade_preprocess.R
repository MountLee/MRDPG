#### local ####
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# # source("../generate_data.R")
# source("../simulation_wrap.R")
# Rcpp::sourceCpp("../cpd_hpca.cpp")
# Rcpp::sourceCpp("../cpd_uase.cpp")
# source("../other_methods.R")


#### main ####

# library(mvtnorm)
# library(dirmult)
library(rTensor)
# library(Rcpp)
library(tictoc)

# library(multiness)
# library(gStream)

# library(egg)
# library(ggplot2)



require(data.table)
data_whole = fread("Trade_DetailedTradeMatrix_E_All_Data.csv", header = T)

a = data_whole[, 70]


dim(data_whole)

colnames(data_whole)[c(1:15, 78, 79)]

data_whole = data_whole[, c(1:9, seq(10, 78, 2))]

colnames(data_whole)

a = data_whole[, 40]


unique(data_whole[, 8])

data_whole_ex = data_whole[which(data_whole[, 8] == "Import Value"),]

country_report = rle(sort(as.matrix(data_whole_ex[,2])))
country_report_100 = country_report$value[order(country_report$lengths, decreasing = T)[1:75]]

country_partner= rle(sort(as.matrix(data_whole_ex[,4])))
country_partner_100 = country_partner$value[order(country_partner$lengths, decreasing = T)[1:75]]

product_sort = rle(sort(as.matrix(data_whole_ex[,6])))
product_sort_4 = product_sort$value[order(product_sort$lengths, decreasing = T)[1:4]]


data_names = c()

for (i in 1:35) {
  data_names = c(data_names, paste("data_year", as.character(i + 1985), sep = "."))
}




dim = c(
  length(country_report_100),
  length(country_partner_100),
  length(product_sort_4),
  length(data_names)
)



data_year_product_tensor = array(rep(NA, length(country_report_100) * length(country_partner_100) * 
                                       length(product_sort_4) * length(data_names)), dim)

dim(data_year_product_tensor)


# u = 1
# year = colnames(data_whole_ex)[10 + (u - 1)]
# 
# 
# ix = c(2, 4, 6, 8, 10 + (u - 1))
# data_year = data_whole[, ..ix]
# a = data_year[, year]
# 
# a = data_year[, 'Item']
# 
# l = 1
# product_sort_4[l]
# 
# i = 1
# j = 1
# 
# country_report_100[i]
# country_partner_100[j]
# x = which(
#   data_year_product[, "Reporter Countries"] == country_report_100[i] &
#     data_year_product[, "Partner Countries"] ==  country_partner_100[j]
# )
# a = data_year_product[x, ]
# a
# 
# 
# 
# country_report_100[i]
# country_partner_100[j]
# x = which(
#   data_year_product[, "Reporter Countries"] == country_partner_100[i] &
#     data_year_product[, "Partner Countries"] ==  country_report_100[j]
# )
# a = data_year_product[x, ]
# a
# data_year_product[x, 5]


n_zeros = rep(0, length(data_names))
for (u in 1:length(data_names)) {
  year = colnames(data_whole_ex)[10 + (u - 1)]

  ix = c(2, 4, 6, 8, 10 + (u - 1))
  data_year = data_whole_ex[, ..ix]
  
  data_year[as.vector(is.na(data_year[, ..year])), year] = 0
  
  for (l in 1:length(product_sort_4)) {
    data_year_product = data_year[which(data_year[, "Item"] == product_sort_4[l]),]
    
    data_year_matrix =  matrix(rep(
      NA,
      length(country_report_100) * length(country_partner_100)
    ), ncol = length(country_partner_100))
    
    for (i in 1:length(country_report_100)) {
      for (j in 1:length(country_partner_100)) {
        x = which(
          data_year_product[, "Reporter Countries"] == country_report_100[i] &
            data_year_product[, "Partner Countries"] ==  country_partner_100[j]
        )
        
        if (length(x) > 0){
          data_year_matrix[i, j] = as.numeric(data_year_product[x, 5])
        }
        else{
          data_year_matrix[i, j] = 0
          n_zeros[u] = n_zeros[u] + 1
        }
        
      }
    }
    
    data_year_product_tensor[, , l, u] = data_year_matrix
  }
  
}


n_zeros
sum(n_zeros)

# data_year_product_tensor[is.na(data_year_product_tensor)] = 0

help(hist)

all_values = as.vector(data_year_product_tensor)
hist(all_values, )

log_values = log1p(all_values)
hist(log_values)



TT = 35
year_list = c()
value_list = c()
for (t in 1:TT){
  data_year = as.numeric(data_year_product_tensor[, , , t])
  pos_year = data_year[data_year > 0]
  value_list = c(value_list, pos_year)
  year_list = c(year_list, rep(1985 + t, length(pos_year)))
}

df = data.frame("values" = value_list, "year" = year_list)
df$log_values = log1p(df$values)

library(dplyr)
# library(hrbrthemes)

p <- df %>%
  ggplot( aes(x = as.factor(year), y = log_values) ) +
  geom_boxplot(fill="#69b3a2") +
  theme(axis.text.x = element_text(angle = 80, vjust = 0.8, hjust=1)) + 
  xlab("year")
p



save(
  data_year_product_tensor,
  country_report_100,
  country_partner_100,
  product_sort_4,
  data_names,
  df,
  file = "trade_tensor.RData"
)
