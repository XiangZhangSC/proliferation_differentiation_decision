library(dplyr)
library(tidyr)
library(stringr)

relative_expression <- read.csv("model_genes_relative_expression.csv")
lin1 <- relative_expression %>% 
  dplyr::filter(external_gene_id == "lin-1") %>% 
  tidyr::pivot_longer(hour00:hour20, names_to = "hour", values_to = "expression") %>% 
  dplyr::mutate(hour = as.numeric(str_sub(hour, 5, 6))) %>% 
  dplyr::select(hour, expression)
  
f_lin1 <- approxfun(lin1, rule = 2)

mls2 <- relative_expression %>% 
  dplyr::filter(external_gene_id == "mls-2") %>% 
  tidyr::pivot_longer(hour00:hour20, names_to = "hour", values_to = "expression") %>% 
  dplyr::mutate(hour = as.numeric(str_sub(hour, 5, 6))) %>% 
  dplyr::select(hour, expression)

f_mls2 <- approxfun(mls2, rule = 2)
