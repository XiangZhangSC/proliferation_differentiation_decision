library(dplyr)
library(tidyr)
library(stringr)

input_lin1 <- function(scenario) {
  if (scenario == "Molly") {
    relative_expression <- read.csv("model_genes_relative_expression.csv")
    lin1 <- relative_expression %>% 
      dplyr::filter(external_gene_id == "lin-1") %>% 
      tidyr::pivot_longer(hour00:hour20, names_to = "hour", values_to = "expression") %>% 
      dplyr::mutate(hour = as.numeric(str_sub(hour, 5, 6))) %>% 
      dplyr::select(hour, expression)
    
  } else if (scenario == "lin-1 ON mls-2 OFF") {
    lin1 <- data.frame(
      hour = seq(from = 0, to = 20, by = 2), 
      expression = c(rep(0, times = 4), rep(1, times = 4), rep(0, times = 3))
    )
  } else if (scenario == "lin-1 OFF mls-2 ON") {
    lin1 <- data.frame(
      hour = seq(from = 0, to = 20, by = 2), 
      expression = c(rep(0, times = 4), rep(0, times = 4), rep(0, times = 3))
    )
  }
  f_lin1 <- approxfun(lin1, rule = 2)
  return(f_lin1)
}

input_mls2 <- function(scenario) {
  if (scenario == "Molly") {
    relative_expression <- read.csv("model_genes_relative_expression.csv")
    mls2 <- relative_expression %>% 
      dplyr::filter(external_gene_id == "mls-2") %>% 
      tidyr::pivot_longer(hour00:hour20, names_to = "hour", values_to = "expression") %>% 
      dplyr::mutate(hour = as.numeric(str_sub(hour, 5, 6))) %>% 
      dplyr::select(hour, expression)
    
  } else if (scenario == "lin-1 ON mls-2 OFF") {
    mls2 <- data.frame(
      hour = seq(from = 0, to = 20, by = 2), 
      expression = c(rep(0, times = 4), rep(0, times = 4), rep(0, times = 3))
    )
  } else if (scenario == "lin-1 OFF mls-2 ON") {
    mls2 <- data.frame(
      hour = seq(from = 0, to = 20, by = 2), 
      expression = c(rep(0, times = 4), rep(1, times = 4), rep(0, times = 3))
    )
  }
  f_mls2 <- approxfun(mls2, rule = 2)
  return(f_mls2)
}
