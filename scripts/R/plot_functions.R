library(dplyr)
library(tidyr)
library(ggrepel)
library(forcats)

chart_dynamics <- function(ode.out, display_dat) {
  ode.out.df <- ode.out %>% 
    as.data.frame() %>% 
    pivot_longer(-time, names_to = "what", values_to = "concentration")
  
  state_panel <- data.frame(
    what = c("MLS2", "LIN1", "HLH1", "FOS1", "CYD1", "CYE1", "CKI1", "LIN35", "E2F", "MEF2", "Proliferation", "Differentiation"), 
    panel = c("panel1", "panel1", "panel2", "panel2", "panel3", "panel3", "panel3", "panel4", "panel5", "panel6", "panel7", "panel7")
  )
  
  state_colors <- c("MLS2" = "black", 
                    "LIN1" = "blue", 
                    "HLH1" = "black", 
                    "FOS1" = "orange", 
                    "CYD1" = "darkgreen", 
                    "CYE1" = "blue", 
                    "CKI1" = "red", 
                    "LIN35" = 'black', 
                    "E2F" = "black", 
                    "MEF2" = "black", 
                    "Proliferation" = "black", 
                    "Differentiation" = "blue")
  
  state_linetypes <- c("MLS2" = "solid", 
                   "LIN1" = "solid", 
                   "HLH1" = "solid", 
                   "FOS1" = "solid", 
                   "CYD1" = "solid", 
                   "CYE1" = "solid", 
                   "CKI1" = "solid", 
                   "LIN35" = "solid", 
                   "E2F" = "solid", 
                   "MEF2" = "solid", 
                   "Proliferation" = "solid", 
                   "Differentiation" = "dashed")
  
  state_peak <- ode.out.df %>% 
    group_by(what) %>% 
    summarize(concentration = max(concentration)) %>% 
    ungroup() %>% 
    left_join(ode.out.df, by = c("what", "concentration")) %>% 
    left_join(state_panel, by = "what") %>% 
    group_by(what, panel) %>% 
    summarize(concentration = concentration, 
              time = min(time)) %>% 
    ungroup() %>% 
    distinct()
  
  p <- ode.out.df %>% 
    left_join(state_panel, by = "what") %>% 
    ggplot(aes(time, concentration)) + 
    geom_line(aes(group = what, color = what, linetype = what), size = 1) + 
    scale_color_manual("", values = state_colors) + 
    scale_linetype_manual("", values = state_linetypes) + 
    geom_label_repel(data = state_peak, aes(label = what, color = what)) + 
    labs(x = "Hour after hatching", y = "Normalized expression level") + 
    facet_wrap(~panel, ncol = 1) + 
    theme_bw() + 
    theme(strip.text = element_blank(), 
          legend.position = "none", 
          axis.title = element_text(size = 14))
  
  if (display_dat == FALSE) {
    p
  } else {
    model_genes_dat <- readr::read_csv("model_genes_relative_expression.csv") %>% 
      dplyr::mutate(external_gene_id = stringr::str_to_upper(external_gene_id), 
                    external_gene_id = stringr::str_replace(external_gene_id, "-", ""), 
                    what = ifelse(external_gene_id == "RNR1", "Proliferation", external_gene_id), 
                    what = ifelse(external_gene_id == "UNC15", "Differentiation", what), 
                    what = ifelse(external_gene_id == "EFL1", "E2F", what)) %>% 
      dplyr::select(-external_gene_id) %>% 
      tidyr::pivot_longer(hour00:hour20, names_to = "time", values_to = "relative_expression") %>% 
      dplyr::mutate(time = stringr::str_sub(time, 5, 6), 
                    time = as.numeric(time)) %>% 
      dplyr::left_join(state_panel, by = "what")
    
    p + geom_point(data = model_genes_dat, aes(time, relative_expression, color = what), shape = 4, size = 3)
  }
}
