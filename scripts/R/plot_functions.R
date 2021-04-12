library(dplyr)
library(tidyr)
library(ggrepel)
library(forcats)

chart_dynamics <- function(ode.out, display_dat) {
  ode.out.df <- ode.out %>% 
    as.data.frame() %>% 
    pivot_longer(-time, names_to = "what", values_to = "concentration")
  
  state_panel <- data.frame(
    what = c("MLS2", "LIN1", "HLH1", "HLH1LIN35", "HLH1CYD1", "FOS1", "CYD1", "CYE1", "CKI1", "LIN35", "E2F", "E2FLIN35", "MEF2", "RNR1", "UNC15"), 
    panel = c("panel1", "panel1", "panel2", "panel2", "panel2", "panel3", "panel4", "panel4", "panel4", "panel5", "panel5", "panel5", "panel6", "panel7", "panel7")
  )
  
  state_colors <- c("MLS2" = "black", 
                    "LIN1" = "coral", 
                    "HLH1" = "black",
                    "HLH1LIN35" = "purple", 
                    "HLH1CYD1" = "deeppink", 
                    "FOS1" = "black", 
                    "CYD1" = "darkgreen", 
                    "CYE1" = "blue", 
                    "CKI1" = "red", 
                    "LIN35" = 'orchid', 
                    "E2F" = "limegreen", 
                    "E2FLIN35" = "firebrick1",
                    "MEF2" = "black", 
                    "RNR1" = "black", 
                    "UNC15" = "blue")
  
  state_linetypes <- c("MLS2" = "solid", 
                   "LIN1" = "solid", 
                   "HLH1" = "solid", 
                   "HLH1LIN35" = "solid", 
                   "HLH1CYD1" = "solid", 
                   "FOS1" = "solid", 
                   "CYD1" = "solid", 
                   "CYE1" = "solid", 
                   "CKI1" = "solid", 
                   "LIN35" = "dashed", 
                   "E2F" = "dashed", 
                   "E2FLIN35" = "solid", 
                   "MEF2" = "solid", 
                   "RNR1" = "solid", 
                   "UNC15" = "dashed")
  
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
          axis.title = element_text(size = 14)) + 
    scale_x_continuous(breaks = seq(from = 0, to = 20, by = 2))
  
  if (display_dat == FALSE) {
    p
  } else {
    model_genes_dat <- readr::read_csv("model_genes_relative_expression.csv") %>% 
      dplyr::mutate(external_gene_id = stringr::str_to_upper(external_gene_id), 
                    external_gene_id = stringr::str_replace(external_gene_id, "-", ""), 
                    what = ifelse(external_gene_id == "EFL1", "E2F", external_gene_id)) %>% 
      dplyr::select(-external_gene_id) %>% 
      tidyr::pivot_longer(hour00:hour20, names_to = "time", values_to = "relative_expression") %>% 
      dplyr::mutate(time = stringr::str_sub(time, 5, 6), 
                    time = as.numeric(time)) %>% 
      dplyr::left_join(state_panel, by = "what")
    
    p + geom_point(data = model_genes_dat, aes(time, relative_expression, color = what), shape = 4, size = 3)
  }
}
