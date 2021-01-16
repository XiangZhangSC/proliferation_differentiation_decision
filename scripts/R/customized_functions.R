# Mean difference plot for two replicates in a group
# x = (log2(A) + log2(B)) / 2
# y = log2(A) - log2(B)
# endogenous and spike-in genes will be highlighted separately
plot_MD <- function(df) {
  ggplot(dplyr::filter(df, what_rna == "Endogenous"), aes(M, D, color = what_rna)) + 
    geom_point(alpha = 0.3) + 
    geom_point(data = dplyr::filter(df, what_rna == "Spike-in")) + 
    scale_color_manual("", values = c("steelblue", "red")) + 
    geom_hline(yintercept = 0, linetype = "dashed") + 
    labs(x = "Average log2-count", y = "Log2-fold-change") + 
    theme_bw() + 
    theme(legend.position = c(0.9,0.1), 
          legend.background = element_rect(fill = "transparent"))
}

# PCA matrix scatter plot
# PC1, PC2 and PC3
sketch_PCA <- function(cnt.mat) {
  Y <- log(cnt.mat + 1)
  # for prcomp, rows must be samples and columns are genes
  res.pca <- prcomp(t(Y), center = TRUE)
  
  # principal components
  res.pca.df <- res.pca$x %>% 
    as.data.frame() %>% 
    rownames_to_column("library_id") %>% 
    as_tibble() %>% 
    separate(library_id, into = c("genotype", "hour", "sample_id"), remove = FALSE)
  
  plotList <- list()
  
  plotList[[1]] <- ggally_text("PC1", color = "black") + theme_void()
  
  plotList[[2]] <- ggplot(res.pca.df, aes(PC2, PC1)) + 
    geom_point(aes(color = hour, shape = genotype), size = 3) + 
    scale_color_brewer("", palette = "Dark2")
  
  plotList[[3]] <- ggplot(res.pca.df, aes(PC3, PC1)) + 
    geom_point(aes(color = hour, shape = genotype), size = 3) + 
    scale_color_brewer("", palette = "Dark2")
  
  plotList[[4]] <- ggplot(res.pca.df, aes(PC1, PC2)) + 
    geom_point(aes(color = hour, shape = genotype), size = 3) + 
    scale_color_brewer("", palette = "Dark2")
  
  plotList[[5]] <- ggally_text("PC2", color = "black") + theme_void()
  
  plotList[[6]] <- ggplot(res.pca.df, aes(PC3, PC2)) + 
    geom_point(aes(color = hour, shape = genotype), size = 3) + 
    scale_color_brewer("", palette = "Dark2")
  
  plotList[[7]] <- ggplot(res.pca.df, aes(PC1, PC3)) + 
    geom_point(aes(color = hour, shape = genotype), size = 3) + 
    scale_color_brewer("", palette = "Dark2")
  
  plotList[[8]] <- ggplot(res.pca.df, aes(PC2, PC3)) + 
    geom_point(aes(color = hour, shape = genotype), size = 3) + 
    scale_color_brewer("", palette = "Dark2")
  
  plotList[[9]] <- ggally_text("PC3", color = "black") + theme_void()
  
  ggmatrix(plotList, nrow = 3, ncol = 3, byrow = TRUE, legend = 2) + theme(legend.position = "bottom")
}

# tidy the outcome of edgeR
# convert gene sequence name into gene name by using the Wormbase
tidy_togTags <- function(lrt, wormbase.df) {
  topTags(lrt, n = Inf, adjust.method = "fdr") %>% 
    as.data.frame() %>% 
    rownames_to_column("wormbase_gseq") %>% 
    as_tibble() %>% 
    left_join(annotation.df, by = "wormbase_gseq")
}
