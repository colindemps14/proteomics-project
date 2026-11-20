library(pcalg) ## install pcalg and related dependencies
library(Rgraphviz)
library(readxl)
library(tidyverse)
library(pheatmap)

protein_data <- read_excel("Wang_GGM_protein_data.xlsx")
protein_numeric <- protein_data[,-1]

head(protein_numeric)
n <- nrow(protein_numeric)
p <- ncol(protein_numeric)

alphas <- c(0.001, 0.01, 0.05, 0.10, 0.5)

pc_results <- list()

for(a in alphas){
  
  suffStat <- list(
    C = cor(protein_numeric, use = "pairwise.complete.obs"),
    n = n
  )
  
  pc_results[[as.character(a)]] <-
    pc(
      suffStat = suffStat,
      indepTest = gaussCItest,
      alpha = a,
      p = p,
      verbose = FALSE
    )
}

adj_mats <- list()

for(a in names(pc_results)){
  mat <- as(pc_results[[a]]@graph, "matrix")
  adj_mats[[a]] <- mat
}

for(a in names(adj_mats)) {
  
  # Get adjacency matrix
  M <- adj_mats[[a]]
  
  # Make undirected (edge exists in either direction)
  M.sym <- pmax(M, t(M))
  diag(M.sym) <- 0
  
  # Save heatmap
  png(
    filename = paste0("heatmap_alpha_clustered", a, ".png"),
    width = 1200,
    height = 1200,
    res = 150
  )
  
  pheatmap(
    M.sym,
    color = c("white", "#0072B2"),  # colorblind-friendly blue
    breaks = c(-0.5, 0.5, 1.5),
    show_rownames = FALSE,
    show_colnames = FALSE,
    border_color = "grey90",
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    legend_breaks = c(0, 1),
    legend_labels = c("No Edge", "Edge"),
    main = paste("Undirected PC Network (alpha =", a, ")")
  )
  
}


