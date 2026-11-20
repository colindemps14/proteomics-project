library(pcalg)
library(readxl)

set.seed(15)

####################################
# Read data
####################################

protein_data <- read_excel("Wang_GGM_protein_data.xlsx")
protein_numeric <- protein_data[, -1]
protein_mat <- as.matrix(protein_numeric)

n <- nrow(protein_mat)
p <- ncol(protein_mat)

####################################
# Sufficient statistics for Gaussian PC
####################################

suffStat <- list(
  C = cor(protein_mat),
  n = n
)

####################################
# Alpha values
####################################

alphas <- c(0.01, 0.05, 0.10)

node_stats_gaussian <- list()

####################################
# Run Gaussian PC algorithm
####################################

for(a in alphas){
  
  cat("\nRunning Gaussian PC algorithm for alpha =", a, "\n")
  
  pc.est <- pc(
    suffStat = suffStat,
    indepTest = gaussCItest,
    alpha = a,
    p = p,
    verbose = FALSE
  )
  
  ####################################
  # Undirected adjacency matrix
  ####################################
  
  M <- as(pc.est@graph, "matrix")
  M.sym <- pmax(M, t(M))
  
  ####################################
  # Degree calculations
  ####################################
  
  degrees <- rowSums(M.sym)
  
  degree_df <- data.frame(
    Protein = colnames(M.sym),
    Degree = degrees
  )
  
  degree_df <- degree_df[order(-degree_df$Degree), ]
  
  node_stats_gaussian[[as.character(a)]] <- degree_df
  
  ####################################
  # Print summary
  ####################################
  
  cat("\n========================\n")
  cat("Alpha =", a, "\n")
  cat("========================\n")
  cat("Average degree:", mean(degrees), "\n")
  cat("SD degree:", sd(degrees), "\n")
  cat("Maximum degree:", max(degrees), "\n\n")
  
  cat("Top 10 connected proteins:\n")
  print(head(degree_df, 10))
  
  ####################################
  # Heatmap
  ####################################
  
  dev.new()
  
  heatmap(
    M.sym,
    scale = "none",
    col = c("white", "#D55E00"),
    Rowv = NA,
    Colv = NA,
    labRow = TRUE,
    labCol = TRUE,
    main = paste("Gaussian PC Network (alpha =", a, ")")
  )
  
  legend(
    "bottomleft",
    legend = c("No Edge", "Edge"),
    fill = c("white", "#D55E00"),
    border = "black",
    bty = "o",
    cex = 0.8
  )
}

####################################
# Save degree tables
####################################

for(a in names(node_stats_gaussian)){
  
  write.csv(
    node_stats_gaussian[[a]],
    paste0("gaussian_node_degrees_alpha_", a, ".csv"),
    row.names = FALSE
  )
  
}