library(weightedGCM)
library(pcalg)
library(readxl)

set.seed(15)

####################################
# Conditional independence test
####################################

wgcmCItest <- function(x, y, S, suffStat){

  dat <- suffStat$data

  X <- dat[, x]
  Y <- dat[, y]

  # Unconditional independence
  if(length(S) == 0){
    return(cor.test(X, Y)$p.value)
  }

  Z <- dat[, S, drop = FALSE]

  wgcm.fix(
    X, Y, Z,
    regr.meth = "gam",
    weight.num = 7,
    weight.meth = "sign"
  )
}

####################################
# Read data
####################################

protein_data <- read_excel("Wang_GGM_protein_data.xlsx")
protein_numeric <- protein_data[, -1]
protein_mat <- as.matrix(protein_numeric)

p <- ncol(protein_mat)

suffStat <- list(
  data = protein_mat
)

####################################
# Alpha values
####################################

alphas <- c(0.001, 0.01, 0.05, 0.10, 0.5)

####################################
# Store degree tables
####################################

node_stats <- list()

####################################
# Run PC algorithm for each alpha
####################################

for(a in alphas){

  cat("\nRunning PC algorithm for alpha =", a, "\n")

  pc.est <- pc(
    suffStat = suffStat,
    indepTest = wgcmCItest,
    alpha = a,
    p = p,
    verbose = FALSE
  )

  ####################################
  # Create undirected adjacency matrix
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

  node_stats[[as.character(a)]] <- degree_df

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
    col = c("white", "#0072B2"),
    Rowv = NA,
    Colv = NA,
    labRow = TRUE,
    labCol = TRUE,
    main = paste("Undirected PC Network (alpha =", a, ")")
  )

  legend(
    "bottomleft",
    legend = c("No Edge", "Edge"),
    fill = c("white", "#0072B2"),
    border = "black",
    bty = "o",
    cex = 0.8
  )

}

####################################
# Save degree tables
####################################

for(a in names(node_stats)){

  write.csv(
    node_stats[[a]],
    file = paste0("node_degrees_alpha_", a, ".csv"),
    row.names = FALSE
  )

}