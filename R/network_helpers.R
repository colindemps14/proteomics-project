# Adapted from weightedGCM_test_colin.R and weightedGCMtest_allie.R.
# Original summer work: Colin Dempsey and Allison Burnett.

make_wgcm_test <- function(conditional_test = NULL, regression = "gam",
                           weight_num = 7L, weight_method = "sign", nsim = 499L) {
  # Keep the original Pearson test when the conditioning set is empty.
  # Allow injection of a conditional test for the small verification script.
  if (is.null(conditional_test)) conditional_test <- weightedGCM::wgcm.fix
  function(x, y, S, suffStat) {
    dat <- suffStat$data
    if (!length(S)) {
      p_value <- stats::cor.test(dat[, x], dat[, y], method = "pearson")$p.value
    } else {
      p_value <- conditional_test(
        dat[, x, drop = FALSE], dat[, y, drop = FALSE],
        dat[, S, drop = FALSE], regr.meth = regression,
        weight.num = weight_num, weight.meth = weight_method, nsim = nsim
      )
    }
    # Stop rather than interpreting a failed test as independence.
    if (!is.numeric(p_value) || length(p_value) != 1L ||
        !is.finite(p_value) || p_value < 0 || p_value > 1) {
      stop("Invalid p-value for proteins ", colnames(dat)[x], " and ",
           colnames(dat)[y], "; conditioning set: ",
           paste(colnames(dat)[S], collapse = ", "))
    }
    as.numeric(p_value)
  }
}

undirected_adjacency <- function(fit, proteins) {
  M <- as(fit@graph, "matrix")
  M <- M[proteins, proteins, drop = FALSE]
  A <- 1L * ((M != 0) | (t(M) != 0))
  diag(A) <- 0L
  stopifnot(identical(A, t(A)), all(A %in% c(0L, 1L)))
  A
}

network_metrics <- function(A) {
  degree <- rowSums(A)
  edges <- sum(A[upper.tri(A)])
  stopifnot(sum(degree) == 2 * edges)
  data.frame(Edges = edges, Density = edges / choose(nrow(A), 2),
             Connected_proteins = sum(degree > 0),
             Average_degree = mean(degree), Maximum_degree = max(degree))
}

jaccard_matrix <- function(networks) {
  stopifnot(length(networks) > 0L)
  reference <- dimnames(networks[[1]])
  stopifnot(all(vapply(networks, function(A) identical(dimnames(A), reference), logical(1))))
  vectors <- lapply(networks, function(A) as.logical(A[upper.tri(A)]))
  J <- matrix(NA_real_, length(networks), length(networks),
              dimnames = list(names(networks), names(networks)))
  for (i in seq_along(networks)) for (j in seq_along(networks)) {
    union_size <- sum(vectors[[i]] | vectors[[j]])
    if (union_size) J[i, j] <- sum(vectors[[i]] & vectors[[j]]) / union_size
  }
  J
}
