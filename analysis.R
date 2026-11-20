# Run from the project folder: source("analysis.R")
# Also called by index.qmd. Results are shared between both entry points.
source("R/network_helpers.R")

run_project <- function() {
  required <- c("readxl", "pcalg", "graph", "weightedGCM", "mgcv")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Run source('setup.R') first. Missing: ", paste(missing, collapse = ", "))

  path <- "data/Wang_GGM_protein_data.xlsx"
  raw <- readxl::read_excel(path, sheet = "ADProteinData_forGGM", .name_repair = "unique")
  ids <- as.character(raw[[1]])
  values <- as.data.frame(raw[-1])
  stopifnot(!anyNA(ids), !anyDuplicated(ids))
  stopifnot(all(vapply(values, is.numeric, logical(1))))
  original <- as.matrix(values)
  stopifnot(nrow(original) == 59L, ncol(original) == 192L)
  stopifnot(all(is.finite(original)), all(apply(original, 2, sd) > 0))
  X <- scale(original)

  config <- list(alphas = c(0.01, 0.05, 0.10), seed = 15L,
                 max_conditioning = nrow(X) - 4L,
                 regression = "gam", weight_num = 7L,
                 weight_method = "sign", nsim = 499L,
                 unconditional = "Pearson", standardize = TRUE)
  versions <- vapply(c(required, "stats"), function(p) as.character(utils::packageVersion(p)), character(1))
  # Invalidate saved fits if inputs, analysis code, settings, or packages change.
  signature <- list(data_md5 = unname(tools::md5sum(path)),
                    code_md5 = unname(tools::md5sum(c("analysis.R", "R/network_helpers.R"))),
                    config = config, packages = versions, R = R.version.string)
  settings <- expand.grid(Alpha = config$alphas,
                           Method = c("Fisher z", "Pearson + WGCM"),
                           stringsAsFactors = FALSE)
  settings$ID <- paste0(ifelse(settings$Method == "Fisher z", "fisher", "pearson_wgcm"),
                         "_", sprintf("%.2f", settings$Alpha))
  settings$Setting <- paste(settings$Method, sprintf("alpha = %.2f", settings$Alpha), sep = " | ")
  settings$Short <- paste(ifelse(settings$Method == "Fisher z", "F", "W"),
                          sprintf("%.2f", settings$Alpha))
  dir.create("results/checkpoints", recursive = TRUE, showWarnings = FALSE)
  dir.create("results/adjacency", recursive = TRUE, showWarnings = FALSE)
  fits <- networks <- setNames(vector("list", nrow(settings)), settings$Setting)
  elapsed <- numeric(nrow(settings))
  warning_tables <- vector("list", nrow(settings))
  wgcm_test <- make_wgcm_test(regression = config$regression,
                              weight_num = config$weight_num,
                              weight_method = config$weight_method, nsim = config$nsim)

  for (i in seq_len(nrow(settings))) {
    checkpoint <- file.path("results/checkpoints", paste0(settings$ID[i], ".rds"))
    request <- list(signature = signature, method = settings$Method[i], alpha = settings$Alpha[i])
    saved <- if (file.exists(checkpoint)) readRDS(checkpoint) else NULL
    if (!is.null(saved) && identical(saved$request, request)) {
      cat("Loading saved fit: ", settings$Setting[i], "\n", sep = "")
    } else {
      cat("Fitting: ", settings$Setting[i], "\n", sep = "")
      is_fisher <- settings$Method[i] == "Fisher z"
      suff <- if (is_fisher) list(C = cor(X), n = nrow(X)) else list(data = X)
      test <- if (is_fisher) pcalg::gaussCItest else wgcm_test
      notices <- character()
      # Reset the seed for each fit so fitting or skipping another setting
      # does not change this fit's simulation-based WGCM p-values.
      set.seed(config$seed)
      timing <- system.time({
        fit <- tryCatch(withCallingHandlers(
          pcalg::pc(suffStat = suff, indepTest = test, alpha = settings$Alpha[i],
                    labels = colnames(X), skel.method = "stable",
                    m.max = config$max_conditioning, NAdelete = FALSE,
                    verbose = FALSE),
          warning = function(w) {
            notices <<- c(notices, conditionMessage(w))
            invokeRestart("muffleWarning")
          }), error = function(e) {
            stop(settings$Setting[i], " failed: ", conditionMessage(e),
                 "\nCompleted fits are saved; send this error for diagnosis.", call. = FALSE)
          })
      })[["elapsed"]]
      saved <- list(request = request, fit = fit,
                    A = undirected_adjacency(fit, colnames(X)),
                    elapsed = timing, warnings = unique(notices))
      saveRDS(saved, checkpoint)
    }
    fits[[i]] <- saved$fit
    networks[[i]] <- saved$A
    elapsed[i] <- saved$elapsed
    if (length(saved$warnings)) {
      warning_tables[[i]] <- data.frame(Setting = settings$Setting[i], Warning = saved$warnings)
      warning(length(saved$warnings), " distinct warning(s) in ", settings$Setting[i],
              "; inspect results/model_warnings.csv before interpretation.", call. = FALSE)
    }
    write.csv(saved$A, file.path("results/adjacency", paste0(settings$ID[i], ".csv")))
  }

  summary_table <- cbind(settings, do.call(rbind, lapply(networks, network_metrics)),
                          Runtime_seconds = elapsed)
  rownames(summary_table) <- NULL
  degrees <- do.call(cbind, lapply(networks, rowSums))
  colnames(degrees) <- settings$Setting
  degree_table <- do.call(rbind, lapply(seq_len(nrow(settings)), function(i) {
    d <- degrees[, i]
    data.frame(Protein = rownames(degrees), Method = settings$Method[i],
               Alpha = settings$Alpha[i], Setting = settings$Setting[i],
               Degree = d, Rank = rank(-d, ties.method = "min"), row.names = NULL)
  }))
  degree_table <- degree_table[order(degree_table$Setting, -degree_table$Degree, degree_table$Protein), ]
  selected <- unique(degree_table$Protein[degree_table$Rank <= 5 & degree_table$Degree > 0])
  if (length(selected)) selected <- selected[order(-rowMeans(degrees[selected, , drop = FALSE]), selected)]
  J <- jaccard_matrix(networks)
  same_alpha <- do.call(rbind, lapply(config$alphas, function(a) {
    indices <- which(settings$Alpha == a)
    v1 <- networks[[indices[1]]][upper.tri(networks[[indices[1]]])]
    v2 <- networks[[indices[2]]][upper.tri(networks[[indices[2]]])]
    data.frame(Alpha = a, Shared_edges = sum(v1 == 1 & v2 == 1),
               Fisher_only = sum(v1 == 1 & v2 == 0),
               WGCM_pipeline_only = sum(v1 == 0 & v2 == 1),
               Jaccard = J[indices[1], indices[2]])
  }))
  model_warnings <- do.call(rbind, warning_tables)
  if (is.null(model_warnings)) model_warnings <- data.frame(Setting = character(), Warning = character())
  write.csv(summary_table, "results/network_summary.csv", row.names = FALSE)
  write.csv(degree_table, "results/protein_degrees.csv", row.names = FALSE)
  write.csv(J, "results/jaccard_similarity.csv")
  write.csv(same_alpha, "results/same_alpha_comparison.csv", row.names = FALSE)
  write.csv(model_warnings, "results/model_warnings.csv", row.names = FALSE)
  writeLines(capture.output(utils::sessionInfo()), "results/session_info.txt")
  result <- list(settings = settings, config = config, fits = fits, networks = networks,
                 summary = summary_table, degrees = degrees, degree_table = degree_table,
                 selected = selected, jaccard = J, same_alpha = same_alpha,
                 model_warnings = model_warnings,
                 data_summary = data.frame(Measure = c("Samples", "Proteins", "Measurements", "Missing"),
                                           Value = c(nrow(X), ncol(X), length(X), sum(is.na(X)))))
  saveRDS(result, "results/project_results.rds")
  cat("Done. All six networks and comparison tables are saved under results/.\n")
  result
}

project_results <- run_project()
