
## Censored nearest-neighbor GP code for 1D and borehole 8D data.
##
## Main functions:
##   simulate_1d_censored_data()
##   simulate_borehole_censored_data()
##   fit_censored_gp_nn()
##   predict_censored_gp_nn_with_censoring()
##
## Composition Sampler: when censored neighbors are present, the interval is based
## on empirical predictive quantiles from the two-stage composition sampler.
##
## Default parameter order for d predictors and intercept-only mean:
##   c(tau_sq, ell_1, ..., ell_d, sigma_sq, beta0)

## Check that packages needed by the selected censoring method and optimizer are installed.
check_required_packages <- function(censor_method = "sov", use_parallel = TRUE) {
  method <- normalize_censor_method(censor_method)
  pkgs <- c("FNN", "mvtnorm")
  
  if (use_parallel) {
    pkgs <- c(pkgs, "optimParallel", "parallel")
  }
  
  if (method == "met") {
    pkgs <- c(pkgs, "TruncatedNormal")
  }
  
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs)) {
    stop("Please install missing package(s): ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
  }
  
  invisible(TRUE)
}

## Convert user-facing censoring method names and aliases to one canonical method name.
normalize_censor_method <- function(method) {
  method <- tolower(method[1])
  aliases <- c(
    truncatednormal = "met",
    tn_pmvnorm = "met",
    tn_mc = "met",
    minimax_tilting = "met",
    botev = "met",
    mvtnorm = "genzbretz",
    mcsov = "mc_sov",
    mc_sov = "mc_sov"
  )
  
  if (method %in% names(aliases)) {
    method <- aliases[[method]]
  }
  
  match.arg(method, c("sov", "mc_sov", "miwa", "genzbretz", "met"))
}

## True 1D signal used for the simulated 1D examples
default_1d_fun <- function(x) {
  sin(2 * (7 * x - 2)) +
    3 * exp(-15^2 * (x - 0.5)^2) +
    sin(2 * (2 * x - 0.5))
}

## Build a fixed 1D prediction grid and attach true/censored response columns
make_1d_test_grid <- function(n_test = 100L,
                              threshold_C = 1.221277,
                              f = default_1d_fun) {
  x1 <- seq(0, 1, length.out = n_test)
  y <- f(x1)
  y_censored <- pmin(y, threshold_C)
  
  data.frame(
    x1 = x1,
    y = y,
    ID = seq_len(n_test),
    y_censored = y_censored,
    censored = as.integer(y > threshold_C)
  )
}

## Simulate replicated noisy 1D data, apply right censoring, and split train/test data
simulate_1d_censored_data <- function(n = 800L,
                                      n_train = NULL,
                                      n_test = NULL,
                                      n_test_grid = 100L,
                                      test_frac = 0.3,
                                      seed = 2024L,
                                      seed_split = 2024L,
                                      seed_train = NULL,
                                      split_method = c("id", "row", "fixed_grid", "none"),
                                      fixed_test_grid = NULL,
                                      threshold_C = 1.221277,
                                      f = default_1d_fun,
                                      noise_sd = 0.2,
                                      reps_min = 1L,
                                      reps_max = 3L) {
  if (!is.null(n_train)) {
    n <- n_train
  }
  if (!is.null(n_test)) {
    n_test_grid <- n_test
  }
  if (!is.null(seed_train)) {
    seed <- seed_train
  }
  if (!is.null(fixed_test_grid)) {
    split_method <- if (isTRUE(fixed_test_grid)) "fixed_grid" else "id"
  }
  
  split_method <- match.arg(split_method)
  
  set.seed(seed)
  x <- runif(n, 0, 1)
  rep_counts <- sample(reps_min:reps_max, n, replace = TRUE)
  x_rep <- rep(x, times = rep_counts)
  f_true <- f(x_rep)
  y_raw <- f_true + rnorm(length(x_rep), mean = 0, sd = noise_sd)
  
  full_data <- data.frame(
    x1 = x_rep,
    y = y_raw,
    ID = rep(seq_len(n), rep_counts)
  )
  full_data$y_censored <- full_data$y
  full_data$y_censored[full_data$y > threshold_C] <- threshold_C
  full_data$censored <- as.integer(full_data$y > threshold_C)
  
  if (split_method == "id") {
    set.seed(seed_split)
    unique_ids <- unique(full_data$ID)
    n_test_ids <- floor(test_frac * length(unique_ids))
    test_ids <- sample(unique_ids, size = n_test_ids)
    test_data <- full_data[full_data$ID %in% test_ids, , drop = FALSE]
    train_data <- full_data[!full_data$ID %in% test_ids, , drop = FALSE]
  } else if (split_method == "row") {
    set.seed(seed_split)
    n_test_rows <- floor(test_frac * nrow(full_data))
    test_idx <- sample(seq_len(nrow(full_data)), size = n_test_rows)
    test_data <- full_data[test_idx, , drop = FALSE]
    train_data <- full_data[-test_idx, , drop = FALSE]
  } else if (split_method == "fixed_grid") {
    train_data <- full_data
    test_data <- make_1d_test_grid(
      n_test = n_test_grid,
      threshold_C = threshold_C,
      f = f
    )
  } else {
    train_data <- full_data
    test_data <- NULL
  }
  
  rownames(train_data) <- NULL
  if (!is.null(test_data)) {
    rownames(test_data) <- NULL
  }
  
  list(
    train_data = train_data,
    test_data = test_data,
    full_data = full_data,
    threshold_C = threshold_C,
    split_method = split_method
  )
}

## Wrap a function defined on physical ranges so it accepts inputs on the unit cube
unit_scale <- function(fun, ranges) {
  function(x) {
    xx <- x * (ranges[, 2] - ranges[, 1]) + ranges[, 1]
    fun(xx)
  }
}

## Evaluate the standard 8-dimensional borehole computer-model function
borehole <- function(xx) {
  rw <- xx[1]
  r <- xx[2]
  Tu <- xx[3]
  Hu <- xx[4]
  Tl <- xx[5]
  Hl <- xx[6]
  L <- xx[7]
  Kw <- xx[8]
  
  frac1 <- 2 * pi * Tu * (Hu - Hl)
  frac2a <- 2 * L * Tu / (log(r / rw) * rw^2 * Kw)
  frac2b <- Tu / Tl
  frac2 <- log(r / rw) * (1 + frac2a + frac2b)
  frac1 / frac2
}

## Return the physical input ranges used by the borehole example
default_borehole_ranges <- function() {
  matrix(c(
    0.05, 0.15,
    100, 50000,
    63070, 115600,
    990, 1110,
    63.1, 116,
    700, 820,
    1120, 1680,
    9855, 12045
  ), ncol = 2, byrow = TRUE)
}

## Expand base design rows into replicated observations and add Gaussian noise
expand_replicates_with_noise <- function(df,
                                         x_cols,
                                         y_true_col = "y_true",
                                         reps_min = 1L,
                                         reps_max = 3L,
                                         seed_reps = NULL,
                                         seed_noise = NULL,
                                         noise_sd = 1) {
  if (!is.null(seed_reps)) {
    set.seed(seed_reps)
  }
  
  reps <- sample(reps_min:reps_max, nrow(df), replace = TRUE)
  idx <- rep(seq_len(nrow(df)), reps)
  out <- df[idx, , drop = FALSE]
  rownames(out) <- NULL
  out$ID <- rep(seq_len(nrow(df)), reps)
  out$rep <- ave(out$ID, out$ID, FUN = seq_along)
  
  if (!is.null(seed_noise)) {
    set.seed(seed_noise)
  }
  out$y_true <- out[[y_true_col]]
  out$y_raw <- out$y_true + rnorm(nrow(out), mean = 0, sd = noise_sd)
  out[, c(x_cols, "ID", "rep", "y_true", "y_raw"), drop = FALSE]
}

## Apply a right-censoring threshold and create y, y_true, y_w, and censored columns
apply_right_censoring <- function(data,
                                  threshold_C,
                                  raw_y_col = "y_raw",
                                  true_y_col = "y_true") {
  raw_y <- data[[raw_y_col]]
  
  if (!is.null(true_y_col) && true_y_col %in% names(data)) {
    data$y_w <- data[[true_y_col]]
  }
  
  data$y_true <- raw_y
  data$y <- raw_y
  data$y[data$y > threshold_C] <- threshold_C
  data$censored <- as.integer(raw_y > threshold_C)
  
  if (raw_y_col %in% names(data) && raw_y_col != "y") {
    data[[raw_y_col]] <- NULL
  }
  
  data
}

## Simulate replicated noisy borehole train/test data with right-censoring
simulate_borehole_censored_data <- function(n_train = 2000L,
                                            n_test = 200L,
                                            d = 8L,
                                            seed_train = 1L,
                                            seed_test = 2L,
                                            seed_train_reps = 777L,
                                            seed_train_noise = 778L,
                                            seed_test_reps = 887L,
                                            seed_test_noise = 888L,
                                            reps_min = 1L,
                                            reps_max = 3L,
                                            noise_sd = 1,
                                            censor_quantile = 0.8,
                                            threshold_C = NULL,
                                            threshold_C_test = NULL,
                                            bore_fun = unit_scale(borehole, default_borehole_ranges())) {
  if (!requireNamespace("lhs", quietly = TRUE)) {
    stop("Please install package 'lhs' to simulate borehole data.", call. = FALSE)
  }
  
  threshold_was_null <- is.null(threshold_C)
  x_cols <- paste0("x", seq_len(d))
  
  set.seed(seed_train)
  X_train <- lhs::randomLHS(n_train, d)
  colnames(X_train) <- x_cols
  y_train_true <- apply(X_train, 1, bore_fun)
  
  set.seed(seed_test)
  X_test <- matrix(runif(n_test * d), nrow = n_test, ncol = d)
  colnames(X_test) <- x_cols
  y_test_true <- apply(X_test, 1, bore_fun)
  
  train_base <- as.data.frame(X_train)
  train_base$y_true <- y_train_true
  test_base <- as.data.frame(X_test)
  test_base$y_true <- y_test_true
  
  train_data <- expand_replicates_with_noise(
    train_base,
    x_cols = x_cols,
    seed_reps = seed_train_reps,
    seed_noise = seed_train_noise,
    reps_min = reps_min,
    reps_max = reps_max,
    noise_sd = noise_sd
  )
  test_data <- expand_replicates_with_noise(
    test_base,
    x_cols = x_cols,
    seed_reps = seed_test_reps,
    seed_noise = seed_test_noise,
    reps_min = reps_min,
    reps_max = reps_max,
    noise_sd = noise_sd
  )
  
  if (is.null(threshold_C)) {
    threshold_C <- as.numeric(stats::quantile(train_data$y_raw, probs = censor_quantile))
  }
  if (is.null(threshold_C_test)) {
    original_test_thresholds <- c(
      "0.95" = 164.9711,
      "0.9" = 143.7556,
      "0.8" = 115.8938,
      "0.7" = 99.69009
    )
    q_key <- as.character(censor_quantile)
    if (threshold_was_null && q_key %in% names(original_test_thresholds)) {
      threshold_C_test <- original_test_thresholds[[q_key]]
    } else {
      threshold_C_test <- threshold_C
    }
  }
  
  train_data <- apply_right_censoring(train_data, threshold_C = threshold_C)
  test_data <- apply_right_censoring(test_data, threshold_C = threshold_C_test)
  
  final_cols <- c(x_cols, "ID", "rep", "y_true", "y", "censored", "y_w")
  train_data <- train_data[, final_cols, drop = FALSE]
  test_data <- test_data[, final_cols, drop = FALSE]
  
  list(train_data = train_data, test_data = test_data, threshold_C = threshold_C, x_cols = x_cols)
}

## Standardize input data into the y, censored, and ID columns used by the fit
prepare_training_data <- function(train_data,
                                  x_cols,
                                  y_col = "y",
                                  censored_col = "censored",
                                  id_col = "ID",
                                  threshold_C = NULL) {
  data <- as.data.frame(train_data)
  
  missing_x <- setdiff(x_cols, names(data))
  if (length(missing_x)) {
    stop("Missing x column(s): ", paste(missing_x, collapse = ", "))
  }
  if (!y_col %in% names(data)) {
    stop("Missing y_col: ", y_col)
  }
  
  data$y <- as.numeric(data[[y_col]])
  
  if (!is.null(censored_col) && censored_col %in% names(data)) {
    data$censored <- as.integer(data[[censored_col]])
  } else {
    if (is.null(threshold_C)) {
      stop("threshold_C is required when censored_col is missing.")
    }
    data$censored <- as.integer(data$y > threshold_C)
    data$y <- pmin(data$y, threshold_C)
  }
  
  if (!is.null(id_col) && id_col %in% names(data)) {
    data$ID <- data[[id_col]]
  } else {
    data$ID <- seq_len(nrow(data))
  }
  
  data
}

## Order rows before the sequential nearest-neighbor likelihood is constructed
order_data <- function(data,
                       x_cols,
                       order_method = c("weighted_sum", "lexicographic", "none"),
                       order_values = NULL) {
  order_method <- match.arg(order_method)
  
  if (!is.null(order_values)) {
    data <- data[order(order_values), , drop = FALSE]
  } else if (order_method == "weighted_sum") {
    X <- as.matrix(data[, x_cols, drop = FALSE])
    denom <- rowSums(X)
    denom[abs(denom) < .Machine$double.eps] <- 1
    weighted_sum <- rowSums((X / denom) * X)
    data <- data[order(weighted_sum), , drop = FALSE]
  } else if (order_method == "lexicographic") {
    data <- data[do.call(order, as.data.frame(data[, x_cols, drop = FALSE])), , drop = FALSE]
  }
  
  rownames(data) <- NULL
  data
}

## Assign marker IDs to unique predictor locations within one observed/censored block
assign_marker_id_block <- function(data, x_cols, offset = 0L) {
  data$marker_id <- NA_integer_
  if (nrow(data) == 0L) {
    return(data)
  }
  
  complete <- complete.cases(data[, x_cols, drop = FALSE])
  if (any(complete)) {
    keys <- do.call(paste, c(data[complete, x_cols, drop = FALSE], sep = "_"))
    data$marker_id[complete] <- offset + as.integer(factor(keys, levels = unique(keys)))
  }
  
  data
}

## Add marker IDs after separating observed and censored rows into likelihood blocks
add_marker_id <- function(data, x_cols) {
  observed_data <- data[data$censored == 0L, , drop = FALSE]
  censored_data <- data[data$censored == 1L, , drop = FALSE]
  
  observed_data <- assign_marker_id_block(observed_data, x_cols, offset = 0L)
  max_observed <- suppressWarnings(max(observed_data$marker_id, na.rm = TRUE))
  if (!is.finite(max_observed)) {
    max_observed <- 0L
  }
  censored_data <- assign_marker_id_block(censored_data, x_cols, offset = max_observed)
  
  out <- rbind(observed_data, censored_data)
  rownames(out) <- NULL
  out
}

## For each marker ID, find previous uncensored nearest-neighbor marker IDs
find_nn_by_unique_id <- function(data, x_cols, k = 20L) {
  row_id <- seq_len(nrow(data))
  
  rep_unc_idx <- tapply(row_id, data$marker_id, function(ix) {
    unc <- ix[data$censored[ix] == 0L]
    if (length(unc)) unc[1] else NA_integer_
  })
  rep_unc_idx <- rep_unc_idx[!is.na(rep_unc_idx)]
  
  rows_by_id <- split(row_id, data$marker_id)
  out_list <- vector("list", length(rows_by_id))
  names(out_list) <- names(rows_by_id)
  
  for (i in seq_along(rows_by_id)) {
    this_id <- as.numeric(names(rows_by_id)[i])
    these_rows <- rows_by_id[[i]]
    first_row_i <- min(these_rows)
    
    cand_rep_idx <- rep_unc_idx[rep_unc_idx < first_row_i]
    
    if (!length(cand_rep_idx)) {
      out_list[[i]] <- data.frame(
        marker_id = this_id,
        neighbor_marker_ids = NA_character_,
        neighbor_indices = NA_character_,
        distances = NA_character_,
        stringsAsFactors = FALSE
      )
      next
    }
    
    k_use <- min(k, length(cand_rep_idx))
    nn <- FNN::get.knnx(
      data[cand_rep_idx, x_cols, drop = FALSE],
      data[first_row_i, x_cols, drop = FALSE],
      k = k_use
    )
    
    loc <- nn$nn.index[1, ]
    dist <- nn$nn.dist[1, ]
    picks <- cand_rep_idx[loc]
    pick_ids <- data$marker_id[picks]
    dist[pick_ids == this_id] <- 0
    
    if (any(dist == 0) && k_use < length(cand_rep_idx)) {
      k_use2 <- min(k + 1L, length(cand_rep_idx))
      nn2 <- FNN::get.knnx(
        data[cand_rep_idx, x_cols, drop = FALSE],
        data[first_row_i, x_cols, drop = FALSE],
        k = k_use2
      )
      loc <- nn2$nn.index[1, ]
      dist <- nn2$nn.dist[1, ]
      picks <- cand_rep_idx[loc]
      pick_ids <- data$marker_id[picks]
      dist[pick_ids == this_id] <- 0
    }
    
    all_pick_rows <- unlist(rows_by_id[as.character(pick_ids)], use.names = FALSE)
    all_pick_ids <- data$marker_id[all_pick_rows]
    rep_counts <- lengths(rows_by_id[as.character(pick_ids)])
    dist_rep <- rep(dist, rep_counts)
    
    out_list[[i]] <- data.frame(
      marker_id = this_id,
      neighbor_marker_ids = paste(all_pick_ids, collapse = ","),
      neighbor_indices = paste(all_pick_rows, collapse = ","),
      distances = paste(dist_rep, collapse = ","),
      stringsAsFactors = FALSE
    )
  }
  
  do.call(rbind, out_list)
}

## Expand marker-level nearest-neighbor information back to observation-level rows
get_nn_per_obs <- function(data, x_cols, k = 20L) {
  nn_unique <- find_nn_by_unique_id(data, x_cols = x_cols, k = k)
  rep_counts <- table(data$marker_id)
  nn_unique$n_rep_query <- as.integer(rep_counts[match(nn_unique$marker_id, names(rep_counts))])
  
  data$query_index <- seq_len(nrow(data))
  res <- merge(
    data[, c("query_index", "marker_id")],
    nn_unique,
    by.x = "marker_id",
    by.y = "marker_id",
    all.x = TRUE,
    sort = FALSE
  )
  
  names(res)[names(res) == "marker_id"] <- "query_marker_id"
  res <- res[, c("query_index", "query_marker_id", "n_rep_query",
                 "neighbor_indices", "neighbor_marker_ids", "distances")]
  res[order(res$query_index), , drop = FALSE]
}

## Keep one nearest-neighbor record per query marker ID
first_nn_row_by_marker <- function(nearest_neighbors) {
  nearest_neighbors[!duplicated(nearest_neighbors$query_marker_id), , drop = FALSE]
}

## Parse comma-separated neighbor ID/index strings into integer vectors
parse_id_string <- function(x) {
  if (length(x) == 0L || is.na(x[1]) || !nzchar(x[1])) {
    return(integer(0))
  }
  as.integer(strsplit(x[1], ",", fixed = TRUE)[[1]])
}

## Compute per-dimension squared-distance arrays for all pairs of rows in X
sqdist_array <- function(X) {
  X <- as.matrix(X)
  n <- nrow(X)
  d <- ncol(X)
  out <- array(0, dim = c(n, n, d))
  
  for (j in seq_len(d)) {
    out[, , j] <- outer(X[, j], X[, j], "-")^2
  }
  
  out
}

## Compute per-dimension squared-distance arrays between rows of two matrices
sqdist_cross_array <- function(X1, X2) {
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  d <- ncol(X1)
  out <- array(0, dim = c(nrow(X1), nrow(X2), d))
  for (j in seq_len(d)) {
    out[, , j] <- outer(X1[, j], X2[, j], "-")^2
  }
  out
}

## Convert squared distances into a squared-exponential GP covariance matrix
cov_from_sqdist <- function(D, tau_sq, ell) {
  scaled <- matrix(0, nrow = dim(D)[1], ncol = dim(D)[2])
  for (j in seq_along(ell)) {
    scaled <- scaled + D[, , j] / (2 * ell[j]^2)
  }
  tau_sq * exp(-scaled)
}

## Build the GP covariance matrix for a design matrix X
gp_cov_matrix <- function(X, tau_sq, ell) {
  cov_from_sqdist(sqdist_array(X), tau_sq = tau_sq, ell = ell)
}

## Precompute local squared-distance arrays needed repeatedly in the likelihood
compute_distance_cache <- function(nearest_neighbors, data, x_cols) {
  nn_by_marker <- first_nn_row_by_marker(nearest_neighbors)
  cache <- list()
  
  for (i in seq_len(nrow(nn_by_marker))) {
    query_marker_id <- nn_by_marker$query_marker_id[i]
    neighbor_ids <- parse_id_string(nn_by_marker$neighbor_marker_ids[i])
    
    if (!length(neighbor_ids)) {
      next
    }
    
    query_rows <- which(data$marker_id == query_marker_id)
    neighbor_rows <- which(data$marker_id %in% neighbor_ids)
    combined_rows <- c(query_rows, neighbor_rows)
    
    cache[[as.character(query_marker_id)]] <- list(
      D = sqdist_array(data[combined_rows, x_cols, drop = FALSE]),
      n_query = length(query_rows),
      neighbor_rows = neighbor_rows
    )
  }
  
  cache
}

## Symmetrize a covariance matrix and add jitter if needed to make it positive definite
make_pd <- function(Sigma, jitter = 1e-6) {
  Sigma <- as.matrix(Sigma)
  Sigma <- (Sigma + t(Sigma)) / 2
  eig <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
  min_eig <- min(eig)
  
  if (!is.finite(min_eig)) {
    stop("Covariance matrix has non-finite eigenvalues.")
  }
  
  if (min_eig <= 0) {
    Sigma <- Sigma + (abs(min_eig) + jitter) * diag(nrow(Sigma))
  }
  
  Sigma
}

## Construct the mean-model design matrix from a formula
mean_design <- function(data, mean_formula) {
  stats::model.matrix(mean_formula, data = as.data.frame(data))
}

## Create parameter names for variance, lengthscale, noise, and mean parameters
param_names <- function(x_cols, mean_names) {
  mean_names <- gsub("^\\(Intercept\\)$", "Intercept", mean_names)
  c("tau_sq", paste0("ell_", x_cols), "sigma_sq", paste0("beta_", make.names(mean_names)))
}

## Build default starting values for GP covariance and mean parameters
default_init_params <- function(data, x_cols, mean_formula) {
  obs <- data$censored == 0L
  if (!any(obs)) {
    obs <- rep(TRUE, nrow(data))
  }
  
  X_mean <- mean_design(data[obs, , drop = FALSE], mean_formula)
  y_obs <- data$y[obs]
  beta_init <- tryCatch(
    as.numeric(solve(crossprod(X_mean) + diag(1e-8, ncol(X_mean)), crossprod(X_mean, y_obs))),
    error = function(e) rep(mean(y_obs), ncol(X_mean))
  )
  
  y_var <- stats::var(y_obs)
  if (!is.finite(y_var) || y_var <= 0) {
    y_var <- 1
  }
  
  out <- c(y_var, rep(0.5, length(x_cols)), max(0.05 * y_var, 1e-4), beta_init)
  names(out) <- param_names(x_cols, colnames(X_mean))
  out
}

## Build default lower and upper bounds for the optimizer
default_bounds <- function(x_cols, mean_names) {
  nm <- param_names(x_cols, mean_names)
  lower <- c(1e-6, rep(1e-6, length(x_cols)), 1e-6, rep(-Inf, length(mean_names)))
  upper <- c(Inf, rep(Inf, length(x_cols)), Inf, rep(Inf, length(mean_names)))
  names(lower) <- nm
  names(upper) <- nm
  list(lower = lower, upper = upper)
}

## Validate user-supplied parameter vectors and align them to the expected order
normalize_params <- function(x, default, expected_names, label) {
  if (is.null(x)) {
    return(default)
  }
  
  x_names <- names(x)
  x <- as.numeric(x)
  names(x) <- x_names
  
  if (!is.null(names(x)) && all(expected_names %in% names(x))) {
    x <- x[expected_names]
  } else if (length(x) == length(expected_names)) {
    names(x) <- expected_names
  } else {
    stop(
      label, " must have length ", length(expected_names), ". Expected order: ",
      paste(expected_names, collapse = ", "),
      call. = FALSE
    )
  }
  
  x
}

## Split one optimizer parameter vector into GP and mean-model components
split_params <- function(params, x_cols, mean_names) {
  d <- length(x_cols)
  p <- length(mean_names)
  
  list(
    tau_sq = params[1],
    ell = as.numeric(params[seq.int(2, d + 1L)]),
    sigma_sq = params[d + 2L],
    beta = as.numeric(params[seq.int(d + 3L, d + 2L + p)])
  )
}

## Deterministic SOV integration for censored tail probabilities with block size 1 to 3
sov_censored_nd <- function(mu,
                            Sigma,
                            C,
                            rel.tol = 1e-4,
                            abs.tol = 0,
                            subdivisions = 1000L) {
  r <- length(mu)
  stopifnot(is.numeric(mu), is.matrix(Sigma), nrow(Sigma) == r, ncol(Sigma) == r, r %in% 1:3)
  
  U <- chol(Sigma)
  L <- t(U)
  qhi <- qnorm(1 - 1e-12)
  clip <- function(x) pmax(pmin(x, qhi), -qhi)
  
  if (r == 1L) {
    a1 <- clip((C - mu[1]) / L[1, 1])
    if (a1 >= qhi) {
      return(0)
    }
    return(integrate(
      dnorm,
      lower = a1,
      upper = qhi,
      rel.tol = rel.tol,
      abs.tol = abs.tol,
      subdivisions = subdivisions,
      stop.on.error = FALSE
    )$value)
  }
  
  if (r == 2L) {
    a1 <- clip((C - mu[1]) / L[1, 1])
    if (a1 >= qhi) {
      return(0)
    }
    a2 <- function(z1) (C - mu[2] - L[2, 1] * z1) / L[2, 2]
    integrand2 <- Vectorize(function(z1) {
      low2 <- clip(a2(z1))
      if (low2 >= qhi) {
        return(0)
      }
      dnorm(z1) * (1 - pnorm(low2))
    })
    return(integrate(
      integrand2,
      lower = a1,
      upper = qhi,
      rel.tol = rel.tol,
      abs.tol = abs.tol,
      subdivisions = subdivisions,
      stop.on.error = FALSE
    )$value)
  }
  
  a1 <- clip((C - mu[1]) / L[1, 1])
  if (a1 >= qhi) {
    return(0)
  }
  a2 <- function(z1) (C - mu[2] - L[2, 1] * z1) / L[2, 2]
  a3 <- function(z1, z2) (C - mu[3] - L[3, 1] * z1 - L[3, 2] * z2) / L[3, 3]
  
  inner_density <- function(z2, z1) {
    low3 <- clip(a3(z1, z2))
    dnorm(z2) * (1 - pnorm(low3))
  }
  
  outer_integrand <- Vectorize(function(z1) {
    low2 <- clip(a2(z1))
    if (low2 >= qhi) {
      return(0)
    }
    val2 <- integrate(
      function(z2) inner_density(z2, z1),
      lower = low2,
      upper = qhi,
      rel.tol = rel.tol,
      abs.tol = abs.tol,
      subdivisions = subdivisions,
      stop.on.error = FALSE
    )$value
    dnorm(z1) * val2
  })
  
  integrate(
    outer_integrand,
    lower = a1,
    upper = qhi,
    rel.tol = rel.tol,
    abs.tol = abs.tol,
    subdivisions = subdivisions,
    stop.on.error = FALSE
  )$value
}

## Monte Carlo SOV/GHK estimator for censored tail probabilities of any block size
mc_sov_censored_nd <- function(mu,
                               Sigma,
                               C,
                               M = 10000L,
                               seed = NULL,
                               eps = 1e-300,
                               return_se = FALSE) {
  mu <- as.numeric(mu)
  Sigma <- as.matrix(Sigma)
  r <- length(mu)
  stopifnot(is.numeric(mu), nrow(Sigma) == r, ncol(Sigma) == r, r >= 1L)
  C_vec <- as.numeric(C)
  if (length(C_vec) == 1L) {
    C_vec <- rep(C_vec, r)
  }
  if (length(C_vec) != r) {
    stop("C must have length 1 or length(mu).", call. = FALSE)
  }
  
  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    } else {
      NULL
    }
    on.exit({
      if (is.null(old_seed)) {
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }
  
  M <- as.integer(M)
  if (M <= 0L) {
    stop("M must be a positive integer.", call. = FALSE)
  }
  
  U <- chol(Sigma)
  L <- t(U)
  z <- matrix(0, nrow = M, ncol = r)
  log_w <- rep(0, M)
  active <- rep(TRUE, M)
  
  for (j in seq_len(r)) {
    if (j == 1L) {
      conditional_shift <- rep(0, M)
    } else {
      conditional_shift <- as.numeric(z[, seq_len(j - 1L), drop = FALSE] %*% L[j, seq_len(j - 1L)])
    }
    
    alpha <- (C_vec[j] - mu[j] - conditional_shift) / L[j, j]
    tail_prob <- stats::pnorm(alpha, lower.tail = FALSE)
    tail_prob[!is.finite(tail_prob)] <- 0
    tail_prob <- pmax(tail_prob, 0)
    
    active <- active & tail_prob > 0
    if (!any(active)) {
      if (return_se) {
        return(list(prob = 0, se = 0, weights = rep(0, M)))
      }
      return(0)
    }
    
    log_w[active] <- log_w[active] + log(pmax(tail_prob[active], eps))
    log_w[!active] <- -Inf
    
    u <- stats::runif(sum(active))
    survival_q <- pmax(u * tail_prob[active], eps)
    z[active, j] <- stats::qnorm(survival_q, lower.tail = FALSE)
  }
  
  max_log_w <- max(log_w)
  if (!is.finite(max_log_w)) {
    if (return_se) {
      return(list(prob = 0, se = 0, weights = rep(0, M)))
    }
    return(0)
  }
  
  weights_scaled <- exp(log_w - max_log_w)
  prob <- exp(max_log_w) * mean(weights_scaled)
  
  if (!return_se) {
    return(as.numeric(prob))
  }
  
  weights <- exp(log_w)
  se <- stats::sd(weights) / sqrt(M)
  list(prob = as.numeric(prob), se = as.numeric(se), weights = weights)
}

## Compute Pr(Y > C) for one censored block using the selected probability method
censored_tail_prob <- function(mu,
                               Sigma,
                               C,
                               method = c("sov", "mc_sov", "miwa", "genzbretz", "met"),
                               p_floor = 1e-300,
                               genzbretz_maxpts = 1e6,
                               genzbretz_abseps = 1e-6,
                               genzbretz_releps = 0,
                               sov_rel_tol = 1e-4,
                               sov_abs_tol = 0,
                               sov_subdivisions = 1000L,
                               mc_sov_M = 10000L,
                               mc_sov_seed = NULL,
                               tn_B = 2e5,
                               tn_type = "mc") {
  method <- normalize_censor_method(method)
  mu <- as.numeric(mu)
  Sigma <- make_pd(Sigma)
  r <- length(mu)
  
  p <- switch(
    method,
    sov = {
      if (r > 3L) {
        stop("SOV is implemented here for replicate block size 1, 2, or 3 only.")
      }
      sov_censored_nd(
        mu = mu,
        Sigma = Sigma,
        C = C,
        rel.tol = sov_rel_tol,
        abs.tol = sov_abs_tol,
        subdivisions = sov_subdivisions
      )
    },
    mc_sov = mc_sov_censored_nd(
      mu = mu,
      Sigma = Sigma,
      C = C,
      M = mc_sov_M,
      seed = mc_sov_seed
    ),
    miwa = mvtnorm::pmvnorm(
      lower = rep(C, r),
      upper = rep(Inf, r),
      mean = mu,
      sigma = Sigma,
      algorithm = mvtnorm::Miwa()
    )[1],
    genzbretz = mvtnorm::pmvnorm(
      lower = rep(C, r),
      upper = rep(Inf, r),
      mean = mu,
      sigma = Sigma,
      algorithm = mvtnorm::GenzBretz(
        maxpts = genzbretz_maxpts,
        abseps = genzbretz_abseps,
        releps = genzbretz_releps
      )
    )[1],
    #'    met = TruncatedNormal::pmvnorm(
    #      lb = rep(C, r),
    #      ub = rep(Inf, r),
    #      mu = mu,
    #      sigma = Sigma,
    #      type = tn_type,
    #      B = tn_B
    #    )[1]'
    met = TruncatedNormal::mvNcdf(
      l = rep(C, r) - as.numeric(mu),
      u = rep(Inf, r) - as.numeric(mu),
      Sig = Sigma#,
      #n = tn_B
    )$prob
  )
  
  p <- as.numeric(p)
  if (!is.finite(p)) {
    p <- 0
  }
  max(p, p_floor)
}

## Compute conditional mean and covariance for one query marker given its neighbors
conditional_params <- function(query_marker_id,
                               data,
                               nn_by_marker,
                               distance_cache,
                               params,
                               x_cols,
                               mean_formula,
                               mean_names) {
  par <- split_params(params, x_cols, mean_names)
  query_rows <- which(data$marker_id == query_marker_id)
  
  nn_row <- nn_by_marker[nn_by_marker$query_marker_id == query_marker_id, , drop = FALSE]
  neighbor_ids <- parse_id_string(nn_row$neighbor_marker_ids)
  
  x_query <- as.matrix(data[query_rows, x_cols, drop = FALSE])
  X_query <- mean_design(data[query_rows, , drop = FALSE], mean_formula)
  mean_query <- as.numeric(X_query %*% par$beta)
  cov_query <- gp_cov_matrix(x_query, par$tau_sq, par$ell) + par$sigma_sq * diag(nrow(x_query))
  
  if (!length(neighbor_ids)) {
    return(list(mean = mean_query, covariance = make_pd(cov_query)))
  }
  
  cache <- distance_cache[[as.character(query_marker_id)]]
  if (is.null(cache)) {
    return(list(mean = mean_query, covariance = make_pd(cov_query)))
  }
  
  neighbor_rows <- cache$neighbor_rows
  y_neighbors <- data$y[neighbor_rows]
  X_neighbors <- mean_design(data[neighbor_rows, , drop = FALSE], mean_formula)
  residual_neighbors <- y_neighbors - as.numeric(X_neighbors %*% par$beta)
  
  n_query <- cache$n_query
  n_total <- dim(cache$D)[1]
  neighbor_idx <- seq.int(n_query + 1L, n_total)
  
  D_qn <- cache$D[seq_len(n_query), neighbor_idx, , drop = FALSE]
  D_nn <- cache$D[neighbor_idx, neighbor_idx, , drop = FALSE]
  
  cov_qn <- cov_from_sqdist(D_qn, tau_sq = par$tau_sq, ell = par$ell)
  cov_nn <- cov_from_sqdist(D_nn, tau_sq = par$tau_sq, ell = par$ell) +
    par$sigma_sq * diag(length(neighbor_rows))
  cov_nn <- make_pd(cov_nn)
  
  inv_cov_nn <- solve(cov_nn)
  cond_mean <- as.numeric(mean_query + cov_qn %*% inv_cov_nn %*% residual_neighbors)
  cond_cov <- cov_query - cov_qn %*% inv_cov_nn %*% t(cov_qn)
  
  list(mean = cond_mean, covariance = make_pd(cond_cov))
}

## Evaluate the nearest-neighbor censored GP negative log-likelihood for optim()
neg_log_likelihood <- function(params,
                               data,
                               nearest_neighbors,
                               distance_cache,
                               threshold_C,
                               x_cols,
                               mean_formula,
                               mean_names,
                               censor_method = "sov",
                               censor_prob_args = list(),
                               p_floor = 1e-300,
                               custom_censor_prob = NULL) {
  par <- split_params(params, x_cols, mean_names)
  if (par$tau_sq <= 0 || par$sigma_sq <= 0 || any(par$ell <= 0)) {
    return(Inf)
  }
  
  nn_by_marker <- first_nn_row_by_marker(nearest_neighbors)
  unique_ids <- unique(nearest_neighbors$query_marker_id)
  last_observed_id <- tail(data$marker_id[data$censored == 0L], 1)
  end_pos <- match(last_observed_id, unique_ids)
  
  observed_ids <- if (is.na(end_pos)) numeric(0) else unique_ids[seq_len(end_pos)]
  censored_ids <- if (is.na(end_pos)) {
    unique_ids
  } else if (end_pos < length(unique_ids)) {
    unique_ids[seq.int(end_pos + 1L, length(unique_ids))]
  } else {
    numeric(0)
  }
  
  log_likelihood <- 0
  
  for (query_marker_id in observed_ids) {
    query_rows <- which(data$marker_id == query_marker_id)
    cond <- conditional_params(
      query_marker_id = query_marker_id,
      data = data,
      nn_by_marker = nn_by_marker,
      distance_cache = distance_cache,
      params = params,
      x_cols = x_cols,
      mean_formula = mean_formula,
      mean_names = mean_names
    )
    
    ll <- mvtnorm::dmvnorm(
      data$y[query_rows],
      mean = cond$mean,
      sigma = cond$covariance,
      log = TRUE
    )
    if (!is.finite(ll)) {
      return(Inf)
    }
    log_likelihood <- log_likelihood + ll
  }
  
  for (query_marker_id in censored_ids) {
    cond <- conditional_params(
      query_marker_id = query_marker_id,
      data = data,
      nn_by_marker = nn_by_marker,
      distance_cache = distance_cache,
      params = params,
      x_cols = x_cols,
      mean_formula = mean_formula,
      mean_names = mean_names
    )
    
    if (is.null(custom_censor_prob)) {
      prob_args <- c(
        list(
          mu = cond$mean,
          Sigma = cond$covariance,
          C = threshold_C,
          method = censor_method,
          p_floor = p_floor
        ),
        censor_prob_args
      )
      p_cens <- do.call(censored_tail_prob, prob_args)
    } else {
      prob_args <- c(
        list(mu = cond$mean, Sigma = cond$covariance, C = threshold_C),
        censor_prob_args
      )
      p_cens <- do.call(custom_censor_prob, prob_args)
      p_cens <- max(as.numeric(p_cens), p_floor)
    }
    
    log_likelihood <- log_likelihood + log(p_cens)
  }
  
  -log_likelihood
}

## Fit the censored nearest-neighbor GP by optimizing the negative log-likelihood
fit_censored_gp_nn <- function(train_data,
                               x_cols,
                               threshold_C,
                               y_col = "y",
                               censored_col = "censored",
                               id_col = "ID",
                               k = 20L,
                               init_params = NULL,
                               lower_bounds = NULL,
                               upper_bounds = NULL,
                               mean_formula = ~ 1,
                               order_method = c("weighted_sum", "lexicographic", "none"),
                               order_values = NULL,
                               censor_method = c("sov", "mc_sov", "miwa", "genzbretz", "met"),
                               censor_prob_args = list(),
                               custom_censor_prob = NULL,
                               p_floor = 1e-300,
                               use_parallel = TRUE,
                               n_cores = max(1L, parallel::detectCores() - 1L),
                               optim_control = list(),
                               return_prepared = TRUE,
                               drop_id_before_fit = TRUE) {
  censor_method <- normalize_censor_method(censor_method)
  check_required_packages(censor_method = censor_method, use_parallel = use_parallel)
  
  data <- prepare_training_data(
    train_data = train_data,
    x_cols = x_cols,
    y_col = y_col,
    censored_col = censored_col,
    id_col = id_col,
    threshold_C = threshold_C
  )
  if (isTRUE(drop_id_before_fit) && "ID" %in% names(data)) {
    data$ID <- NULL
  }
  data <- order_data(data, x_cols = x_cols, order_method = order_method, order_values = order_values)
  data <- add_marker_id(data, x_cols = x_cols)
  
  mean_template <- mean_design(data[seq_len(1L), , drop = FALSE], mean_formula)
  mean_names <- colnames(mean_template)
  expected_names <- param_names(x_cols, mean_names)
  defaults <- default_init_params(data, x_cols, mean_formula)
  bounds <- default_bounds(x_cols, mean_names)
  
  init_params <- normalize_params(init_params, defaults, expected_names, "init_params")
  lower_bounds <- normalize_params(lower_bounds, bounds$lower, expected_names, "lower_bounds")
  upper_bounds <- normalize_params(upper_bounds, bounds$upper, expected_names, "upper_bounds")
  
  nearest_neighbors <- get_nn_per_obs(data, x_cols = x_cols, k = k)
  distance_cache <- compute_distance_cache(nearest_neighbors, data = data, x_cols = x_cols)
  
  full_start_time <- Sys.time()
  ran_parallel <- FALSE
  
  if (use_parallel && n_cores > 1L) {
    cl <- NULL
    parallel_fit <- tryCatch(
      {
        cl <- parallel::makeCluster(n_cores)
        parallel::setDefaultCluster(cl)
        
        parallel::clusterExport(
          cl,
          c(
            "normalize_censor_method", "make_pd", "sqdist_array", "sqdist_cross_array",
            "cov_from_sqdist", "gp_cov_matrix", "mean_design", "split_params",
            "sov_censored_nd", "mc_sov_censored_nd", "censored_tail_prob", "conditional_params",
            "neg_log_likelihood", "first_nn_row_by_marker", "parse_id_string"
          ),
          envir = environment()
        )
        parallel::clusterEvalQ(cl, {
          library(mvtnorm)
          NULL
        })
        
        optim_time <- system.time({
          result <- optimParallel::optimParallel(
            par = init_params,
            fn = neg_log_likelihood,
            method = "L-BFGS-B",
            lower = lower_bounds,
            upper = upper_bounds,
            control = optim_control,
            parallel = list(cl = cl, forward = FALSE, loginfo = FALSE),
            data = data,
            nearest_neighbors = nearest_neighbors,
            distance_cache = distance_cache,
            threshold_C = threshold_C,
            x_cols = x_cols,
            mean_formula = mean_formula,
            mean_names = mean_names,
            censor_method = censor_method,
            censor_prob_args = censor_prob_args,
            p_floor = p_floor,
            custom_censor_prob = custom_censor_prob
          )
        })
        
        list(result = result, optim_time = optim_time)
      },
      error = function(e) {
        warning(
          "Parallel optimization failed; falling back to serial optim(): ",
          conditionMessage(e),
          call. = FALSE
        )
        NULL
      },
      finally = {
        if (!is.null(cl)) {
          try(parallel::stopCluster(cl), silent = TRUE)
        }
      }
    )
    
    if (!is.null(parallel_fit)) {
      result <- parallel_fit$result
      optim_time <- parallel_fit$optim_time
      ran_parallel <- TRUE
    }
  }
  
  if (!ran_parallel) {
    optim_time <- system.time({
      result <- stats::optim(
        par = init_params,
        fn = neg_log_likelihood,
        method = "L-BFGS-B",
        lower = lower_bounds,
        upper = upper_bounds,
        control = optim_control,
        data = data,
        nearest_neighbors = nearest_neighbors,
        distance_cache = distance_cache,
        threshold_C = threshold_C,
        x_cols = x_cols,
        mean_formula = mean_formula,
        mean_names = mean_names,
        censor_method = censor_method,
        censor_prob_args = censor_prob_args,
        p_floor = p_floor,
        custom_censor_prob = custom_censor_prob
      )
    })
  }
  
  full_end_time <- Sys.time()
  
  out <- list(
    result = result,
    estimated_params = result$par,
    final_log_likelihood = -result$value,
    x_cols = x_cols,
    threshold_C = threshold_C,
    k = k,
    mean_formula = mean_formula,
    mean_names = mean_names,
    censor_method = censor_method,
    init_params = init_params,
    lower_bounds = lower_bounds,
    upper_bounds = upper_bounds,
    used_parallel = ran_parallel,
    optim_time = optim_time,
    total_seconds = as.numeric(full_end_time - full_start_time, units = "secs"),
    train_data = train_data
  )
  
  if (return_prepared) {
    out$prepared_data <- data
    out$nearest_neighbors <- nearest_neighbors
    out$distance_cache <- distance_cache
  }
  
  class(out) <- "censored_gp_nn_fit"
  out
}

## Predict from the fitted GP without conditioning on censored-neighbor latent values
predict_censored_gp_nn <- function(fit,
                                   test_data,
                                   max_marker_neighbors = fit$k,
                                   interval_level = 0.95,
                                   include_noise = TRUE) {
  if (is.null(fit$prepared_data)) {
    stop("fit must be created with return_prepared = TRUE to predict.")
  }
  
  par <- split_params(fit$estimated_params, fit$x_cols, fit$mean_names)
  train_obs <- fit$prepared_data[fit$prepared_data$censored == 0L, , drop = FALSE]
  if (nrow(train_obs) == 0L) {
    stop("Prediction requires at least one uncensored training row.")
  }
  
  results <- data.frame(
    test_index = seq_len(nrow(test_data)),
    predicted_mean = NA_real_,
    predicted_variance = NA_real_
  )
  
  for (i in seq_len(nrow(test_data))) {
    x_star <- as.matrix(test_data[i, fit$x_cols, drop = FALSE])
    x_train <- as.matrix(train_obs[, fit$x_cols, drop = FALSE])
    dist <- rowSums((t(t(x_train) - as.numeric(x_star)))^2)
    
    group_dist <- tapply(dist, train_obs$marker_id, min)
    selected_groups <- as.integer(names(sort(group_dist))[seq_len(min(max_marker_neighbors, length(group_dist)))])
    neighbor_rows <- which(train_obs$marker_id %in% selected_groups)
    
    x_neighbors <- as.matrix(train_obs[neighbor_rows, fit$x_cols, drop = FALSE])
    y_neighbors <- train_obs$y[neighbor_rows]
    
    cov_nn <- gp_cov_matrix(x_neighbors, par$tau_sq, par$ell) +
      par$sigma_sq * diag(nrow(x_neighbors))
    cov_nn <- make_pd(cov_nn)
    
    D_star_neighbors <- sqdist_cross_array(x_star, x_neighbors)
    cov_star_neighbors <- cov_from_sqdist(D_star_neighbors, tau_sq = par$tau_sq, ell = par$ell)
    
    X_star <- mean_design(test_data[i, , drop = FALSE], fit$mean_formula)
    X_neighbors <- mean_design(train_obs[neighbor_rows, , drop = FALSE], fit$mean_formula)
    residual_neighbors <- y_neighbors - as.numeric(X_neighbors %*% par$beta)
    
    inv_cov_nn <- solve(cov_nn)
    pred_mean <- as.numeric(X_star %*% par$beta + cov_star_neighbors %*% inv_cov_nn %*% residual_neighbors)
    
    marginal_var <- par$tau_sq
    if (include_noise) {
      marginal_var <- marginal_var + par$sigma_sq
    }
    pred_var <- as.numeric(marginal_var - cov_star_neighbors %*% inv_cov_nn %*% t(cov_star_neighbors))
    
    results$predicted_mean[i] <- pred_mean
    results$predicted_variance[i] <- max(pred_var, 0)
  }
  
  alpha <- 1 - interval_level
  sd_pred <- sqrt(results$predicted_variance)
  results$CI_lower <- stats::qnorm(alpha / 2, mean = results$predicted_mean, sd = sd_pred)
  results$CI_upper <- stats::qnorm(1 - alpha / 2, mean = results$predicted_mean, sd = sd_pred)
  results$interval_length <- results$CI_upper - results$CI_lower
  results
}

## Make prediction covariance matrices numerically positive definite
make_posdef_prediction <- function(mat, jitter = 1e-8) {
  mat <- (mat + t(mat)) / 2
  eigvals <- eigen(mat, symmetric = TRUE, only.values = TRUE)$values
  if (any(eigvals <= 0)) {
    mat <- mat + diag(abs(min(eigvals)) + jitter, nrow(mat))
  } else {
    mat <- mat + diag(jitter, nrow(mat))
  }
  mat
}

## Solve linear systems by Cholesky factorization with prediction-time jitter.
chol_solve_prediction <- function(S, b = diag(nrow(S)), jitter = 1e-8) {
  S <- make_posdef_prediction(S, jitter = jitter)
  R <- chol(S)
  backsolve(R, forwardsolve(t(R), b))
}

## Row-bind two data frames while filling columns that appear in only one input.
bind_rows_fill <- function(a, b) {
  a <- as.data.frame(a)
  b <- as.data.frame(b)
  all_names <- union(names(a), names(b))
  
  for (nm in setdiff(all_names, names(a))) {
    a[[nm]] <- NA
  }
  for (nm in setdiff(all_names, names(b))) {
    b[[nm]] <- NA
  }
  
  rbind(a[, all_names, drop = FALSE], b[, all_names, drop = FALSE])
}

## Compute moments of a multivariate normal truncated to a rectangular region.
##
## For more than one censored response, method = "auto" first uses
## tmvtnorm::mtmvnorm(). If that is unavailable or fails, it estimates the
## moments from exact TruncatedNormal::rtmvnorm() draws.
truncated_mvn_moments <- function(mu,
                                  Sigma,
                                  lower,
                                  upper,
                                  method = c("auto", "tmvtnorm", "TruncatedNormal", "mc", "univariate"),
                                  mc_samples = 20000L,
                                  mc_burn_in = 2000L,
                                  mc_thinning = 5L) {
  method <- match.arg(method)
  
  mu <- as.numeric(mu)
  Sigma <- make_posdef_prediction(as.matrix(Sigma))
  p <- length(mu)
  
  expand_bound <- function(x, name) {
    x <- as.numeric(x)
    if (length(x) == 1L) {
      return(rep(x, p))
    }
    if (length(x) != p) {
      stop(sprintf("%s must have length 1 or length(mu).", name))
    }
    x
  }
  
  lower <- expand_bound(lower, "lower")
  upper <- expand_bound(upper, "upper")
  
  if (!all(dim(Sigma) == c(p, p))) {
    stop("Sigma must be a square matrix with dimension length(mu).")
  }
  if (any(lower >= upper)) {
    stop("Every lower truncation bound must be smaller than its upper bound.")
  }
  
  ## Exact analytical moments for the one-sided univariate case used here.
  if (p == 1L &&
      method %in% c("auto", "univariate") &&
      is.infinite(upper[1L]) && upper[1L] > 0) {
    sd_value <- sqrt(max(Sigma[1L, 1L], 1e-12))
    alpha <- (lower[1L] - mu[1L]) / sd_value
    
    log_tail <- stats::pnorm(alpha, lower.tail = FALSE, log.p = TRUE)
    lambda <- exp(stats::dnorm(alpha, log = TRUE) - log_tail)
    
    truncated_mean <- mu[1L] + sd_value * lambda
    truncated_var <- sd_value^2 * (1 + alpha * lambda - lambda^2)
    truncated_var <- max(truncated_var, 1e-10)
    
    return(list(
      mean = truncated_mean,
      var = matrix(truncated_var, nrow = 1L, ncol = 1L),
      method_used = "univariate exact"
    ))
  }
  
  ## Deterministic Tallis-moment calculation.
  if (method %in% c("auto", "tmvtnorm") &&
      requireNamespace("tmvtnorm", quietly = TRUE)) {
    mt <- tryCatch(
      tmvtnorm::mtmvnorm(
        mean = mu,
        sigma = Sigma,
        lower = lower,
        upper = upper
      ),
      error = function(e) NULL
    )
    
    if (!is.null(mt)) {
      mt_mean <- as.numeric(mt$tmean)
      mt_var <- as.matrix(mt$tvar)
      if (p == 1L) {
        mt_var <- matrix(mt_var, nrow = 1L, ncol = 1L)
      }
      if (all(is.finite(mt_mean)) && all(is.finite(mt_var))) {
        return(list(
          mean = mt_mean,
          var = mt_var,
          method_used = "tmvtnorm::mtmvnorm"
        ))
      }
    }
    
    if (method == "tmvtnorm") {
      stop("tmvtnorm::mtmvnorm() did not return finite truncated moments.")
    }
  }
  
  ## Monte Carlo moment calculation using exact minimax-tilting draws.
  if (method %in% c("auto", "TruncatedNormal", "mc") &&
      requireNamespace("TruncatedNormal", quietly = TRUE)) {
    samples <- tryCatch(
      TruncatedNormal::rtmvnorm(
        n = as.integer(mc_samples),
        mu = mu,
        sigma = Sigma,
        lb = lower,
        ub = upper
      ),
      error = function(e) NULL
    )
    
    if (!is.null(samples)) {
      samples <- as.matrix(samples)
      if (ncol(samples) != p && nrow(samples) == p) {
        samples <- t(samples)
      }
      if (p == 1L && ncol(samples) != 1L) {
        samples <- matrix(as.numeric(samples), ncol = 1L)
      }
      
      mt_mean <- colMeans(samples)
      mt_var <- stats::cov(samples)
      if (p == 1L) {
        mt_var <- matrix(mt_var, nrow = 1L, ncol = 1L)
      }
      
      if (all(is.finite(mt_mean)) && all(is.finite(mt_var))) {
        return(list(
          mean = as.numeric(mt_mean),
          var = as.matrix(mt_var),
          method_used = "TruncatedNormal Monte Carlo"
        ))
      }
    }
    
    if (method == "TruncatedNormal") {
      stop("TruncatedNormal::rtmvnorm() did not return usable draws.")
    }
  }
  
  ## Last-resort MCMC moment calculation.
  if (method %in% c("auto", "mc") &&
      requireNamespace("tmvtnorm", quietly = TRUE)) {
    samples <- tryCatch(
      tmvtnorm::rtmvnorm(
        n = as.integer(mc_samples),
        mean = mu,
        sigma = Sigma,
        lower = lower,
        upper = upper,
        algorithm = "gibbs",
        burn.in.samples = as.integer(mc_burn_in),
        thinning = as.integer(mc_thinning)
      ),
      error = function(e) NULL
    )
    
    if (!is.null(samples)) {
      samples <- as.matrix(samples)
      if (ncol(samples) != p && nrow(samples) == p) {
        samples <- t(samples)
      }
      if (p == 1L && ncol(samples) != 1L) {
        samples <- matrix(as.numeric(samples), ncol = 1L)
      }
      
      mt_mean <- colMeans(samples)
      mt_var <- stats::cov(samples)
      if (p == 1L) {
        mt_var <- matrix(mt_var, nrow = 1L, ncol = 1L)
      }
      
      if (all(is.finite(mt_mean)) && all(is.finite(mt_var))) {
        return(list(
          mean = as.numeric(mt_mean),
          var = as.matrix(mt_var),
          method_used = "tmvtnorm Gibbs Monte Carlo"
        ))
      }
    }
  }
  
  if (method == "univariate" && p > 1L) {
    stop("method = 'univariate' is valid only for one censored response.")
  }
  
  stop(
    paste0(
      "Could not compute finite truncated-normal moments. ",
      "Install tmvtnorm and/or TruncatedNormal, or try fewer censored neighbors."
    )
  )
}

## Draw from a multivariate normal distribution truncated to a rectangle.
## The default uses the minimax-tilting sampler in TruncatedNormal.
sample_truncated_mvn <- function(n,
                                 mu,
                                 Sigma,
                                 lower,
                                 upper,
                                 method = c("TruncatedNormal", "tmvtnorm_gibbs", "tmvtnorm_rejection"),
                                 burn_in = 2000L,
                                 thinning = 5L) {
  method <- match.arg(method)
  
  n <- as.integer(n)
  if (length(n) != 1L || is.na(n) || n < 2L) {
    stop("n must be a single integer of at least 2.")
  }
  
  mu <- as.numeric(mu)
  Sigma <- make_posdef_prediction(as.matrix(Sigma))
  p <- length(mu)
  
  expand_bound <- function(x, name) {
    x <- as.numeric(x)
    if (length(x) == 1L) {
      return(rep(x, p))
    }
    if (length(x) != p) {
      stop(sprintf("%s must have length 1 or length(mu).", name))
    }
    x
  }
  
  lower <- expand_bound(lower, "lower")
  upper <- expand_bound(upper, "upper")
  
  if (!all(dim(Sigma) == c(p, p))) {
    stop("Sigma must be a square matrix with dimension length(mu).")
  }
  if (any(lower >= upper)) {
    stop("Every lower truncation bound must be smaller than its upper bound.")
  }
  
  if (method == "TruncatedNormal") {
    if (!requireNamespace("TruncatedNormal", quietly = TRUE)) {
      stop(
        paste0(
          "sample_method = 'TruncatedNormal' requires the TruncatedNormal package. ",
          "Install it with install.packages('TruncatedNormal')."
        )
      )
    }
    
    draws <- TruncatedNormal::rtmvnorm(
      n = n,
      mu = mu,
      sigma = Sigma,
      lb = lower,
      ub = upper
    )
  } else {
    if (!requireNamespace("tmvtnorm", quietly = TRUE)) {
      stop(
        paste0(
          "This sample_method requires the tmvtnorm package. ",
          "Install it with install.packages('tmvtnorm')."
        )
      )
    }
    
    if (method == "tmvtnorm_gibbs") {
      draws <- tmvtnorm::rtmvnorm(
        n = n,
        mean = mu,
        sigma = Sigma,
        lower = lower,
        upper = upper,
        algorithm = "gibbs",
        burn.in.samples = as.integer(burn_in),
        thinning = as.integer(thinning)
      )
    } 
  }
  
  draws <- as.matrix(draws)
  if (ncol(draws) != p && nrow(draws) == p) {
    draws <- t(draws)
  }
  if (p == 1L && ncol(draws) != 1L) {
    draws <- matrix(as.numeric(draws), ncol = 1L)
  }
  
  if (nrow(draws) != n || ncol(draws) != p) {
    stop("The truncated-normal sampler returned an unexpected dimension.")
  }
  if (any(!is.finite(draws))) {
    stop("The truncated-normal sampler returned non-finite values.")
  }
  
  tolerance <- 1e-8
  below_lower <- sweep(draws, 2, lower, FUN = "-") < -tolerance
  above_upper <- sweep(draws, 2, upper, FUN = "-") > tolerance
  if (any(below_lower) || any(above_upper)) {
    stop("Some truncated-normal draws violate the requested bounds.")
  }
  
  draws
}

## Predict at one test point while integrating over censored neighbors.
##
## The analytical mean and variance use the truncated moments m_C and V_C.
## The interval is based on composition samples and is generally asymmetric.
gp_predict_with_censoring <- function(x_star,
                                      x_neighbors,
                                      y_neighbors = NULL,
                                      params = NULL,
                                      tau_sq = NULL,
                                      ell = NULL,
                                      sigma_sq = NULL,
                                      beta0 = NULL,
                                      threshold_C,
                                      x_cols = NULL,
                                      predict_y = TRUE,
                                      moment_method = c("auto", "tmvtnorm", "TruncatedNormal", "mc", "univariate"),
                                      moment_mc_samples = 20000L,
                                      n_pred_samples = 5000L,
                                      sample_method = c("TruncatedNormal", "tmvtnorm_gibbs"),
                                      gibbs_burn_in = 2000L,
                                      gibbs_thinning = 5L,
                                      interval_level = 0.95,
                                      prediction_seed = NULL,
                                      return_draws = FALSE,
                                      jitter = 1e-8) {
  moment_method <- match.arg(moment_method)
  sample_method <- match.arg(sample_method)
  
  if (!is.numeric(interval_level) || length(interval_level) != 1L ||
      is.na(interval_level) || interval_level <= 0 || interval_level >= 1) {
    stop("interval_level must be a number strictly between 0 and 1.")
  }
  n_pred_samples <- as.integer(n_pred_samples)
  if (length(n_pred_samples) != 1L || is.na(n_pred_samples) || n_pred_samples < 2L) {
    stop("n_pred_samples must be a single integer of at least 2.")
  }
  
  if (!is.null(prediction_seed)) {
    set.seed(as.integer(prediction_seed))
  }
  
  if (is.null(x_cols)) {
    x_cols <- grep("^x[0-9]+$", names(x_neighbors), value = TRUE)
  }
  if (!length(x_cols)) {
    stop("Please provide x_cols.")
  }
  
  d <- length(x_cols)
  if (!is.null(params)) {
    tau_sq <- params[1]
    ell <- as.numeric(params[seq.int(2, d + 1L)])
    sigma_sq <- params[d + 2L]
    beta0 <- params[d + 3L]
  }
  
  if (is.null(tau_sq) || is.null(ell) || is.null(sigma_sq) || is.null(beta0)) {
    stop("Provide either params or tau_sq, ell, sigma_sq, and beta0.")
  }
  if (length(ell) != d) {
    stop("Length of ell must equal length(x_cols).")
  }
  
  if (is.null(y_neighbors)) {
    y_neighbors <- x_neighbors
  }
  
  X_nn <- as.matrix(x_neighbors[, x_cols, drop = FALSE])
  x_star <- as.matrix(x_star[, x_cols, drop = FALSE])
  y_nn <- as.numeric(y_neighbors$y)
  censored <- as.integer(y_neighbors$censored)
  n <- nrow(X_nn)
  
  if (length(y_nn) != n || length(censored) != n) {
    stop("x_neighbors and y_neighbors must have the same number of rows.")
  }
  
  idx_O <- which(censored == 0L)
  idx_C <- which(censored == 1L)
  m <- length(idx_C)
  
  ## Training responses include nugget variance; the new latent response does not.
  K_all <- gp_cov_matrix(rbind(X_nn, x_star), tau_sq = tau_sq, ell = ell)
  K_all <- make_posdef_prediction(
    K_all + diag(c(rep(sigma_sq, n), 0)),
    jitter = jitter
  )
  
  K_OO <- if (length(idx_O) > 0L) {
    as.matrix(K_all[idx_O, idx_O, drop = FALSE])
  } else {
    NULL
  }
  K_CO <- if (m > 0L && length(idx_O) > 0L) {
    as.matrix(K_all[idx_C, idx_O, drop = FALSE])
  } else {
    NULL
  }
  K_CC <- if (m > 0L) {
    as.matrix(K_all[idx_C, idx_C, drop = FALSE])
  } else {
    NULL
  }
  k_starO <- if (length(idx_O) > 0L) {
    matrix(K_all[n + 1L, idx_O, drop = FALSE], nrow = 1L)
  } else {
    NULL
  }
  k_starC <- if (m > 0L) {
    matrix(K_all[n + 1L, idx_C, drop = FALSE], nrow = 1L)
  } else {
    NULL
  }
  k_starstar <- as.numeric(K_all[n + 1L, n + 1L])
  
  mu_new <- as.numeric(beta0)
  mu_O <- rep(mu_new, length(idx_O))
  mu_C <- rep(mu_new, length(idx_C))
  
  alpha_interval <- 1 - interval_level
  
  ## Case 1: no censored neighbors. The predictive distribution is Gaussian.
  if (m == 0L) {
    if (length(idx_O) > 0L) {
      SOO_y <- chol_solve_prediction(
        K_OO,
        y_nn[idx_O] - mu_O,
        jitter = jitter
      )
      mean_pred <- as.numeric(mu_new + k_starO %*% SOO_y)
      
      SOO_k <- chol_solve_prediction(
        K_OO,
        t(k_starO),
        jitter = jitter
      )
      latent_var <- as.numeric(k_starstar - k_starO %*% SOO_k)
    } else {
      mean_pred <- mu_new
      latent_var <- k_starstar
    }
    
    latent_var <- max(latent_var, 1e-10)
    var_pred <- latent_var + if (predict_y) sigma_sq else 0
    var_pred <- max(as.numeric(var_pred), 1e-10)
    
    ci_lower <- stats::qnorm(
      alpha_interval / 2,
      mean = mean_pred,
      sd = sqrt(var_pred)
    )
    ci_upper <- stats::qnorm(
      1 - alpha_interval / 2,
      mean = mean_pred,
      sd = sqrt(var_pred)
    )
    
    normal_draws <- NULL
    if (return_draws) {
      normal_draws <- stats::rnorm(
        n = n_pred_samples,
        mean = mean_pred,
        sd = sqrt(var_pred)
      )
    }
    
    return(list(
      mean = mean_pred,
      var = var_pred,
      CI_lower = as.numeric(ci_lower),
      CI_upper = as.numeric(ci_upper),
      interval_length = as.numeric(ci_upper - ci_lower),
      sample_mean = if (is.null(normal_draws)) mean_pred else mean(normal_draws),
      sample_var = if (is.null(normal_draws)) var_pred else stats::var(normal_draws),
      interval_method = "Gaussian quantiles",
      moment_method = "not needed",
      n_pred_samples = if (is.null(normal_draws)) 0L else n_pred_samples,
      draws = normal_draws
    ))
  }
  
  ## Condition jointly on the uncensored neighbors.
  if (length(idx_O) > 0L) {
    S_OO_inv <- chol_solve_prediction(K_OO, jitter = jitter)
    
    mean_new_given_O <- as.numeric(
      mu_new + k_starO %*% S_OO_inv %*% (y_nn[idx_O] - mu_O)
    )
    var_new_given_O <- as.numeric(
      k_starstar - k_starO %*% S_OO_inv %*% t(k_starO)
    )
    
    mean_C_given_O <- as.numeric(
      mu_C + K_CO %*% S_OO_inv %*% (y_nn[idx_O] - mu_O)
    )
    cov_C_given_O <- make_posdef_prediction(
      K_CC - K_CO %*% S_OO_inv %*% t(K_CO),
      jitter = jitter
    )
    cov_new_C_given_O <- matrix(
      k_starC - k_starO %*% S_OO_inv %*% t(K_CO),
      nrow = 1L
    )
  } else {
    mean_new_given_O <- mu_new
    var_new_given_O <- k_starstar
    mean_C_given_O <- mu_C
    cov_C_given_O <- make_posdef_prediction(K_CC, jitter = jitter)
    cov_new_C_given_O <- matrix(k_starC, nrow = 1L)
  }
  
  ## Allow one common threshold, one threshold per neighbor, or one per
  ## censored neighbor.
  threshold_C <- as.numeric(threshold_C)
  if (length(threshold_C) == 1L) {
    c_vec <- rep(threshold_C, m)
  } else if (length(threshold_C) == n) {
    c_vec <- threshold_C[idx_C]
  } else if (length(threshold_C) == m) {
    c_vec <- threshold_C
  } else {
    stop(
      paste0(
        "threshold_C must have length 1, the number of all neighbors, ",
        "or the number of censored neighbors."
      )
    )
  }
  
  ## Step 1 of the composition sampler:
  ## draw the unknown censored-neighbor responses from their TMVN law.
  y_C_samples <- sample_truncated_mvn(
    n = n_pred_samples,
    mu = mean_C_given_O,
    Sigma = cov_C_given_O,
    lower = c_vec,
    upper = rep(Inf, m),
    method = sample_method,
    burn_in = gibbs_burn_in,
    thinning = gibbs_thinning
  )
  
  ## Compute m_C and V_C for the analytical predictive mean and variance.
  ## If the requested moment routine fails, use the same valid TMVN draws.
  mt <- tryCatch(
    truncated_mvn_moments(
      mu = mean_C_given_O,
      Sigma = cov_C_given_O,
      lower = c_vec,
      upper = rep(Inf, m),
      method = moment_method,
      mc_samples = moment_mc_samples,
      mc_burn_in = gibbs_burn_in,
      mc_thinning = gibbs_thinning
    ),
    error = function(e) NULL
  )
  
  if (is.null(mt)) {
    m_C <- colMeans(y_C_samples)
    V_C <- stats::cov(y_C_samples)
    if (m == 1L) {
      V_C <- matrix(V_C, nrow = 1L, ncol = 1L)
    }
    moment_method_used <- "empirical moments from predictive TMVN draws"
  } else {
    m_C <- as.numeric(mt$mean)
    V_C <- as.matrix(mt$var)
    if (m == 1L) {
      V_C <- matrix(V_C, nrow = 1L, ncol = 1L)
    }
    moment_method_used <- mt$method_used
  }
  
  S_CC_O_inv <- chol_solve_prediction(cov_C_given_O, jitter = jitter)
  
  ## B = Sigma_new,c|o Sigma_cc|o^{-1}
  B <- cov_new_C_given_O %*% S_CC_O_inv
  
  ## Residual variance of y_new conditional on exact y_C values.
  latent_residual_var <- as.numeric(
    var_new_given_O - B %*% t(cov_new_C_given_O)
  )
  latent_residual_var <- max(latent_residual_var, 1e-10)
  
  ## Exact mean and variance formulas, given m_C and V_C.
  mean_pred <- as.numeric(
    mean_new_given_O +
      B %*% matrix(m_C - mean_C_given_O, ncol = 1L)
  )
  var_pred <- as.numeric(
    latent_residual_var + B %*% V_C %*% t(B)
  )
  
  if (predict_y) {
    var_pred <- var_pred + sigma_sq
  }
  var_pred <- max(var_pred, 1e-10)
  
  ## Step 2 of the composition sampler:
  ## for every sampled censored vector, draw one new response.
  centered_C_samples <- sweep(
    y_C_samples,
    MARGIN = 2,
    STATS = mean_C_given_O,
    FUN = "-"
  )
  conditional_means <- as.numeric(
    mean_new_given_O + centered_C_samples %*% t(B)
  )
  
  conditional_residual_var <- latent_residual_var +
    if (predict_y) sigma_sq else 0
  conditional_residual_var <- max(conditional_residual_var, 1e-10)
  
  y_new_samples <- stats::rnorm(
    n = n_pred_samples,
    mean = conditional_means,
    sd = sqrt(conditional_residual_var)
  )
  
  empirical_interval <- stats::quantile(
    y_new_samples,
    probs = c(alpha_interval / 2, 1 - alpha_interval / 2),
    names = FALSE,
    na.rm = TRUE,
    type = 7
  )
  
  list(
    mean = mean_pred,
    var = var_pred,
    CI_lower = as.numeric(empirical_interval[1L]),
    CI_upper = as.numeric(empirical_interval[2L]),
    interval_length = as.numeric(empirical_interval[2L] - empirical_interval[1L]),
    sample_mean = mean(y_new_samples),
    sample_var = stats::var(y_new_samples),
    interval_method = "empirical predictive quantiles",
    moment_method = moment_method_used,
    n_pred_samples = n_pred_samples,
    draws = if (return_draws) y_new_samples else NULL
  )
}

## Compute distances from training rows to one test point for prediction neighbor search.
compute_prediction_distance <- function(training_subset,
                                        x_star,
                                        x_cols,
                                        method = c("euclidean", "correlation", "mahalanobis"),
                                        ell = NULL) {
  method <- match.arg(method)
  X <- as.matrix(training_subset[, x_cols, drop = FALSE])
  x_star <- as.numeric(x_star)
  
  if (method == "euclidean") {
    return(rowSums((t(t(X) - x_star))^2))
  }
  
  if (method == "correlation") {
    if (is.null(ell)) {
      stop("distance_method = 'correlation' requires ell.")
    }
    scaled_diff <- sweep(t(t(X) - x_star), 2, ell, "/")
    return(1 - exp(-0.5 * rowSums(scaled_diff^2)))
  }
  
  Sigma <- stats::cov(X) + diag(1e-6, ncol(X))
  Sigma_inv <- solve(Sigma)
  diff <- sweep(X, 2, x_star, "-")
  quad <- rowSums((diff %*% Sigma_inv) * diff)
  1 - exp(-0.5 * quad)
}

## Select prediction neighbors while limiting censored IDs and keeping all reps per ID.
select_prediction_neighbors_with_censoring <- function(train_data,
                                                       test_row,
                                                       x_cols,
                                                       k = 20L,
                                                       max_censored_ids = 5L,
                                                       min_uncensored_ids = NULL,
                                                       id_col = "ID",
                                                       distance_method = c("euclidean", "correlation", "mahalanobis"),
                                                       ell = NULL) {
  distance_method <- match.arg(distance_method)
  train_work <- as.data.frame(train_data)
  test_work <- as.data.frame(test_row)
  
  if (!id_col %in% names(train_work)) {
    train_work[[id_col]] <- seq_len(nrow(train_work))
  }
  if (!id_col %in% names(test_work)) {
    test_work[[id_col]] <- -seq_len(nrow(test_work))
  }
  
  train_work$data_type <- "train"
  test_work$data_type <- "test"
  combined_data <- bind_rows_fill(train_work, test_work)
  combined_data <- add_marker_id(combined_data, x_cols)
  
  training_subset <- combined_data[combined_data$data_type == "train", , drop = FALSE]
  x_star <- as.matrix(test_work[1, x_cols, drop = FALSE])
  training_subset$distance <- compute_prediction_distance(
    training_subset = training_subset,
    x_star = x_star,
    x_cols = x_cols,
    method = distance_method,
    ell = ell
  )
  
  group_distance <- stats::aggregate(
    distance ~ marker_id,
    data = training_subset,
    FUN = min,
    na.rm = TRUE
  )
  group_distance <- group_distance[order(group_distance$distance), , drop = FALSE]
  
  selected_groups <- head(group_distance$marker_id, min(k, nrow(group_distance)))
  selected_rows <- training_subset[training_subset$marker_id %in% selected_groups, , drop = FALSE]
  
  min_by_id <- stats::aggregate(
    distance ~ ID,
    data = selected_rows,
    FUN = min,
    na.rm = TRUE
  )
  cens_by_id <- stats::aggregate(
    censored ~ ID,
    data = selected_rows,
    FUN = function(z) as.integer(all(z == 1L))
  )
  sum_selected <- merge(min_by_id, cens_by_id, by = "ID", all = TRUE)
  names(sum_selected)[names(sum_selected) == "censored"] <- "is_cens_id"
  sum_selected <- sum_selected[order(sum_selected$distance), , drop = FALSE]
  
  outside_unc <- training_subset[
    training_subset$censored == 0L & !(training_subset$ID %in% sum_selected$ID),
    ,
    drop = FALSE
  ]
  
  if (nrow(outside_unc) > 0L) {
    sum_train_unc_outside <- stats::aggregate(
      distance ~ ID,
      data = outside_unc,
      FUN = min,
      na.rm = TRUE
    )
    sum_train_unc_outside$is_cens_id <- 0L
    sum_train_unc_outside <- sum_train_unc_outside[order(sum_train_unc_outside$distance), , drop = FALSE]
  } else {
    sum_train_unc_outside <- data.frame(ID = integer(0), distance = numeric(0), is_cens_id = integer(0))
  }
  
  candidates <- head(sum_selected, k)
  if (is.null(min_uncensored_ids)) {
    min_uncensored_ids <- max(0L, k - max_censored_ids)
  }
  
  unc_now <- sum(candidates$is_cens_id == 0L)
  if (unc_now < min_uncensored_ids) {
    need_more_unc <- min_uncensored_ids - unc_now
    cens_candidates <- candidates[candidates$is_cens_id == 1L, , drop = FALSE]
    cens_candidates <- cens_candidates[order(cens_candidates$distance, decreasing = TRUE), , drop = FALSE]
    drop_cens <- head(cens_candidates$ID, need_more_unc)
    
    add_unc <- sum_train_unc_outside[
      !(sum_train_unc_outside$ID %in% candidates$ID),
      ,
      drop = FALSE
    ]
    add_unc <- head(add_unc, need_more_unc)
    
    candidates <- candidates[!(candidates$ID %in% drop_cens), , drop = FALSE]
    candidates <- rbind(candidates, add_unc[, names(candidates), drop = FALSE])
    candidates <- candidates[order(candidates$distance), , drop = FALSE]
    candidates <- head(candidates, k)
  }
  
  final_ids <- candidates$ID
  neighbor_rows <- training_subset[training_subset$ID %in% final_ids, , drop = FALSE]
  neighbor_rows$.ord <- match(neighbor_rows$ID, final_ids)
  neighbor_rows <- neighbor_rows[order(neighbor_rows$.ord, neighbor_rows$distance), , drop = FALSE]
  
  x_neighbors <- neighbor_rows[, c("ID", "censored", x_cols), drop = FALSE]
  y_neighbors <- neighbor_rows[, c("ID", "censored", "y"), drop = FALSE]
  
  n_total <- length(unique(x_neighbors$ID))
  n_unc_id <- length(unique(x_neighbors$ID[x_neighbors$censored == 0L]))
  n_cens_id <- n_total - n_unc_id
  
  list(
    x_neighbors = x_neighbors,
    y_neighbors = y_neighbors,
    final_ids = final_ids,
    n_total = n_total,
    n_uncens = n_unc_id,
    n_cens = n_cens_id,
    training_subset = training_subset
  )
}

## Predict over test data using censoring-aware neighbor selection.
## Intervals are empirical predictive quantiles whenever censored neighbors occur.
predict_censored_gp_nn_with_censoring <- function(fit = NULL,
                                                  train_data = NULL,
                                                  test_data,
                                                  params = NULL,
                                                  x_cols = NULL,
                                                  threshold_C = NULL,
                                                  k = NULL,
                                                  max_censored_ids = 5L,
                                                  min_uncensored_ids = NULL,
                                                  id_col = "ID",
                                                  y_col = "y",
                                                  censored_col = "censored",
                                                  distance_method = c("euclidean", "correlation", "mahalanobis"),
                                                  moment_method = c("auto", "tmvtnorm", "TruncatedNormal", "mc", "univariate"),
                                                  moment_mc_samples = 20000L,
                                                  n_pred_samples = 5000L,
                                                  sample_method = c("TruncatedNormal", "tmvtnorm_gibbs"),
                                                  gibbs_burn_in = 2000L,
                                                  gibbs_thinning = 5L,
                                                  predict_y = TRUE,
                                                  interval_level = 0.95,
                                                  prediction_seed = 123L,
                                                  return_draws = FALSE,
                                                  unique_test_inputs = TRUE,
                                                  verbose = TRUE) {
  distance_method <- match.arg(distance_method)
  moment_method <- match.arg(moment_method)
  sample_method <- match.arg(sample_method)
  
  if (!is.null(fit)) {
    if (is.null(params)) {
      params <- fit$estimated_params
    }
    if (is.null(x_cols)) {
      x_cols <- fit$x_cols
    }
    if (is.null(threshold_C)) {
      threshold_C <- fit$threshold_C
    }
    if (is.null(k)) {
      k <- fit$k
    }
    if (is.null(train_data)) {
      if (!is.null(fit$train_data)) {
        train_data <- fit$train_data
      } else {
        train_data <- fit$prepared_data
      }
    }
  }
  
  if (is.null(train_data) || is.null(params) ||
      is.null(x_cols) || is.null(threshold_C)) {
    stop("Provide fit, or provide train_data, params, x_cols, and threshold_C.")
  }
  if (is.null(k)) {
    k <- 20L
  }
  
  train_work <- as.data.frame(train_data)
  test_work <- as.data.frame(test_data)
  train_work$y <- train_work[[y_col]]
  train_work$censored <- as.integer(train_work[[censored_col]])
  
  d <- length(x_cols)
  ell <- as.numeric(params[seq.int(2, d + 1L)])
  
  if (unique_test_inputs) {
    keep <- !duplicated(test_work[, x_cols, drop = FALSE])
    test_unique <- test_work[keep, , drop = FALSE]
  } else {
    test_unique <- test_work
  }
  rownames(test_unique) <- NULL
  
  if (verbose) {
    cat(sprintf(
      "Reduced from %d to %d unique test inputs.\n",
      nrow(test_work),
      nrow(test_unique)
    ))
  }
  
  results_list <- vector("list", nrow(test_unique))
  predictive_draws <- if (return_draws) {
    vector("list", nrow(test_unique))
  } else {
    NULL
  }
  
  for (i in seq_len(nrow(test_unique))) {
    if (verbose) {
      cat(sprintf("Processing test row %d / %d...\n", i, nrow(test_unique)))
    }
    
    test_row <- test_unique[i, , drop = FALSE]
    selected <- select_prediction_neighbors_with_censoring(
      train_data = train_work,
      test_row = test_row,
      x_cols = x_cols,
      k = k,
      max_censored_ids = max_censored_ids,
      min_uncensored_ids = min_uncensored_ids,
      id_col = id_col,
      distance_method = distance_method,
      ell = ell
    )
    
    seed_i <- if (is.null(prediction_seed)) {
      NULL
    } else {
      as.integer(prediction_seed) + i - 1L
    }
    
    out <- tryCatch(
      gp_predict_with_censoring(
        x_star = test_row[, x_cols, drop = FALSE],
        x_neighbors = selected$x_neighbors,
        y_neighbors = selected$y_neighbors,
        params = params,
        threshold_C = threshold_C,
        x_cols = x_cols,
        predict_y = predict_y,
        moment_method = moment_method,
        moment_mc_samples = moment_mc_samples,
        n_pred_samples = n_pred_samples,
        sample_method = sample_method,
        gibbs_burn_in = gibbs_burn_in,
        gibbs_thinning = gibbs_thinning,
        interval_level = interval_level,
        prediction_seed = seed_i,
        return_draws = return_draws
      ),
      error = function(e) {
        if (verbose) {
          message(sprintf("Error at test row %d: %s", i, e$message))
        }
        list(
          mean = NA_real_,
          var = NA_real_,
          CI_lower = NA_real_,
          CI_upper = NA_real_,
          interval_length = NA_real_,
          sample_mean = NA_real_,
          sample_var = NA_real_,
          interval_method = NA_character_,
          moment_method = NA_character_,
          n_pred_samples = NA_integer_,
          draws = NULL
        )
      }
    )
    
    results_list[[i]] <- data.frame(
      ID_test = if (id_col %in% names(test_row)) test_row[[id_col]] else i,
      mean_pred = out$mean,
      var_pred = out$var,
      sample_mean = out$sample_mean,
      sample_var = out$sample_var,
      CI_lower = out$CI_lower,
      CI_upper = out$CI_upper,
      interval_length = out$interval_length,
      interval_method = out$interval_method,
      moment_method = out$moment_method,
      n_pred_samples = out$n_pred_samples,
      n_total = selected$n_total,
      n_cens = selected$n_cens,
      n_uncens = selected$n_uncens,
      stringsAsFactors = FALSE
    )
    
    if (return_draws) {
      predictive_draws[[i]] <- out$draws
    }
  }
  
  results_df <- do.call(rbind, results_list)
  average_interval <- mean(results_df$interval_length, na.rm = TRUE)
  
  merge_columns <- c(
    "ID_test",
    "mean_pred",
    "var_pred",
    "sample_mean",
    "sample_var",
    "CI_lower",
    "CI_upper",
    "interval_length",
    "interval_method",
    "moment_method",
    "n_pred_samples"
  )
  
  if (id_col %in% names(test_work)) {
    merged_df <- merge(
      test_work,
      results_df[, merge_columns, drop = FALSE],
      by.x = id_col,
      by.y = "ID_test",
      all.x = TRUE,
      sort = FALSE
    )
  } else {
    merged_df <- test_work
  }
  
  if (verbose) {
    cat("Prediction loop finished.\n")
    cat(sprintf(
      "Average %.1f%% interval length = %.4f\n",
      interval_level * 100,
      average_interval
    ))
  }
  
  list(
    results_df = results_df,
    merged_df = merged_df,
    average_interval = average_interval,
    unique_test_data = test_unique,
    predictive_draws = predictive_draws
  )
}

## Print a compact summary of a censored GP nearest-neighbor fit object.
print.censored_gp_nn_fit <- function(x, ...) {
  cat("Censored GP nearest-neighbor fit\n")
  cat("Dimension:", length(x$x_cols), "\n")
  cat("Predictors:", paste(x$x_cols, collapse = ", "), "\n")
  cat("Censored probability method:", x$censor_method, "\n")
  cat("Threshold C:", x$threshold_C, "\n")
  cat("Neighbors k:", x$k, "\n")
  cat("Used parallel optimization:", x$used_parallel, "\n")
  cat("Estimated parameters:\n")
  print(x$estimated_params)
  cat("Final log-likelihood:", x$final_log_likelihood, "\n")
  cat("Counts:\n")
  print(x$result$counts)
  cat("Optim elapsed seconds:", x$optim_time["elapsed"], "\n")
  invisible(x)
}

## ------------------------------------------------------------------
## Minimal usage examples
## ------------------------------------------------------------------
##
## ## 1D with real train/test split by original replicate ID:
## dat1 <- simulate_1d_censored_data(
##   n = 800,
##   split_method = "id",
##   test_frac = 0.3,
##   threshold_C = 1.221277
## )
## fit1 <- fit_censored_gp_nn(
##   train_data = dat1$train_data,
##   x_cols = "x1",
##   threshold_C = dat1$threshold_C,
##   k = 20,
##   censor_method = "sov",
##   init_params = c(0.3462413, 1.5, 0.02911056, 0.01410213),
##   lower_bounds = c(1e-6, 1e-6, 1e-6, -Inf),
##   upper_bounds = c(Inf, Inf, Inf, Inf)
## )
## pred1 <- predict_censored_gp_nn(fit1, dat1$test_data)
## pred1_cens <- predict_censored_gp_nn_with_censoring(
##  fit = fit1,
##  test_data = test_data,
##  k = 20,
##  max_censored_ids = 5,
##  distance_method = "euclidean",
##  n_pred_samples = 20000,
##  sample_method = "TruncatedNormal",
##  prediction_seed = 123,
##  predict_y = TRUE,
##  verbose = TRUE
##)
##
## ## 1D with fixed 100-point test grid from 0 to 1:
## dat1_grid <- simulate_1d_censored_data(
##   n = 800,
##   split_method = "fixed_grid",
##   n_test_grid = 100,
##   threshold_C = 1.221277
## )
##
## ## 8D borehole:
## dat8 <- simulate_borehole_censored_data(censor_quantile = 0.8)
## init8 <- c(28.48859, 0.9, 0.5, 0.2, 0.46, 0.2, 0.35, 0.2, 0.7, 1.027467, 59.97952)
## fit8 <- fit_censored_gp_nn(
##   train_data = dat8$train_data,
##   x_cols = dat8$x_cols,
##   threshold_C = dat8$threshold_C,
##   k = 20,
##   censor_method = "sov",
##   init_params = init8,
##   lower_bounds = c(1e-6, rep(1e-6, 8), 1e-6, -Inf),
##   upper_bounds = c(Inf, rep(Inf, 8), Inf, Inf)
## )
## pred8 <- predict_censored_gp_nn(fit8, dat8$test_data)
## pred8_cens <- predict_censored_gp_nn_with_censoring(
##  fit = fit8,
##  test_data = dat8$test_data,
##  k = 20,
##  max_censored_ids = 5,
##  distance_method = "euclidean",
##  n_pred_samples = 20000,
##  sample_method = "TruncatedNormal",
##  prediction_seed = 123,
##  predict_y = TRUE,
##  verbose = TRUE
## )
