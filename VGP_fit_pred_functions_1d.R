## =====================================================
## VGP_compare.R
## Compare three 1D VGP data treatments using the same fitting
## and fixed-grid prediction workflow.
##
## Three training-data cases:
##   1. "raw_all":         use original noisy y as all-observed train rows
##   2. "remove_censored": replace y by y_censored, then remove censored rows
##   3. "threshold_all":   replace y by y_censored and keep all rows
##
## Prediction test_data is always:
##   data.frame(x1 = seq(0, 1, length.out = 100), ...)
## =====================================================

required_vgp_packages <- function(use_parallel = TRUE) {
  pkgs <- c("FNN", "mvtnorm", "ggplot2")
  if (use_parallel) {
    pkgs <- c(pkgs, "parallel", "optimParallel")
  }
  
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs)) {
    stop("Please install missing package(s): ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
  }
  
  invisible(TRUE)
}

true_f_1d <- function(x) {
  sin(2 * (7 * x - 2)) +
    3 * exp(-15^2 * (x - 0.5)^2) +
    sin(2 * (2 * x - 0.5))
}


make_fixed_test_grid_1d <- function(n_test = 100L,
                                    threshold_C = 1.221277,
                                    f = true_f_1d) {
  x1 <- seq(0, 1, length.out = n_test)
  y <- f(x1)
  has_threshold <- !is.null(threshold_C) &&
    length(threshold_C) == 1L &&
    is.finite(threshold_C)
  y_censored <- if (has_threshold) pmin(y, threshold_C) else y
  
  data.frame(
    x1 = x1,
    y = y,
    ID = seq_len(n_test),
    y_censored = y_censored,
    censored = if (has_threshold) as.integer(y > threshold_C) else rep(0L, n_test)
  )
}



simulate_original_1d_base <- function(n = 800L,
                                      seed = 2024L,
                                      reps_min = 1L,
                                      reps_max = 3L,
                                      noise_sd = 0.2,
                                      f = true_f_1d) {
  set.seed(seed)
  x <- runif(n, 0, 1)
  rep_counts <- sample(reps_min:reps_max, n, replace = TRUE)
  x_rep <- rep(x, times = rep_counts)
  y <- f(x_rep) + rnorm(length(x_rep), mean = 0, sd = noise_sd)
  
  data <- data.frame(x = x_rep, y = y)
  data$ID <- rep(seq_len(n), rep_counts)
  data
}

add_censor_columns_original <- function(data,
                                        threshold_C = NULL,
                                        censor_quantile = NULL) {
  if (is.null(threshold_C)) {
    if (is.null(censor_quantile)) {
      stop("Provide threshold_C or censor_quantile.")
    }
    threshold_C <- as.numeric(stats::quantile(data$y, censor_quantile))
  }
  
  data$y_censored <- data$y
  censor_candidates <- which(data$y > threshold_C)
  data$y_censored[censor_candidates] <- threshold_C
  data$censored <- 0L
  data$censored[censor_candidates] <- 1L
  
  list(data = data, threshold_C = threshold_C)
}

split_original_by_id <- function(data,
                                 test_frac = 0.3,
                                 seed_split = 2024L) {
  set.seed(seed_split)
  unique_ids <- unique(data$ID)
  test_ids <- sample(unique_ids, size = floor(test_frac * length(unique_ids)))
  
  list(
    train_data = data[!data$ID %in% test_ids, , drop = FALSE],
    split_test_data = data[data$ID %in% test_ids, , drop = FALSE],
    test_ids = test_ids
  )
}

prepare_vgp_case <- function(case = c("raw_all", "remove_censored", "threshold_all"),
                             n = 800L,
                             seed = 2024L,
                             seed_split = 2024L,
                             test_frac = 0.3,
                             threshold_C = NULL,
                             censor_quantile = NULL,
                             n_test_grid = 100L,
                             f = true_f_1d) {
  case <- match.arg(case)
  
  if (case %in% c("remove_censored", "threshold_all") &&
      is.null(threshold_C) && is.null(censor_quantile)) {
    threshold_C <- 1.221277
  }
  
  data <- simulate_original_1d_base(n = n, seed = seed, f = f)
  
  if (case == "raw_all") {
    threshold_C <- NA_real_
  } else {
    cens <- add_censor_columns_original(
      data = data,
      threshold_C = threshold_C,
      censor_quantile = censor_quantile
    )
    data <- cens$data
    threshold_C <- cens$threshold_C
  }
  
  split <- split_original_by_id(
    data = data,
    test_frac = test_frac,
    seed_split = seed_split
  )
  train_data <- split$train_data
  
  if (case == "raw_all") {
    names(train_data)[names(train_data) == "x"] <- "x1"
  } else {
    train_data$y_true <- train_data$y
    train_data$y <- NULL
    names(train_data)[names(train_data) == "y_censored"] <- "y"
    names(train_data)[names(train_data) == "x"] <- "x1"
    
    if (case == "remove_censored") {
      train_data <- train_data[train_data$censored == 0L, , drop = FALSE]
    }
  }
  
  test_data <- make_fixed_test_grid_1d(n_test = n_test_grid)
  rownames(train_data) <- NULL
  rownames(test_data) <- NULL
  
  list(
    case = case,
    train_data = train_data,
    test_data = test_data,
    split_test_data = split$split_test_data,
    threshold_C = threshold_C,
    x_cols = "x1",
    f = f
  )
}

make_positive_definite <- function(Sigma, jitter = 1e-6) {
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

cov_from_sqdist <- function(D, tau_sq, ell) {
  scaled <- matrix(0, nrow = dim(D)[1], ncol = dim(D)[2])
  
  for (j in seq_along(ell)) {
    scaled <- scaled + D[, , j] / (2 * ell[j]^2)
  }
  
  tau_sq * exp(-scaled)
}

gp_cov_matrix_vgp <- function(X, tau_sq, ell) {
  cov_from_sqdist(sqdist_array(X), tau_sq = tau_sq, ell = ell)
}

order_vgp_data <- function(train_data,
                           x_cols,
                           order_method = c("weighted_sum", "lexicographic", "none")) {
  order_method <- match.arg(order_method)
  data <- train_data
  if ("ID" %in% names(data)) {
    data$ID <- NULL
  }
  
  if (order_method == "weighted_sum") {
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

add_marker_id_vgp <- function(data, x_cols) {
  marker_id <- rep(NA_integer_, nrow(data))
  complete <- complete.cases(data[, x_cols, drop = FALSE])
  
  if (any(complete)) {
    keys <- do.call(paste, c(data[complete, x_cols, drop = FALSE], sep = "_"))
    marker_id[complete] <- as.integer(factor(keys, levels = unique(keys)))
  }
  
  data$marker_id <- marker_id
  data
}

find_marker_id_neighbors_vgp <- function(data, x_cols, max_marker_neighbors = 20L) {
  data <- data[order(data$marker_id), , drop = FALSE]
  results <- vector("list", nrow(data))
  
  for (query_idx in seq_len(nrow(data))) {
    query_marker_id <- data$marker_id[query_idx]
    valid_marker_ids <- unique(data$marker_id[data$marker_id < query_marker_id])
    valid_marker_ids <- tail(valid_marker_ids, max_marker_neighbors)
    ref_data <- data[data$marker_id %in% valid_marker_ids, , drop = FALSE]
    
    if (nrow(ref_data) > 0L) {
      nn_result <- FNN::get.knnx(
        ref_data[, x_cols, drop = FALSE],
        data[query_idx, x_cols, drop = FALSE],
        k = nrow(ref_data)
      )
      matched_indices <- nn_result$nn.index[1, ]
      matched_distances <- nn_result$nn.dist[1, ]
      
      results[[query_idx]] <- data.frame(
        query_index = query_idx,
        query_marker_id = query_marker_id,
        neighbor_indices = paste(rownames(ref_data)[matched_indices], collapse = ","),
        neighbor_marker_ids = paste(ref_data$marker_id[matched_indices], collapse = ","),
        distances = paste(matched_distances, collapse = ",")
      )
    } else {
      results[[query_idx]] <- data.frame(
        query_index = query_idx,
        query_marker_id = query_marker_id,
        neighbor_indices = NA_character_,
        neighbor_marker_ids = NA_character_,
        distances = NA_character_
      )
    }
  }
  
  do.call(rbind, results)
}

parse_id_string <- function(x) {
  if (length(x) == 0L || is.na(x[1]) || !nzchar(x[1])) {
    return(integer(0))
  }
  as.integer(strsplit(x[1], ",", fixed = TRUE)[[1]])
}

compute_distance_cache_vgp <- function(nearest_neighbors, data, x_cols) {
  cache <- list()
  nn_by_marker <- nearest_neighbors[!duplicated(nearest_neighbors$query_marker_id), , drop = FALSE]
  
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

split_vgp_params <- function(params, x_cols) {
  d <- length(x_cols)
  list(
    tau_sq = params[1],
    ell = as.numeric(params[seq.int(2, d + 1L)]),
    sigma_sq = params[d + 2L],
    beta0 = params[d + 3L]
  )
}

default_vgp_init <- function(train_data, x_cols) {
  y_var <- stats::var(train_data$y)
  if (!is.finite(y_var) || y_var <= 0) {
    y_var <- 1
  }
  c(y_var, rep(2, length(x_cols)), max(0.05 * y_var, 1e-4), mean(train_data$y))
}

vgp_case_initial_params <- function(case = c("raw_all", "remove_censored", "threshold_all")) {
  case <- match.arg(case)
  
  switch(
    case,
    raw_all = c(
      tau_sq = 1.612907,
      ell1 = 2,
      sigma_sq = 0.034895,
      beta0 = 0.8218356
    ),
    remove_censored = c(
      tau_sq = 0.3462413,
      ell1 = 2,
      sigma_sq = 0.02911056,
      beta0 = -0.06377002
    ),
    threshold_all = c(
      tau_sq = 0.4797904,
      ell1 = 2,
      sigma_sq = 0.01978341,
      beta0 = 0.4348702
    )
  )
}

vgp_case_bounds <- function(case = c("raw_all", "remove_censored", "threshold_all")) {
  match.arg(case)
  list(
    lower = c(tau_sq = 1e-6, ell1 = 1e-6, sigma_sq = 1e-6, beta0 = -Inf),
    upper = c(tau_sq = Inf, ell1 = Inf, sigma_sq = Inf, beta0 = Inf)
  )
}

default_vgp_bounds <- function(x_cols) {
  list(
    lower = c(1e-6, rep(1e-6, length(x_cols)), 1e-6, -Inf),
    upper = c(Inf, rep(Inf, length(x_cols)), Inf, Inf)
  )
}

conditional_vgp_params <- function(query_marker_id,
                                   data,
                                   nn_by_marker,
                                   distance_cache,
                                   params,
                                   x_cols) {
  par <- split_vgp_params(params, x_cols)
  query_rows <- which(data$marker_id == query_marker_id)
  x_query <- as.matrix(data[query_rows, x_cols, drop = FALSE])
  mean_query <- rep(par$beta0, length(query_rows))
  cov_query <- gp_cov_matrix_vgp(x_query, par$tau_sq, par$ell) +
    par$sigma_sq * diag(length(query_rows))
  
  nn_row <- nn_by_marker[nn_by_marker$query_marker_id == query_marker_id, , drop = FALSE]
  neighbor_ids <- parse_id_string(nn_row$neighbor_marker_ids)
  
  if (!length(neighbor_ids)) {
    return(list(mean = mean_query, covariance = make_positive_definite(cov_query)))
  }
  
  cache <- distance_cache[[as.character(query_marker_id)]]
  if (is.null(cache)) {
    return(list(mean = mean_query, covariance = make_positive_definite(cov_query)))
  }
  
  neighbor_rows <- cache$neighbor_rows
  y_neighbors <- data$y[neighbor_rows]
  n_query <- cache$n_query
  n_total <- dim(cache$D)[1]
  neighbor_idx <- seq.int(n_query + 1L, n_total)
  
  D_qn <- cache$D[seq_len(n_query), neighbor_idx, , drop = FALSE]
  D_nn <- cache$D[neighbor_idx, neighbor_idx, , drop = FALSE]
  
  cov_qn <- cov_from_sqdist(D_qn, tau_sq = par$tau_sq, ell = par$ell)
  cov_nn <- cov_from_sqdist(D_nn, tau_sq = par$tau_sq, ell = par$ell) +
    par$sigma_sq * diag(length(neighbor_rows))
  cov_nn <- make_positive_definite(cov_nn)
  
  inv_cov_nn <- solve(cov_nn)
  residual_neighbors <- y_neighbors - par$beta0
  cond_mean <- as.numeric(mean_query + cov_qn %*% inv_cov_nn %*% residual_neighbors)
  cond_cov <- cov_query - cov_qn %*% inv_cov_nn %*% t(cov_qn)
  
  list(mean = cond_mean, covariance = make_positive_definite(cond_cov))
}

neg_log_likelihood_vgp <- function(params,
                                   simulated_data,
                                   nearest_neighbors,
                                   distance_cache,
                                   x_cols) {
  par <- split_vgp_params(params, x_cols)
  if (par$tau_sq <= 0 || par$sigma_sq <= 0 || any(par$ell <= 0)) {
    return(Inf)
  }
  
  nn_by_marker <- nearest_neighbors[!duplicated(nearest_neighbors$query_marker_id), , drop = FALSE]
  unique_marker_ids <- unique(nearest_neighbors$query_marker_id)
  log_likelihood <- 0
  
  for (query_marker_id in unique_marker_ids) {
    query_rows <- which(simulated_data$marker_id == query_marker_id)
    cond <- conditional_vgp_params(
      query_marker_id = query_marker_id,
      data = simulated_data,
      nn_by_marker = nn_by_marker,
      distance_cache = distance_cache,
      params = params,
      x_cols = x_cols
    )
    
    ll <- mvtnorm::dmvnorm(
      simulated_data$y[query_rows],
      mean = cond$mean,
      sigma = cond$covariance,
      log = TRUE
    )
    
    if (!is.finite(ll)) {
      return(Inf)
    }
    log_likelihood <- log_likelihood + ll
  }
  
  -log_likelihood
}

prepare_vgp_fit_data <- function(train_data,
                                 x_cols = "x1",
                                 max_marker_neighbors = 20L,
                                 order_method = c("weighted_sum", "lexicographic", "none")) {
  order_method <- match.arg(order_method)
  data <- order_vgp_data(train_data, x_cols = x_cols, order_method = order_method)
  simulated_data <- add_marker_id_vgp(data, x_cols = x_cols)
  nearest_neighbors <- find_marker_id_neighbors_vgp(
    simulated_data,
    x_cols = x_cols,
    max_marker_neighbors = max_marker_neighbors
  )
  distance_cache <- compute_distance_cache_vgp(
    nearest_neighbors = nearest_neighbors,
    data = simulated_data,
    x_cols = x_cols
  )
  
  list(
    simulated_data = simulated_data,
    nearest_neighbors = nearest_neighbors,
    distance_cache = distance_cache
  )
}

fit_vgp_model <- function(train_data,
                          x_cols = "x1",
                          init_params = NULL,
                          lower_bounds = NULL,
                          upper_bounds = NULL,
                          max_marker_neighbors = 20L,
                          order_method = c("weighted_sum", "lexicographic", "none"),
                          use_parallel = TRUE,
                          n_cores = max(1L, parallel::detectCores() - 1L),
                          optim_control = list()) {
  required_vgp_packages(use_parallel = use_parallel)
  order_method <- match.arg(order_method)
  
  prepared <- prepare_vgp_fit_data(
    train_data = train_data,
    x_cols = x_cols,
    max_marker_neighbors = max_marker_neighbors,
    order_method = order_method
  )
  
  if (is.null(init_params)) {
    init_params <- default_vgp_init(train_data, x_cols)
  }
  bounds <- default_vgp_bounds(x_cols)
  if (is.null(lower_bounds)) {
    lower_bounds <- bounds$lower
  }
  if (is.null(upper_bounds)) {
    upper_bounds <- bounds$upper
  }
  
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
            "make_positive_definite",
            "sqdist_array",
            "cov_from_sqdist",
            "gp_cov_matrix_vgp",
            "parse_id_string",
            "split_vgp_params",
            "conditional_vgp_params",
            "neg_log_likelihood_vgp"
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
            fn = neg_log_likelihood_vgp,
            method = "L-BFGS-B",
            lower = lower_bounds,
            upper = upper_bounds,
            control = optim_control,
            parallel = list(cl = cl, forward = FALSE, loginfo = FALSE),
            simulated_data = prepared$simulated_data,
            nearest_neighbors = prepared$nearest_neighbors,
            distance_cache = prepared$distance_cache,
            x_cols = x_cols
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
        fn = neg_log_likelihood_vgp,
        method = "L-BFGS-B",
        lower = lower_bounds,
        upper = upper_bounds,
        control = optim_control,
        simulated_data = prepared$simulated_data,
        nearest_neighbors = prepared$nearest_neighbors,
        distance_cache = prepared$distance_cache,
        x_cols = x_cols
      )
    })
  }
  
  full_end_time <- Sys.time()
  
  out <- list(
    result = result,
    estimated_params = result$par,
    final_log_likelihood = -result$value,
    x_cols = x_cols,
    train_data = train_data,
    max_marker_neighbors = max_marker_neighbors,
    init_params = init_params,
    lower_bounds = lower_bounds,
    upper_bounds = upper_bounds,
    prepared = prepared,
    used_parallel = ran_parallel,
    optim_time = optim_time,
    total_seconds = as.numeric(full_end_time - full_start_time, units = "secs")
  )
  
  class(out) <- "vgp_fit"
  out
}

predict_vgp <- function(fit,
                        test_data,
                        train_data = fit$train_data,
                        x_cols = fit$x_cols,
                        params = fit$estimated_params,
                        max_marker_neighbors = fit$max_marker_neighbors,
                        interval_level = 0.95) {
  par <- split_vgp_params(params, x_cols)
  
  predictive_results <- data.frame(
    test_index = seq_len(nrow(test_data)),
    predicted_mean = NA_real_,
    predicted_variance = NA_real_
  )
  
  bind_cols <- unique(c(x_cols, "y"))
  missing_train_cols <- setdiff(bind_cols, names(train_data))
  missing_test_cols <- setdiff(bind_cols, names(test_data))
  if (length(missing_train_cols)) {
    stop("train_data is missing column(s): ", paste(missing_train_cols, collapse = ", "), call. = FALSE)
  }
  if (length(missing_test_cols)) {
    stop("test_data is missing column(s): ", paste(missing_test_cols, collapse = ", "), call. = FALSE)
  }
  
  train_work <- train_data[, bind_cols, drop = FALSE]
  for (i in seq_len(nrow(test_data))) {
    test_row <- test_data[i, bind_cols, drop = FALSE]
    x_star <- as.matrix(test_row[, x_cols, drop = FALSE])
    
    test_row$data_type <- "test"
    train_work$data_type <- "train"
    combined_data <- rbind(train_work, test_row)
    combined_data <- add_marker_id_vgp(combined_data, x_cols)
    training_subset <- combined_data[combined_data$data_type == "train", , drop = FALSE]
    
    X_train <- as.matrix(training_subset[, x_cols, drop = FALSE])
    distances <- rowSums((t(t(X_train) - as.numeric(x_star)))^2)
    training_subset$distance <- distances
    
    min_dist <- stats::aggregate(distance ~ marker_id, data = training_subset, FUN = min)
    min_dist <- min_dist[order(min_dist$distance), , drop = FALSE]
    selected_groups <- head(min_dist$marker_id, min(max_marker_neighbors, nrow(min_dist)))
    selected_indices <- which(training_subset$marker_id %in% selected_groups)
    
    x_neighbors <- as.matrix(training_subset[selected_indices, x_cols, drop = FALSE])
    y_neighbors <- training_subset$y[selected_indices]
    
    cov_neighbors <- gp_cov_matrix_vgp(x_neighbors, par$tau_sq, par$ell)
    cov_neighbors_y <- cov_neighbors + par$sigma_sq * diag(nrow(x_neighbors))
    cov_neighbors_y <- make_positive_definite(cov_neighbors_y)
    
    D_star_neighbors <- sqdist_cross_array(x_star, x_neighbors)
    cov_star_neighbors <- cov_from_sqdist(D_star_neighbors, tau_sq = par$tau_sq, ell = par$ell)
    cov_star <- par$tau_sq + par$sigma_sq
    
    inv_cov_neighbors <- solve(cov_neighbors_y)
    u_star <- as.numeric(
      par$beta0 +
        cov_star_neighbors %*%
        inv_cov_neighbors %*%
        (y_neighbors - par$beta0)
    )
    sigma_star <- as.numeric(
      cov_star -
        cov_star_neighbors %*%
        inv_cov_neighbors %*%
        t(cov_star_neighbors)
    )
    
    predictive_results$predicted_mean[i] <- u_star
    predictive_results$predicted_variance[i] <- max(sigma_star, 1e-10)
  }
  
  alpha <- 1 - interval_level
  pred_sd <- sqrt(predictive_results$predicted_variance)
  predictive_results$CI_lower <- stats::qnorm(alpha / 2, predictive_results$predicted_mean, pred_sd)
  predictive_results$CI_upper <- stats::qnorm(1 - alpha / 2, predictive_results$predicted_mean, pred_sd)
  predictive_results$interval_length <- predictive_results$CI_upper - predictive_results$CI_lower
  
  predictive_results
}

plot_vgp_1d <- function(test_data,
                        predictions,
                        threshold_C,
                        train_data = NULL,
                        f = true_f_1d,
                        title = "Fitted GP Outcome Based on Training Data",
                        model_label = "vGP") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Please install ggplot2 to plot.", call. = FALSE)
  }
  
  vgp <- test_data
  vgp$y_pred <- predictions$predicted_mean
  vgp$pred_var <- predictions$predicted_variance
  vgp$CI_lower <- predictions$CI_lower
  vgp$CI_upper <- predictions$CI_upper
  vgp$interval_length <- predictions$interval_length
  vgp$true_f <- f(vgp$x1)
  
  ci_label <- paste(model_label, "95% CI")
  show_threshold <- !is.null(threshold_C) &&
    length(threshold_C) == 1L &&
    is.finite(threshold_C)
  
  p <- ggplot2::ggplot(vgp, ggplot2::aes(x = x1)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = CI_lower, ymax = CI_upper, fill = ci_label),
      alpha = 0.18
    )
  
  if (!is.null(train_data)) {
    p <- p +
      ggplot2::geom_point(
        data = train_data,
        ggplot2::aes(x = x1, y = y),
        color = "grey60",
        alpha = 0.5,
        size = 1.5,
        inherit.aes = FALSE
      )
  }
  
  p +
    ggplot2::geom_line(ggplot2::aes(y = true_f, color = "true f"), linewidth = 1) +
    ggplot2::geom_line(ggplot2::aes(y = y_pred, color = model_label), linewidth = 1) +
    ggplot2::labs(
      title = title,
      #subtitle = if (show_threshold) sprintf("C = %.3f", threshold_C) else NULL,
      x = "x",
      y = "y",
      color = NULL,
      fill = NULL
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::scale_color_manual(
      values = stats::setNames(c("#F8766D", "#00BFC4"), c("true f", model_label)),
      breaks = c("true f", model_label)
    ) +
    ggplot2::scale_fill_manual(
      values = stats::setNames("#F8766D33", ci_label)
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(order = 1),
      fill = ggplot2::guide_legend(order = 2)
    ) +
    {
      if (show_threshold) {
        ggplot2::geom_hline(yintercept = threshold_C, linetype = "dashed")
      }
    }
}

run_vgp_case <- function(case = c("raw_all", "remove_censored", "threshold_all"),
                         n = 800L,
                         threshold_C = NULL,
                         censor_quantile = NULL,
                         init_params = NULL,
                         lower_bounds = NULL,
                         upper_bounds = NULL,
                         max_marker_neighbors = 20L,
                         n_test_grid = 100L,
                         use_parallel = TRUE,
                         n_cores = max(1L, parallel::detectCores() - 1L),
                         optim_control = list(),
                         make_plot = TRUE) {
  case <- match.arg(case)
  dat <- prepare_vgp_case(
    case = case,
    n = n,
    threshold_C = threshold_C,
    censor_quantile = censor_quantile,
    n_test_grid = n_test_grid
  )
  
  if (is.null(init_params)) {
    init_params <- vgp_case_initial_params(case)
  }
  if (is.null(lower_bounds) || is.null(upper_bounds)) {
    bounds <- vgp_case_bounds(case)
    if (is.null(lower_bounds)) {
      lower_bounds <- bounds$lower
    }
    if (is.null(upper_bounds)) {
      upper_bounds <- bounds$upper
    }
  }
  
  fit <- fit_vgp_model(
    train_data = dat$train_data,
    x_cols = dat$x_cols,
    init_params = init_params,
    lower_bounds = lower_bounds,
    upper_bounds = upper_bounds,
    max_marker_neighbors = max_marker_neighbors,
    use_parallel = use_parallel,
    n_cores = n_cores,
    optim_control = optim_control
  )
  
  predictions <- predict_vgp(
    fit = fit,
    test_data = dat$test_data,
    max_marker_neighbors = max_marker_neighbors
  )
  
  p <- NULL
  if (make_plot) {
    p <- plot_vgp_1d(
      test_data = dat$test_data,
      predictions = predictions,
      threshold_C = dat$threshold_C,
      train_data = dat$train_data,
      f = dat$f,
      title = paste("VGP:", case),
      model_label = "vGP"
    )
  }
  
  list(
    case = case,
    data = dat,
    fit = fit,
    predictions = predictions,
    plot = p
  )
}

fit_vgp_raw_all <- function(init_params = vgp_case_initial_params("raw_all"),
                            n = 800L,
                            threshold_C = NULL,
                            censor_quantile = NULL,
                            lower_bounds = vgp_case_bounds("raw_all")$lower,
                            upper_bounds = vgp_case_bounds("raw_all")$upper,
                            max_marker_neighbors = 20L,
                            n_test_grid = 100L,
                            use_parallel = TRUE,
                            n_cores = max(1L, parallel::detectCores() - 1L),
                            optim_control = list(maxit = 200),
                            make_plot = TRUE) {
  run_vgp_case(
    case = "raw_all",
    n = n,
    threshold_C = threshold_C,
    censor_quantile = censor_quantile,
    init_params = init_params,
    lower_bounds = lower_bounds,
    upper_bounds = upper_bounds,
    max_marker_neighbors = max_marker_neighbors,
    n_test_grid = n_test_grid,
    use_parallel = use_parallel,
    n_cores = n_cores,
    optim_control = optim_control,
    make_plot = make_plot
  )
}

fit_vgp_remove_censored <- function(init_params = vgp_case_initial_params("remove_censored"),
                                    n = 800L,
                                    threshold_C = NULL,
                                    censor_quantile = NULL,
                                    lower_bounds = vgp_case_bounds("remove_censored")$lower,
                                    upper_bounds = vgp_case_bounds("remove_censored")$upper,
                                    max_marker_neighbors = 20L,
                                    n_test_grid = 100L,
                                    use_parallel = TRUE,
                                    n_cores = max(1L, parallel::detectCores() - 1L),
                                    optim_control = list(maxit = 200),
                                    make_plot = TRUE) {
  run_vgp_case(
    case = "remove_censored",
    n = n,
    threshold_C = threshold_C,
    censor_quantile = censor_quantile,
    init_params = init_params,
    lower_bounds = lower_bounds,
    upper_bounds = upper_bounds,
    max_marker_neighbors = max_marker_neighbors,
    n_test_grid = n_test_grid,
    use_parallel = use_parallel,
    n_cores = n_cores,
    optim_control = optim_control,
    make_plot = make_plot
  )
}

fit_vgp_threshold_all <- function(init_params = vgp_case_initial_params("threshold_all"),
                                  n = 800L,
                                  threshold_C = NULL,
                                  censor_quantile = NULL,
                                  lower_bounds = vgp_case_bounds("threshold_all")$lower,
                                  upper_bounds = vgp_case_bounds("threshold_all")$upper,
                                  max_marker_neighbors = 20L,
                                  n_test_grid = 100L,
                                  use_parallel = TRUE,
                                  n_cores = max(1L, parallel::detectCores() - 1L),
                                  optim_control = list(maxit = 200),
                                  make_plot = TRUE) {
  run_vgp_case(
    case = "threshold_all",
    n = n,
    threshold_C = threshold_C,
    censor_quantile = censor_quantile,
    init_params = init_params,
    lower_bounds = lower_bounds,
    upper_bounds = upper_bounds,
    max_marker_neighbors = max_marker_neighbors,
    n_test_grid = n_test_grid,
    use_parallel = use_parallel,
    n_cores = n_cores,
    optim_control = optim_control,
    make_plot = make_plot
  )
}

run_vgp_comparison <- function(cases = c("raw_all", "remove_censored", "threshold_all"),
                               ...) {
  out <- lapply(cases, function(case) run_vgp_case(case = case, ...))
  names(out) <- cases
  out
}

is_vgp_result <- function(x) {
  is.list(x) && all(c("case", "data", "fit", "predictions") %in% names(x))
}

get_vgp_plot <- function(vgp_result,
                         model_label = "vGP",
                         rebuild = TRUE) {
  if (!is_vgp_result(vgp_result)) {
    stop("vgp_result must be returned by run_vgp_case() or one of the fit_vgp_*() functions.", call. = FALSE)
  }
  
  if (!rebuild && !is.null(vgp_result$plot)) {
    return(vgp_result$plot)
  }
  
  plot_vgp_1d(
    test_data = vgp_result$data$test_data,
    predictions = vgp_result$predictions,
    threshold_C = vgp_result$data$threshold_C,
    train_data = vgp_result$data$train_data,
    f = vgp_result$data$f,
    title = paste("VGP:", vgp_result$case),
    model_label = model_label
  )
}

print_vgp_plots <- function(results) {
  if (is_vgp_result(results)) {
    results <- list(results)
  }
  
  invisible(lapply(results, function(vgp_result) {
    print(get_vgp_plot(vgp_result))
  }))
}

save_vgp_plots <- function(results,
                           out_dir = "vgp_plots",
                           width = 7,
                           height = 5,
                           dpi = 300) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Please install ggplot2 to save plots.", call. = FALSE)
  }
  if (is_vgp_result(results)) {
    results <- list(results)
    names(results) <- results[[1]]$case
  }
  if (is.null(names(results)) || any(names(results) == "")) {
    names(results) <- vapply(results, function(x) x$case, character(1))
  }
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  out_files <- file.path(out_dir, paste0(names(results), ".png"))
  
  for (i in seq_along(results)) {
    ggplot2::ggsave(
      filename = out_files[i],
      plot = get_vgp_plot(results[[i]]),
      width = width,
      height = height,
      dpi = dpi
    )
  }
  
  out_files
}

print.vgp_fit <- function(x, ...) {
  cat("VGP fit\n")
  cat("Dimension:", length(x$x_cols), "\n")
  cat("Predictors:", paste(x$x_cols, collapse = ", "), "\n")
  cat("Estimated parameters:\n")
  print(x$estimated_params)
  cat("Final log-likelihood:", x$final_log_likelihood, "\n")
  cat("Counts:\n")
  print(x$result$counts)
  cat("Used parallel:", x$used_parallel, "\n")
  cat("Optim elapsed seconds:", x$optim_time["elapsed"], "\n")
  invisible(x)
}

## =====================================================
## Example usage
## =====================================================
##
## source("VGP_compare.R")
##
## ## Check the three different starting values.
## vgp_case_initial_params("raw_all")
## vgp_case_initial_params("remove_censored")
## vgp_case_initial_params("threshold_all")
##
## ## Fit each case separately, using its own default initial values.
## ## Parallel optimization is the intended full-run setting.
## raw_fit <- fit_vgp_raw_all(
##   init_params = c(1.612907, 2, 0.034895, 0.8218356),
##   lower_bounds = c(1e-6, 1e-6, 1e-6, -Inf),
##   upper_bounds = c(Inf, Inf, Inf, Inf),
##   max_marker_neighbors = 20,
##   use_parallel = TRUE,
##   optim_control = list(maxit = 200)
## )
## raw_fit$fit$estimated_params
## print_vgp_plots(raw_fit)
##
## remove_fit <- fit_vgp_remove_censored(
##   init_params = c(0.3462413, 2, 0.02911056, -0.06377002),
##   lower_bounds = c(1e-6, 1e-6, 1e-6, -Inf),
##   upper_bounds = c(Inf, Inf, Inf, Inf),
##   max_marker_neighbors = 20,
##   use_parallel = TRUE,
##   optim_control = list(maxit = 200)
## )
## remove_fit$fit$estimated_params
##
## threshold_fit <- fit_vgp_threshold_all(
##   init_params = c(0.4797904, 2, 0.01978341, 0.4348702),
##   lower_bounds = c(1e-6, 1e-6, 1e-6, -Inf),
##   upper_bounds = c(Inf, Inf, Inf, Inf),
##   max_marker_neighbors = 20,
##   use_parallel = TRUE,
##   optim_control = list(maxit = 200)
## )
## threshold_fit$fit$estimated_params
##
## ## Or run all three cases. Each case still uses its own default init_params.
## res <- run_vgp_comparison(
##   use_parallel = TRUE,
##   optim_control = list(maxit = 200)
## )
##
## ## Access results:
## res$raw_all$fit
## res$raw_all$predictions
## print_vgp_plots(res$raw_all)
## print_vgp_plots(res)
## save_vgp_plots(res, out_dir = "vgp_plots")
##
## res$remove_censored$fit
## res$threshold_all$fit
##
## ## Fast test:
## small <- run_vgp_comparison(
##   n = 40,
##   ## Use FALSE only for quick debugging or restricted environments.
##   use_parallel = FALSE,
##   optim_control = list(maxit = 1),
##   make_plot = FALSE
## )
