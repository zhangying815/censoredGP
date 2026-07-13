## =====================================================
## C++ MC-SOV usage helper
##
## Run this after source("cenGP_fit_pred_functions.R").
## It compiles mc_sov_censored_nd_cpp() and creates a wrapper
## that can be passed to fit_censored_gp_nn(custom_censor_prob = ...).
## =====================================================
setwd("~/Desktop/censoredGP")

if (!requireNamespace("Rcpp", quietly = TRUE)) {
  stop("Please install Rcpp to use the C++ MC-SOV implementation.", call. = FALSE)
}

cpp_file <- "mc_sov_censored_nd.cpp"
if (!file.exists(cpp_file)) {
  cpp_file <- "mc_sov_censored_nd.cpp"
}

## On some macOS/R installations, sourceCpp() may fail with
## "fatal error: 'cmath' file not found" because R does not pass the
## Command Line Tools SDK include path to clang++.
sdk_path <- tryCatch(system2("xcrun", "--show-sdk-path", stdout = TRUE, stderr = FALSE), error = function(e) "")
restore_pkg_cppflags <- FALSE
old_pkg_cppflags <- Sys.getenv("PKG_CPPFLAGS", unset = "")
if (length(sdk_path) > 0L && nzchar(sdk_path[1])) {
  sdk_path <- sdk_path[1]
  cxx_include <- file.path(sdk_path, "usr", "include", "c++", "v1")
  sdk_flags <- sprintf("-isysroot %s -I%s", shQuote(sdk_path), shQuote(cxx_include))
  Sys.setenv(PKG_CPPFLAGS = paste(old_pkg_cppflags, sdk_flags))
  restore_pkg_cppflags <- TRUE
}

if (restore_pkg_cppflags) {
  tryCatch(
    Rcpp::sourceCpp(cpp_file),
    finally = Sys.setenv(PKG_CPPFLAGS = old_pkg_cppflags)
  )
} else {
  Rcpp::sourceCpp(cpp_file)
}

mc_sov_censor_prob_cpp <- function(mu,
                                   Sigma,
                                   C,
                                   mc_sov_M = 5000L,
                                   mc_sov_seed = NULL,
                                   ...) {
  if (!is.null(mc_sov_seed)) {
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
    set.seed(mc_sov_seed)
  }
  
  mc_sov_censored_nd_cpp(
    mu = as.numeric(mu),
    Sigma = as.matrix(Sigma),
    C = C,
    M = as.integer(mc_sov_M)
  )
}

## Example:
##
## source("flexible_censored_gp_nn.R")
## source("mc_sov_cpp_usage.R")
##
## fit8_mc_sov_cpp <- fit_censored_gp_nn(
##   train_data = dat8$train_data,
##   x_cols = dat8$x_cols,
##   threshold_C = dat8$threshold_C,
##   k = 20,
##   custom_censor_prob = mc_sov_censor_prob_cpp,
##   censor_prob_args = list(
##     mc_sov_M = 5000,
##     mc_sov_seed = 123
##   ),
##   init_params = init8,
##   lower_bounds = c(1e-6, rep(1e-6, 8), 1e-6, -Inf),
##   upper_bounds = c(Inf, rep(Inf, 8), Inf, Inf),
##   use_parallel = FALSE
## )
