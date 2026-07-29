
library(Rcpp)
library(dplyr)
library(tidyr)
library(scoringRules)

setwd("~/Desktop/censoredGP")

sdk_path <- system2("xcrun", "--show-sdk-path", stdout = TRUE)[1]
cxx_include <- file.path(sdk_path, "usr", "include", "c++", "v1")
Sys.setenv(
  PKG_CPPFLAGS = paste(
    Sys.getenv("PKG_CPPFLAGS"),
    sprintf("-isysroot %s -I%s", shQuote(sdk_path), shQuote(cxx_include))
  )
)
Rcpp::sourceCpp("mc_sov_censored_nd.cpp")
source("cenGP_fit_pred_functions.R")
source("mc_sov_cpp_usage.R")
source("VGP_fit_pred_functions_1d.R")




## 1D data
dat1 <- simulate_1d_censored_data(split_method = "id")
train_data <- dat1$train_data
test_data  <- dat1$test_data
C <- dat1$threshold_C
## fit cvGP model
fit1 <- fit_censored_gp_nn(
  train_data = dat1$train_data,
  x_cols = "x1",
  threshold_C = dat1$threshold_C,
  k = 20,
  censor_method = "sov", # sov miwa
  init_params = c(0.3462413, 1.5, 0.02911056, 0.01410213),
  lower_bounds = c(1e-6, 1e-6, 1e-6, -Inf),
  upper_bounds = c(Inf, Inf, Inf, Inf)
)
###### avoid the data points is not split in average, let us do set x from 0 to 1 with seq 100
test_data <- make_1d_test_grid(
  n_test = 100,
  threshold_C = C,
  f = default_1d_fun)
## prediction
pred1_cens <- predict_censored_gp_nn_with_censoring(
  fit = fit1,
  test_data = test_data,
  k = 20,
  max_censored_ids = 5,
  distance_method = "euclidean",
  n_pred_samples = 20000,
  sample_method = "TruncatedNormal",
  prediction_seed = 123,
  predict_y = TRUE,
  verbose = TRUE
)
results_df <- pred1_cens$results_df
merged_df  <- pred1_cens$merged_df
p4 <- ggplot(merged_df, aes(x = x1)) +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper, fill = "cvGP 95% PI"),
              alpha = 0.18) +
  geom_point(data = train_data, aes(x = x1, y = y),
             color = "grey60", alpha = 0.6, size = 1.8)+
  # ---- correct mappings for legend labels/colors ----
geom_line(aes(y = true_f_1d(x1),  color = "true f"), linewidth = 1) +   # true f line -> "true f"
  geom_line(aes(y = mean_pred, color = "cvGP"),  linewidth = 1) +    # model line -> "cvGP"
  geom_hline(yintercept = C, linetype = "dashed") +
  labs(
    title = "cvGP",
    subtitle = sprintf("C = %.3f", C),
    x = "x", y = "y", color = NULL, fill = NULL
  ) +
  theme_minimal(base_size = 13) +
scale_color_manual(values = c("true f" = "#F8766D",  
                              "cvGP"   = "#00BFC4"),  
                   breaks = c("true f", "cvGP")) +
  scale_fill_manual(values = c("cvGP 95% PI" = "#F8766D33")) +  # light red ribbon
  guides(
    color = guide_legend(order = 1),  # lines first
    fill  = guide_legend(order = 2)   # ribbon below
  )
print(p4)






## 8D borehole data
dat8 <- simulate_borehole_censored_data(censor_quantile = 0.8)
init8 <- c(28.48859,
           0.9, 0.5, 0.2, 0.46, 0.2, 0.35, 0.2, 0.7,
           1.027467,
           59.97952)
# fit cvGP
fit8_sov <- fit_censored_gp_nn(
  train_data = dat8$train_data,
  x_cols = dat8$x_cols,
  threshold_C = dat8$threshold_C,
  k = 20,
  censor_method = "sov",
  init_params = init8,
  lower_bounds = c(1e-6, rep(1e-6, 8), 1e-6, -Inf),
  upper_bounds = c(Inf, rep(Inf, 8), Inf, Inf),
  use_parallel = TRUE
)
fit8_sov

## do prediction
pred8_cens <- predict_censored_gp_nn_with_censoring(
  fit = fit8_sov,
  test_data = dat8$test_data,
  k = 20,
  max_censored_ids = 5,
  distance_method = "euclidean",
  n_pred_samples = 20000,
  sample_method = "TruncatedNormal",
  prediction_seed = 123,
  predict_y = TRUE,
  verbose = TRUE
)
results_df_sov <- pred8_cens$results_df
merged_df_sov  <- pred8_cens$merged_df





## using mc-sov
fit8_mc_sov_cpp <- fit_censored_gp_nn(
  train_data = dat8$train_data,
  x_cols = dat8$x_cols,
  threshold_C = dat8$threshold_C,
  k = 20,
  censor_method = "mc_sov",
  custom_censor_prob = mc_sov_censor_prob_cpp,
  censor_prob_args = list(
    mc_sov_M = 5000,
    mc_sov_seed = 123
  ),
  init_params = init8,
  lower_bounds = c(1e-6, rep(1e-6, 8), 1e-6, -Inf),
  upper_bounds = c(Inf, rep(Inf, 8), Inf, Inf),
  use_parallel = FALSE
)
fit8_mc_sov_cpp

## do prediction
pred8_cens_mc_sov <- predict_censored_gp_nn_with_censoring(
  fit = fit8_mc_sov_cpp,
  test_data = dat8$test_data,
  k = 20,
  max_censored_ids = 5,
  distance_method = "euclidean",
  n_pred_samples = 20000,
  sample_method = "TruncatedNormal",
  prediction_seed = 123,
  predict_y = TRUE,
  verbose = TRUE
)
results_df_mv_sov <- pred8_cens_mc_sov$results_df
merged_df_mc_sov  <- pred8_cens_mc_sov$merged_df



## using MET
fit8_met <- fit_censored_gp_nn(
  train_data = dat8$train_data,
  x_cols = dat8$x_cols,
  threshold_C = dat8$threshold_C,
  k = 20,
  censor_method = "met",
  init_params = init8,
  lower_bounds = c(1e-6, rep(1e-6, 8), 1e-6, -Inf),
  upper_bounds = c(Inf, rep(Inf, 8), Inf, Inf),
  use_parallel = TRUE
)
fit8_met

## do prediction
pred8_cens_met <- predict_censored_gp_nn_with_censoring(
  fit = fit8_met,
  test_data = dat8$test_data,
  k = 20,
  max_censored_ids = 5,
  distance_method = "euclidean",
  n_pred_samples = 20000,
  sample_method = "TruncatedNormal",
  prediction_seed = 123,
  predict_y = TRUE,
  verbose = TRUE
)
results_df_met <- pred8_cens_met$results_df
merged_df_met  <- pred8_cens_met$merged_df






############################################################################
############################### Performance ################################
############################################################################
calc_one_group_metrics <- function(df,
                                   group_name,
                                   threshold_C,
                                   truth_col = "y_w",
                                   pred_col = "mean_pred",
                                   var_col = "var_pred",
                                   lower_col = "CI_lower",
                                   upper_col = "CI_upper",
                                   interval_col = "interval_length") {
  y_true <- df[[truth_col]]
  y_pred <- df[[pred_col]]
  var_pred <- df[[var_col]]
  
  mspe <- mean((y_true - y_pred)^2, na.rm = TRUE)
  mae  <- mean(abs(y_true - y_pred), na.rm = TRUE)
  mape <- mean(abs(y_true - y_pred) / abs(y_true), na.rm = TRUE) * 100
  crps_values <- scoringRules::crps_norm(
    y = y_true,
    mean = y_pred,
    sd = sqrt(var_pred)
  )
  mean_crps <- mean(crps_values, na.rm = TRUE)
  avg_interval_length <- mean(df[[interval_col]], na.rm = TRUE)
  coverage <- mean(
    y_true >= df[[lower_col]] & y_true <= df[[upper_col]],
    na.rm = TRUE
  )
  pct_above_threshold <- if (group_name == "Censored") {
    mean(y_pred > threshold_C, na.rm = TRUE) * 100
  } else {
    NA_real_
  }
  data.frame(
    Group = group_name,
    N = nrow(df),
    MSPE = mspe,
    MAE = mae,
    MAPE = mape,
    Mean_CRPS = mean_crps,
    Avg_Interval_Length = avg_interval_length,
    Coverage = coverage,
    Pct_Above_Censored_Threshold = pct_above_threshold
  )
}

calc_cengp_performance <- function(merged_df,
                                   threshold_C,
                                   model_name = "Model",
                                   truth_col = "y_w",
                                   censored_col = "censored",
                                   pred_col = "mean_pred",
                                   var_col = "var_pred",
                                   lower_col = "CI_lower",
                                   upper_col = "CI_upper",
                                   interval_col = "interval_length",
                                   digits = 4) {
  
  test_data0 <- merged_df %>%
    dplyr::filter(.data[[censored_col]] == 0)
  test_data1 <- merged_df %>%
    dplyr::filter(.data[[censored_col]] == 1)
  
  uncensored_metrics <- calc_one_group_metrics(
    df = test_data0,
    group_name = "Uncensored",
    threshold_C = threshold_C,
    truth_col = truth_col,
    pred_col = pred_col,
    var_col = var_col,
    lower_col = lower_col,
    upper_col = upper_col,
    interval_col = interval_col
  )
  
  censored_metrics <- calc_one_group_metrics(
    df = test_data1,
    group_name = "Censored",
    threshold_C = threshold_C,
    truth_col = truth_col,
    pred_col = pred_col,
    var_col = var_col,
    lower_col = lower_col,
    upper_col = upper_col,
    interval_col = interval_col
  )
  summary_long <- bind_rows(uncensored_metrics, censored_metrics) %>%
    mutate(Model = model_name, .before = 1) %>%
    mutate(across(where(is.numeric), ~ round(.x, digits)))
  summary_wide <- summary_long %>%
    select(-Model, -N) %>%
    pivot_longer(
      cols = -Group,
      names_to = "Metric",
      values_to = "Value"
    ) %>%
    pivot_wider(
      names_from = Group,
      values_from = Value
    )
  list(
    model = model_name,
    test_data0 = test_data0,
    test_data1 = test_data1,
    summary_long = summary_long,
    summary_wide = summary_wide
  )
}

## performance of SOV
perf_sov <- calc_cengp_performance(
  merged_df = merged_df_sov,
  threshold_C = dat8$threshold_C,
  model_name = "SOV",
  truth_col = "y_true"
)
perf_sov$summary_wide
perf_sov$summary_long

## performance of MC-SOV
perf_mc_sov <- calc_cengp_performance(
  merged_df = merged_df_mc_sov,
  threshold_C = dat8$threshold_C,
  model_name = "MC-SOV",
  truth_col = "y_true"
)
perf_mc_sov$summary_wide
perf_mc_sov$summary_long

## performance of MET
perf_met <- calc_cengp_performance(
  merged_df = merged_df_met,
  threshold_C = dat8$threshold_C,
  model_name = "MET",
  truth_col = "y_true"
)
perf_met$summary_wide
perf_met$summary_long

## combine results
all_perf <- bind_rows(perf_sov$summary_long, perf_mc_sov$summary_long, perf_met$summary_long)
all_perf
