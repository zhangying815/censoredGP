
library(patchwork)
library(ggplot2)
library(grid)

setwd("~/Desktop/censoredGP")
source("VGP_fit_pred_functions_1d.R")
source("cenGP_fit_pred_functions.R")

## vGP-Oracle 
raw_fit <- fit_vgp_raw_all(
  init_params = c(0.7950512,5, 0.034895,0.5512484),
  lower_bounds = c(1e-6, 1e-6, 1e-6, -Inf),
  upper_bounds = c(Inf, Inf, Inf, Inf),
  max_marker_neighbors = 20,
  use_parallel = TRUE,
  optim_control = list(maxit = 200)
)
raw_fit$fit$estimated_params
raw_fit$fit$optim_time
raw_fit$fit$total_seconds
print_vgp_plots(raw_fit)



## vGP-Fixed 
threshold_fit <- fit_vgp_threshold_all(
  init_params = c(0.4797904, 2, 0.01978341, 0.4348702),
  lower_bounds = c(1e-6, 1e-6, 1e-6, -Inf),
  upper_bounds = c(Inf, Inf, Inf, Inf),
  max_marker_neighbors = 20,
  use_parallel = TRUE,
  optim_control = list(maxit = 200)
)
threshold_fit$fit$estimated_params
threshold_fit$fit$optim_time
threshold_fit$fit$total_seconds
print_vgp_plots(threshold_fit)



## vGP-Removed 
remove_fit <- fit_vgp_remove_censored(
  init_params = c(0.3462413, 2, 0.02911056, -0.06377002),
  lower_bounds = c(1e-6, 1e-6, 1e-6, -Inf),
  upper_bounds = c(Inf, Inf, Inf, Inf),
  max_marker_neighbors = 20,
  use_parallel = TRUE,
  optim_control = list(maxit = 200)
)
print(remove_fit$fit)
print_vgp_plots(remove_fit)




## cvGP-SOV
## 1D censored data generate
dat1 <- simulate_1d_censored_data(split_method = "id")
train_data <- dat1$train_data
#test_data  <- dat1$test_data
C <- dat1$threshold_C
## fit cvGP model
fit1 <- fit_censored_gp_nn(
  train_data = dat1$train_data,
  x_cols = "x1",
  threshold_C = dat1$threshold_C,
  k = 20,
  censor_method = "sov", 
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
    #title = "cvGP",
    #subtitle = sprintf("C = %.3f", C),
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






## dashboard plot
common_theme <- theme(
  plot.title = element_text(size = 18),
  legend.position = "right",
  legend.title = element_blank(),
  legend.text = element_text(size = 11),
  legend.key.size = unit(0.65, "cm"),
  legend.box.margin = margin(0, 0, 0, 0),
  plot.margin = margin(5.5, 5.5, 5.5, 5.5)
)

p_raw <- get_vgp_plot(raw_fit) +
  labs(title = "vGP-Oracle") +
  coord_cartesian(xlim = c(0, 1), ylim = c(-2, 5)) +
  common_theme
p_fixed <- get_vgp_plot(threshold_fit) +
  labs(title = "vGP-Fixed") +
  coord_cartesian(xlim = c(0, 1), ylim = c(-2, 5)) +
  common_theme
p_removed <- get_vgp_plot(remove_fit) +
  labs(title = "vGP-Removed") +
  coord_cartesian(xlim = c(0, 1), ylim = c(-2, 5)) +
  common_theme
p_cvgp <- p4 +
  labs(title = "cvGP") +
  coord_cartesian(xlim = c(0, 1), ylim = c(-2, 5)) +
  common_theme

dashboard <- (p_raw | p_fixed) /
  (p_removed | p_cvgp) +
  plot_layout(
    guides = "keep",
    widths = c(1, 1),
    heights = c(1, 1)
  )
print(dashboard)


