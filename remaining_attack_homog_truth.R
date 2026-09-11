## Standalone complement to Fig_3_alt_remaining_attack_rate
## where the data-generating outbreak is a homogeneous SIR epidemic (nu = 0)

library(cmdstanr)
library(posterior)

source("R/functions.R")
source("../HSIR_fitting/R/VOI.R")

stan_sir  <- "../HSIR_fitting/Stan/VOI_SIR.stan"
stan_hsir <- "../HSIR_fitting/Stan/VOI_HSIR.stan"
out_pdf   <- "outputs/SI_Fig_remaining_attack_homog_truth.pdf"

truth_beta  <- 1.2
truth_gamma <- 0.4
truth_cv    <- 0

P             <- 1e4
i0            <- 1e-4
seed          <- 123
max_days      <- 365
end_threshold <- 1

chains          <- 4
parallel_chains <- 4
iter_warmup     <- 2000
iter_sampling   <- 2000
refresh         <- 100
thin_draws      <- NULL

regimes      <- c("2GI", "1GI", "0GI", "full")
regime_theta <- c("2GI" = 1L, "1GI" = 2L, "0GI" = 3L, "full" = 4L)

theta_labels <- c(
  "Theta_1" = "2 GI pre-peak",
  "Theta_2" = "1 GI pre-peak",
  "Theta_3" = "At peak",
  "Theta_4" = "Complete"
)

mod_sir  <- cmdstan_model(stan_sir)
mod_hsir <- cmdstan_model(stan_hsir)

fit_one <- function(mod, regime) {
  params <- list(
    beta = truth_beta, gamma = truth_gamma, cv = truth_cv,
    regime = regime, P = P, i0 = i0, seed = seed,
    max_days = max_days, end_threshold = end_threshold
  )
  sim <- voi_simulate_truth(params)
  stan_data <- voi_build_stan_data(sim, params)
  fit <- mod$sample(
    data = stan_data,
    seed = seed,
    chains = chains,
    parallel_chains = parallel_chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    refresh = refresh
  )
  s <- fit$summary()
  d <- fit$diagnostic_summary(quiet = TRUE)
  message(sprintf(
    "  %-4s | max rhat %.3f | min ess_bulk %.0f | divergences %d",
    regime,
    max(s$rhat, na.rm = TRUE),
    min(s$ess_bulk, na.rm = TRUE),
    sum(d$num_divergent)
  ))
  as_draws_df(fit$draws())
}

message("Fitting HSIR model to homogeneous-truth data:")
draws_hs <- lapply(regimes, function(r) fit_one(mod_hsir, r))
message("Fitting SIR model to homogeneous-truth data:")
draws_homog <- lapply(regimes, function(r) fit_one(mod_sir, r))
names(draws_hs) <- names(draws_homog) <- regimes

keys_hs    <- paste0("HS_", unname(regime_theta[regimes]))
keys_homog <- paste0("Homog_", unname(regime_theta[regimes]))

spec_ht <- rbind(
  data.frame(arm = "HS",    theta = unname(regime_theta[regimes]), key = keys_hs),
  data.frame(arm = "Homog", theta = unname(regime_theta[regimes]), key = keys_homog)
)

y_array <- function(d) as_draws_array(subset_draws(d, variable = "y"))
states_ht <- c(
  stats::setNames(lapply(draws_hs, y_array), keys_hs),
  stats::setNames(lapply(draws_homog, y_array), keys_homog)
)

param_df <- function(d, arm, theta) {
  data.frame(
    arm = arm,
    theta = theta,
    beta = d$beta,
    gamma = d$gamma,
    cv = if ("cv" %in% names(d)) d$cv else NA_real_,
    R0 = d$R0
  )
}
params_ht <- rbind(
  do.call(rbind, Map(param_df, draws_hs, "HS", unname(regime_theta[regimes]))),
  do.call(rbind, Map(param_df, draws_homog, "Homog", unname(regime_theta[regimes])))
)

remaining_HS_ht <- summarize_remaining_attack_by_model(
  states = states_ht,
  params = params_ht,
  spec = spec_ht,
  arm = "HS",
  model_type = "Heterogeneous",
  probs = c(0.025, 0.975),
  t_extend = 5000,
  thin = thin_draws
)

remaining_Homog_ht <- summarize_remaining_attack_by_model(
  states = states_ht,
  params = params_ht,
  spec = spec_ht,
  arm = "Homog",
  model_type = "Homogeneous",
  probs = c(0.025, 0.975),
  t_extend = 5000,
  thin = thin_draws
)

fig3_remaining_df_ht <- make_fig3_remaining_df(
  remaining_HS = remaining_HS_ht,
  remaining_Homog = remaining_Homog_ht
)

true_params_ht <- data.frame(
  Theta_num = unname(regime_theta[regimes]),
  beta = truth_beta,
  gamma = truth_gamma,
  cv = truth_cv
)

remaining_attack_true_ht <- compute_true_remaining_attack(
  true_params = true_params_ht,
  remaining_df = fig3_remaining_df_ht,
  i0 = i0
)

fig_3_remaining_ht <- plot_fig3_remaining_attack(
  fig3_df = fig3_remaining_df_ht,
  theta_labels_named = theta_labels,
  true_df = remaining_attack_true_ht
)

save_ggplot_pdf(
  plot = fig_3_remaining_ht,
  path = out_pdf,
  width = 12,
  height = 9
)

message("Figure written to: ", out_pdf)

dir.create("outputs", showWarnings = FALSE, recursive = TRUE)
write.csv(
  fig3_remaining_df_ht,
  file.path("outputs", "source_data_SuppFig3_homog_truth.csv"),
  row.names = FALSE
)
