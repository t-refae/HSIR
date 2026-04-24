library(targets)
library(cmdstanr)
library(stantargets)

tar_option_set(
  packages = c(
    "deSolve", "data.table", "dplyr", "tidyr", "ggplot2", "patchwork", "posterior",
    "cmdstanr", "stantargets"
  ),
  format = "rds"
)

tar_source()

list(
  #### SIR vs. HSIR model comparisons ####
  
  ## shared params
  tar_target(shared_tps, seq(1, 60)),
  tar_target(proportion_initially_infected, 1e-4),
  tar_target(population_size, 1e4),
  
  tar_target(
    shared_init_state,
    c(
      S = 1 - proportion_initially_infected,
      I = proportion_initially_infected,
      R = 0,
      C = 0
    )
  ),
  
  tar_target(homog_params, c(beta = 0.6, gamma = 0.2)),
  tar_target(HS_params, c(beta = 0.6, gamma = 0.2, cv = 1)),
  
  # simulate trajectories
  tar_target(
    homog_out,
    simulate_sir(
      init_state = shared_init_state,
      times = shared_tps,
      params = homog_params,
      model = "homogeneous"
    )
  ),
  
  tar_target(
    HS_out,
    simulate_sir(
      init_state = shared_init_state,
      times = shared_tps,
      params = HS_params,
      model = "heterogeneous"
    )
  ),
  
  tar_target(
    combined_data,
    combine_trajectory_outputs(homog_out, HS_out)
  ),
  
  tar_target(
    combined_long,
    make_combined_long(combined_data)
  ),
  
  tar_target(
    delta_df,
    make_delta_df(homog_out, HS_out)
  ),
  
  ## Fig. 1 trajectory panels
  tar_target(
    fig_1_a,
    plot_trajectories_faceted(combined_long)
  ),
  
  tar_target(
    fig_1_b,
    plot_absolute_differences(delta_df)
  ),
  

  ## generate model fits

  ## MCMC params
  tar_target(ncores, 4),
  tar_target(niter_in, 1000),
  tar_target(warmup_iter_in, round(niter_in / 10)),
  tar_target(n_chains_in, ncores),

  tar_target(stan_t0, 0.9999999999),

  tar_target(
    stan_ts,
    seq(1, length(shared_tps), 1)
  ),

  tar_target(
    my_ndays,
    length(shared_tps)
  ),

  # fitting the HS model to homogeneous-generated data.
  tar_target(
    homog_cases_stan,
    make_stan_cases_from_cumulative(
      out = homog_out,
      population_size = population_size
    )
  ),

  tar_target(
    stan_data_hs_fit_to_homog,
    make_stan_data(
      n_days = my_ndays,
      y0 = shared_init_state,
      t0 = stan_t0,
      ts = stan_ts,
      population_size = population_size,
      cases = homog_cases_stan
    )
  ),
  
  tar_stan_mcmc(
    name = hs_fit_to_homog,
    stan_files = c(hs = "stan/HetSus.stan"),
    data = stan_data_hs_fit_to_homog,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 100,
    iter_sampling = 900,
    seed = 2,
    refresh = 100,
    return_draws = TRUE,
    return_summary = TRUE,
    return_diagnostics = TRUE
  ),

  ## fitting the homogeneous model to HS-generated data
  tar_target(
    HS_cases_stan,
    make_stan_cases_from_cumulative(
      out = HS_out,
      population_size = population_size
    )
  ),

  tar_target(
    stan_data_homog_fit_to_hs,
    make_stan_data(
      n_days = my_ndays,
      y0 = shared_init_state,
      t0 = stan_t0,
      ts = stan_ts,
      population_size = population_size,
      cases = HS_cases_stan
    )
  ),

  tar_stan_mcmc(
    name = homog_fit_to_hs,
    stan_files = c(homog = "stan/SIR_homog.stan"),
    data = stan_data_homog_fit_to_hs,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 100,
    iter_sampling = 900,
    seed = 2,
    refresh = 100,
    return_draws = TRUE,
    return_summary = TRUE,
    return_diagnostics = TRUE
  ),

  ## posterior predictive checks

  # HS model fit to homogeneous data
  tar_target(
    posterior_dt_1,
    extract_posterior_draws(
      fit = hs_fit_to_homog_mcmc_hs,
      variables = c("beta", "gamma", "cv"),
      n_draws = 500,
      seed = 2
    )
  ),

  tar_target(
    posterior_incidence_1,
    simulate_posterior_incidence(
      posterior_dt = posterior_dt_1,
      init_state = shared_init_state,
      times = shared_tps,
      population_size = population_size,
      model = "heterogeneous"
    )
  ),

  tar_target(
    incidence_summary_1,
    summarize_posterior_incidence(posterior_incidence_1)
  ),

  tar_target(
    plot_fit_1,
    plot_fit_to_data(
      incidence_summary = incidence_summary_1,
      observed_cases = data.frame(Day = shared_tps[-1],
                                  Data=homog_cases_stan),
      ribbon_fill = "#ff7f0e",
      line_color = "#ff7f0e",
      point_color = "#1f77b4",
      line_label = "Heterogeneous Model Fit",
      point_label = "Data (Homogeneous Model)"
    )
  ),

  ## homogeneous model fit to HS data
  tar_target(
    posterior_dt_2,
    extract_posterior_draws(
      fit = homog_fit_to_hs_mcmc_homog,
      variables = c("beta", "gamma"),
      n_draws = 500,
      seed = 2
    )
  ),

  tar_target(
    posterior_incidence_2,
    simulate_posterior_incidence(
      posterior_dt = posterior_dt_2,
      init_state = shared_init_state,
      times = shared_tps,
      population_size = population_size,
      model = "homogeneous"
    )
  ),

  tar_target(
    incidence_summary_2,
    summarize_posterior_incidence(posterior_incidence_2)
  ),

  tar_target(
    plot_fit_2,
    plot_fit_to_data(
      incidence_summary = incidence_summary_2,
      observed_cases = data.frame(Day = shared_tps[-1],
                                  Data=HS_cases_stan),
      ribbon_fill = "#1f77b4",
      line_color = "#1f77b4",
      point_color = "#ff7f0e",
      line_label = "Homogeneous Model Fit",
      point_label = "Data (Heterogeneous Model)"
    )
  ),


  ## Fig. 1
  tar_target(
    fig_1_panel,
    assemble_fig_1_panel(
      fig_1_a = fig_1_a,
      fig_1_b = fig_1_b,
      plot_fit_1 = plot_fit_1,
      plot_fit_2 = plot_fit_2
    )
  ),

  tar_target(
    manuscript_figures_dir,
    {
      dir.create("outputs", showWarnings = FALSE, recursive = TRUE)
      "outputs"
    }
  ),

  tar_target(
    fig_1_pdf,
    save_ggplot_pdf(
      plot = fig_1_panel,
      path = file.path(manuscript_figures_dir, "Fig_1_model_comparisons.pdf"),
      width = 15,
      height = 15
    ),
    format = "file"
  )
)