library(targets)
library(cmdstanr)
library(stantargets)

tar_option_set(
  packages = c(
    "deSolve", "data.table", "dplyr", "tidyr", "ggplot2", "patchwork", "posterior",
    "cmdstanr", "stantargets", "purrr", "tibble", "ggh4x", "ggridges", "scales",
    "viridisLite"
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
  ),
  
  #### Inference-based identifiability analysis of cv as a function of speed of dynamics ####
  tar_target(
    dynamic_speed_comp_parent_dir,
    "data/processed/dynamic_speed_comp"
  ),
  
  
  tar_target(
    dynamic_speed_comp_dt,
    make_dynamic_speed_comp_dt()
  ),
  
  tar_target(
    dynamic_speed_comp_model,
    "HS"
  ),
  
  tar_target(
    dynamic_speed_comp_LL,
    "pois"
  ),

  tar_target(
    dynamic_speed_comp_TS_type,
    "PTS"
  ),

  tar_target(
    dynamic_speed_comp_theta,
    1:12
  ),

  tar_target(
    dynamic_speed_window_levels,
    c("-2 GI pre-peak", "-1 GI pre-peak", "At peak", "Complete")
  ),

  tar_target(
    dynamic_speed_setting_levels,
    c("Fast", "Reference", "Slow")
  ),

  tar_target(
    dynamic_speed_desired_thetas,
    paste(
      "Theta",
      dynamic_speed_comp_theta,
      dynamic_speed_comp_LL,
      dynamic_speed_comp_TS_type,
      sep = "_"
    )
  ),

  tar_target(
    dynamic_speed_matched_folders,
    find_dynamic_speed_matched_folders(
      parent_dir = dynamic_speed_comp_parent_dir,
      desired_thetas = dynamic_speed_desired_thetas
    )
  ),

  tar_target(
    dynamic_speed_rdata_files,
    list_dynamic_speed_rdata_files(
      parent_dir = dynamic_speed_comp_parent_dir,
      matched_folders = dynamic_speed_matched_folders
    ),
    format = "file"
  ),

  tar_target(
    dynamic_speed_loaded_objects,
    load_dynamic_speed_rdata_objects(
      rdata_files = dynamic_speed_rdata_files,
      desired_thetas = dynamic_speed_desired_thetas
    )
  ),

  tar_target(
    dynamic_speed_true_params,
    make_dynamic_speed_true_params(
      dt = dynamic_speed_comp_dt,
      theta = dynamic_speed_comp_theta
    )
  ),

  tar_target(
    dynamic_speed_cv_draws_df,
    make_dynamic_speed_cv_draws_df(
      loaded_objects = dynamic_speed_loaded_objects,
      desired_thetas = dynamic_speed_desired_thetas,
      window_levels = dynamic_speed_window_levels,
      setting_levels = dynamic_speed_setting_levels
    )
  ),

  tar_target(
    dynamic_speed_true_cv_df,
    make_dynamic_speed_true_cv_df(
      true_params = dynamic_speed_true_params,
      window_levels = dynamic_speed_window_levels,
      setting_levels = dynamic_speed_setting_levels
    )
  ),

  tar_target(
    dynamic_speed_p_cv,
    plot_dynamic_speed_cv_density_matrix(
      cv_draws_df = dynamic_speed_cv_draws_df,
      true_cv_df = dynamic_speed_true_cv_df
    )
  ),

  tar_target(
    dynamic_speed_cv_ridge_df,
    make_dynamic_speed_cv_ridge_df(
      cv_draws_df = dynamic_speed_cv_draws_df,
      setting_levels = dynamic_speed_setting_levels,
      window_levels = dynamic_speed_window_levels
    )
  ),

  tar_target(
    dynamic_speed_ridge_plot_cv_4rows,
    plot_dynamic_speed_cv_ridges_by_window(
      cv_ridge_df = dynamic_speed_cv_ridge_df
    )
  ),

  tar_target(
    dynamic_speed_p_cv_pdf,
    save_ggplot_pdf(
      plot = dynamic_speed_p_cv,
      path = file.path(
        "outputs",
        "dynamic_speed_comp_cv_plot_matrix.pdf"
      ),
      width = 12,
      height = 9
    ),
    format = "file"
  ),

  tar_target(
    dynamic_speed_ridge_plot_cv_4rows_pdf,
    save_ggplot_pdf(
      plot = dynamic_speed_ridge_plot_cv_4rows,
      path = file.path(
        "outputs",
        "Fig_2_dynamic_speed_comp_ridge_plot.pdf"
      ),
      width = 12,
      height = 9
    ),
    format = "file"
  )

  
)