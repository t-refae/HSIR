library(targets)
library(cmdstanr)
library(stantargets)

tar_option_set(
  packages = c(
    "deSolve", "data.table", "dplyr", "tidyr", "ggplot2", "patchwork", "posterior",
    "cmdstanr", "stantargets", "purrr", "tibble", "ggh4x", "ggridges", "scales",
    "viridisLite", "latex2exp", "scico"
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
  ),
  
  #### S/HIT differences and VOI posterior summaries ####
  
  tar_target(
    hit_differences_parent_dir,
    "data/processed/HIT_differences"
  ),
  
  tar_target(
    hit_differences_hs_dir,
    file.path(hit_differences_parent_dir, "HS")
  ),
  
  tar_target(
    hit_differences_homog_dir,
    file.path(hit_differences_parent_dir, "Homog")
  ),
  
  tar_target(
    VOI_thetas,
    1:4
  ),
  
  tar_target(
    VOI_desired_thetas,
    paste0("Theta_", VOI_thetas)
  ),
  
  tar_target(
    VOI_true_params,
    make_voi_true_params(
      true_params = dynamic_speed_comp_dt,
      voi_thetas = VOI_thetas
    )
  ),
  
  tar_target(
    VOI_hs_rdata_files,
    list_voi_rdata_files(
      parent_dir = hit_differences_hs_dir,
      desired_thetas = VOI_desired_thetas
    ),
    format = "file"
  ),
  
  tar_target(
    VOI_homog_rdata_files,
    list_voi_rdata_files(
      parent_dir = hit_differences_homog_dir,
      desired_thetas = VOI_desired_thetas
    ),
    format = "file"
  ),
  
  tar_target(
    VOI_hs_loaded_objects,
    load_voi_rdata_objects(
      rdata_files = VOI_hs_rdata_files,
      desired_thetas = VOI_desired_thetas,
      prefix = "VOI",
      suffix = NULL
    )
  ),
  
  tar_target(
    VOI_homog_loaded_objects,
    load_voi_rdata_objects(
      rdata_files = VOI_homog_rdata_files,
      desired_thetas = VOI_desired_thetas,
      prefix = NULL,
      suffix = "homog"
    )
  ),
  
  tar_target(
    VOI_hs_vars,
    c("beta", "gamma", "cv", "R0")
  ),
  
  tar_target(
    VOI_hs_all_draws,
    make_voi_parameter_draws_df(
      loaded_objects = VOI_hs_loaded_objects,
      desired_thetas = VOI_desired_thetas,
      voi_thetas = VOI_thetas,
      vars = VOI_hs_vars
    )
  ),
  
  tar_target(
    VOI_true_vals_long,
    make_voi_true_values_long(
      true_params = VOI_true_params,
      voi_thetas = VOI_thetas,
      vars = VOI_hs_vars
    )
  ),
  
  tar_target(
    VOI_hs_all_draws_trim,
    trim_voi_parameter_draws(
      draws_df = VOI_hs_all_draws,
      r0_upper = 6.5
    )
  ),
  
  tar_target(
    VOI_ridge_plot,
    plot_voi_parameter_ridges(
      draws_df = VOI_hs_all_draws_trim,
      true_values_df = VOI_true_vals_long
    )
  ),
  
  tar_target(
    VOI_ridge_plot_pdf,
    save_ggplot_pdf_plain(
      plot = VOI_ridge_plot,
      path = file.path(manuscript_figures_dir, "VOI_ridge_plot.pdf"),
      width = 12,
      height = 16
    ),
    format = "file"
  ),
  
  tar_target(
    HS_x_stars,
    summarize_hs_x_stars(
      loaded_objects = VOI_hs_loaded_objects,
      voi_thetas = VOI_thetas
    )
  ),
  
  tar_target(
    Homog_x_stars,
    summarize_homog_x_stars(
      loaded_objects = VOI_homog_loaded_objects,
      voi_thetas = VOI_thetas
    )
  ),
  
  tar_target(
    true_HS_x_star,
    calculate_true_hs_x_star(VOI_true_params)
  ),
  
  tar_target(
    S_HS,
    summarize_S_by_model(
      loaded_objects = VOI_hs_loaded_objects,
      voi_thetas = VOI_thetas,
      model_type = "Heterogeneous",
      object_prefix = "VOI",
      object_suffix = NULL
    )
  ),
  
  tar_target(
    S_Homog,
    summarize_S_by_model(
      loaded_objects = VOI_homog_loaded_objects,
      voi_thetas = VOI_thetas,
      model_type = "Homogeneous",
      object_prefix = NULL,
      object_suffix = "homog"
    )
  ),
  
  tar_target(
    S_HS_diff,
    make_S_minus_xstar_df(
      S_summary = S_HS,
      x_stars = HS_x_stars
    )
  ),
  
  tar_target(
    S_Homog_diff,
    make_S_minus_xstar_df(
      S_summary = S_Homog,
      x_stars = Homog_x_stars
    )
  ),
  
  tar_target(
    fig3_df,
    make_fig3_df(
      S_HS_diff = S_HS_diff,
      S_Homog_diff = S_Homog_diff
    )
  ),
  
  tar_target(
    theta_labels_named_fig3,
    c(
      "Theta_1" = "2 GI pre-peak",
      "Theta_2" = "1 GI pre-peak",
      "Theta_3" = "At peak",
      "Theta_4" = "Complete"
    )
  ),
  
  tar_target(
    fig_3,
    plot_fig3_S_minus_xstar(
      fig3_df = fig3_df,
      theta_labels_named = theta_labels_named_fig3
    )
  ),
  
  tar_target(
    fig_3_pdf,
    save_ggplot_pdf(
      plot = fig_3,
      path = file.path(manuscript_figures_dir, "Fig_3_HIT_differences.pdf"),
      width = 12,
      height = 9
    ),
    format = "file"
  ),
  
  #### NPI policy simulations using posterior median fitted parameters ####
  
  tar_target(
    npi_complete_theta,
    4
  ),
  
  tar_target(
    npi_t_max,
    230
  ),
  
  tar_target(
    npi_population_size,
    1e6
  ),
  
  tar_target(
    npi_i0_prop,
    1 / npi_population_size
  ),
  
  tar_target(
    VOI_HS_medians,
    get_voi_hs_parameter_medians(
      loaded_objects = VOI_hs_loaded_objects,
      theta = npi_complete_theta
    )
  ),
  
  tar_target(
    VOI_homog_medians,
    get_voi_homog_parameter_medians(
      loaded_objects = VOI_homog_loaded_objects,
      theta = npi_complete_theta
    )
  ),
  
  tar_target(
    npi_eff_seq,
    seq(0, 0.9, by = 0.05)
  ),
  
  tar_target(
    npi_start_seq,
    seq(0, 50, by = 5)
  ),
  
  tar_target(
    npi_dur_seq,
    seq(30, 180, by = 30)
  ),
  
  tar_target(
    npi_policy_grid,
    make_npi_policy_grid(
      eff_seq = npi_eff_seq,
      start_seq = npi_start_seq,
      dur_seq = npi_dur_seq
    )
  ),
  
  tar_target(
    npi_results_HS,
    simulate_npi_policy_grid(
      policy_grid = npi_policy_grid,
      model = "HS",
      medians = VOI_HS_medians,
      t_max = npi_t_max,
      population_size = npi_population_size,
      i0_prop = npi_i0_prop
    )
  ),
  
  tar_target(
    npi_results_homog,
    simulate_npi_policy_grid(
      policy_grid = npi_policy_grid,
      model = "homog",
      medians = VOI_homog_medians,
      t_max = npi_t_max,
      population_size = npi_population_size,
      i0_prop = npi_i0_prop
    )
  ),
  
  tar_target(
    npi_results_all,
    dplyr::bind_rows(npi_results_HS, npi_results_homog)
  ),
  
  tar_target(
    npi_results_wide,
    make_npi_results_wide(npi_results_all)
  ),
  
  tar_target(
    npi_attack_midpoint,
    calculate_npi_attack_midpoint(npi_results_wide)
  ),
  
  tar_target(
    hm_attack,
    plot_npi_heatmap(
      df = npi_results_wide,
      z = "attack_rate_diff",
      title = "Delta Final Attack Rate (Homog median - HS median)",
      fill_lab = "Delta Attack Rate",
      midpoint = npi_attack_midpoint,
      label_fun = scales::percent_format(accuracy = 1)
    )
  ),
  
  tar_target(
    hm_peak_inc,
    plot_npi_heatmap(
      df = npi_results_wide,
      z = "peak_incidence_diff",
      title = "Delta Peak Incidence (Homog median - HS median)",
      fill_lab = "Delta Peak incidence",
      midpoint = 0,
      label_fun = scales::label_number(big.mark = ",", accuracy = 1)
    )
  ),
  
  tar_target(
    hm_propS,
    plot_npi_heatmap(
      df = npi_results_wide,
      z = "prop_S_end_diff",
      title = "Delta Prop Susceptible at NPI end (Homog median - HS median)",
      fill_lab = "Delta Prop S",
      midpoint = 0,
      label_fun = scales::percent_format(accuracy = 1)
    )
  ),
  
  tar_target(
    hm_attack_pdf,
    save_ggplot_pdf(
      plot = hm_attack,
      path = file.path(
        manuscript_figures_dir,
        "Fig_4_a_attack_rate_heat_map.pdf"
      ),
      width = 9,
      height = 6
    ),
    format = "file"
  ),
  
  tar_target(
    hm_peak_inc_pdf,
    save_ggplot_pdf(
      plot = hm_peak_inc,
      path = file.path(
        manuscript_figures_dir,
        "Fig_4_b_peak_incidence_heat_map.pdf"
      ),
      width = 9,
      height = 6
    ),
    format = "file"
  ),
  
  tar_target(
    hm_propS_pdf,
    save_ggplot_pdf(
      plot = hm_propS,
      path = file.path(
        manuscript_figures_dir,
        "Fig_4_c_prop_susceptible_heat_map.pdf"
      ),
      width = 9,
      height = 6
    ),
    format = "file"
  ),
  
  tar_target(
    npi_tradeoff_start_keep,
    c(5, 10, 15, 20)
  ),
  
  tar_target(
    npi_tradeoff_dur_keep,
    c(30, 60, 90)
  ),
  
  tar_target(
    npi_tradeoff_grid_df,
    make_npi_tradeoff_grid_df(
      results_all = npi_results_all,
      start_keep = npi_tradeoff_start_keep,
      dur_keep = npi_tradeoff_dur_keep
    )
  ),
  
  tar_target(
    p_tradeoff_grid,
    plot_npi_tradeoff_grid(npi_tradeoff_grid_df)
  ),
  
  tar_target(
    p_tradeoff_grid_pdf,
    save_ggplot_pdf(
      plot = p_tradeoff_grid,
      path = file.path(
        manuscript_figures_dir,
        "Fig_4_comparison_attack_vs_susceptible.pdf"
      ),
      width = 9,
      height = 9
    ),
    format = "file"
  ),
  
  #### Post-NPI information gain: CV posteriors by observed generation interval ####
  
  tar_target(
    npi_informing_parent_dir,
    "data/processed/NPI_inform_release"
  ),
  
  tar_target(
    npi_informing_cv_files,
    list_npi_informing_cv_files(
      parent_dir = npi_informing_parent_dir
    ),
    format = "file"
  ),
  
  tar_target(
    npi_informing_cv_draws_df,
    make_npi_informing_cv_draws_df(
      cv_files = npi_informing_cv_files
    )
  ),
  
  tar_target(
    NPI_informing_plot,
    plot_npi_informing_cv_posteriors(
      cv_draws_df = npi_informing_cv_draws_df
    )
  ),
  
  tar_target(
    NPI_informing_plot_pdf,
    save_ggplot_pdf(
      plot = NPI_informing_plot,
      path = file.path(
        manuscript_figures_dir,
        "CV_posteriors_post_NPI_to_inform_release.pdf"
      ),
      width = 6,
      height = 9
    ),
    format = "file"
  )

  
)