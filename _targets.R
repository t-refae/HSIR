library(targets)
library(cmdstanr)
library(stantargets)

ggplot2::theme_set(
  ggplot2::theme_minimal(base_size = 8, base_family = "Helvetica")
)
ggplot2::update_geom_defaults("text", list(family = "Helvetica"))
ggplot2::update_geom_defaults("label", list(family = "Helvetica"))

tar_option_set(
  packages = c(
    "deSolve", "data.table", "dplyr", "tidyr", "ggplot2", "patchwork", "posterior",
    "cmdstanr", "stantargets", "purrr", "tibble", "ggh4x", "ggridges", "scales",
    "viridisLite", "latex2exp", "scico", "ggridges", "yaml", "openxlsx"
  ),
  format = "rds"
)

tar_source()

list(
  #### Susceptibility contour plot (introductory figure) ####
  tar_target(fig_0_dist_df, make_susceptibility_dist_df()),
  tar_target(fig_0_mean_sus_df, make_mean_sus_df()),
  tar_target(
    fig_0_contour,
    plot_susceptibility_and_hit(fig_0_dist_df, fig_0_mean_sus_df)
  ),
  tar_target(
    fig_0_contour_pdf,
    save_ggplot_pdf(
      plot = fig_0_contour,
      path = file.path(manuscript_figures_dir, "Fig_0_susceptibility_contour.pdf"),
      width = 11, height = 4.8
    ),
    format = "file"
  ),
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
    plot_trajectory_panel(combined_long, "Cumulative Incidence", "Cumulative incidence")
  ),
  
  tar_target(
    fig_1_b,
    plot_trajectory_panel(combined_long, "Prevalence", "Prevalence")
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

  # fitting the HS model to homogeneous-generated data
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
      observed_cases    = data.frame(Day = shared_tps[-1], Data = homog_cases_stan),
      fit_label  = "Heterogeneous",
      data_label = "Homogeneous"
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
      observed_cases    = data.frame(Day = shared_tps[-1], Data = HS_cases_stan),
      fit_label  = "Homogeneous",
      data_label = "Heterogeneous"
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
      width = 7,
      height = 6.5
    ),
    format = "file"
  ),
  
  #### Inference-based identifiability analysis of cv as a function of speed of dynamics ####
  
  tar_target(pts_fitting_root, "../HSIR_fitting"),
  
  tar_target(pts_spec, pts_fit_spec(pts_fitting_root, cv = 1)),
  
  tar_target(pts_bundle_files, pts_fit_files(pts_spec), format = "file"),
  
  tar_target(pts_cv_draws_df, pts_cv_draws(pts_spec, pts_bundle_files)),
  
  tar_target(pts_true_cv_df, pts_true_cv(pts_spec)),

  tar_target(
    dynamic_speed_p_cv,
    plot_dynamic_speed_cv_density_matrix(
      cv_draws_df = pts_cv_draws_df,
      true_cv_df = pts_true_cv_df
    )
  ),
  
  tar_target(
    dynamic_speed_cv_ridge_df,
    make_dynamic_speed_cv_ridge_df(cv_draws_df = pts_cv_draws_df)
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
  
  tar_target(voi_fitting_root, "../HSIR_fitting"),
  tar_target(voi_spec, voi_fit_spec(voi_fitting_root)),
  tar_target(voi_bundle_files, voi_fit_files(voi_spec), format = "file"),
  tar_target(voi_params, voi_param_draws(voi_spec, voi_bundle_files)),
  tar_target(voi_states, voi_state_draws(voi_spec, voi_bundle_files)),
  
  tar_target(VOI_thetas, sort(unique(voi_spec$theta))),
  tar_target(VOI_true_params, voi_true_params(voi_spec)),
  tar_target(VOI_hs_vars, c("beta", "gamma", "cv", "R0")),
  tar_target(VOI_hs_all_draws, voi_parameter_draws_df(voi_params, VOI_hs_vars)),
  
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
    true_HS_x_star,
    calculate_true_hs_x_star(VOI_true_params)
  ),
  
  tar_target(HS_x_stars,    voi_x_stars(voi_params, "HS")),
  tar_target(Homog_x_stars, voi_x_stars(voi_params, "Homog")),
  
  tar_target(S_HS,    voi_S_summary(voi_states, voi_spec, "HS", "Heterogeneous")),
  tar_target(S_Homog, voi_S_summary(voi_states, voi_spec, "Homog", "Homogeneous")),
  
  tar_target(VOI_HS_medians,    voi_medians(voi_params, "HS", npi_complete_theta)),
  tar_target(VOI_homog_medians, voi_medians(voi_params, "Homog", npi_complete_theta)),
  
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
  
  #### Alternative Fig. 3: remaining attack rate (unmitigated remaining burden) ####
  
  tar_target(
    remaining_attack_HS,
    summarize_remaining_attack_by_model(
      states = voi_states,
      params = voi_params,
      spec = voi_spec,
      arm = "HS",
      model_type = "Heterogeneous",
      probs = c(0.025, 0.975),
      t_extend = 5000,
      thin = NULL
    )
  ),
  
  tar_target(
    remaining_attack_Homog,
    summarize_remaining_attack_by_model(
      states = voi_states,
      params = voi_params,
      spec = voi_spec,
      arm = "Homog",
      model_type = "Homogeneous",
      probs = c(0.025, 0.975),
      t_extend = 5000,
      thin = NULL
    )
  ),
  
  tar_target(
    remaining_attack_true,
    compute_true_remaining_attack(
      true_params = VOI_true_params,
      remaining_df = fig3_remaining_df,
      i0 = 1e-4
    )
  ),
  
  tar_target(
    fig3_remaining_df,
    make_fig3_remaining_df(
      remaining_HS = remaining_attack_HS,
      remaining_Homog = remaining_attack_Homog
    )
  ),
  
  tar_target(
    fig_3_remaining,
    plot_fig3_remaining_attack(
      fig3_df = fig3_remaining_df,
      theta_labels_named = theta_labels_named_fig3,
      true_df = remaining_attack_true
    )
  ),
  
  tar_target(
    fig_3_remaining_pdf,
    save_ggplot_pdf(
      plot = fig_3_remaining,
      path = file.path(manuscript_figures_dir, "Fig_3_alt_remaining_attack_rate.pdf"),
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
    1080
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
    calculate_npi_midpoint(npi_results_wide, z = "attack_rate_diff")
  ),
  
  tar_target(
    hm_attack,
    plot_npi_heatmap(
      df = npi_results_wide,
      z = "attack_rate_diff",
      title = NULL,
      fill_lab = "Final attack rate\n(SIR - hSIR)",
      midpoint = npi_attack_midpoint,
      label_fun = scales::percent_format(accuracy = 1)
    )
  ),
  
  tar_target(
    hm_peak_inc,
    plot_npi_heatmap(
      df = npi_results_wide,
      z = "peak_incidence_diff",
      title = NULL,
      fill_lab = "Peak incidence\n(SIR - hSIR)",
      midpoint = 0,
      label_fun = scales::percent_format(accuracy = 1)
    )
  ),
  
  tar_target(
    hm_propS,
    plot_npi_heatmap(
      df = npi_results_wide,
      z = "prop_S_end_diff",
      title = NULL,
      fill_lab = "Proportion susceptible at NPI end\n(SIR - hSIR)",
      midpoint = 0,
      label_fun = scales::percent_format(accuracy = 1)
    )
  ),
  
  tar_target(
    hm_attack_pdf,
    save_ggplot_pdf(
      plot = hm_attack,
      path = file.path(manuscript_figures_dir, "Fig_4_a_attack_rate_heat_map.pdf"),
      width = 9, height = 6
    ),
    format = "file"
  ),
  
  tar_target(
    hm_peak_inc_pdf,
    save_ggplot_pdf(
      plot = hm_peak_inc,
      path = file.path(manuscript_figures_dir, "Fig_4_b_peak_incidence_heat_map.pdf"),
      width = 9, height = 6
    ),
    format = "file"
  ),
  
  tar_target(
    hm_propS_pdf,
    save_ggplot_pdf(
      plot = hm_propS,
      path = file.path(manuscript_figures_dir, "Fig_4_c_prop_susceptible_heat_map.pdf"),
      width = 9, height = 6
    ),
    format = "file"
  ),
  
  tar_target(
    npi_tradeoff_start_keep,
    c(5, 10, 15)
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
  
  
  #### Alternative Fig. 4: continuous delta in remaining attack rate ####
  
  tar_target(
    npi_remaining_eff_seq,
    c(0.3, 0.5, 0.8)
  ),
  
  tar_target(
    npi_remaining_policy_grid,
    make_npi_policy_grid(
      eff_seq = npi_remaining_eff_seq,
      start_seq = npi_tradeoff_start_keep,
      dur_seq = npi_tradeoff_dur_keep
    )
  ),
  
  tar_target(
    npi_remaining_grid,
    simulate_npi_remaining_grid(
      policy_grid = npi_remaining_policy_grid,
      medians_HS = VOI_HS_medians,
      medians_homog = VOI_homog_medians,
      t_max = npi_t_max,
      population_size = npi_population_size,
      i0_prop = npi_i0_prop
    )
  ),
  
  tar_target(
    fig_4_remaining_by_model,
    plot_npi_remaining_by_model_grid(
      remaining_grid = npi_remaining_grid,
      start_keep = npi_tradeoff_start_keep,
      dur_keep = npi_tradeoff_dur_keep,
      eff_keep = npi_remaining_eff_seq,
      x_max_display = 200
    )
  ),
  
  tar_target(
    fig_4_remaining_by_model_pdf,
    save_ggplot_pdf(
      plot = fig_4_remaining_by_model,
      path = file.path(manuscript_figures_dir, "Fig_4_alt_remaining_attack_by_model.pdf"),
      width = 11,
      height = 9
    ),
    format = "file"
  ),
  
  #### NPI efficiency: marginal value of an extra NPI day ####
  
  tar_target(
    npi_mv_eff_seq,
    c(0.3, 0.5, 0.8)
  ),
  
  tar_target(
    npi_mv_dur_seq,
    seq(0, 150, by = 2)   # fine duration grid; use by = 1 for exact per-day
  ),
  
  tar_target(
    npi_mv_policy_grid,
    make_npi_policy_grid(
      eff_seq = npi_mv_eff_seq,
      start_seq = npi_tradeoff_start_keep,
      dur_seq = npi_mv_dur_seq
    )
  ),
  
  tar_target(
    npi_mv_attack_curve,
    simulate_npi_attack_curve(
      policy_grid = npi_mv_policy_grid,
      medians_HS = VOI_HS_medians,
      medians_homog = VOI_homog_medians,
      t_max = npi_t_max,
      population_size = npi_population_size,
      i0_prop = npi_i0_prop
    )
  ),
  
  tar_target(
    npi_mv_df,
    make_npi_marginal_value_df(
      attack_curve = npi_mv_attack_curve,
      population_size = npi_population_size
    )
  ),
  
  tar_target(
    fig_npi_marginal_value,
    plot_npi_marginal_value(
      marginal_df = npi_mv_df,
      start_keep = npi_tradeoff_start_keep,
      dur_max_display = 150
    )
  ),
  
  tar_target(
    fig_npi_marginal_value_pdf,
    save_ggplot_pdf(
      plot = fig_npi_marginal_value,
      path = file.path(manuscript_figures_dir, "Fig_NPI_marginal_value_per_day.pdf"),
      width = 11,
      height = 8
    ),
    format = "file"
  ),
  
  #### Post-NPI information gain: CV posteriors by observed generation interval ####
  
  tar_target(npi_fitting_root, "../HSIR_fitting"),
  
  tar_target(npi_spec, npi_fit_spec(npi_fitting_root)),
  
  tar_target(npi_bundle_files, npi_fit_files(npi_spec), format = "file"),
  
  tar_target(npi_informing_cv_draws_df, npi_cv_draws(npi_spec, npi_bundle_files)),
  
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
  ),

  #### Supplementary figs ####
  
  # infectious dynamics across settings 
  tar_target(
    I_dynamics_df,
    make_I_dynamics_df(
      grid_csv = FTS_grid_csv,
      population_size = 1e6,
      t_max = 800
    )
  ),
  
  tar_target(
    I_dynamics_plot,
    plot_I_dynamics(I_dynamics_df)
  ),
  
  tar_target(
    I_dynamics_pdf,
    save_ggplot_pdf(
      plot = I_dynamics_plot,
      path = file.path(manuscript_figures_dir, "SI_Fig_I_dynamics.pdf"),
      width = 11,
      height = 8
    ),
    format = "file"
  ),
  
  # full time series (FTS) fits
  tar_target(
    FTS_fit_files, list.files("../HSIR_fitting/outputs/FTS", "\\.rds$", full.names = TRUE),
    format = "file"
  ),
  
  tar_target(FTS_fit_bundles, lapply(FTS_fit_files, readRDS)),
  
  tar_target(FTS_cv_ridge,    param_ridge_plot(FTS_fit_bundles, FTS_fit_files, "cv")),
  tar_target(FTS_beta_ridge,  param_ridge_plot(FTS_fit_bundles, FTS_fit_files, "beta")),
  tar_target(FTS_gamma_ridge, param_ridge_plot(FTS_fit_bundles, FTS_fit_files, "gamma")),
  
  tar_target(
    FTS_ridge_files,
    {
      dir.create("outputs", showWarnings = FALSE)
      out <- c(
        cv    = file.path("outputs", "FTS_cv_ridge.png"),
        beta  = file.path("outputs", "FTS_beta_ridge.png"),
        gamma = file.path("outputs", "FTS_gamma_ridge.png")
      )
      ggplot2::ggsave(out[["cv"]],    FTS_cv_ridge,    width = 7, height = 8, dpi = 300)
      ggplot2::ggsave(out[["beta"]],  FTS_beta_ridge,  width = 7, height = 8, dpi = 300)
      ggplot2::ggsave(out[["gamma"]], FTS_gamma_ridge, width = 7, height = 8, dpi = 300)
      unname(out)
    },
    format = "file"
  ),
  
  tar_target(
    FTS_grid_csv,
    {
      cfg <- yaml::read_yaml("../HSIR_fitting/config.yml")
      file.path("../HSIR_fitting", cfg$scenarios_csv)
    },
    format = "file"
  ),
  
  tar_target(
    FTS_param_ridge_df,
    fts_param_ridge_df(
      bundles = FTS_fit_bundles,
      files = FTS_fit_files,
      grid_csv = FTS_grid_csv,
      params = c("beta", "gamma", "R0", "cv"),
      r0_keep = 3
    )
  ),
  
  tar_target(
    FTS_param_ridges,
    fts_param_ridge_plot(FTS_param_ridge_df)
  ),
  
  tar_target(
    FTS_ridge_pdf,
    save_ggplot_pdf(
      plot = FTS_param_ridges,
      path = file.path(manuscript_figures_dir, "FTS_ridge_plots.pdf"),
      width = 11,
      height = 8
    ),
    format = "file"
  ),
  
  #### Supplementary tables ####
  tar_target(
    convergence_pts_df,
    convergence_df_pts(pts_spec, pts_bundle_files)
  ),
  
  tar_target(
    convergence_voi_df,
    convergence_df_voi(voi_spec, voi_bundle_files)
  ),
  
  tar_target(
    convergence_npi_df,
    convergence_df_npi(npi_spec, npi_bundle_files)
  ),
  
  tar_target(
    convergence_fts_df,
    convergence_df_fts(
      files = FTS_fit_files,
      grid_csv = FTS_grid_csv,
      r0_keep = 3
    )
  ),
  
  tar_target(
    convergence_workbook,
    save_convergence_workbook(
      dfs = list(
        PTS = convergence_pts_df,
        VOI = convergence_voi_df,
        NPI = convergence_npi_df,
        FTS = convergence_fts_df
      ),
      path = file.path("outputs", "tables", "Supplementary_Data_1_convergence.xlsx")
    ),
    format = "file"
  ),
  
  tar_target(
    supp_fig3_source_csv,
    "outputs/source_data_SuppFig3_homog_truth.csv",
    format = "file"
  ),
  
  tar_target(
    supp_fig3_source_df,
    utils::read.csv(supp_fig3_source_csv, stringsAsFactors = FALSE)
  ),
  
  tar_target(
    source_data_workbook,
    save_source_data_workbook(
      sheets = list(
        "Fig1a_susceptibility"     = fig_0_dist_df,
        "Fig1b_mean_susceptibility"= fig_0_mean_sus_df,
        "Fig2ab_trajectories"      = combined_long,
        "Fig2c_hSIR_fit_to_SIR"    = incidence_summary_1,
        "Fig2d_SIR_fit_to_hSIR"    = incidence_summary_2,
        "Fig3_nu_posteriors"       = dynamic_speed_cv_ridge_df,
        "Fig4_remaining_attack"    = fig3_remaining_df,
        "Fig4_truth"               = remaining_attack_true,
        "Fig5_tradeoff"            = npi_tradeoff_grid_df,
        "Fig6_remaining_by_model"  = npi_remaining_grid,
        "Fig7_marginal_value"      = npi_mv_df,
        "Fig8_post_NPI_nu"         = npi_informing_cv_draws_df,
        "SuppFig1_I_dynamics"      = I_dynamics_df,
        "SuppFig2_FTS_ridges"      = FTS_param_ridge_df,
        "SuppFig3_homog_truth"     = supp_fig3_source_df,
        "SuppFig4_6_heatmaps"      = npi_results_wide
      ),
      path = file.path("outputs", "Source_Data.xlsx")
    ),
    format = "file"
  )
  
)