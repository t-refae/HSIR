library(targets)

tar_option_set(
  packages = c(
    "deSolve",
    "data.table",
    "dplyr",
    "tidyr",
    "ggplot2",
    "patchwork",
    "posterior",
    "pracma"
  ),
  format = "rds"
)

source("R/functions_pipeline.R")

list(
  # ---------------------------------------------------------------------------
  # Global constants
  # ---------------------------------------------------------------------------
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
  
  # ---------------------------------------------------------------------------
  # Simulated trajectories
  # ---------------------------------------------------------------------------
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
  
  tar_target(
    auc_differences,
    calculate_auc_differences(homog_out, HS_out)
  ),
  
  # ---------------------------------------------------------------------------
  # Figure 1 trajectory panels
  # ---------------------------------------------------------------------------
  tar_target(
    fig_1_c,
    plot_trajectories_faceted(combined_long)
  ),
  
  tar_target(
    fig_1_d,
    plot_absolute_differences(delta_df)
  ),
  
  # ---------------------------------------------------------------------------
  # Existing fitted model files
  # ---------------------------------------------------------------------------
  tar_target(fit_sir_rds, "fit_sir.RDS", format = "file"),
  tar_target(fit_sir2_rds, "fit_sir2.RDS", format = "file"),
  
  tar_target(fit_sir, readRDS(fit_sir_rds)),
  tar_target(fit_sir2, readRDS(fit_sir2_rds)),
  
  # ---------------------------------------------------------------------------
  # Posterior predictive incidence: HS model fit to homogeneous data
  # ---------------------------------------------------------------------------
  tar_target(
    posterior_dt_1,
    extract_posterior_draws(
      fit = fit_sir,
      variables = c("beta", "gamma", "cv"),
      n_draws = 500,
      seed = 2
    )
  ),
  
  tar_target(
    homog_cases,
    make_cases_from_cumulative(
      out = homog_out,
      times = shared_tps,
      population_size = population_size
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
      observed_cases = homog_cases,
      ribbon_fill = "#ff7f0e",
      line_color = "#ff7f0e",
      point_color = "#1f77b4",
      line_label = "Heterogeneous Model Fit",
      point_label = "Data (Homogeneous Model)"
    )
  ),
  
  # ---------------------------------------------------------------------------
  # Posterior predictive incidence: homogeneous model fit to HS data
  # ---------------------------------------------------------------------------
  tar_target(
    posterior_dt_2,
    extract_posterior_draws(
      fit = fit_sir2,
      variables = c("beta", "gamma"),
      n_draws = 500,
      seed = 2
    )
  ),
  
  tar_target(
    HS_cases,
    make_cases_from_cumulative(
      out = HS_out,
      times = shared_tps,
      population_size = population_size
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
      observed_cases = HS_cases,
      ribbon_fill = "#1f77b4",
      line_color = "#1f77b4",
      point_color = "#ff7f0e",
      line_label = "Homogeneous Model Fit",
      point_label = "Data (Heterogeneous Model)"
    )
  ),
  
  # ---------------------------------------------------------------------------
  # Final manuscript figure
  # ---------------------------------------------------------------------------
  tar_target(
    fig_1_panel,
    assemble_fig_1_panel(
      fig_1_c = fig_1_c,
      fig_1_d = fig_1_d,
      plot_fit_1 = plot_fit_1,
      plot_fit_2 = plot_fit_2
    )
  ),
  
  tar_target(
    manuscript_figures_dir,
    {
      dir.create("manuscript_figures", showWarnings = FALSE, recursive = TRUE)
      "manuscript_figures"
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