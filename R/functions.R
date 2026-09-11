#### ODE models ####

homog_SIR <- function(t, state, params) {
  with(as.list(c(state, params)), {
    dS <- -beta * S * I
    dI <- beta * S * I - gamma * I
    dR <- gamma * I
    dC <- beta * S * I
    
    list(c(dS, dI, dR, dC))
  })
}

HS_SIR <- function(t, state, params) {
  with(as.list(c(state, params)), {
    dS <- -beta * I * S^(1 + cv^2)
    dI <- beta * I * S^(1 + cv^2) - gamma * I
    dR <- gamma * I
    dC <- beta * I * S^(1 + cv^2)
    
    list(c(dS, dI, dR, dC))
  })
}

get_model_function <- function(model) {
  model <- match.arg(model, c("homogeneous", "heterogeneous"))
  
  if (model == "homogeneous") {
    homog_SIR
  } else {
    HS_SIR
  }
}

simulate_sir <- function(init_state, times, params, model) {
  model_fun <- get_model_function(model)
  
  out <- deSolve::ode(
    y = init_state,
    times = times,
    func = model_fun,
    parms = params,
    method = "lsoda"
  )
  
  data.table::as.data.table(out)
}

#### data prep ####

combine_trajectory_outputs <- function(homog_out, HS_out) {
  homog_out <- data.table::copy(homog_out)
  HS_out <- data.table::copy(HS_out)
  
  homog_out[, susceptibility := "Homogeneous"]
  HS_out[, susceptibility := "Heterogeneous"]
  
  dplyr::bind_rows(homog_out, HS_out) |>
    dplyr::mutate(
      susceptibility = factor(
        susceptibility,
        levels = c("Homogeneous", "Heterogeneous")
      )
    )
}

make_combined_long <- function(combined_data) {
  combined_data |>
    dplyr::select(time, susceptibility, I, C) |>
    tidyr::pivot_longer(
      cols = c(I, C),
      names_to = "Compartment",
      values_to = "Value"
    ) |>
    dplyr::mutate(
      Compartment = factor(
        Compartment,
        levels = c("C", "I"),
        labels = c("Cumulative Incidence", "Prevalence")
      )
    )
}

make_delta_df <- function(homog_out, HS_out) {
  data.frame(
    time = homog_out$time,
    delta_infected = abs(homog_out$I - HS_out$I),
    delta_cumulative_incidence = abs(homog_out$C - HS_out$C)
  ) |>
    tidyr::pivot_longer(
      cols = -time,
      names_to = "Metric",
      values_to = "Absolute_Difference"
    )
}

make_stan_cases_from_cumulative <- function(out, population_size) {
  round(diff(out$C) * population_size)
}

make_stan_data <- function(
    n_days,
    y0,
    t0,
    ts,
    population_size,
    cases
) {
  list(
    n_days = n_days,
    y0 = y0,
    t0 = t0,
    ts = ts,
    N = population_size,
    cases = cases[seq_len(n_days - 1)]
  )
}

#### posterior processing ####

extract_posterior_draws <- function(fit, variables, n_draws = 500, seed = 2) {
  posterior_draws <- posterior::as_draws_df(fit$draws(variables))
  
  posterior_dt <- data.table::as.data.table(posterior_draws)[, ..variables]
  
  set.seed(seed)
  posterior_dt <- posterior_dt[sample(.N, min(.N, n_draws))]
  
  posterior_dt[, draw := .I]
  
  posterior_dt
}

simulate_posterior_incidence <- function(
    posterior_dt,
    init_state,
    times,
    population_size,
    model
) {
  model <- match.arg(model, c("homogeneous", "heterogeneous"))
  model_fun <- get_model_function(model)
  
  incidence_list <- lapply(seq_len(nrow(posterior_dt)), function(i) {
    params <- posterior_dt[i]
    
    if (model == "homogeneous") {
      ode_params <- c(
        beta = params$beta,
        gamma = params$gamma
      )
    } else {
      ode_params <- c(
        beta = params$beta,
        gamma = params$gamma,
        cv = params$cv
      )
    }
    
    sim <- deSolve::ode(
      y = init_state,
      times = times,
      func = model_fun,
      parms = ode_params,
      method = "lsoda"
    )
    
    inc <- c(1, round(diff(sim[, "C"]) * population_size))
    
    data.table::data.table(
      Day = times,
      Incidence = inc,
      draw = i
    )
  })
  
  data.table::rbindlist(incidence_list)
}

summarize_posterior_incidence <- function(posterior_incidence) {
  posterior_incidence[
    ,
    .(
      Lower = stats::quantile(Incidence, 0.025),
      Upper = stats::quantile(Incidence, 0.975),
      Median = stats::median(Incidence)
    ),
    by = Day
  ]
}

#### dynamic speed / partial time-series cv inference ####

plot_dynamic_speed_cv_density_matrix <- function(cv_draws_df, true_cv_df) {
  ggplot2::ggplot(
    cv_draws_df,
    ggplot2::aes(x = value, fill = window)
  ) +
    ggplot2::geom_density(alpha = 0.6, color = NA) +
    ggh4x::facet_grid2(
      rows = ggplot2::vars(setting),
      cols = ggplot2::vars(window),
      scales = "free",
      independent = "all",
      strip = ggh4x::strip_nested(
        text_x = list(ggplot2::element_text(angle = 0)),
        text_y = list(ggplot2::element_text(angle = 0))
      )
    ) +
    ggplot2::geom_vline(
      data = true_cv_df,
      ggplot2::aes(xintercept = true_value),
      inherit.aes = FALSE,
      linetype = "dashed",
      color = "black",
      linewidth = 1
    ) +
    ggplot2::scale_fill_viridis_d(
      option = "plasma",
      end = 0.85,
      name = NULL
    ) +
    ggplot2::labs(
      x = "CV Value",
      y = "Density"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      strip.text = ggplot2::element_text(size = 14, face = "bold"),
      legend.position = "none",
      axis.title.x = ggplot2::element_text(size = 14, face = "bold"),
      axis.title.y = ggplot2::element_text(size = 14, face = "bold"),
      panel.spacing = grid::unit(1, "lines")
    ) +
    ggplot2::scale_x_continuous(
      breaks = scales::pretty_breaks(n = 3)
    ) +
    ggplot2::scale_y_continuous(
      breaks = scales::pretty_breaks(n = 3)
    )
}

make_dynamic_speed_cv_ridge_df <- function(
    cv_draws_df,
    setting_levels = c("Fast", "Reference", "Slow"),
    window_levels = c("-2 GI pre-peak", "-1 GI pre-peak", "At peak", "Complete")
) {
  ridge_levels <- as.vector(
    outer(setting_levels, window_levels, paste, sep = " — ")
  )
  
  cv_draws_df |>
    dplyr::mutate(
      setting = factor(as.character(setting), levels = setting_levels),
      window = factor(as.character(window), levels = window_levels),
      ridge_label = factor(
        paste(setting, window, sep = " — "),
        levels = ridge_levels
      )
    )
}

plot_dynamic_speed_cv_ridges_by_window <- function(cv_ridge_df) {
  ggplot2::ggplot(
    cv_ridge_df,
    ggplot2::aes(
      x = value,
      y = setting,
      fill = ggplot2::after_stat(x)
    )
  ) +
    ggridges::geom_density_ridges_gradient(
      scale = 3,
      rel_min_height = 0.01,
      color = "white"
    ) +
    ggplot2::geom_vline(
      xintercept = 1,
      linetype = "dashed",
      color = "black",
      linewidth = 0.8
    ) +
    ggplot2::facet_grid(
      window ~ .,
      scales = "free_y",
      space = "free_y"
    ) +
    ggplot2::scale_fill_viridis_c(
      option = "plasma",
      end = 0.85,
      name = NULL
    ) +
    ggplot2::labs(
      x = "CV value",
      y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      strip.text.y = ggplot2::element_text(face = "bold"),
      axis.text.y = ggplot2::element_text(size = 11),
      legend.position = "none"
    )
}



#### plotting ####

susceptibility_palette <- c(
  "Homogeneous"   = "#1f77b4",
  "Heterogeneous" = "#ff7f0e"
)

theme_fig1 <- function(base_size = 11) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      legend.position  = "bottom",
      legend.box       = "horizontal",
      legend.title     = ggplot2::element_text(face = "bold"),
      plot.title       = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin      = ggplot2::margin(t = 8, r = 8, b = 2, l = 6)
    )
}

plot_trajectory_panel <- function(combined_long, compartment, y_label,
                                  palette = susceptibility_palette) {
  ggplot2::ggplot(
    dplyr::filter(combined_long, Compartment == compartment),
    ggplot2::aes(x = time, y = Value, colour = susceptibility)
  ) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::scale_colour_manual(
      values = palette,
      limits = names(palette),
      name   = "Susceptibility"
    ) +
    ggplot2::guides(colour = ggplot2::guide_legend(order = 1)) +
    ggplot2::labs(x = "Time (days)", y = y_label) +
    theme_fig1()
}

plot_fit_to_data <- function(incidence_summary, observed_cases,
                             fit_label, data_label,
                             palette = susceptibility_palette) {
  glyph_levels <- c("data", "model fit")
  
  ggplot2::ggplot() +
    ggplot2::geom_ribbon(
      data = incidence_summary,
      ggplot2::aes(x = Day, ymin = Lower, ymax = Upper, fill = fit_label),
      alpha = 0.25
    ) +
    ggplot2::geom_line(
      data = incidence_summary,
      ggplot2::aes(x = Day, y = Median, colour = fit_label, linetype = "model fit"),
      linewidth = 1
    ) +
    ggplot2::geom_point(
      data = observed_cases,
      ggplot2::aes(x = Day, y = Data, colour = data_label, shape = "data"),
      size = 2
    ) +
    ggplot2::scale_colour_manual(
      values = palette,
      limits = names(palette),
      guide  = "none"
    ) +
    ggplot2::scale_fill_manual(
      values = palette,
      limits = names(palette),
      guide  = "none"
    ) +
    ggplot2::scale_shape_manual(
      name   = "Data type",
      values = c("data" = 16, "model fit" = NA),
      limits = glyph_levels
    ) +
    ggplot2::scale_linetype_manual(
      name   = "Data type",
      values = c("data" = "blank", "model fit" = "solid"),
      limits = glyph_levels
    ) +
    ggplot2::guides(
      shape    = ggplot2::guide_legend(
        order = 2,
        override.aes = list(colour = "black")
      ),
      linetype = ggplot2::guide_legend(
        order = 2,
        override.aes = list(colour = "black")
      )
    ) +
    ggplot2::labs(x = "Time (days)", y = "Incidence") +
    theme_fig1()
}

assemble_fig_1_panel <- function(fig_1_a, fig_1_b, plot_fit_1, plot_fit_2) {
  (
    (fig_1_a + fig_1_b) /
      (plot_fit_1 + plot_fit_2)
  ) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(
      tag_levels = "a",
      tag_prefix = "(",
      tag_suffix = ")"
    ) &
    ggplot2::theme(
      legend.position = "bottom",
      legend.box      = "horizontal",
      plot.tag        = ggplot2::element_text(face = "bold")
    )
}

save_ggplot_pdf <- function(plot, path, width, height, dpi = 300,
                            device = grDevices::cairo_pdf) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  
  ggplot2::ggsave(
    filename = path,
    plot = plot,
    device = grDevices::cairo_pdf,
    width = width,
    height = height,
    dpi = dpi
  )
  
  path
}

#### Fig. 3: S/HIT differences and VOI-probability summaries ####

make_voi_true_values_long <- function(
    true_params,
    voi_thetas,
    vars = c("beta", "gamma", "cv", "R0")
) {
  theta_levels <- rev(paste0("Theta[", voi_thetas, "]"))
  
  param_labels <- c(
    beta = "beta",
    gamma = "gamma",
    cv = "nu",
    R0 = "R[0]"
  )
  
  true_params |>
    dplyr::mutate(
      theta_label = factor(
        paste0("Theta[", Theta_num, "]"),
        levels = theta_levels
      )
    ) |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(vars),
      names_to = "param",
      values_to = "true_value"
    ) |>
    dplyr::mutate(
      param_label = factor(
        param,
        levels = names(param_labels),
        labels = param_labels
      )
    )
}

trim_voi_parameter_draws <- function(draws_df, r0_upper = 6.5) {
  draws_df |>
    dplyr::filter(param != "R0" | value < r0_upper)
}

plot_voi_parameter_ridges <- function(draws_df, true_values_df) {
  ggplot2::ggplot(
    draws_df,
    ggplot2::aes(
      x = value,
      y = theta_label,
      fill = ggplot2::after_stat(x)
    )
  ) +
    ggridges::geom_density_ridges_gradient(
      scale = 3,
      rel_min_height = 0.01,
      color = "white"
    ) +
    ggplot2::facet_wrap(
      ~ param_label,
      scales = "free_x",
      labeller = ggplot2::label_parsed
    ) +
    ggplot2::geom_vline(
      data = true_values_df,
      ggplot2::aes(
        xintercept = true_value,
        group = interaction(param_label, theta_label)
      ),
      inherit.aes = FALSE,
      linetype = "dashed",
      color = "black",
      linewidth = 1
    ) +
    ggplot2::scale_fill_viridis_c(
      name = "Value",
      option = "plasma",
      end = 0.85
    ) +
    ggplot2::scale_y_discrete(
      labels = function(x) parse(text = x)
    ) +
    ggplot2::labs(
      x = "Value",
      y = "Parameter Set"
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold", size = 14),
      axis.text.y = ggplot2::element_text(size = 12),
      legend.position = "none"
    )
}

calculate_true_hs_x_star <- function(VOI_true_params) {
  R0_vals <- VOI_true_params$beta / VOI_true_params$gamma
  cv_vals <- VOI_true_params$cv
  
  (1 / R0_vals)^(1 / (1 + cv_vals^2))
}

extract_SIR_summary <- function(y_draws, theta_label = "Theta_x") {
  var_names <- dimnames(y_draws)[[3]]
  
  if (is.null(var_names)) {
    stop("y_draws must have variable names in dimnames(y_draws)[[3]].")
  }
  
  get_state_draws <- function(state_index) {
    pattern <- paste0("^y\\[([0-9]+),", state_index, "\\]$")
    vars <- grep(pattern, var_names, value = TRUE)
    
    if (length(vars) == 0) {
      stop("No y_draws variables found for state index ", state_index)
    }
    
    time_indices <- as.integer(
      sub("y\\[([0-9]+),[0-9]+\\]", "\\1", vars)
    )
    
    vars_ordered <- vars[order(time_indices)]
    state_draws <- y_draws[, , vars_ordered, drop = FALSE]
    
    samples <- apply(state_draws, 3, function(x) as.vector(x))
    samples <- t(samples)
    
    summary_df <- apply(samples, 1, function(x) {
      c(
        median = stats::median(x, na.rm = TRUE),
        lower = stats::quantile(x, 0.025, na.rm = TRUE),
        upper = stats::quantile(x, 0.975, na.rm = TRUE)
      )
    }) |>
      t() |>
      as.data.frame()
    
    summary_df$time <- seq_len(nrow(summary_df))
    summary_df
  }
  
  S <- get_state_draws(1) |>
    dplyr::mutate(state = "S")
  
  I <- get_state_draws(2) |>
    dplyr::mutate(state = "I")
  
  R <- get_state_draws(3) |>
    dplyr::mutate(state = "R")
  
  dplyr::bind_rows(S, I, R) |>
    dplyr::mutate(theta = theta_label)
}

make_S_minus_xstar_df <- function(S_summary, x_stars) {
  S_summary |>
    dplyr::left_join(
      x_stars |>
        dplyr::rename(
          Theta_num = Theta,
          xstar_med = median,
          xstar_low = lower,
          xstar_up = upper
        ),
      by = "Theta_num"
    ) |>
    dplyr::mutate(
      diff_median = median - xstar_med,
      diff_lower = lower - xstar_up,
      diff_upper = upper - xstar_low
    ) |>
    dplyr::select(
      theta,
      time,
      model_type,
      diff_median,
      diff_lower,
      diff_upper
    )
}

make_fig3_df <- function(S_HS_diff, S_Homog_diff) {
  dplyr::bind_rows(S_HS_diff, S_Homog_diff) |>
    dplyr::mutate(
      model_type = factor(
        model_type,
        levels = c("Homogeneous", "Heterogeneous")
      )
    )
}

plot_fig3_S_minus_xstar <- function(fig3_df, theta_labels_named) {
  fig3_df <- fig3_df |>
    dplyr::mutate(
      theta = factor(theta, levels = names(theta_labels_named))
    )
  
  ggplot2::ggplot(
    fig3_df,
    ggplot2::aes(
      x = time,
      y = diff_median,
      color = model_type,
      fill = model_type
    )
  ) +
    ggplot2::geom_line(linewidth = 1.3) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = diff_lower, ymax = diff_upper),
      alpha = 0.3,
      color = NA
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      color = "black",
      linetype = "dashed"
    ) +
    ggplot2::facet_wrap(
      ~ theta,
      ncol = 2,
      labeller = ggplot2::as_labeller(theta_labels_named),
      scales = "free_y"
    ) +
    ggplot2::scale_color_manual(
      values = c(
        "Homogeneous" = "#1f77b4",
        "Heterogeneous" = "#ff7f0e"
      )
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        "Homogeneous" = "#1f77b4",
        "Heterogeneous" = "#ff7f0e"
      )
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::percent_format(accuracy = 1)
    ) +
    ggplot2::labs(
      x = "Time (days)",
      y = latex2exp::TeX("$S(t) - x^*$"),
      color = "Susceptibility",
      fill = "Susceptibility"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position = "top",
      legend.title = ggplot2::element_text(face = "bold"),
      legend.box.background = ggplot2::element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.6
      ),
      legend.background = ggplot2::element_blank(),
      legend.key.height = grid::unit(0.5, "cm"),
      axis.title = ggplot2::element_text(size = 12),
      axis.text = ggplot2::element_text(size = 10),
      strip.text = ggplot2::element_text(size = 12, face = "bold")
    )
}

save_ggplot_pdf_plain <- function(plot, path, width, height) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  
  ggplot2::ggsave(
    filename = path,
    plot = plot,
    device = "pdf",
    width = width,
    height = height,
    units = "in",
    useDingbats = FALSE
  )
  
  path
}

#### Fig. 4: NPI policy simulations ####

HS_NPI <- function(t, state, params) {
  with(as.list(c(state, params)), {
    if (t >= NPI_start && t <= (NPI_start + NPI_dur)) {
      beta0 <- (1 - eff) * beta
    } else {
      beta0 <- beta
    }
    
    dS <- -beta0 * S^(1 + cv^2) * I
    dI <- beta0 * S^(1 + cv^2) * I - gamma * I
    dR <- gamma * I
    dC <- beta0 * S^(1 + cv^2) * I
    
    list(c(dS, dI, dR, dC))
  })
}

homog_NPI <- function(t, state, params) {
  with(as.list(c(state, params)), {
    if (t >= NPI_start && t <= (NPI_start + NPI_dur)) {
      beta0 <- (1 - eff) * beta
    } else {
      beta0 <- beta
    }
    
    dS <- -beta0 * S * I
    dI <- beta0 * S * I - gamma * I
    dR <- gamma * I
    dC <- beta0 * S * I
    
    list(c(dS, dI, dR, dC))
  })
}

make_npi_policy_grid <- function(eff_seq, start_seq, dur_seq) {
  expand.grid(
    eff = eff_seq,
    NPI_start = start_seq,
    NPI_dur = dur_seq,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  ) |>
    tibble::as_tibble()
}

simulate_policy_metrics <- function(
    model = c("HS", "homog"),
    eff,
    NPI_start,
    NPI_dur,
    beta,
    gamma,
    cv = NULL,
    t_max = 230,
    population_size = 1e6,
    i0_prop = 1 / population_size
) {
  model <- match.arg(model)
  
  y0 <- c(
    S = 1 - i0_prop,
    I = i0_prop,
    R = 0,
    C = 0
  )
  
  times <- seq(1, t_max)
  
  parms <- c(
    beta = beta,
    gamma = gamma,
    eff = eff,
    NPI_start = NPI_start,
    NPI_dur = NPI_dur
  )
  
  if (model == "HS") {
    if (is.null(cv)) {
      stop("cv must be supplied when model = 'HS'.")
    }
    
    parms <- c(parms, cv = cv)
    ode_fun <- HS_NPI
  } else {
    ode_fun <- homog_NPI
  }
  
  traj <- deSolve::ode(
    y = y0,
    times = times,
    func = ode_fun,
    parms = parms,
    method = "lsoda"
  ) |>
    as.data.frame()
  
  I_count <- traj$I * population_size
  peak_prev <- max(I_count, na.rm = TRUE)
  time_to_peak <- traj$time[which.max(I_count)]
  
  inc <- pmax(0, diff(traj$C) * population_size)
  peak_incidence <- max(inc, na.rm = TRUE)
  
  attack_rate <- traj$C[traj$time == t_max][1]
  
  end_day <- NPI_start + NPI_dur
  prop_S_end <- if (end_day <= t_max) {
    traj$S[traj$time == end_day][1]
  } else {
    NA_real_
  }
  
  tibble::tibble(
    eff = eff,
    NPI_start = NPI_start,
    NPI_dur = NPI_dur,
    attack_rate = attack_rate,
    time_to_peak = time_to_peak,
    prop_S_end = prop_S_end,
    peak_prev = peak_prev,
    peak_incidence = peak_incidence
  )
}

simulate_npi_policy_grid <- function(
    policy_grid,
    model = c("HS", "homog"),
    medians,
    t_max = 230,
    population_size = 1e6,
    i0_prop = 1 / population_size
) {
  model <- match.arg(model)
  
  beta <- medians$beta[[1]]
  gamma <- medians$gamma[[1]]
  cv <- if ("cv" %in% names(medians)) medians$cv[[1]] else NULL
  
  purrr::pmap_dfr(
    policy_grid,
    function(eff, NPI_start, NPI_dur) {
      simulate_policy_metrics(
        model = model,
        eff = eff,
        NPI_start = NPI_start,
        NPI_dur = NPI_dur,
        beta = beta,
        gamma = gamma,
        cv = cv,
        t_max = t_max,
        population_size = population_size,
        i0_prop = i0_prop
      ) |>
        dplyr::mutate(model = model)
    }
  )
}

make_npi_results_wide <- function(results_all) {
  results_all |>
    tidyr::pivot_wider(
      id_cols = c(eff, NPI_start, NPI_dur),
      names_from = model,
      values_from = c(
        attack_rate,
        prop_S_end,
        peak_incidence,
        time_to_peak,
        peak_prev
      ),
      names_glue = "{.value}_{model}"
    ) |>
    dplyr::mutate(
      attack_rate_diff = attack_rate_homog - attack_rate_HS,
      prop_S_end_diff = prop_S_end_homog - prop_S_end_HS,
      peak_incidence_diff = peak_incidence_homog - peak_incidence_HS,
      time_to_peak_diff = time_to_peak_homog - time_to_peak_HS,
      peak_prev_diff = peak_prev_homog - peak_prev_HS
    )
}

calculate_npi_midpoint <- function(results_wide, z) {
  midpoint <- results_wide |>
    dplyr::filter(eff == 0, NPI_start == 0, NPI_dur == 30) |>
    dplyr::summarise(mid = mean(.data[[z]], na.rm = TRUE)) |>
    dplyr::pull(mid)
  if (length(midpoint) == 0 || is.na(midpoint)) midpoint <- 0
  midpoint
}

theme_npi_heat <- function() {
  ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold"),
      legend.title = ggplot2::element_text(face = "bold")
    )
}

plot_npi_heatmap <- function(
    df,
    z,
    title,
    fill_lab,
    midpoint = 0,
    label_fun = ggplot2::waiver()
) {
  div_cols <- scico::scico(3, palette = "vik")
  
  ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = NPI_start,
      y = eff,
      fill = .data[[z]]
    )
  ) +
    ggplot2::geom_tile() +
    ggplot2::facet_wrap(
      ~ NPI_dur,
      labeller = ggplot2::label_both
    ) +
    ggplot2::scale_fill_gradient2(
      low = div_cols[1],
      mid = div_cols[2],
      high = div_cols[3],
      midpoint = midpoint,
      labels = label_fun
    ) +
    ggplot2::labs(
      title = title,
      x = "NPI Start Day",
      y = "NPI Effectiveness",
      fill = fill_lab
    ) +
    theme_npi_heat()
}

make_npi_tradeoff_grid_df <- function(
    results_all,
    start_keep = c(5, 10, 15),
    dur_keep = c(30, 60, 90)
) {
  results_all |>
    dplyr::filter(
      NPI_start %in% start_keep,
      NPI_dur %in% dur_keep
    ) |>
    dplyr::mutate(
      model = factor(
        model,
        levels = c("homog", "HS"),
        labels = c("Homogeneous", "Heterogeneous")
      ),
      NPI_dur = factor(
        NPI_dur,
        levels = dur_keep,
        labels = c(paste0("Duration (days): ", dur_keep[1]), as.character(dur_keep[-1]))
      ),
      NPI_start = factor(
        NPI_start,
        levels = start_keep,
        labels = c(paste0("Start day: ", start_keep[1]), as.character(start_keep[-1]))
      )
    )
}

plot_npi_tradeoff_grid <- function(df_tradeoff_grid) {
  ggplot2::ggplot(
    df_tradeoff_grid,
    ggplot2::aes(
      x = prop_S_end,
      y = attack_rate,
      color = eff,
      shape = model
    )
  ) +
    ggplot2::geom_point(size = 2.7, alpha = 0.95) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(NPI_start),
      cols = ggplot2::vars(NPI_dur)
    ) +
    ggplot2::scale_color_viridis_c(
      option = "plasma",
      end = 0.85,
      limits = c(0, 1),
      labels = scales::percent_format(accuracy = 1)
    ) +
    ggplot2::scale_shape_manual(
      values = c(
        "Homogeneous" = 16,
        "Heterogeneous" = 17
      )
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.25),
      labels = scales::percent_format(accuracy = 1)
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0,1),
      labels = scales::percent_format(accuracy = 1)
    ) +
    ggplot2::labs(
      x = "Proportion susceptible at NPI end",
      y = "Final attack rate",
      color = "NPI effectiveness",
      shape = "Susceptibility type"
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      strip.placement = "outside",
      strip.text.x = ggplot2::element_text(face = "bold"),
      strip.text.y.right = ggplot2::element_text(face = "bold", angle = 270),
      strip.text.y.left = ggplot2::element_blank(),
      strip.background.y.left = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(
        color = "grey40",
        fill = NA,
        linewidth = 0.6
      ),
      panel.spacing = grid::unit(0.6, "lines"),
      panel.grid.major.y = ggplot2::element_line(
        color = "grey85",
        linewidth = 0.4
      ),
      panel.grid.major.x = ggplot2::element_line(
        color = "grey90",
        linewidth = 0.3
      ),
      legend.position = "top",
      legend.box = "horizontal",
      legend.direction = "horizontal",
      legend.title = ggplot2::element_text(face = "bold"),
      legend.margin = ggplot2::margin(b = 6),
      legend.spacing.x = grid::unit(14, "pt")
    ) +
    ggplot2::guides(
      color = ggplot2::guide_colorbar(
        order = 1,
        direction = "horizontal",
        title.position = "top",
        label.position = "bottom",
        barwidth = grid::unit(3.5, "in"),
        barheight = grid::unit(0.18, "in")
      ),
      shape = ggplot2::guide_legend(
        order = 2,
        direction = "horizontal",
        title.position = "top",
        label.position = "right",
        nrow = 1,
        byrow = TRUE,
        override.aes = list(size = 3, alpha = 1)
      )
    )
}

#### Post-NPI information gain: CV posterior distributions ####

list_npi_informing_cv_files <- function(parent_dir) {
  if (!dir.exists(parent_dir)) {
    stop("Directory does not exist: ", parent_dir)
  }
  
  cv_files <- list.files(
    parent_dir,
    pattern = "_cv_draws\\.RData$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  cv_files <- cv_files[grepl("NPI_HS", cv_files)]
  
  if (length(cv_files) == 0) {
    stop(
      "No NPI_HS *_cv_draws.RData files found under: ",
      parent_dir
    )
  }
  
  sort(cv_files)
}

parse_npi_informing_theta_info <- function(path) {
  fname <- basename(path)
  
  theta_chr <- sub(
    pattern = ".*Theta_([0-9]+).*",
    replacement = "\\1",
    x = fname
  )
  
  theta <- suppressWarnings(as.integer(theta_chr))
  
  if (is.na(theta)) {
    stop("Could not parse Theta number from file name: ", fname)
  }
  
  tibble::tibble(
    theta = theta,
    eff = dplyr::if_else(theta <= 4, 0.4, 0.8),
    GI = (theta - 1) %% 4
  )
}

load_single_npi_cv_draw_file <- function(path) {
  temp_env <- new.env(parent = emptyenv())
  loaded_names <- load(path, envir = temp_env)
  
  if (!"cv_draws" %in% loaded_names) {
    stop("File does not contain an object named cv_draws: ", path)
  }
  
  cv_draws <- get("cv_draws", envir = temp_env)
  
  info <- parse_npi_informing_theta_info(path)
  
  cv_df <- posterior::as_draws_df(cv_draws)
  
  if (!"cv" %in% names(cv_df)) {
    stop("Object cv_draws does not contain a column named cv in: ", path)
  }
  
  cv_df |>
    dplyr::select(cv) |>
    dplyr::mutate(
      theta = info$theta,
      eff = info$eff,
      GI = info$GI,
      source_file = path
    )
}

make_npi_informing_cv_draws_df <- function(cv_files) {
  purrr::map_dfr(
    cv_files,
    load_single_npi_cv_draw_file
  )
}

plot_npi_informing_cv_posteriors <- function(cv_draws_df) {
  cv_draws_df |>
    dplyr::mutate(
      GI = factor(
        GI,
        levels = 0:3,
        labels = c("+0 GI", "+1 GI", "+2 GI", "+3 GI")
      ),
      eff = factor(eff)
    ) |>
    ggplot2::ggplot(
      ggplot2::aes(
        x = cv,
        y = GI,
        fill = eff
      )
    ) +
    ggridges::geom_density_ridges(
      alpha = 0.6,
      scale = 1.1,
      rel_min_height = 0.01,
      panel_scaling = FALSE
    ) +
    ggplot2::geom_vline(
      xintercept = 1,
      linetype = "dashed",
      linewidth = 0.8,
      color = "black"
    ) +
    ggplot2::facet_wrap(
      ~ eff,
      nrow = 1,
      labeller = ggplot2::as_labeller(
        c(
          `0.4` = "NPI eff = 0.4",
          `0.8` = "NPI eff = 0.8"
        )
      )
    ) +
    ggplot2::coord_cartesian(xlim = c(0, 2)) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::labs(
      x = "Posterior CV",
      y = "Observed window after NPI"
    ) +
    ggplot2::theme(
      legend.position = "none",
      strip.text = ggplot2::element_text(face = "bold")
    )
}

#### Supplementary: full time series fits ####
param_ridge_data <- function(bundles, files, param) {
  ids <- regmatches(basename(files), regexpr("id[0-9]+", basename(files)))
  dplyr::bind_rows(Map(function(b, id) {
    stopifnot(param %in% colnames(b$draws))
    data.frame(id = id, value = as.numeric(b$draws[[param]]))
  }, bundles, ids))
}

param_ridge_plot <- function(bundles, files, param) {
  d  <- param_ridge_data(bundles, files, param)
  lv <- unique(d$id)
  lv <- lv[order(as.integer(sub("id", "", lv)), decreasing = TRUE)]
  d$id <- factor(d$id, levels = lv)
  ggplot2::ggplot(d, ggplot2::aes(x = value, y = id, fill = id)) +
    ggridges::geom_density_ridges(
      scale = 1.1, alpha = 0.8, colour = "grey30", rel_min_height = 0.01
    ) +
    ggplot2::scale_fill_viridis_d(guide = "none") +
    ggplot2::labs(x = param, y = NULL) +
    ggplot2::theme_minimal(base_size = 11)
}

#### Alternative Fig. 3: remaining attack rate (unmitigated remaining burden) ####

extract_state_draws_matrix <- function(y_draws, state_index) {
  var_names <- dimnames(y_draws)[[3]]
  if (is.null(var_names)) {
    stop("y_draws must have variable names in dimnames(y_draws)[[3]].")
  }
  pattern <- paste0("^y\\[([0-9]+),", state_index, "\\]$")
  vars <- grep(pattern, var_names, value = TRUE)
  if (length(vars) == 0) {
    stop("No y_draws variables found for state index ", state_index)
  }
  time_indices <- as.integer(sub("y\\[([0-9]+),[0-9]+\\]", "\\1", vars))
  vars <- vars[order(time_indices)]
  vapply(
    vars,
    function(v) as.vector(y_draws[, , v]),
    numeric(prod(dim(y_draws)[1:2])),
    USE.NAMES = FALSE
  )
}

final_size_from_state <- function(
    S_end, I_end, R_end, beta, gamma, cv = NULL,
    model = c("heterogeneous", "homogeneous"),
    t_extend = 5000, i_tol = 1e-9
) {
  model <- match.arg(model)
  if (is.na(I_end) || I_end <= i_tol) return(unname(S_end))
  
  parms <- c(beta = beta, gamma = gamma)
  if (model == "heterogeneous") {
    if (is.null(cv) || is.na(cv)) stop("cv required for heterogeneous final size.")
    parms <- c(parms, cv = cv)
    ode_fun <- HS_SIR
  } else {
    ode_fun <- homog_SIR
  }
  
  y0 <- c(S = unname(S_end), I = unname(I_end), R = unname(R_end), C = 0)
  out <- deSolve::ode(
    y = y0, times = c(0, t_extend),
    func = ode_fun, parms = parms, method = "lsoda"
  )
  unname(out[nrow(out), "S"])
}

summarize_remaining_attack_by_model <- function(
    states,
    params,
    spec,
    arm,
    model_type,
    probs = c(0.025, 0.975),
    t_extend = 5000,
    thin = NULL,
    seed = 2
) {
  model <- if (model_type == "Heterogeneous") "heterogeneous" else "homogeneous"
  s <- spec[spec$arm == arm, ]
  
  purrr::map_dfr(seq_len(nrow(s)), function(i) {
    theta_val <- s$theta[i]
    y_draws <- states[[s$key[i]]]
    
    S_mat <- extract_state_draws_matrix(y_draws, 1)
    I_mat <- extract_state_draws_matrix(y_draws, 2)
    R_mat <- extract_state_draws_matrix(y_draws, 3)
    n_draws <- nrow(S_mat); n_time <- ncol(S_mat)
    
    p <- params[params$arm == arm & params$theta == theta_val, ]
    beta_v  <- p$beta
    gamma_v <- p$gamma
    cv_v <- if (model == "heterogeneous") p$cv else rep(NA_real_, n_draws)
    if (length(beta_v) != n_draws || length(gamma_v) != n_draws) {
      stop("Parameter draws not aligned with trajectory draws for Theta_", theta_val)
    }
    
    idx <- seq_len(n_draws)
    if (!is.null(thin) && thin < n_draws) {
      set.seed(seed); idx <- sort(sample.int(n_draws, thin))
    }
    
    S_inf <- vapply(idx, function(k) {
      final_size_from_state(
        S_end = S_mat[k, n_time], I_end = I_mat[k, n_time], R_end = R_mat[k, n_time],
        beta = beta_v[k], gamma = gamma_v[k],
        cv = if (model == "heterogeneous") cv_v[k] else NULL,
        model = model, t_extend = t_extend
      )
    }, numeric(1))
    
    S_sub <- S_mat[idx, , drop = FALSE]
    rem_mat <- S_sub - S_inf
    
    summ <- t(apply(rem_mat, 2, function(col) {
      c(median = stats::median(col, na.rm = TRUE),
        lower  = stats::quantile(col, probs[1], na.rm = TRUE, names = FALSE),
        upper  = stats::quantile(col, probs[2], na.rm = TRUE, names = FALSE))
    }))
    
    tibble::tibble(
      time = seq_len(n_time),
      rem_median = summ[, "median"],
      rem_lower  = summ[, "lower"],
      rem_upper  = summ[, "upper"],
      theta = paste0("Theta_", theta_val),
      model_type = model_type,
      Theta_num = theta_val
    )
  })
}

compute_true_remaining_attack <- function(
    true_params,
    remaining_df,
    i0 = 1e-4,
    t_extend = 5000
) {
  max_times <- remaining_df |>
    dplyr::group_by(Theta_num) |>
    dplyr::summarise(t_end = max(time), .groups = "drop")
  
  purrr::map_dfr(seq_len(nrow(max_times)), function(i) {
    theta_val <- max_times$Theta_num[i]
    t_end <- max_times$t_end[i]
    row <- true_params[true_params$Theta_num == theta_val, ]
    
    y0 <- c(S = 1 - i0, I = i0, R = 0, C = 0)
    traj <- as.data.frame(deSolve::ode(
      y = y0,
      times = c(seq_len(t_end + 1L), t_extend),
      func = HS_SIR,
      parms = c(beta = row$beta, gamma = row$gamma, cv = row$cv),
      method = "lsoda"
    ))
    S_inf <- traj$S[nrow(traj)]
    
    tibble::tibble(
      time = seq_len(t_end),
      rem_true = traj$S[1L + seq_len(t_end)] - S_inf,
      theta = paste0("Theta_", theta_val)
    )
  })
}

make_fig3_remaining_df <- function(remaining_HS, remaining_Homog) {
  dplyr::bind_rows(remaining_HS, remaining_Homog) |>
    dplyr::mutate(
      model_type = factor(
        model_type,
        levels = c("Homogeneous", "Heterogeneous")
      )
    )
}

plot_fig3_remaining_attack <- function(
    fig3_df,
    theta_labels_named,
    true_df = NULL
) {
  fig3_df <- fig3_df |>
    dplyr::mutate(theta = factor(theta, levels = names(theta_labels_named)))
  
  p <- ggplot2::ggplot(
    fig3_df,
    ggplot2::aes(x = time, y = rem_median, color = model_type, fill = model_type)
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = rem_lower, ymax = rem_upper),
      alpha = 0.3, color = NA
    ) +
    ggplot2::geom_line(linewidth = 1.3) +
    ggplot2::facet_wrap(
      ~ theta, ncol = 2,
      labeller = ggplot2::as_labeller(theta_labels_named)
    ) +
    ggplot2::scale_color_manual(
      values = c("Homogeneous" = "#1f77b4", "Heterogeneous" = "#ff7f0e")
    ) +
    ggplot2::scale_fill_manual(
      values = c("Homogeneous" = "#1f77b4", "Heterogeneous" = "#ff7f0e")
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::percent_format(accuracy = 1)
    ) +
    ggplot2::labs(
      x = "Time (days)",
      y = "Remaining attack rate (% of population)",
      color = "Susceptibility",
      fill = "Susceptibility"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position = "top",
      legend.title = ggplot2::element_text(face = "bold"),
      legend.box.background = ggplot2::element_rect(
        color = "black", fill = NA, linewidth = 0.6
      ),
      legend.background = ggplot2::element_blank(),
      legend.key.height = grid::unit(0.5, "cm"),
      axis.title = ggplot2::element_text(size = 12),
      axis.text = ggplot2::element_text(size = 10),
      strip.text = ggplot2::element_text(size = 12, face = "bold")
    )
  
  if (!is.null(true_df)) {
    true_df <- true_df |>
      dplyr::mutate(theta = factor(theta, levels = names(theta_labels_named)))
    p <- p +
      ggplot2::geom_line(
        data = true_df,
        ggplot2::aes(x = time, y = rem_true),
        inherit.aes = FALSE,
        color = "black",
        linetype = "dashed",
        linewidth = 0.9
      )
  }
  
  p
}

#### Alternative Fig. 4: continuous delta in remaining attack rate over time ####

simulate_policy_remaining <- function(
    model = c("HS", "homog"),
    eff, NPI_start, NPI_dur,
    beta, gamma, cv = NULL,
    t_max = 1080,
    population_size = 1e6,
    i0_prop = 1 / population_size,
    denom_floor = 1e-8
) {
  model <- match.arg(model)
  
  y0 <- c(S = 1 - i0_prop, I = i0_prop, R = 0, C = 0)
  times <- seq(1, t_max)
  parms <- c(
    beta = beta, gamma = gamma, eff = eff,
    NPI_start = NPI_start, NPI_dur = NPI_dur
  )
  
  if (model == "HS") {
    if (is.null(cv)) stop("cv must be supplied when model = 'HS'.")
    parms <- c(parms, cv = cv)
    ode_fun <- HS_NPI
  } else {
    ode_fun <- homog_NPI
  }
  
  traj <- deSolve::ode(
    y = y0, times = times, func = ode_fun,
    parms = parms, method = "lsoda"
  ) |>
    as.data.frame()
  
  S_inf <- traj$S[nrow(traj)]        # final size under THIS policy
  denom <- traj$S[1] - S_inf         # eventual total epidemic size for this policy
  if (denom <= denom_floor) denom <- NA_real_
  
  tibble::tibble(
    time = traj$time,
    remaining = (traj$S - S_inf) / denom,  # relative: 1 at t0, decays to 0
    eff = eff, NPI_start = NPI_start, NPI_dur = NPI_dur, model = model
  )
}

simulate_npi_remaining_grid <- function(
    policy_grid, medians_HS, medians_homog,
    t_max = 1080, population_size = 1e6, i0_prop = 1 / population_size
) {
  run_one <- function(model, medians) {
    beta <- medians$beta[[1]]
    gamma <- medians$gamma[[1]]
    cv <- if ("cv" %in% names(medians)) medians$cv[[1]] else NULL
    purrr::pmap_dfr(policy_grid, function(eff, NPI_start, NPI_dur) {
      simulate_policy_remaining(
        model = model, eff = eff, NPI_start = NPI_start, NPI_dur = NPI_dur,
        beta = beta, gamma = gamma, cv = cv,
        t_max = t_max, population_size = population_size, i0_prop = i0_prop
      )
    })
  }
  dplyr::bind_rows(
    run_one("HS", medians_HS),
    run_one("homog", medians_homog)
  )
}

plot_npi_remaining_by_model_grid <- function(
    remaining_grid,
    start_keep = c(5, 10, 15, 20),
    dur_keep = c(30, 60, 90),
    eff_keep = c(0.25, 0.5, 0.75),
    x_max_display = 200
) {
  df <- remaining_grid |>
    dplyr::filter(
      NPI_start %in% start_keep,
      NPI_dur %in% dur_keep,
      eff %in% eff_keep
    ) |>
    dplyr::mutate(
      model = factor(
        model, levels = c("homog", "HS"),
        labels = c("Homogeneous", "Heterogeneous")
      ),
      eff_f = factor(
        eff, levels = sort(eff_keep),
        labels = paste0(round(sort(eff_keep) * 100), "%")
      ),
      NPI_start_f = factor(NPI_start, levels = start_keep,
                           labels = paste0("Start: day ", start_keep)),
      NPI_dur_f = factor(NPI_dur, levels = dur_keep,
                         labels = paste0("Duration: ", dur_keep, "d"))
    )
  
  windows <- df |>
    dplyr::distinct(NPI_start, NPI_dur, NPI_start_f, NPI_dur_f) |>
    dplyr::mutate(xmin = NPI_start, xmax = NPI_start + NPI_dur)
  
  ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = time, y = remaining,
      color = eff_f, linetype = model,
      group = interaction(eff_f, model)
    )
  ) +
    ggplot2::geom_rect(
      data = windows,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE, fill = "grey88"
    ) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(NPI_start_f),
      cols = ggplot2::vars(NPI_dur_f)
    ) +
    ggplot2::coord_cartesian(xlim = c(0, x_max_display)) +
    ggplot2::scale_color_viridis_d(option = "plasma", end = 0.85) +
    ggplot2::scale_linetype_manual(
      values = c("Homogeneous" = "solid", "Heterogeneous" = "dashed")
    ) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::labs(
      x = "Time (days)",
      y = "Remaining attack rate (relative to eventual total)",
      color = "NPI effectiveness",
      linetype = "Susceptibility type"
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      strip.text.x = ggplot2::element_text(face = "bold"),
      strip.text.y = ggplot2::element_text(face = "bold", angle = 270),
      panel.border = ggplot2::element_rect(color = "grey40", fill = NA, linewidth = 0.6),
      panel.spacing = grid::unit(0.6, "lines"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "top",
      legend.title = ggplot2::element_text(face = "bold"),
      legend.key.width = grid::unit(2.2, "cm")
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(
        order = 1, title.position = "top",
        override.aes = list(linewidth = 1.1, linetype = "solid")
      ),
      linetype = ggplot2::guide_legend(
        order = 2, title.position = "top",
        keywidth = grid::unit(2.2, "cm"),
        override.aes = list(
          linewidth = 1.1,
          linetype = c("solid", "dashed"),
          colour = "grey20"
        )
      )
    )
}

#### NPI efficiency: marginal value of extending the NPI by one day ####

simulate_npi_attack_curve <- function(
    policy_grid, medians_HS, medians_homog,
    t_max = 1080, population_size = 1e6, i0_prop = 1 / population_size
) {
  run_one <- function(model, medians) {
    beta <- medians$beta[[1]]
    gamma <- medians$gamma[[1]]
    cv <- if ("cv" %in% names(medians)) medians$cv[[1]] else NULL
    purrr::pmap_dfr(policy_grid, function(eff, NPI_start, NPI_dur) {
      m <- simulate_policy_metrics(
        model = model, eff = eff, NPI_start = NPI_start, NPI_dur = NPI_dur,
        beta = beta, gamma = gamma, cv = cv,
        t_max = t_max, population_size = population_size, i0_prop = i0_prop
      )
      tibble::tibble(
        eff = eff, NPI_start = NPI_start, NPI_dur = NPI_dur,
        attack_rate = m$attack_rate, model = model
      )
    })
  }
  dplyr::bind_rows(
    run_one("HS", medians_HS),
    run_one("homog", medians_homog)
  )
}

make_npi_marginal_value_df <- function(attack_curve, population_size = 1e6) {
  attack_curve |>
    dplyr::arrange(model, eff, NPI_start, NPI_dur) |>
    dplyr::group_by(model, eff, NPI_start) |>
    dplyr::mutate(
      d_dur = NPI_dur - dplyr::lag(NPI_dur),
      # infections averted by the marginal day: -dAR/dD * N
      marginal_value_per_day =
        -(attack_rate - dplyr::lag(attack_rate)) / d_dur * population_size
    ) |>
    dplyr::ungroup() |>
    dplyr::filter(!is.na(marginal_value_per_day))
}

plot_npi_marginal_value <- function(
    marginal_df,
    start_keep = c(5, 10, 15, 20),
    dur_max_display = 150
) {
  df <- marginal_df |>
    dplyr::filter(NPI_start %in% start_keep) |>
    dplyr::mutate(
      model = factor(model, levels = c("homog", "HS"),
                     labels = c("Homogeneous", "Heterogeneous")),
      NPI_start_f = factor(NPI_start, levels = start_keep,
                           labels = paste0("Start: day ", start_keep)),
      eff_f = factor(
        eff, levels = sort(unique(eff)),
        labels = paste0(round(sort(unique(eff)) * 100), "%")
      )
    )
  
  ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = NPI_dur, y = marginal_value_per_day,
      color = eff_f, linetype = model,
      group = interaction(eff_f, model)
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 0, color = "black", linetype = "dashed", linewidth = 0.4
    ) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::facet_wrap(~ NPI_start_f) +
    ggplot2::coord_cartesian(xlim = c(0, dur_max_display)) +
    ggplot2::scale_color_viridis_d(option = "plasma", end = 0.85) +
    ggplot2::scale_linetype_manual(
      values = c("Homogeneous" = "solid", "Heterogeneous" = "dashed")
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::label_number(big.mark = ",", accuracy = 1)
    ) +
    ggplot2::labs(
      x = "NPI duration (days)",
      y = "Infections averted per additional NPI day",
      color = "NPI effectiveness",
      linetype = "Susceptibility type"
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold"),
      panel.border = ggplot2::element_rect(color = "grey40", fill = NA, linewidth = 0.6),
      panel.spacing = grid::unit(0.6, "lines"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "top",
      legend.title = ggplot2::element_text(face = "bold"),
      legend.key.width = grid::unit(2.2, "cm")
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(
        order = 1, title.position = "top",
        override.aes = list(linewidth = 1.1, linetype = "solid")
      ),
      linetype = ggplot2::guide_legend(
        order = 2, title.position = "top",
        keywidth = grid::unit(2.2, "cm"),
        override.aes = list(
          linewidth = 1.1,
          linetype = c("solid", "dashed"),
          colour = "grey20"
        )
      )
    )
}

#### Supplementary plots ####

# FTS ridge plot

fts_param_ridge_df <- function(
    bundles,
    files,
    grid_csv,
    params = c("beta", "gamma", "R0", "cv"),
    r0_keep = 3
) {
  grid <- utils::read.csv(grid_csv, stringsAsFactors = FALSE)
  ids <- as.integer(sub(
    "^id", "",
    regmatches(basename(files), regexpr("id[0-9]+", basename(files)))
  ))
  meta <- grid[match(ids, grid$id), , drop = FALSE]
  
  ok <- !is.na(meta$id)
  if (!is.null(r0_keep)) ok <- ok & meta$R0 %in% r0_keep
  bundles <- bundles[ok]
  meta <- meta[ok, , drop = FALSE]
  if (nrow(meta) == 0) stop("No FTS bundles matched the parameter grid.")
  
  meta$setting <- factor(
    meta$param_combo_label,
    levels = c("slow", "reference", "fast"),
    labels = c("Slow", "Reference", "Fast")
  )
  cv_levels <- sort(unique(meta$cv))
  meta$cv_lab <- factor(
    sprintf("cv = %g", meta$cv),
    levels = rev(sprintf("cv = %g", cv_levels))
  )
  
  disp <- c(beta = "beta", gamma = "gamma", R0 = "R[0]", cv = "nu")
  
  d <- purrr::map_dfr(seq_len(nrow(meta)), function(i) {
    dr <- as.data.frame(bundles[[i]]$draws)
    keep_p <- intersect(params, names(dr))
    if (length(keep_p) == 0) {
      stop("None of the requested parameters found in bundle ", i)
    }
    purrr::map_dfr(keep_p, function(p) {
      data.frame(
        setting = meta$setting[i],
        cv_lab = meta$cv_lab[i],
        parameter = disp[[p]],
        value = as.numeric(dr[[p]]),
        true = switch(
          p,
          beta = meta$beta[i],
          gamma = meta$gamma[i],
          R0 = meta$R0[i],
          cv = meta$cv[i],
          NA_real_
        )
      )
    })
  })
  
  d$parameter <- factor(d$parameter, levels = unname(disp[params]))
  d
}

fts_param_ridge_plot <- function(ridge_df) {
  tv <- unique(ridge_df[, c("setting", "cv_lab", "parameter", "true")])
  
  ggplot2::ggplot(ridge_df, ggplot2::aes(x = value, y = cv_lab, fill = cv_lab)) +
    ggridges::geom_density_ridges(
      scale = 1.05, alpha = 0.8, colour = "grey30", rel_min_height = 0.01
    ) +
    ggplot2::geom_point(
      data = tv, inherit.aes = FALSE,
      ggplot2::aes(x = true, y = cv_lab),
      shape = 124, size = 3
    ) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(setting),
      cols = ggplot2::vars(parameter),
      scales = "free_x",
      labeller = ggplot2::labeller(parameter = ggplot2::label_parsed)
    ) +
    ggplot2::scale_fill_viridis_d(guide = "none") +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold"),
      strip.text.y.right = ggplot2::element_text(angle = 270)
    )
}


#### Supplementary tables ####

# convergence diagnostics
.spec_truths <- function(spec) {
  lapply(seq_len(nrow(spec)), function(i) {
    tr <- c(
      beta  = if ("beta"  %in% names(spec)) spec$beta[i]  else NA_real_,
      gamma = if ("gamma" %in% names(spec)) spec$gamma[i] else NA_real_,
      cv    = if ("cv"    %in% names(spec)) spec$cv[i]    else NA_real_
    )
    tr <- c(tr, R0 = unname(tr[["beta"]] / tr[["gamma"]]))
    tr[!is.na(tr)]
  })
}

bundle_convergence_summary <- function(
    path,
    label,
    params = c("beta", "gamma", "cv", "R0"),
    truth = NULL
) {
  b <- readRDS(path)
  dr <- posterior::as_draws_df(as.data.frame(b$draws))
  keep <- intersect(params, posterior::variables(dr))
  if (length(keep) == 0) {
    stop("None of the requested parameters found in ", basename(path))
  }
  
  s <- posterior::summarise_draws(
    posterior::subset_draws(dr, variable = keep),
    "median",
    ~posterior::quantile2(.x, probs = c(0.025, 0.975)),
    posterior::default_convergence_measures()
  )
  
  
  true_vals <- if (is.null(truth)) {
    rep(NA_real_, nrow(s))
  } else {
    unname(truth[s$variable])
  }
  
  s |>
    tibble::as_tibble() |>
    dplyr::mutate(
      fit = label,
      true = true_vals,
      .before = 1
    )
}

make_convergence_df <- function(
    paths,
    labels,
    truths = NULL,
    params = c("beta", "gamma", "cv", "R0")
) {
  if (is.null(truths)) truths <- vector("list", length(paths))
  purrr::pmap_dfr(
    list(paths, labels, truths),
    function(p, l, tr) {
      bundle_convergence_summary(p, l, params = params, truth = tr)
    }
  )
}

convergence_df_pts <- function(spec, paths) {
  make_convergence_df(
    paths,
    sprintf("%s \u2014 %s", spec$setting, spec$window),
    truths = .spec_truths(spec)
  )
}

convergence_df_voi <- function(spec, paths) {
  arm_lab <- ifelse(spec$arm == "HS", "hSIR", "SIR")
  win <- if ("window" %in% names(spec)) {
    as.character(spec$window)
  } else {
    paste("Theta", spec$theta)
  }
  make_convergence_df(
    paths,
    sprintf("%s \u2014 %s", arm_lab, win),
    truths = .spec_truths(spec)
  )
}

convergence_df_npi <- function(spec, paths) {
  lab <- sprintf(
    "NPI \u2014 %g GI post-implementation, eff = %g%%",
    spec$n_intervals, 100 * spec$eff
  )
  if ("cv" %in% names(spec)) lab <- sprintf("%s, cv = %g", lab, spec$cv)
  make_convergence_df(paths, lab, truths = .spec_truths(spec))
}

convergence_df_fts <- function(files, grid_csv, r0_keep = 3) {
  grid <- utils::read.csv(grid_csv, stringsAsFactors = FALSE)
  ids <- as.integer(sub(
    "^id", "",
    regmatches(basename(files), regexpr("id[0-9]+", basename(files)))
  ))
  meta <- grid[match(ids, grid$id), , drop = FALSE]
  ok <- !is.na(meta$id)
  if (!is.null(r0_keep)) ok <- ok & meta$R0 %in% r0_keep
  if (!any(ok)) stop("No FTS bundles matched the parameter grid.")
  
  speed_rank <- match(meta$param_combo_label[ok], c("slow", "reference", "fast"))
  ord <- order(speed_rank, meta$cv[ok])
  paths <- files[ok][ord]
  meta_ok <- meta[ok, , drop = FALSE][ord, , drop = FALSE]
  
  make_convergence_df(
    paths,
    sprintf(
      "%s (R0 = %g), cv = %g",
      meta_ok$param_combo_label, meta_ok$R0, meta_ok$cv
    ),
    truths = .spec_truths(meta_ok)
  )
}

make_convergence_gt <- function(df, title) {
  parm_lab <- c(
    beta = "\u03b2", gamma = "\u03b3", cv = "\u03bd",
    R0 = "R\u2080", D = "D"
  )
  
  df |>
    dplyr::mutate(
      parameter = dplyr::coalesce(parm_lab[variable], variable)
    ) |>
    dplyr::select(
      fit, parameter, true, median, q2.5, q97.5,
      rhat, ess_bulk, ess_tail
    ) |>
    gt::gt(groupname_col = "fit") |>
    gt::cols_label(
      parameter = "Parameter",
      true = "True value",
      median = "Median",
      q2.5 = "2.5%",
      q97.5 = "97.5%",
      rhat = gt::html("R&#770;"),
      ess_bulk = "Bulk ESS",
      ess_tail = "Tail ESS"
    ) |>
    gt::tab_spanner(label = "95% CrI", columns = c(q2.5, q97.5)) |>
    gt::fmt_number(
      columns = c(true, median, q2.5, q97.5, rhat),
      decimals = 3
    ) |>
    gt::fmt_number(
      columns = c(ess_bulk, ess_tail),
      decimals = 0, use_seps = TRUE
    ) |>
    gt::sub_missing(missing_text = "\u2014") |>
    gt::tab_header(title = title) |>
    gtExtras::gt_theme_nytimes()
}

save_gt_table <- function(tbl, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  gt::gtsave(tbl, path)
  path
}
