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

make_cases_from_cumulative <- function(out, times, population_size) {
  data.frame(
    Day = times,
    Data = c(1, round(diff(out$C) * population_size))
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

#### plotting ####

plot_trajectories_faceted <- function(combined_long) {
  ggplot2::ggplot(
    combined_long,
    ggplot2::aes(x = time, y = Value, color = susceptibility)
  ) +
    ggplot2::geom_line(linewidth = 1.3) +
    ggplot2::facet_wrap(~Compartment, scales = "free_y") +
    ggplot2::scale_color_manual(
      values = c(
        "Homogeneous" = "#1f77b4",
        "Heterogeneous" = "#ff7f0e"
      )
    ) +
    ggplot2::labs(
      x = "Time (Days)",
      y = "Proportion",
      color = "Susceptibility"
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      legend.position = "top",
      plot.title = ggplot2::element_blank(),
      legend.background = ggplot2::element_blank(),
      legend.box.background = ggplot2::element_rect(color = "black")
    )
}

plot_absolute_differences <- function(delta_df) {
  ggplot2::ggplot(
    delta_df,
    ggplot2::aes(
      x = time,
      y = Absolute_Difference,
      linetype = Metric
    )
  ) +
    ggplot2::geom_line(linewidth = 1.3) +
    ggplot2::scale_linetype_manual(
      values = c(
        delta_infected = "solid",
        delta_cumulative_incidence = "dotted"
      ),
      labels = c(
        delta_infected = "Prevalence",
        delta_cumulative_incidence = "Cumulative Incidence"
      )
    ) +
    ggplot2::labs(
      title = expression(paste("|", Homogeneous - Heterogeneous, "|")),
      x = "Time (Days)",
      y = "Absolute Difference",
      linetype = "Metric"
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      legend.position = c(0.1, 0.98),
      legend.justification = c(0, 1),
      legend.background = ggplot2::element_rect(
        fill = "white",
        color = "black"
      ),
      plot.title = ggplot2::element_text(
        face = "bold",
        hjust = 0.5
      ),
      legend.box.margin = ggplot2::margin(6, 6, 6, 6)
    )
}

plot_fit_to_data <- function(
    incidence_summary,
    observed_cases,
    ribbon_fill,
    line_color,
    point_color,
    line_label,
    point_label
) {
  ggplot2::ggplot() +
    ggplot2::geom_ribbon(
      data = incidence_summary,
      ggplot2::aes(x = Day, ymin = Lower, ymax = Upper),
      fill = ribbon_fill,
      alpha = 0.3
    ) +
    ggplot2::geom_line(
      data = incidence_summary,
      ggplot2::aes(
        x = Day,
        y = Median,
        linetype = line_label
      ),
      linewidth = 1.3,
      color = line_color
    ) +
    ggplot2::geom_point(
      data = observed_cases,
      ggplot2::aes(
        x = Day,
        y = Data,
        shape = point_label
      ),
      size = 3,
      color = point_color
    ) +
    ggplot2::scale_shape_manual(values = setNames(1, point_label)) +
    ggplot2::scale_linetype_manual(values = setNames("solid", line_label)) +
    ggplot2::labs(
      x = "Time (Days)",
      y = "Incidence",
      shape = "",
      linetype = ""
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      legend.position = "top",
      legend.box.background = ggplot2::element_rect(color = "black"),
      legend.background = ggplot2::element_blank()
    )
}

assemble_fig_1_panel <- function(fig_1_a, fig_1_b, plot_fit_1, plot_fit_2) {
  (
    (fig_1_a + fig_1_b) /
      (plot_fit_1 + plot_fit_2)
  ) +
    patchwork::plot_layout(guides = "keep") +
    patchwork::plot_annotation(
      tag_levels = "a",
      tag_prefix = "(",
      tag_suffix = ")",
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(
          hjust = 0.5,
          face = "bold"
        )
      )
    ) &
    ggplot2::theme(legend.position = "top")
}

save_ggplot_pdf <- function(plot, path, width, height, dpi = 300) {
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