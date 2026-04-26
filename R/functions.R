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

#### dynamic speed / partial time-series cv inference ####
make_dynamic_speed_comp_dt <- function() {
  # fast, ref, slow (in order of beta)
  
  beta_vals <- c(1.2, 0.9, 0.6)
  gamma_vals <- c(0.4, 0.3, 0.2)
  R0_val <- 3 # fixed here
  n_intervals <- c(-2, -1, 0, 1000) # partial time-series cropping (last is full TS)
  cv_val <- 1
  theta_vals <- seq(1, length(beta_vals)*length(n_intervals))
  
  out <- data.table(Theta = paste0("Theta_{", theta_vals, "}"),
                    beta = rep(beta_vals, each=length(n_intervals)),
                    gamma = rep(gamma_vals, each=length(n_intervals)),
                    cv = rep(cv_val, length(theta_vals)),
                    n_intervals = rep(n_intervals, times = length(beta_vals)),
                    R0 = rep(R0_val, length(n_intervals))
  )
  
  out
}

theta_to_dynamic_speed_setting <- function(theta_val) {
  dplyr::case_when(
    theta_val %in% 1:4 ~ "Fast",
    theta_val %in% 5:8 ~ "Reference",
    theta_val %in% 9:12 ~ "Slow",
    TRUE ~ NA_character_
  )
}

find_dynamic_speed_matched_folders <- function(parent_dir, desired_thetas) {
  all_files <- list.files(parent_dir)
  
  matched_folders <- all_files[
    grepl(
      paste(desired_thetas, collapse = "|"),
      all_files
    )
  ]
  
  if (length(matched_folders) == 0) {
    stop(
      "No dynamic-speed folders matched the requested theta patterns in: ",
      parent_dir
    )
  }
  
  matched_folders
}

list_dynamic_speed_rdata_files <- function(parent_dir, matched_folders) {
  rdata_files <- list.files(
    file.path(parent_dir, matched_folders),
    pattern = "\\.RData$",
    full.names = TRUE
  )
  
  if (length(rdata_files) == 0) {
    stop(
      "No .RData files found under matched dynamic-speed folders in: ",
      parent_dir
    )
  }
  
  rdata_files
}

load_dynamic_speed_rdata_objects <- function(rdata_files, desired_thetas) {
  loaded_objects <- list()
  
  for (theta_name in desired_thetas) {
    scoped_files <- rdata_files[grep(theta_name, rdata_files)]
    
    if (length(scoped_files) == 0) {
      warning("No .RData files found for ", theta_name)
      next
    }
    
    for (file in scoped_files) {
      temp_env <- new.env(parent = emptyenv())
      loaded_names <- load(file, envir = temp_env)
      
      for (loaded_name in loaded_names) {
        new_name <- paste0(theta_name, "_", loaded_name)
        loaded_objects[[new_name]] <- get(loaded_name, envir = temp_env)
      }
    }
  }
  
  if (length(loaded_objects) == 0) {
    stop("No dynamic-speed RData objects were loaded.")
  }
  
  loaded_objects
}

make_dynamic_speed_true_params <- function(dt, theta) {
  true_params <- as.data.frame(dt)
  true_params$theta_raw <- theta
  true_params
}

make_dynamic_speed_cv_draws_df <- function(
    loaded_objects,
    desired_thetas,
    window_levels,
    setting_levels = c("Fast", "Reference", "Slow")
) {
  purrr::map_dfr(seq_along(desired_thetas), function(theta_val) {
    obj_name <- paste0(desired_thetas[theta_val], "_cv_draws")
    
    if (!obj_name %in% names(loaded_objects)) {
      warning("Missing object: ", obj_name)
      return(NULL)
    }
    
    samples <- as.vector(loaded_objects[[obj_name]])
    window_idx <- ((theta_val - 1) %% length(window_levels)) + 1
    
    tibble::tibble(
      value = samples,
      Theta = theta_val,
      setting = factor(
        theta_to_dynamic_speed_setting(theta_val),
        levels = setting_levels
      ),
      window = factor(
        window_levels[window_idx],
        levels = window_levels
      )
    )
  })
}

make_dynamic_speed_true_cv_df <- function(
    true_params,
    window_levels,
    setting_levels = c("Fast", "Reference", "Slow")
) {
  true_params |>
    dplyr::select(Theta, true_value = cv, theta_raw) |>
    dplyr::mutate(
      setting = factor(
        theta_to_dynamic_speed_setting(theta_raw),
        levels = setting_levels
      ),
      window_idx = ((theta_raw - 1) %% length(window_levels)) + 1,
      window = factor(
        window_levels[window_idx],
        levels = window_levels
      )
    )
}

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

#### Fig. 3: S/HIT differences and VOI-probability summaries ####

make_voi_true_params <- function(true_params, voi_thetas = 1:4) {
  true_params <- as.data.frame(true_params)
  
  if ("Theta" %in% names(true_params)) {
    theta_chr <- as.character(true_params$Theta)
    
    theta_num <- suppressWarnings(as.integer(theta_chr))
    
    if (all(is.na(theta_num))) {
      theta_num <- suppressWarnings(
        as.integer(gsub("[^0-9]", "", theta_chr))
      )
    }
    
    true_params$Theta_num <- theta_num
  } else {
    true_params$Theta_num <- seq_len(nrow(true_params))
  }
  
  out <- true_params |>
    dplyr::filter(Theta_num %in% voi_thetas) |>
    dplyr::arrange(Theta_num)
  
  required_cols <- c("beta", "gamma", "cv", "R0")
  missing_cols <- setdiff(required_cols, names(out))
  
  if (length(missing_cols) > 0) {
    stop(
      "VOI true parameter table is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  out |>
    dplyr::select(Theta_num, beta, gamma, cv, R0)
}

list_voi_rdata_files <- function(parent_dir, desired_thetas) {
  if (!dir.exists(parent_dir)) {
    stop("Directory does not exist: ", parent_dir)
  }
  
  all_entries <- list.files(parent_dir, full.names = FALSE)
  
  matched_folders <- all_entries[
    grepl(paste(desired_thetas, collapse = "|"), all_entries)
  ]
  
  if (length(matched_folders) == 0) {
    stop(
      "No folders matched ",
      paste(desired_thetas, collapse = ", "),
      " under ",
      parent_dir
    )
  }
  
  rdata_files <- list.files(
    file.path(parent_dir, matched_folders),
    pattern = "\\.RData$",
    full.names = TRUE,
    recursive = TRUE
  )
  
  if (length(rdata_files) == 0) {
    stop("No .RData files found under matched folders in: ", parent_dir)
  }
  
  sort(rdata_files)
}

load_voi_rdata_objects <- function(
    rdata_files,
    desired_thetas,
    prefix = NULL,
    suffix = NULL
) {
  loaded_objects <- list()
  
  for (theta_name in desired_thetas) {
    scoped_files <- rdata_files[grep(theta_name, rdata_files)]
    
    if (length(scoped_files) == 0) {
      warning("No .RData files found for ", theta_name)
      next
    }
    
    for (file in scoped_files) {
      temp_env <- new.env(parent = emptyenv())
      loaded_names <- load(file, envir = temp_env)
      
      for (loaded_name in loaded_names) {
        name_parts <- c(prefix, theta_name, suffix, loaded_name)
        name_parts <- name_parts[!is.null(name_parts) & !is.na(name_parts)]
        new_name <- paste(name_parts, collapse = "_")
        
        loaded_objects[[new_name]] <- get(loaded_name, envir = temp_env)
      }
    }
  }
  
  if (length(loaded_objects) == 0) {
    stop("No VOI .RData objects were loaded.")
  }
  
  loaded_objects
}

make_voi_parameter_draws_df <- function(
    loaded_objects,
    desired_thetas,
    voi_thetas,
    vars = c("beta", "gamma", "cv", "R0")
) {
  theta_levels <- rev(paste0("Theta[", voi_thetas, "]"))
  
  out <- purrr::map_dfr(seq_along(desired_thetas), function(theta_idx) {
    theta_name <- desired_thetas[theta_idx]
    theta_val <- voi_thetas[theta_idx]
    
    purrr::map_dfr(vars, function(param) {
      obj_name <- paste0("VOI_", theta_name, "_", param, "_draws")
      
      if (!obj_name %in% names(loaded_objects)) {
        warning("Missing object: ", obj_name)
        return(NULL)
      }
      
      tibble::tibble(
        value = as.vector(loaded_objects[[obj_name]]),
        theta_label = factor(
          paste0("Theta[", theta_val, "]"),
          levels = theta_levels
        ),
        param = param
      )
    })
  })
  
  param_labels <- c(
    beta = "beta",
    gamma = "gamma",
    cv = "nu",
    R0 = "R[0]"
  )
  
  out |>
    dplyr::mutate(
      param_label = factor(
        param,
        levels = names(param_labels),
        labels = param_labels
      )
    )
}

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

summarize_hs_x_stars <- function(loaded_objects, voi_thetas = 1:4) {
  purrr::map_dfr(voi_thetas, function(theta_val) {
    R0_name <- paste0("VOI_Theta_", theta_val, "_R0_draws")
    cv_name <- paste0("VOI_Theta_", theta_val, "_cv_draws")
    
    if (!R0_name %in% names(loaded_objects)) {
      stop("Missing object: ", R0_name)
    }
    
    if (!cv_name %in% names(loaded_objects)) {
      stop("Missing object: ", cv_name)
    }
    
    R0_temp_draws <- as.vector(loaded_objects[[R0_name]])
    cv_temp_draws <- as.vector(loaded_objects[[cv_name]])
    
    x_stars <- (1 / R0_temp_draws)^(1 / (1 + cv_temp_draws^2))
    bounds <- stats::quantile(x_stars, probs = c(0.05, 0.95), na.rm = TRUE)
    
    tibble::tibble(
      Theta = theta_val,
      median = stats::median(x_stars, na.rm = TRUE),
      lower = unname(bounds[[1]]),
      upper = unname(bounds[[2]])
    )
  })
}

summarize_homog_x_stars <- function(loaded_objects, voi_thetas = 1:4) {
  purrr::map_dfr(voi_thetas, function(theta_val) {
    R0_name <- paste0("Theta_", theta_val, "_homog_R0_draws")
    
    if (!R0_name %in% names(loaded_objects)) {
      stop("Missing object: ", R0_name)
    }
    
    R0_temp_draws <- as.vector(loaded_objects[[R0_name]])
    
    x_stars <- 1 / R0_temp_draws
    bounds <- stats::quantile(x_stars, probs = c(0.05, 0.95), na.rm = TRUE)
    
    tibble::tibble(
      Theta = theta_val,
      median = stats::median(x_stars, na.rm = TRUE),
      lower = unname(bounds[[1]]),
      upper = unname(bounds[[2]])
    )
  })
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

summarize_S_by_model <- function(
    loaded_objects,
    voi_thetas = 1:4,
    model_type,
    object_prefix = NULL,
    object_suffix = NULL
) {
  purrr::map_dfr(voi_thetas, function(theta_val) {
    name_parts <- c(
      object_prefix,
      paste0("Theta_", theta_val),
      object_suffix,
      "y_draws"
    )
    name_parts <- name_parts[!is.null(name_parts) & !is.na(name_parts)]
    obj_name <- paste(name_parts, collapse = "_")
    
    if (!obj_name %in% names(loaded_objects)) {
      stop("Missing object: ", obj_name)
    }
    
    SIR_sum <- extract_SIR_summary(
      y_draws = loaded_objects[[obj_name]],
      theta_label = paste0("Theta_", theta_val)
    )
    
    names(SIR_sum)[1:3] <- c("median", "lower", "upper")
    
    SIR_sum |>
      dplyr::filter(state == "S") |>
      dplyr::select(time, median, lower, upper) |>
      dplyr::mutate(
        theta = paste0("Theta_", theta_val),
        model_type = model_type,
        Theta_num = theta_val
      )
  })
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
      x = "Time (Days)",
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