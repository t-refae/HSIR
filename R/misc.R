nu_palette <- function(
    nu_levels,
    palette_fun = function(n) viridisLite::mako(n, begin = 0.3, end = 0.7),
    zero_colour = "grey40"
) {
  nu_levels <- as.character(nu_levels)
  het <- nu_levels[nu_levels != "0"]
  cols <- stats::setNames(as.character(palette_fun(length(het))), het)
  if ("0" %in% nu_levels) cols <- c("0" = zero_colour, cols)
  cols[nu_levels]
}

#### susceptibility contour plot ####
hsir_hit <- function(R0, nu) {
  1 - R0^(-1 / (1 + nu^2))
}

make_susceptibility_dist_df <- function(
    nu_values = c(0.5, 1, 2),
    x_max = 4,
    n = 400
) {
  purrr::map_dfr(nu_values, function(nu) {
    shape <- 1 / nu^2
    x <- seq(0.001, x_max, length.out = n)
    tibble::tibble(
      nu = nu,
      x = x,
      density = stats::dgamma(x, shape = shape, rate = shape)
    )
  }) |>
    dplyr::mutate(nu = factor(nu, levels = nu_values))
}

mean_susceptibility_remaining <- function(q, nu) {
  (1 - q)^(nu^2)
}

make_mean_sus_df <- function(
    q_seq = seq(0, 0.80, by = 0.002),
    nu_seq = seq(0, 3, by = 0.01)
) {
  tidyr::expand_grid(consumption = q_seq, nu = nu_seq) |>
    dplyr::mutate(
      mean_sus  = mean_susceptibility_remaining(consumption, nu),
      rt_factor = (1 - consumption)^(1 + nu^2)
    )
}

make_iso_label_df <- function(
    R0_vals = c(1.5, 3, 4.5),
    nu_at = c(1.5, 1.0, 0.6),
    dx = c(0.055, 0, 0),
    dy = c(0, -0.40, 0.40),
    q_max = 0.80,
    nu_max = 3
) {
  stopifnot(length(R0_vals) == length(nu_at))
  dx <- rep_len(dx, length(R0_vals))
  dy <- rep_len(dy, length(R0_vals))
  purrr::pmap_dfr(
    list(R0_vals, nu_at, dx, dy),
    function(R0, nu, ox, oy) {
      tibble::tibble(
        consumption = pmin(pmax(hsir_hit(R0, nu) + ox, 0.03), q_max - 0.03),
        nu = pmin(pmax(nu + oy, 0.08), nu_max - 0.08),
        label = sprintf("R[0] * {phantom() == phantom()} * %s", format(R0, trim = TRUE))
      )
    }
  )
}

make_binned_fill <- function(
    breaks = seq(0, 1, by = 0.1),
    palette_fun = function(n) paletteer::paletteer_c("grDevices::Lajolla", n),
    direction = 1,
    na_colour = NULL,
    limits = range(breaks),
    bar_width = grid::unit(0.35, "cm"),
    bar_height = grid::unit(7, "cm"),
    show_limits = TRUE,
    aesthetics = "fill"
) {
  cols_for <- function(n) {
    cols <- as.character(palette_fun(n))
    if (direction == -1) cols <- rev(cols)
    cols
  }
  n_bins <- max(length(breaks) - 1L, 1L)
  fallback <- if (is.null(na_colour)) cols_for(n_bins)[n_bins] else na_colour
  
  ggplot2::binned_scale(
    aesthetics = aesthetics,
    palette = function(x) cols_for(length(x)),
    breaks = breaks,
    limits = limits,
    oob = scales::squish,
    na.value = fallback,
    guide = ggplot2::guide_colorsteps(
      barwidth = bar_width,
      barheight = bar_height,
      show.limits = show_limits
    )
  )
}

plot_susceptibility_and_hit <- function(
    dist_df,
    mean_sus_df,
    R0_vals = c(1.5, 3, 4.5),
    nu_at = c(1.5, 1.0, 0.6),
    label_dx = c(0.055, 0, 0),
    label_dy = c(0, -0.40, 0.40),
    breaks = seq(0, 1, by = 0.1),
    palette_fun = function(n) paletteer::paletteer_c("grDevices::ag_Sunset", n),
    palette_direction = -1,
    na_colour = NULL,
    fill_lab = "Mean\nsusceptibility",
    contour_colour = "black",
    label_fill = "white",
    label_alpha = 0.85,
    base_size = 12
) {
  q_max <- max(mean_sus_df$consumption)
  nu_max <- max(mean_sus_df$nu)
  iso_labs <- make_iso_label_df(
    R0_vals, nu_at, label_dx, label_dy, q_max = q_max, nu_max = nu_max
  )
  
  p_a <- ggplot2::ggplot(dist_df, ggplot2::aes(x = x, y = density, colour = nu)) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
    ggplot2::geom_line(linewidth = 1.1) +
    ggplot2::annotate(
      "text", x = 1.05, y = 1.9,
      label = "\u03bd = 0 (SIR)", hjust = 0, colour = "grey40", size = 3.5
    ) +
    ggplot2::coord_cartesian(ylim = c(0, 2)) +
    ggplot2::scale_colour_manual(
      values = nu_palette(levels(dist_df$nu)),
      labels = function(l) paste0("\u03bd = ", l)
    ) +
    ggplot2::guides(colour = ggplot2::guide_legend(reverse = TRUE)) +
    ggplot2::labs(
      x = "Relative susceptibility to infection",
      y = "Density", colour = NULL, tag = "a"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position = "inside",
      legend.position.inside = c(0.98, 0.98),
      legend.justification = c(1, 1),
      legend.background = ggplot2::element_rect(
        fill = ggplot2::alpha("white", 0.75), colour = NA
      ),
      legend.key.height = grid::unit(0.9, "lines"),
      legend.margin = ggplot2::margin(2, 4, 2, 4)
    )
  
  fill_scale <- make_binned_fill(
    breaks = breaks,
    palette_fun = palette_fun,
    direction = palette_direction,
    na_colour = na_colour
  )
  
  p_b <- ggplot2::ggplot(mean_sus_df, ggplot2::aes(x = consumption, y = nu)) +
    metR::geom_contour_fill(ggplot2::aes(z = mean_sus), breaks = breaks) +
    fill_scale +
    ggplot2::geom_contour(
      ggplot2::aes(z = rt_factor), breaks = sort(1 / R0_vals),
      colour = "black", linewidth = 0.8
    ) +
    ggplot2::geom_label(
      data = iso_labs,
      ggplot2::aes(x = consumption, y = nu, label = label),
      colour = "black", size = 3.4, fontface = "bold", parse = TRUE,
      fill = "white", alpha = 0.85, label.size = 0,
      label.padding = grid::unit(0.15, "lines"),
      inherit.aes = FALSE
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq(0, q_max, by = 0.1),
      labels = scales::percent
    ) +
    ggplot2::labs(
      x = "Proportion of the population infected",
      y = expression("Coefficient of variation," ~ nu),
      fill = "Mean\nsusceptibility",
      tag = "b"
    ) +
    ggplot2::coord_cartesian(
      xlim = c(0, q_max),
      ylim = range(mean_sus_df$nu),
      expand = FALSE
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
  
  patchwork::wrap_plots(p_a, p_b, ncol = 2, widths = c(1, 1.15))
}
#### infectious dynamics figure across settings ####

simulate_hsir_prevalence <- function(
    beta, gamma, cv,
    population_size = 1e6,
    t_max = 800,
    dt = 1
) {
  i0 <- 1 / population_size
  derivs <- function(t, y, p) {
    S <- y[["S"]]; I <- y[["I"]]
    force <- beta * I * S^(1 + cv^2)
    list(c(-force, force - gamma * I, gamma * I, force))
  }
  out <- deSolve::ode(
    y = c(S = 1 - i0, I = i0, R = 0, C = 0),
    times = seq(0, t_max, by = dt),
    func = derivs,
    method = "lsoda"
  )
  as.data.frame(out)
}

make_I_dynamics_df <- function(
    grid_csv,
    population_size = 1e6,
    t_max = 800,
    tail_quantile = 0.995,
    tail_pad = 14
) {
  grid <- utils::read.csv(grid_csv, stringsAsFactors = FALSE)
  
  traj <- purrr::map_dfr(seq_len(nrow(grid)), function(i) {
    p <- grid[i, ]
    d <- simulate_hsir_prevalence(
      p$beta, p$gamma, p$cv,
      population_size = population_size, t_max = t_max
    )
    cum_inc <- d$S[1] - d$S
    hit <- which(cum_inc >= tail_quantile * cum_inc[length(cum_inc)])[1]
    idx <- if (is.na(hit)) nrow(d) else min(hit + tail_pad, nrow(d))
    data.frame(
      time = d$time[seq_len(idx)],
      I = d$I[seq_len(idx)],
      beta = p$beta,
      gamma = p$gamma,
      R0 = p$R0,
      cv = p$cv
    )
  })
  
  traj$cv <- factor(traj$cv, levels = sort(unique(grid$cv)))
  traj
}

plot_I_dynamics <- function(traj, base_size = 9) {
  ggplot2::ggplot(traj, ggplot2::aes(time, I, colour = cv)) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::facet_grid(
      R0 ~ gamma,
      scales = "free",
      labeller = ggplot2::label_bquote(
        rows = R[0] == .(R0),
        cols = gamma == .(gamma)
      )
    ) +
    ggplot2::scale_colour_viridis_d(name = expression(nu), end = 0.9) +
    ggplot2::labs(x = "Time (days)", y = expression(I(t))) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey92"),
      legend.position = "right",
      axis.title.y = ggplot2::element_text(angle = 0, vjust = 0.5)
    )
}

#### Supplementary figure: strength of selection by infection ####

simulate_hsir_selection <- function(
    cv_values = c(0, 0.5, 1, 2),
    beta = 1.2,
    gamma = 0.4,
    i0 = 1e-4,
    times = seq(0, 400, by = 0.1)
) {
  R0 <- beta / gamma
  
  rhs <- function(t, y, p) {
    with(as.list(c(y, p)), {
      foi <- beta * I * S^(cv^2)
      list(c(S = -foi * S, I = foi * S - gamma * I))
    })
  }
  
  purrr::map_dfr(cv_values, function(cv) {
    out <- as.data.frame(deSolve::ode(
      y = c(S = 1 - i0, I = i0),
      times = times,
      func = rhs,
      parms = c(beta = beta, gamma = gamma, cv = cv),
      method = "lsoda"
    ))
    out$cv <- cv
    out$Re <- R0 * out$S^(1 + cv^2)
    out$dRe_dS <- R0 * (1 + cv^2) * out$S^(cv^2)
    out$amp <- (1 + cv^2) * out$S^(cv^2)
    out$Var <- cv^2 * out$S^(2 * cv^2)
    out
  }) |>
    dplyr::mutate(cvf = factor(cv, levels = sort(unique(cv_values))))
}

plot_hsir_selection <- function(
    sim,
    beta = 1.2,
    gamma = 0.4,
    x_max_display = 50,
    base_size = 9
) {
  R0 <- beta / gamma
  het <- dplyr::filter(sim, cv > 0)
  
  het_levels <- levels(droplevels(het$cvf))
  het_levels <- levels(droplevels(het$cvf))
  palf <- nu_palette(levels(sim$cvf))
  pal <- palf[het_levels]
  
  ltys <- stats::setNames(
    c("dashed", rep("solid", length(het_levels))),
    levels(sim$cvf)
  )
  
  ends <- het |>
    dplyr::group_by(cv, cvf) |>
    dplyr::slice_tail(n = 1) |>
    dplyr::ungroup()
  
  sel_theme <- list(
    ggplot2::scale_x_reverse(),
    ggplot2::scale_colour_manual(values = pal, name = expression(nu)),
    ggplot2::theme_minimal(base_size = base_size),
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "none"
    )
  )
  end_pts <- ggplot2::geom_point(
    data = ends, shape = 21, fill = "white", size = 2.2, stroke = 0.8
  )
  
  p_S <- ggplot2::ggplot(sim, ggplot2::aes(time, S, colour = cvf, linetype = cvf)) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::scale_colour_manual(values = palf, name = expression(nu)) +
    ggplot2::scale_linetype_manual(values = ltys, name = expression(nu)) +
    ggplot2::coord_cartesian(xlim = c(0, x_max_display), ylim = c(0, 1)) +
    ggplot2::labs(x = "Time (days)", y = "S(t)", tag = "a") +
    ggplot2::guides(
      colour = ggplot2::guide_legend(reverse = TRUE),
      linetype = ggplot2::guide_legend(reverse = TRUE)
    ) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
  
  p_A <- ggplot2::ggplot(het, ggplot2::aes(S, dRe_dS, colour = cvf)) +
    ggplot2::geom_hline(yintercept = R0, linetype = "dashed", linewidth = 0.5) +
    ggplot2::geom_line(linewidth = 0.9) + end_pts + sel_theme +
    ggplot2::labs(
      x = "S (epidemic progress \u2192)",
      y = expression(d*R[e]/d*S),
      tag = "b"
    )
  
  p_B <- ggplot2::ggplot(het, ggplot2::aes(S, amp, colour = cvf)) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.5) +
    ggplot2::geom_line(linewidth = 0.9) + end_pts + sel_theme +
    ggplot2::labs(
      x = "S (epidemic progress \u2192)",
      y = expression((1 + nu^2) * S^{nu^2}),
      tag = "c"
    )
  
  p_V <- ggplot2::ggplot(het, ggplot2::aes(S, Var, colour = cvf)) +
    ggplot2::geom_line(linewidth = 0.9) + end_pts + sel_theme +
    ggplot2::labs(
      x = "S (epidemic progress \u2192)",
      y = expression(nu^2 * S^{2 * nu^2}),
      tag = "d"
    )
  
  top <- patchwork::plot_spacer() + p_S + patchwork::plot_spacer() +
    patchwork::plot_layout(widths = c(1, 5, 1))
  
  top / (p_A | p_B | p_V) + patchwork::plot_layout(heights = c(1, 1.18))
}
#### Supporting tables (Excel export) ####

save_convergence_workbook <- function(dfs, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  
  tidy_one <- function(df) {
    parm_lab <- c(beta = "beta", gamma = "gamma", cv = "nu", R0 = "R0", D = "D")
    df |>
      dplyr::mutate(parameter = dplyr::coalesce(parm_lab[variable], variable)) |>
      dplyr::select(
        Fit = fit,
        Parameter = parameter,
        `True value` = true,
        Median = median,
        `CrI 2.5%` = q2.5,
        `CrI 97.5%` = q97.5,
        Rhat = rhat,
        `Bulk ESS` = ess_bulk,
        `Tail ESS` = ess_tail
      ) |>
      dplyr::mutate(dplyr::across(
        c(`True value`, Median, `CrI 2.5%`, `CrI 97.5%`, Rhat),
        ~round(.x, 3)
      ))
  }
  
  wb <- openxlsx::createWorkbook()
  hdr <- openxlsx::createStyle(textDecoration = "bold", halign = "center")
  
  for (nm in names(dfs)) {
    openxlsx::addWorksheet(wb, nm)
    tidy <- tidy_one(dfs[[nm]])
    openxlsx::writeData(wb, nm, tidy, headerStyle = hdr)
    openxlsx::freezePane(wb, nm, firstRow = TRUE)
    openxlsx::setColWidths(wb, nm, cols = seq_len(ncol(tidy)), widths = "auto")
  }
  
  openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
  path
}

#### Source data export ####

save_source_data_workbook <- function(sheets, path, max_rows = 1e6) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  wb <- openxlsx::createWorkbook()
  hdr <- openxlsx::createStyle(textDecoration = "bold")
  
  for (nm in names(sheets)) {
    d <- as.data.frame(sheets[[nm]])
    if (nrow(d) > max_rows) {
      stop("Sheet '", nm, "' has ", nrow(d), " rows, exceeding the sheet limit.")
    }
    sheet <- substr(gsub("[\\\\/*?:\\[\\]]", "_", nm), 1, 31)
    openxlsx::addWorksheet(wb, sheet)
    openxlsx::writeData(wb, sheet, d, headerStyle = hdr)
    openxlsx::freezePane(wb, sheet, firstRow = TRUE)
  }
  
  openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
  path
}