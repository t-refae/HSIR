voi_fit_spec <- function(fitting_root = "../HSIR_fitting") {
  cfg <- yaml::read_yaml(file.path(fitting_root, "config.yml"))
  sc <- utils::read.csv(file.path(fitting_root, cfg$voi_scenarios_csv), stringsAsFactors = FALSE)
  
  thetas <- data.frame(
    regime = c("2GI", "1GI", "0GI", "full"),
    theta  = 1:4,
    window = c("-2 GI pre-peak", "-1 GI pre-peak", "At peak", "Complete"),
    stringsAsFactors = FALSE
  )
  
  arms <- data.frame(
    arm   = c("HS", "Homog"),
    model = c(tools::file_path_sans_ext(basename(cfg$voi_stan_hsir)),
              tools::file_path_sans_ext(basename(cfg$voi_stan_sir))),
    stringsAsFactors = FALSE
  )
  
  miss <- setdiff(thetas$regime, sc$regime)
  if (length(miss)) {
    stop("regimes absent from ", cfg$voi_scenarios_csv, ": ",
         paste(miss, collapse = ", "), call. = FALSE)
  }
  
  spec <- merge(arms, merge(sc, thetas, by = "regime"), by = character(0))
  spec$file <- sprintf("voi_fit_%s_%s.rds", spec$model, spec$id)
  spec$path <- file.path(fitting_root, "outputs", "VOI", spec$file)
  spec$key  <- paste(spec$arm, spec$theta, sep = "_")
  spec$arm    <- factor(spec$arm, levels = arms$arm)
  spec$window <- factor(spec$window, levels = thetas$window)
  spec <- spec[order(spec$arm, spec$theta), ]
  rownames(spec) <- NULL
  spec
}

voi_fit_files <- function(spec) {
  miss <- !file.exists(spec$path)
  if (any(miss)) {
    stop("missing bundles:\n  ", paste(spec$file[miss], collapse = "\n  "), call. = FALSE)
  }
  spec$path
}

voi_param_draws <- function(spec, paths = spec$path) {
  out <- lapply(seq_along(paths), function(i) {
    d <- as.data.frame(readRDS(paths[i])$draws)
    if (!"gamma" %in% names(d) && "D" %in% names(d)) d$gamma <- 1 / d$D
    if (!"R0" %in% names(d) && all(c("beta", "D") %in% names(d))) d$R0 <- d$beta * d$D
    if (!"cv" %in% names(d)) d$cv <- NA_real_
    data.frame(
      arm = spec$arm[i], theta = spec$theta[i], window = spec$window[i],
      .draw = seq_len(nrow(d)),
      beta = d$beta, gamma = d$gamma, cv = d$cv, R0 = d$R0
    )
  })
  do.call(rbind, out)
}

voi_state_draws <- function(spec, paths = spec$path) {
  out <- lapply(paths, function(p) {
    d  <- posterior::as_draws_df(as.data.frame(readRDS(p)$draws))
    yv <- grep("^y\\[", posterior::variables(d), value = TRUE)
    if (!length(yv)) stop("no y[...] columns in ", basename(p), call. = FALSE)
    posterior::as_draws_array(posterior::subset_draws(d, variable = yv))
  })
  names(out) <- spec$key
  out
}

voi_true_params <- function(spec) {
  s <- spec[spec$arm == "HS", ]
  data.frame(
    Theta_num = s$theta,
    beta = s$beta, gamma = s$gamma, cv = s$cv,
    R0 = s$beta / s$gamma
  )
}

voi_parameter_draws_df <- function(params, vars = c("beta", "gamma", "cv", "R0")) {
  d <- params[params$arm == "HS", ]
  theta_levels <- rev(paste0("Theta[", sort(unique(d$theta)), "]"))
  param_labels <- c(beta = "beta", gamma = "gamma", cv = "nu", R0 = "R[0]")
  out <- do.call(rbind, lapply(vars, function(p) {
    data.frame(
      value = d[[p]],
      theta_label = factor(paste0("Theta[", d$theta, "]"), levels = theta_levels),
      param = p
    )
  }))
  out$param_label <- factor(out$param, levels = names(param_labels), labels = param_labels)
  out
}

voi_x_stars <- function(params, arm) {
  d <- params[params$arm == arm, ]
  do.call(rbind, lapply(sort(unique(d$theta)), function(th) {
    r <- d[d$theta == th, ]
    x <- if (arm == "HS") (1 / r$R0)^(1 / (1 + r$cv^2)) else 1 / r$R0
    b <- stats::quantile(x, probs = c(0.05, 0.95), na.rm = TRUE)
    data.frame(Theta = th, median = stats::median(x, na.rm = TRUE),
               lower = unname(b[[1]]), upper = unname(b[[2]]))
  }))
}

voi_medians <- function(params, arm, theta = 4) {
  d <- params[params$arm == arm & params$theta == theta, ]
  if (!nrow(d)) stop("no VOI draws for arm=", arm, " theta=", theta, call. = FALSE)
  out <- data.frame(beta = stats::median(d$beta, na.rm = TRUE),
                    gamma = stats::median(d$gamma, na.rm = TRUE))
  if (arm == "HS") out$cv <- stats::median(d$cv, na.rm = TRUE)
  out
}

voi_S_summary <- function(states, spec, arm, model_type) {
  s <- spec[spec$arm == arm, ]
  do.call(rbind, lapply(seq_len(nrow(s)), function(i) {
    SIR_sum <- extract_SIR_summary(states[[s$key[i]]],
                                   theta_label = paste0("Theta_", s$theta[i]))
    names(SIR_sum)[1:3] <- c("median", "lower", "upper")
    SIR_sum |>
      dplyr::filter(state == "S") |>
      dplyr::select(time, median, lower, upper) |>
      dplyr::mutate(theta = paste0("Theta_", s$theta[i]),
                    model_type = model_type,
                    Theta_num = s$theta[i])
  }))
}
