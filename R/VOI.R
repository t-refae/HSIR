voi_bundles_to_loaded_objects <- function(bundles, files, prefix = NULL, suffix = NULL) {
  regime_to_theta <- c("2GI" = 1L, "1GI" = 2L, "0GI" = 3L, "full" = 4L)
  regimes <- sub("^.*_pts_([^.]+)\\.rds$", "\\1", basename(files))
  stopifnot(all(regimes %in% names(regime_to_theta)))
  thetas <- unname(regime_to_theta[regimes])
  
  nm <- vapply(thetas, function(k) {
    core <- paste0("Theta_", k)
    if (!is.null(prefix)) core <- paste0(prefix, "_", core)
    if (!is.null(suffix)) core <- paste0(core, "_", suffix)
    core
  }, character(1))
  
  objs <- lapply(bundles, function(b) {
    d <- b$draws
    if (!inherits(d, "draws_df")) d <- posterior::as_draws_df(d)
    d
  })
  
  ord <- order(thetas)
  stats::setNames(objs[ord], nm[ord])
}

make_voi_true_params_from_csv <- function(csv_path, voi_thetas = 1:4) {
  regime_to_theta <- c("2GI" = 1L, "1GI" = 2L, "0GI" = 3L, "full" = 4L)
  grid <- utils::read.csv(csv_path, stringsAsFactors = FALSE)
  grid$theta <- unname(regime_to_theta[trimws(as.character(grid$regime))])
  stopifnot(!anyNA(grid$theta), setequal(grid$theta, voi_thetas))
  grid <- grid[order(grid$theta), , drop = FALSE]
  data.frame(
    theta = grid$theta,
    beta  = grid$beta,
    gamma = grid$gamma,
    cv    = grid$cv,
    R0    = grid$beta / grid$gamma
  )
}