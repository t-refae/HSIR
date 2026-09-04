npi_fit_spec <- function(fitting_root = "../HSIR_fitting") {
  cfg <- yaml::read_yaml(file.path(fitting_root, "config.yml"))
  sc <- utils::read.csv(
    file.path(fitting_root, cfg$npi_scenarios_csv),
    stringsAsFactors = FALSE, strip.white = TRUE
  )
  names(sc) <- tolower(trimws(names(sc)))
  model <- tools::file_path_sans_ext(basename(cfg$npi_stan))
  
  sc <- sc[order(sc$eff, sc$n_intervals), , drop = FALSE]
  sc$theta <- seq_len(nrow(sc))
  sc$file  <- sprintf("npi_fit_%s_%s.rds", model, sc$id)
  sc$path  <- file.path(fitting_root, "outputs", "NPI", sc$file)
  rownames(sc) <- NULL
  sc
}

npi_fit_files <- function(spec) {
  miss <- !file.exists(spec$path)
  if (any(miss)) {
    stop("missing bundles:\n  ", paste(spec$file[miss], collapse = "\n  "), call. = FALSE)
  }
  spec$path
}

npi_cv_draws <- function(spec, paths = spec$path) {
  out <- lapply(seq_along(paths), function(i) {
    d <- as.data.frame(readRDS(paths[i])$draws)
    if (!"cv" %in% names(d)) stop("no cv column in ", basename(paths[i]), call. = FALSE)
    data.frame(
      cv    = as.numeric(d$cv),
      theta = spec$theta[i],
      eff   = spec$eff[i],
      GI    = spec$n_intervals[i]
    )
  })
  do.call(rbind, out)
}
