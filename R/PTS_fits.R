pts_fit_spec <- function(fitting_root = "../HSIR_fitting", cv = 1) {
  cfg <- yaml::read_yaml(file.path(fitting_root, "config.yml"))
  model <- cfg$model_type
  
  sc <- utils::read.csv(file.path(fitting_root, cfg$scenarios_csv), stringsAsFactors = FALSE)
  names(sc) <- tolower(trimws(names(sc)))
  
  settings <- data.frame(
    setting = c("Fast", "Reference", "Slow"),
    beta    = c(1.2, 0.9, 0.6),
    gamma   = c(0.4, 0.3, 0.2),
    stringsAsFactors = FALSE
  )
  
  windows <- data.frame(
    window = c("-2 GI pre-peak", "-1 GI pre-peak", "At peak", "Complete"),
    stage  = c("PTS", "PTS", "PTS", "FTS"),
    regime = c("m2GI", "m1GI", "peak", NA_character_),
    stringsAsFactors = FALSE
  )
  
  spec <- merge(settings, windows, by = character(0))
  spec$cv <- cv
  
  k  <- sprintf("%.6f_%.6f_%.6f", spec$beta, spec$gamma, spec$cv)
  ks <- sprintf("%.6f_%.6f_%.6f", sc$beta, sc$gamma, sc$cv)
  i  <- match(k, ks)
  if (anyNA(i)) {
    stop("no scenario in ", cfg$scenarios_csv, " for:\n  ",
         paste(unique(k[is.na(i)]), collapse = "\n  "), call. = FALSE)
  }
  
  id <- gsub("[^A-Za-z0-9_]", "_", as.character(sc$id[i]))
  spec$id <- ifelse(grepl("^[A-Za-z]", id), id, paste0("id", id))
  
  spec$file <- ifelse(
    spec$stage == "PTS",
    sprintf("pts_fit_%s_%s_%s.rds", model, spec$id, spec$regime),
    sprintf("fit_%s_%s.rds", model, spec$id)
  )
  spec$path <- file.path(fitting_root, "outputs", spec$stage, spec$file)
  
  spec$setting <- factor(spec$setting, levels = settings$setting)
  spec$window  <- factor(spec$window,  levels = windows$window)
  spec <- spec[order(spec$setting, spec$window), ]
  rownames(spec) <- NULL
  spec
}

pts_fit_files <- function(spec) {
  miss <- !file.exists(spec$path)
  if (any(miss)) {
    stop("missing bundles:\n  ", paste(spec$file[miss], collapse = "\n  "), call. = FALSE)
  }
  spec$path
}

pts_cv_draws <- function(spec, paths = spec$path) {
  out <- lapply(seq_along(paths), function(i) {
    d <- as.data.frame(readRDS(paths[i])$draws)
    if (!"cv" %in% names(d)) stop("no cv column in ", basename(paths[i]), call. = FALSE)
    data.frame(
      value   = as.numeric(d$cv),
      setting = spec$setting[i],
      window  = spec$window[i]
    )
  })
  do.call(rbind, out)
}

pts_true_cv <- function(spec) {
  data.frame(
    setting    = spec$setting,
    window     = spec$window,
    true_value = spec$cv
  )
}
