#' Resolve and validate the summary selection before scheduling or acquisition
#' @param sim Simulation with bird_summary parameters.
#' @return Validated configuration list.
#' @keywords internal
birdSummaryConfig <- function(sim) {
  p <- params(sim)$bird_summary
  if (!is.character(p$mode) || length(p$mode) != 1L || is.na(p$mode) ||
      !p$mode %in% c("single", "multi")) {
    stop("bird_summary: mode must be single or multi.")
  }
  config <- list(mode = p$mode)
  for (field in c("scenario", "climateScenario", ".studyAreaName")) {
    value <- p[[field]]
    if (!is.character(value) || length(value) != 1L || is.na(value) ||
        !grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$", value) || value %in% c(".", "..")) {
      stop("bird_summary: supply a safe, nonempty ", field, " label.")
    }
    config[[if (field == ".studyAreaName") "studyAreaName" else field]] <- value
  }
  predictor <- params(sim)$bird_modelPredict
  years <- p$years
  if (p$mode == "single" && length(years) == 1L && is.na(years)) {
    if (is.null(predictor)) stop("bird_summary: supply years without a predictor calendar.")
    calendar <- c(predictor$predictStartYear, predictor$predictEndYear, predictor$predictInterval)
    if (length(calendar) != 3L || !is.numeric(calendar) || anyNA(calendar) ||
        any(!is.finite(calendar)) || calendar[3] <= 0 || calendar[1] > calendar[2]) {
      stop("bird_summary: invalid predictor calendar; supply explicit years.")
    }
    years <- seq(calendar[1], calendar[2], by = calendar[3])
  }
  if (!is.numeric(years) || !length(years) || anyNA(years) ||
      any(!is.finite(years)) || any(years != floor(years) | years < 1000 | years > 9999) ||
      anyDuplicated(years)) {
    stop("bird_summary: years must be explicit, unique four-digit simulation years.")
  }
  config$years <- sort(as.integer(years))
  if (p$mode == "single" && any(config$years < start(sim) | config$years > end(sim))) {
    stop("bird_summary: save years must lie inside the simulation calendar.")
  }
  species <- p$species
  if (p$mode == "single" && is.character(species) && !length(species)) {
    index <- sim$modelIndex
    if (is.null(index) || !"species" %in% names(index)) {
      stop("bird_summary: supply species or a modelIndex with species codes.")
    }
    species <- unique(index$species)
    requested <- params(sim)$bird_modelLoad$species
    if (length(requested)) species <- intersect(species, toupper(trimws(requested)))
  }
  if (!is.character(species) || !length(species) || anyNA(species) ||
      any(!grepl("^[A-Za-z][A-Za-z0-9]*$", species)) ||
      anyDuplicated(toupper(species))) {
    stop("bird_summary: species must be explicit unique species codes.")
  }
  config$species <- sort(toupper(species))
  if (p$mode == "single") {
    repID <- p$repID
    if (length(repID) == 1L && is.na(repID)) repID <- predictor$repID
    if (!is.character(repID) || length(repID) != 1L || is.na(repID) ||
        !grepl("^rep[1-9][0-9]*$", repID)) {
      stop("bird_summary: repID must be a canonical label such as rep1.")
    }
    bootstrap <- p$bootstrapID
    if (length(bootstrap) == 1L && is.na(bootstrap)) bootstrap <- predictor$bootstrapIDs
    if (!is.numeric(bootstrap) || length(bootstrap) != 1L || is.na(bootstrap) ||
        !is.finite(bootstrap) || bootstrap < 1 || bootstrap != floor(bootstrap)) {
      stop("bird_summary: exactly one positive whole-number bootstrapID is required.")
    }
    if (!is.null(predictor) &&
        (!identical(as.character(predictor$repID), repID) ||
         !identical(as.numeric(predictor$bootstrapIDs), as.numeric(bootstrap)))) {
      stop("bird_summary: explicit repID/bootstrapID must match the predictor pairing.")
    }
    config$repID <- repID
    config$bootstrapID <- bootstrap
  } else {
    reps <- p$reps
    if (!is.numeric(reps) || !length(reps) || anyNA(reps) ||
        any(!is.finite(reps)) || any(reps < 1 | reps != floor(reps)) || anyDuplicated(reps)) {
      stop("bird_summary: reps must be unique positive whole-number replicate IDs.")
    }
    config$reps <- paste0("rep", sort(reps))
    config$simOutputPath <- p$simOutputPath
  }
  config
}
