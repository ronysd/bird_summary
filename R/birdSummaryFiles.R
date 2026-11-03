#' Construct a scenario-separated summary destination
#' @param sim Simulation with configured outputPath.
#' @param config Validated summary configuration.
#' @param ... Replicate/year or multi subdirectories.
#' @return Existing absolute directory path.
#' @keywords internal
birdSummaryDirectory <- function(sim, config, ...) {
  root <- outputPath(sim)
  if (!is.character(root) || length(root) != 1L || is.na(root) || !nzchar(root)) {
    stop("bird_summary: outputPath must be one nonempty directory.")
  }
  path <- file.path(root, "bird_summary", config$scenario, config$climateScenario,
                    config$studyAreaName, ...)
  if (!dir.exists(path) && !dir.create(path, recursive = TRUE, showWarnings = FALSE)) {
    stop("bird_summary: could not create directory: ", path)
  }
  normalizePath(path, mustWork = TRUE)
}

#' Discover registered pair tables, or scan disk when no registry is provided
#' @param sim Simulation with optional combined outputsDF.
#' @param config Validated multi-mode configuration.
#' @return Unique existing pair-table paths.
#' @keywords internal
birdSummaryFiles <- function(sim, config) {
  registry <- sim$outputsDF
  if (!is.null(registry) && NROW(registry) > 0L) {
    if (!is.data.frame(registry) || !"file" %in% names(registry) ||
        !is.character(registry$file) || anyNA(registry$file)) {
      stop("bird_summary: outputsDF must have a nonmissing character file column.")
    }
    files <- registry$file
  } else {
    root <- config$simOutputPath
    if (!is.character(root) || length(root) != 1L || is.na(root) || !dir.exists(root)) {
      stop("bird_summary: simOutputPath must be an existing directory without outputsDF.")
    }
    files <- list.files(root, recursive = TRUE, full.names = TRUE)
  }
  files <- unique(files[grepl("^birdPairSummary_year[0-9]{4}_rep[1-9][0-9]*_boot-[0-9]+\\.csv$",
                              basename(files))])
  if (!length(files)) stop("bird_summary: no registered or on-disk pair summary tables found.")
  missing <- files[!file.exists(files)]
  if (length(missing)) stop("bird_summary: registered pair files are missing: ", paste(missing, collapse = ", "))
  normalizePath(files, mustWork = TRUE)
}

#' Read, filter and validate a complete species/year/replicate selection
#' @param files Pair-table paths.
#' @param config Validated multi-mode configuration.
#' @return Complete selected pair records, including their source file.
#' @keywords internal
birdReadPairTables <- function(files, config) {
  tables <- vector("list", length(files))
  required <- c("species", "year", "repID", "bootstrapID", "scenario", "climateScenario",
                "studyAreaName", "meanDensity", "medianDensity", "totalAbundance",
                "nValidCells", "validAreaKm2", "densityUnits", "abundanceUnits", "rasterPath")
  for (i in seq_along(files)) {
    table <- data.table::fread(files[i])
    if (!all(required %in% names(table)) || !nrow(table)) {
      stop("bird_summary: invalid or empty pair table: ", files[i])
    }
    table[, sourceFile := files[i]]
    tables[[i]] <- table
  }
  pairs <- data.table::rbindlist(tables, use.names = TRUE)
  selected <- pairs$scenario == config$scenario &
    pairs$climateScenario == config$climateScenario &
    pairs$studyAreaName == config$studyAreaName &
    pairs$species %in% config$species & pairs$year %in% config$years & pairs$repID %in% config$reps
  if (anyNA(selected)) stop("bird_summary: pair-table selection metadata contains missing values.")
  pairs <- pairs[selected, ]
  if (anyNA(pairs[, ..required])) {
    stop("bird_summary: selected pair records contain missing values.")
  }
  if (anyDuplicated(pairs, by = c("species", "year", "repID"))) {
    stop("bird_summary: duplicate species/year/replicate pair records.")
  }
  expected <- data.table::CJ(species = config$species, year = config$years, repID = config$reps)
  missing <- expected[!pairs, on = .(species, year, repID)]
  if (nrow(missing)) {
    stop("bird_summary: missing expected species/year/replicate records: ",
         paste(paste(missing$species, missing$year, missing$repID, sep = "/"), collapse = ", "))
  }
  numeric_fields <- c("bootstrapID", "meanDensity", "medianDensity", "totalAbundance",
                      "nValidCells", "validAreaKm2")
  for (field in numeric_fields) {
    if (!is.numeric(pairs[[field]]) || anyNA(pairs[[field]]) ||
        any(!is.finite(pairs[[field]])) || any(pairs[[field]] < 0)) {
      stop("bird_summary: invalid numeric pair field: ", field)
    }
    data.table::set(pairs, j = field, value = as.numeric(pairs[[field]]))
  }
  if (any(pairs$bootstrapID < 1 | pairs$bootstrapID != floor(pairs$bootstrapID)) ||
      any(pairs$nValidCells < 1 | pairs$nValidCells != floor(pairs$nValidCells)) ||
      any(pairs$validAreaKm2 <= 0) ||
      any(pairs$densityUnits != "males/km2" | pairs$abundanceUnits != "males")) {
    stop("bird_summary: invalid pair identities, area, cell count or response units.")
  }
  mapping <- unique(pairs[, .(repID, bootstrapID)])
  if (anyDuplicated(mapping$repID)) stop("bird_summary: each replicate must have exactly one bootstrapID.")
  data.table::setorder(pairs, species, year, repID)
  pairs
}

#' Verify raster identity, geometry, reporting support and persisted metric values
#' @param pairs Validated pair metrics table.
#' @return Invisible TRUE.
#' @keywords internal
birdValidatePairRasters <- function(pairs) {
  templates <- list()
  for (i in seq_len(nrow(pairs))) {
    row <- pairs[i]
    file <- row$rasterPath
    if (!is.character(file) || length(file) != 1L || is.na(file) || !file.exists(file)) {
      stop("bird_summary: missing pair raster: ", file)
    }
    expected <- paste0(row$year, "_", row$species, "_mosaic_", row$repID,
                        "_boot-", row$bootstrapID, ".tif")
    if (!identical(basename(file), expected)) stop("bird_summary: raster filename conflicts with pair identity: ", file)
    raster <- terra::rast(file)
    if (terra::nlyr(raster) != 1L || !identical(names(raster), paste0("b", row$bootstrapID, "_", row$year))) {
      stop("bird_summary: raster layer conflicts with pair identity: ", file)
    }
    template <- templates[[row$species]]
    if (!is.null(template)) {
      if (!terra::compareGeom(template, raster, stopOnError = FALSE)) {
        stop("bird_summary: incompatible raster geometry for ", row$species)
      }
      mask_difference <- terra::global(is.na(template) != is.na(raster), "sum", na.rm = TRUE)[1, 1]
      if (mask_difference != 0) stop("bird_summary: reporting support differs for ", row$species)
    } else {
      templates[[row$species]] <- raster
    }
    metrics <- birdRasterMetrics(raster)
    for (field in c("meanDensity", "medianDensity", "totalAbundance", "nValidCells", "validAreaKm2")) {
      if (!isTRUE(all.equal(metrics[[field]], row[[field]], tolerance = 1e-6, check.attributes = FALSE))) {
        stop("bird_summary: persisted metrics disagree with raster for ", field, ": ", file)
      }
    }
  }
  invisible(TRUE)
}
