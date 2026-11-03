#' Validate selected single-bootstrap mosaics
#' @param predictions Species/year raster list.
#' @param config Validated summary configuration.
#' @param year Simulation year.
#' @return Named list of selected rasters.
#' @keywords internal
birdPairRasters <- function(predictions, config, year) {
  if (!is.list(predictions) || !length(predictions)) {
    stop("bird_summary: mosaicedPredictions must be populated before save_single.")
  }
  rasters <- vector("list", length(config$species))
  names(rasters) <- config$species
  for (species in config$species) {
    raster <- predictions[[species]][[as.character(year)]]
    if (!inherits(raster, "SpatRaster") || terra::nlyr(raster) != 1L ||
        !identical(names(raster), paste0("b", config$bootstrapID, "_", year))) {
      stop("bird_summary: missing or mismatched single-bootstrap mosaic for ", species,
           " year ", year, " bootstrap ", config$bootstrapID, ".")
    }
    rasters[[species]] <- raster
  }
  rasters
}

#' Summarize a density raster on its nonmissing reporting support
#' @param raster One-layer SpatRaster in males/km2.
#' @return One-row data.table of area-weighted mean density, pixel median density,
#'   integrated abundance, valid cell count, area and units.
#' @keywords internal
birdRasterMetrics <- function(raster) {
  if (!inherits(raster, "SpatRaster") || terra::nlyr(raster) != 1L ||
      !nzchar(terra::crs(raster))) {
    stop("bird_summary: metrics require a one-layer raster with a defined CRS.")
  }
  valid <- !is.na(raster)
  count <- terra::global(valid, "sum", na.rm = TRUE)[1, 1]
  if (!is.finite(count) || count == 0) stop("bird_summary: raster has no valid reporting cells.")
  bounds <- terra::global(raster, c("min", "max"), na.rm = TRUE)
  if (any(!is.finite(as.matrix(bounds))) || bounds[1, "min"] < 0) {
    stop("bird_summary: predicted males/km2 must be finite and nonnegative on valid cells.")
  }
  # Model response is males/km2. A future males/ha conversion divides density by
  # 100; integrating with hectares instead must leave total predicted males unchanged.
  area <- terra::cellSize(raster, unit = "km")
  supported_area <- terra::global(terra::ifel(valid, area, NA), "sum", na.rm = TRUE)[1, 1]
  abundance <- terra::global(raster * area, "sum", na.rm = TRUE)[1, 1]
  if (!is.finite(supported_area) || supported_area <= 0 || !is.finite(abundance)) {
    stop("bird_summary: invalid reporting area or integrated abundance.")
  }
  data.table::data.table(
    meanDensity = abundance / supported_area,
    medianDensity = terra::global(raster, stats::median, na.rm = TRUE)[1, 1],
    totalAbundance = abundance, nValidCells = count, validAreaKm2 = supported_area,
    densityUnits = "males/km2", abundanceUnits = "males"
  )
}

#' Calculate annual descriptive mean and median across complete replicate pairs
#' @param pairs Validated pair/year metrics table.
#' @return Long-format annual metrics table, without SD or NRV bounds.
#' @keywords internal
birdAcrossPairMetrics <- function(pairs) {
  long <- data.table::melt(
    pairs, id.vars = c("species", "year", "repID", "bootstrapID", "scenario",
                       "climateScenario", "studyAreaName"),
    measure.vars = c("meanDensity", "medianDensity", "totalAbundance"),
    variable.name = "metric", value.name = "value"
  )
  annual <- long[, .(nPairs = .N, mean = mean(value), median = stats::median(value)),
                 by = .(species, year, scenario, climateScenario, studyAreaName, metric)]
  annual[, units := data.table::fifelse(metric == "totalAbundance", "males", "males/km2")]
  annual
}
