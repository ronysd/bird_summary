#' Read selected saved single-bootstrap mosaics without copying them
#' @param summaries Predictor species/year list with a meanPath per entry.
#' @param config Validated summary configuration.
#' @param year Simulation year.
#' @return Selected rasters and their absolute existing file paths.
#' @keywords internal
birdPairRasters <- function(summaries, config, year) {
  if (!is.list(summaries) || !length(summaries)) {
    stop("bird_summary: summaries with saved meanPath files must be populated before save_single.")
  }
  rasters <- vector("list", length(config$species))
  names(rasters) <- config$species
  paths <- character(length(rasters))
  for (species in config$species) {
    entry <- summaries[[species]][[as.character(year)]]
    if (!is.list(entry)) {
      stop("bird_summary: missing saved summary for ", species, " year ", year, ".")
    }
    file <- entry$meanPath
    rasters[[species]] <- birdReadPairRaster(file, species, year, config$repID, config$bootstrapID)
    paths[match(species, config$species)] <- normalizePath(file, mustWork = TRUE)
  }
  list(rasters = rasters, paths = paths)
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

#' Plot replicate metric trajectories for each species
#' @param pairs Validated pair/year metrics table.
#' @param destination Existing directory for multi-mode outputs.
#' @param metric One supported metric name.
#' @return Character vector of saved PNG paths.
#' @keywords internal
birdPlotAcrossReplicateMetrics <- function(pairs, destination, metric) {
  label <- switch(metric,
    meanDensity = list(stem = "Density", title = "area-weighted mean density",
                       axis = expression("Area-weighted mean density (males/" * km^2 * ")")),
    totalAbundance = list(stem = "TotalAbundance", title = "total predicted abundance",
                          axis = "Total predicted abundance (males)"),
    medianDensity = list(stem = "MedianDensity", title = "median pixel density",
                         axis = expression("Median pixel density (males/" * km^2 * ")")),
    stop("bird_summary: unsupported plot metric: ", metric)
  )
  species_codes <- sort(unique(pairs$species))
  files <- character(length(species_codes))
  for (i in seq_along(species_codes)) {
    selected_species <- species_codes[i]
    selected_rows <- pairs$species == selected_species
    species_pairs <- pairs[selected_rows, ]
    years <- sort(unique(species_pairs$year))
    replicates <- sort(unique(species_pairs$repID))
    data.table::setorder(species_pairs, repID, year)
    annual_mean <- species_pairs[, .(value = mean(get(metric))), by = year]
    file <- file.path(destination, paste0("birdAcrossReplicate", label$stem, "_",
                                         selected_species, ".png"))
    grDevices::png(filename = file, width = 1200, height = 800, res = 120)
    device <- grDevices::dev.cur()
    on.exit(if (grDevices::dev.cur() == device) grDevices::dev.off(), add = TRUE)
    graphics::par(mar = c(5.1, 6.1, 4.1, 2.1))

    x_limits <- range(years)
    if (x_limits[1] == x_limits[2]) x_limits <- x_limits + c(-1, 1)
    max_value <- max(species_pairs[[metric]])
    y_limits <- c(0, if (max_value > 0) max_value * 1.08 else 1)
    graphics::plot.default(NA, xlim = x_limits, ylim = y_limits, xlab = "Simulation year",
                   ylab = label$axis,
                   main = paste(selected_species, "-", label$title),
                   xaxt = "n")
    ticks <- if (length(years) < 8L) years else pretty(years, n = 7)
    ticks <- ticks[ticks >= min(years) & ticks <= max(years)]
    if (!length(ticks)) ticks <- years
    graphics::axis(1, at = ticks, labels = round(ticks), las = 2)
    graphics::abline(h = pretty(y_limits, n = 5), col = "grey90", lty = "dotted")
    for (replicate in replicates) {
      replicate_pairs <- species_pairs[species_pairs$repID == replicate, ]
      graphics::lines(replicate_pairs$year, replicate_pairs[[metric]],
                      col = "grey65", lty = 1)
      graphics::points(replicate_pairs$year, replicate_pairs[[metric]],
                       col = "grey45", pch = 16)
    }
    graphics::lines(annual_mean$year, annual_mean$value, col = "#1F4E79", lwd = 3)
    graphics::points(annual_mean$year, annual_mean$value,
                     col = "#1F4E79", pch = 18, cex = 1.2)
    graphics::legend("topleft", legend = c("Individual replicate", "Across-replicate mean"),
                     col = c("grey55", "#1F4E79"), lty = 1, lwd = c(1, 3),
                     pch = c(16, 18), bty = "n")
    graphics::mtext("Descriptive replicate comparison; not an NRV range",
                    side = 3, line = 0.2, adj = 1, cex = 0.8, col = "grey35")
    grDevices::dev.off()
    files[i] <- normalizePath(file, mustWork = TRUE)
  }
  files
}
