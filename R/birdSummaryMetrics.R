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

#' Calculate sample min-max bounds from exact FireSense reference snapshots
#' @param pairs Validated complete pair/year metrics table.
#' @param config Validated multi-mode selection.
#' @return Per-species/metric bounds, units, context and reference sample counts.
#' @keywords internal
birdNRVBounds <- function(pairs, config) {
  reference_years <- c(2720L, 2820L, 2920L, 3020L)
  selected <- pairs$year %in% reference_years
  reference <- pairs[selected, ]
  expected <- data.table::CJ(species = config$species, year = reference_years, repID = config$reps)
  if (anyDuplicated(reference, by = c("species", "year", "repID")) ||
      nrow(reference) != nrow(expected) ||
      nrow(expected[!reference, on = .(species, year, repID)])) {
    stop("bird_summary: NRV requires one record for every species/reference year/selected replicate.")
  }
  long <- data.table::melt(
    reference, id.vars = c("species", "year", "repID", "scenario", "climateScenario", "studyAreaName"),
    measure.vars = c("meanDensity", "medianDensity", "totalAbundance"),
    variable.name = "metric", value.name = "value"
  )
  if (anyNA(long$value) || any(!is.finite(long$value)) || any(long$value < 0)) {
    stop("bird_summary: NRV reference metrics must be finite, nonmissing and nonnegative.")
  }
  bounds <- long[, .(lower = min(value), upper = max(value), nObservations = .N,
                     nPairs = data.table::uniqueN(repID), nYears = data.table::uniqueN(year)),
                 by = .(species, scenario, climateScenario, studyAreaName, metric)]
  bounds[, `:=`(
    referenceYears = paste(reference_years, collapse = ";"),
    method = "sample_min_max",
    units = data.table::fifelse(metric == "totalAbundance", "males", "males/km2")
  )]
  data.table::setorder(bounds, species, metric)
  bounds
}

#' Plot replicate metric trajectories for each species
#' @param pairs Validated pair/year metrics table.
#' @param destination Existing directory for multi-mode outputs.
#' @param metric One supported metric name.
#' @param bounds Optional per-species/metric table from birdNRVBounds.
#' @return Character vector of saved PNG paths.
#' @keywords internal
birdPlotAcrossReplicateMetrics <- function(pairs, destination, metric, bounds = NULL) {
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
    annual_median <- species_pairs[, .(value = stats::median(get(metric))), by = year]
    species_bounds <- NULL
    if (!is.null(bounds)) {
      selected_bounds <- bounds$species == selected_species & bounds$metric == metric
      species_bounds <- bounds[selected_bounds, ]
      if (nrow(species_bounds) != 1L || anyNA(species_bounds[, .(lower, upper)]) ||
          any(!is.finite(c(species_bounds$lower, species_bounds$upper))) ||
          species_bounds$lower < 0 || species_bounds$lower > species_bounds$upper) {
        stop("bird_summary: missing or invalid NRV plot bounds for ", selected_species, "/", metric)
      }
    }
    file <- file.path(destination, paste0("birdAcrossReplicate", label$stem, "_",
                                         selected_species, ".png"))
    grDevices::png(filename = file, width = 1200, height = 950, res = 120)
    device <- grDevices::dev.cur()
    on.exit(if (grDevices::dev.cur() == device) grDevices::dev.off(), add = TRUE)
    graphics::layout(matrix(1:2, ncol = 1), heights = c(1, 0.20))
    graphics::par(mar = c(5.1, 6.1, 4.1, 2.1))

    x_limits <- range(years)
    if (x_limits[1] == x_limits[2]) x_limits <- x_limits + c(-1, 1)
    max_value <- max(species_pairs[[metric]])
    if (!is.null(species_bounds)) max_value <- max(max_value, species_bounds$upper)
    y_limits <- c(0, if (max_value > 0) max_value * 1.08 else 1)
    graphics::plot.default(NA, xlim = x_limits, ylim = y_limits, xlab = "Simulation year",
                   ylab = label$axis,
                   main = paste(selected_species, "-", label$title),
                   xaxt = "n")
    reference_limits <- c(max(2720, x_limits[1]), min(3020, x_limits[2]))
    show_reference <- reference_limits[1] < reference_limits[2]
    reference_colour <- "#EEE9D5"
    if (show_reference) {
      graphics::rect(reference_limits[1], y_limits[1], reference_limits[2], y_limits[2],
                     col = reference_colour, border = NA)
    }
    ticks <- if (length(years) < 8L) years else pretty(years, n = 7)
    ticks <- ticks[ticks >= min(years) & ticks <= max(years)]
    if (!length(ticks)) ticks <- years
    graphics::axis(1, at = ticks, labels = round(ticks), las = 2)
    graphics::abline(h = pretty(y_limits, n = 5), col = "grey90", lty = "dotted")
    graphics::abline(v = years, col = "grey90", lty = "dotted")
    for (replicate in replicates) {
      replicate_pairs <- species_pairs[species_pairs$repID == replicate, ]
      graphics::lines(replicate_pairs$year, replicate_pairs[[metric]],
                      col = "grey65", lty = 1)
    }
    graphics::lines(annual_median$year, annual_median$value, col = "#1F4E79", lwd = 3)
    if (!is.null(species_bounds)) {
      graphics::abline(h = c(species_bounds$lower, species_bounds$upper),
                       col = "#9A4B25", lty = 2, lwd = 2)
    }
    legend_labels <- c("Individual replicate", "Across-replicate median")
    legend_colours <- c("grey55", "#1F4E79")
    legend_lines <- c(1, 1)
    legend_widths <- c(1, 3)
    legend_points <- c(NA, NA)
    legend_fills <- c(NA, NA)
    if (show_reference) {
      legend_labels <- c(legend_labels, "FireSense NRV reference years")
      legend_colours <- c(legend_colours, "#B5AA7A")
      legend_lines <- c(legend_lines, NA)
      legend_widths <- c(legend_widths, NA)
      legend_points <- c(legend_points, 22)
      legend_fills <- c(legend_fills, reference_colour)
    }
    if (!is.null(species_bounds)) {
      legend_labels <- c(legend_labels, "NRV sample min-max bounds")
      legend_colours <- c(legend_colours, "#9A4B25")
      legend_lines <- c(legend_lines, 2)
      legend_widths <- c(legend_widths, 2)
      legend_points <- c(legend_points, NA)
      legend_fills <- c(legend_fills, NA)
    }
    subtitle <- if (is.null(species_bounds)) "Descriptive replicate comparison; not an NRV range" else
      paste0("NRV sample envelope: ", species_bounds$nPairs, " replicates, ",
             species_bounds$nObservations, " observations; not confidence limits")
    graphics::mtext(subtitle,
                    side = 3, line = 0.2, adj = 1, cex = 0.8, col = "grey35")
    graphics::par(mar = c(0, 6.1, 0, 2.1))
    graphics::plot.new()
    graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1))
    graphics::legend("center", legend = legend_labels,
                     col = legend_colours, lty = legend_lines, lwd = legend_widths,
                     pch = legend_points, pt.bg = legend_fills,
                     ncol = 2, cex = 0.85, bty = "n")
    grDevices::dev.off()
    files[i] <- normalizePath(file, mustWork = TRUE)
  }
  files
}
