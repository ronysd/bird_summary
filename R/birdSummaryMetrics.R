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
#' @param raster One-layer density SpatRaster.
#' @param densityUnits Raster response units: males/ha or males/km2.
#' @return One-row data.table of area-weighted mean density, pixel median density,
#'   integrated abundance, valid cell count, area and units.
#' @keywords internal
birdRasterMetrics <- function(raster, densityUnits = "males/ha") {
  if (!is.character(densityUnits) || length(densityUnits) != 1L ||
      is.na(densityUnits) || !densityUnits %in% c("males/ha", "males/km2")) {
    stop("bird_summary: densityUnits must be males/ha or males/km2.")
  }
  if (!inherits(raster, "SpatRaster") || terra::nlyr(raster) != 1L ||
      !nzchar(terra::crs(raster))) {
    stop("bird_summary: metrics require a one-layer raster with a defined CRS.")
  }
  if (densityUnits == "males/km2") raster <- raster / 100
  valid <- !is.na(raster)
  count <- terra::global(valid, "sum", na.rm = TRUE)[1, 1]
  if (!is.finite(count) || count == 0) stop("bird_summary: raster has no valid reporting cells.")
  bounds <- terra::global(raster, c("min", "max"), na.rm = TRUE)
  if (any(!is.finite(as.matrix(bounds))) || bounds[1, "min"] < 0) {
    stop("bird_summary: predicted males/ha must be finite and nonnegative on valid cells.")
  }
  # Resolution does not determine density units; normalize before integrating area.
  area <- terra::cellSize(raster, unit = "km")
  supported_area <- terra::global(terra::ifel(valid, area, NA), "sum", na.rm = TRUE)[1, 1]
  abundance <- terra::global(raster * area * 100, "sum", na.rm = TRUE)[1, 1]
  if (!is.finite(supported_area) || supported_area <= 0 || !is.finite(abundance)) {
    stop("bird_summary: invalid reporting area or integrated abundance.")
  }
  data.table::data.table(
    meanDensity = abundance / (supported_area * 100),
    medianDensity = terra::global(raster, stats::median, na.rm = TRUE)[1, 1],
    totalAbundance = abundance, nValidCells = count, validAreaKm2 = supported_area,
    densityUnits = "males/ha", abundanceUnits = "males"
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
  annual[, units := data.table::fifelse(metric == "totalAbundance", "males", "males/ha")]
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
    units = data.table::fifelse(metric == "totalAbundance", "males", "males/ha")
  )]
  data.table::setorder(bounds, species, metric)
  bounds
}

#' Resolve labels shared by trajectory and benchmark comparison plots
#' @param metric Supported metric name.
#' @return File stem, title and axis label.
#' @keywords internal
birdMetricPlotLabel <- function(metric) {
  switch(metric,
    meanDensity = list(stem = "Density", title = "area-weighted mean density",
                       axis = "Area-weighted mean density (males/ha)"),
    totalAbundance = list(stem = "TotalAbundance", title = "total predicted abundance",
                          axis = "Total predicted abundance (males)"),
    medianDensity = list(stem = "MedianDensity", title = "median pixel density",
                         axis = "Median pixel density (males/ha)"),
    stop("bird_summary: unsupported plot metric: ", metric)
  )
}

#' Plot replicate metric trajectories for each species
#' @param pairs Validated pair/year metrics table.
#' @param destination Existing directory for multi-mode outputs.
#' @param metric One supported metric name.
#' @param bounds Optional per-species/metric table from birdNRVBounds.
#' @param benchmark Optional prepared BAM 2020 metrics on the same reporting support.
#' @return Character vector of saved PNG paths.
#' @keywords internal
birdPlotAcrossReplicateMetrics <- function(pairs, destination, metric, bounds = NULL,
                                           benchmark = NULL, currentCondition = "both") {
  label <- birdMetricPlotLabel(metric)
  species_codes <- sort(unique(pairs$species))
  if (!2020L %in% pairs$year && currentCondition %in% c("SIM", "both")) {
    message("bird_summary: 2020 is not selected; plots omit the simulated current-condition line.")
  }
  files <- character(length(species_codes))
  for (i in seq_along(species_codes)) {
    selected_species <- species_codes[i]
    selected_rows <- pairs$species == selected_species
    species_pairs <- pairs[selected_rows, ]
    years <- sort(unique(species_pairs$year))
    replicates <- sort(unique(species_pairs$repID))
    data.table::setorder(species_pairs, repID, year)
    annual_median <- species_pairs[, .(value = stats::median(get(metric))), by = year]
    benchmark_value <- NULL
    if (!is.null(benchmark)) {
      benchmark_value <- benchmark[benchmark$species == selected_species][[metric]]
      if (length(benchmark_value) != 1L || !is.finite(benchmark_value) || benchmark_value < 0) {
        stop("bird_summary: missing or invalid BAM benchmark for ", selected_species, "/", metric)
      }
    }
    current_condition <- NULL
    if (2020L %in% years && currentCondition %in% c("SIM", "both")) {
      baseline <- species_pairs[species_pairs$year == 2020L, ]
      if (anyDuplicated(baseline$repID) ||
          !setequal(baseline$repID, replicates) ||
          any(!is.finite(baseline[[metric]])) || any(baseline[[metric]] < 0)) {
        stop("bird_summary: current condition requires one valid 2020 value per replicate for ",
             selected_species, "/", metric, ".")
      }
      current_condition <- stats::median(baseline[[metric]])
    }
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
    if (!is.null(benchmark_value)) max_value <- max(max_value, benchmark_value)
    if (!is.null(species_bounds)) max_value <- max(max_value, species_bounds$upper)
    y_ticks <- pretty(c(0, if (max_value > 0) max_value * 1.08 else 1), n = 5)
    y_limits <- range(y_ticks)
    graphics::plot.default(NA, xlim = x_limits, ylim = y_limits, xlab = "Simulation year",
                   ylab = label$axis,
                   main = paste(selected_species, "-", label$title),
                   xaxt = "n", yaxt = "n", yaxs = "i")
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
    graphics::axis(2, at = y_ticks, las = 1)
    graphics::abline(h = y_ticks, col = "grey90", lty = "dotted")
    graphics::abline(v = years, col = "grey90", lty = "dotted")
    for (replicate in replicates) {
      replicate_pairs <- species_pairs[species_pairs$repID == replicate, ]
      graphics::lines(replicate_pairs$year, replicate_pairs[[metric]],
                      col = "grey65", lty = 1)
    }
    graphics::lines(annual_median$year, annual_median$value, col = "#1F4E79", lwd = 3)
    if (!is.null(current_condition)) {
      graphics::abline(h = current_condition, col = "red", lty = 2, lwd = 2)
    }
    if (!is.null(benchmark_value)) {
      graphics::abline(h = benchmark_value, col = "#18724B", lty = 4, lwd = 2)
    }
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
    if (!is.null(current_condition)) {
      legend_labels <- c(legend_labels, "Current condition: 2020 median")
      legend_colours <- c(legend_colours, "red")
      legend_lines <- c(legend_lines, 2)
      legend_widths <- c(legend_widths, 2)
      legend_points <- c(legend_points, NA)
      legend_fills <- c(legend_fills, NA)
    }
    if (show_reference) {
      legend_labels <- c(legend_labels, "FireSense NRV reference years")
      legend_colours <- c(legend_colours, "#B5AA7A")
      legend_lines <- c(legend_lines, NA)
      legend_widths <- c(legend_widths, NA)
      legend_points <- c(legend_points, 22)
      legend_fills <- c(legend_fills, reference_colour)
    }
    if (!is.null(benchmark_value)) {
      legend_labels <- c(legend_labels, "BAM 2020: mean; fixed BAM-supported area")
      legend_colours <- c(legend_colours, "#18724B")
      legend_lines <- c(legend_lines, 4)
      legend_widths <- c(legend_widths, 2)
      legend_points <- c(legend_points, NA)
      legend_fills <- c(legend_fills, NA)
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
