#' Obtain official BAM products for the selected species
#' @param sim Simulation providing configured input and scratch paths.
#' @param manifest Optional caller-supplied benchmark manifest.
#' @param config Validated multi-mode configuration.
#' @return Benchmark manifest with persistent local raster paths.
#' @keywords internal
birdPrepareBAMBenchmark <- function(sim, manifest, config) {
  if (!is.null(manifest)) return(manifest)
  destination <- file.path(inputPath(sim), "bird_summary", "BAM_v5_2020", config$bamRegion)
  if (!dir.exists(destination) &&
      !dir.create(destination, recursive = TRUE, showWarnings = FALSE)) {
    stop("bird_summary: cannot create BAM input directory: ", destination)
  }
  paths <- character(length(config$species))
  for (i in seq_along(config$species)) {
    species <- config$species[i]
    file <- file.path(destination, paste0(species, "_", config$bamRegion, "_2020.tif"))
    if (!file.exists(file)) {
      if (!requireNamespace("BAMexploreR", quietly = TRUE)) {
        stop("bird_summary: BAMexploreR is required to download BAM benchmarks.")
      }
      birdDownloadBAMBenchmark(species, config$bamRegion, file, scratchPath(sim))
    } else {
      message("bird_summary: reusing BAM product: ", file)
    }
    raster <- reproducible::prepInputs(
      targetFile = basename(file), destinationPath = destination,
      fun = terra::rast, useCache = FALSE
    )
    if (!inherits(raster, "SpatRaster") ||
        sum(names(raster) == "mean") != 1L || !nzchar(terra::crs(raster))) {
      stop("bird_summary: invalid BAM product at ", file,
           "; expected a raster with a CRS and one mean band.")
    }
    paths[i] <- normalizePath(file, mustWork = TRUE)
  }
  data.table::data.table(
    species = config$species, year = 2020L, rasterPath = paths,
    densityUnits = "males/ha", modelVersion = "v5",
    sourceURL = "https://doi.org/10.5281/zenodo.21632211",
    nBootstraps = 32L, ensembleStatistic = "mean"
  )
}

#' Download a BAM product without overwriting existing inputs
#' @param species Species code.
#' @param region BAM product region.
#' @param file Destination raster path.
#' @param scratchRoot Configured scratch directory.
#' @return Invisibly, the destination path.
#' @keywords internal
birdDownloadBAMBenchmark <- function(species, region, file, scratchRoot) {
  scratch <- tempfile("bird_bam_", tmpdir = scratchRoot)
  if (!dir.create(scratch, recursive = TRUE, showWarnings = FALSE)) {
    stop("bird_summary: cannot create BAM download directory: ", scratch)
  }
  on.exit(unlink(scratch, recursive = TRUE), add = TRUE)
  message("bird_summary: downloading BAM 2020 ", region, " product for ", species)
  downloaded <- BAMexploreR::bam_get_layer(
    version = "v5", spList = species, destfolder = scratch,
    crop_ext = NULL, bcrNM = region, year = "2020"
  )
  staged <- file.path(scratch, basename(file))
  if (!is.list(downloaded) || !identical(names(downloaded), species) ||
      !inherits(downloaded[[1]], "SpatRaster") || !file.exists(staged) ||
      sum(names(downloaded[[1]]) == "mean") != 1L ||
      !nzchar(terra::crs(downloaded[[1]]))) {
    stop("bird_summary: BAM product unavailable or incomplete for ", species,
         " in ", region, ".")
  }
  if (!file.copy(staged, file, overwrite = FALSE)) {
    stop("bird_summary: could not save BAM product without overwriting: ", file)
  }
  invisible(file)
}

#' Recalculate BAM and all simulation records on fixed BAM-supported reporting areas
#' @param manifest Optional data.table of local official BAM raster paths and provenance.
#' @param pairs Validated simulation pair metrics, including existing raster paths.
#' @param config Validated multi-mode configuration.
#' @return List of pair metrics, optional benchmark metrics and per-species masks.
#' @keywords internal
birdBAMComparison <- function(manifest, pairs, config) {
  if (is.null(manifest)) return(list(pairs = pairs, benchmark = NULL, masks = list()))
  fields <- c("species", "year", "rasterPath", "densityUnits", "modelVersion",
              "sourceURL", "nBootstraps", "ensembleStatistic")
  if (!is.data.frame(manifest) || !all(fields %in% names(manifest)) || !nrow(manifest)) {
    stop("bird_summary: bamBenchmarkFiles must be a nonempty manifest with: ",
         paste(fields, collapse = ", "), ".")
  }
  manifest <- data.table::as.data.table(data.table::copy(manifest))
  for (field in setdiff(fields, c("year", "nBootstraps"))) {
    value <- manifest[[field]]
    if (!is.character(value) || anyNA(value) || any(!nzchar(value))) {
      stop("bird_summary: benchmark ", field, " must contain nonempty character values.")
    }
  }
  manifest <- manifest[manifest$species %in% config$species]
  if (anyDuplicated(manifest$species) || !setequal(manifest$species, config$species)) {
    stop("bird_summary: require exactly one BAM benchmark per selected species.")
  }
  if (!is.numeric(manifest$year) || anyNA(manifest$year) || any(manifest$year != 2020) ||
      any(manifest$densityUnits != "males/ha" | manifest$modelVersion != "v5" |
          manifest$ensembleStatistic != "mean") ||
      !is.numeric(manifest$nBootstraps) || anyNA(manifest$nBootstraps) ||
      any(!is.finite(manifest$nBootstraps)) ||
      any(manifest$nBootstraps < 1 | manifest$nBootstraps != floor(manifest$nBootstraps))) {
    stop("bird_summary: benchmark must identify V5 2020 bootstrap-mean surfaces in males/ha ",
         "with a positive whole-number nBootstraps.")
  }
  rows <- vector("list", length(config$species))
  masks <- stats::setNames(vector("list", length(config$species)), config$species)
  updated <- data.table::copy(pairs)
  if ("reportingMaskPath" %in% names(updated)) updated[, reportingMaskPath := NULL]
  for (i in seq_along(config$species)) {
    selected_species <- config$species[i]
    record <- manifest[manifest$species == selected_species]
    if (!file.exists(record$rasterPath) || dir.exists(record$rasterPath)) {
      stop("bird_summary: BAM benchmark raster is missing: ", record$rasterPath)
    }
    raster <- terra::rast(record$rasterPath)
    if (sum(names(raster) == "mean") != 1L || !nzchar(terra::crs(raster))) {
      stop("bird_summary: BAM benchmark requires a defined CRS and one named mean band: ",
           record$rasterPath)
    }
    raster <- raster[["mean"]]
    valid <- !is.na(raster)
    pair <- pairs[pairs$species == selected_species][1]
    template <- birdReadPairRaster(pair$rasterPath, selected_species, pair$year, pair$repID, pair$bootstrapID)
    alignment <- "native_grid"
    if (!terra::compareGeom(raster, template, stopOnError = FALSE)) {
      message("bird_summary: bilinearly aligning BAM 2020 density for ", selected_species,
              " to the simulation grid; preserving its native valid-data footprint with nearest-neighbour alignment.")
      if (terra::same.crs(raster, template)) {
        valid <- terra::resample(valid, template, method = "near")
        raster <- terra::resample(raster, template, method = "bilinear")
      } else {
        valid <- terra::project(valid, template, method = "near")
        raster <- terra::project(raster, template, method = "bilinear")
      }
      alignment <- "bilinear_to_simulation_grid"
    }
    reporting_mask <- terra::ifel(!is.na(template) & !is.na(valid) & valid == 1, 1, NA)
    retained <- terra::global(!is.na(reporting_mask), "sum", na.rm = TRUE)[1, 1]
    if (!is.finite(retained) || retained == 0) {
      stop("bird_summary: no BAM-supported reporting cells for ", selected_species, ".")
    }
    original <- terra::global(!is.na(template), "sum", na.rm = TRUE)[1, 1]
    message("bird_summary: fixed BAM-supported mask for ", selected_species, " retains ", retained,
            " of ", original, " reporting cells; the same mask applies to every selected year and replicate. Valid zeros are retained.")
    missing <- terra::global(!is.na(reporting_mask) & is.na(raster), "sum", na.rm = TRUE)[1, 1]
    if (missing > 0) {
      stop("bird_summary: aligned BAM density is missing within its fixed valid-data mask for ",
           selected_species, "; reporting support was not reduced further.")
    }
    masks[[selected_species]] <- reporting_mask
    selected_rows <- which(updated$species == selected_species)
    for (j in selected_rows) {
      row <- updated[j]
      prediction <- birdReadPairRaster(row$rasterPath, selected_species, row$year, row$repID, row$bootstrapID)
      if (!terra::compareGeom(template, prediction, stopOnError = FALSE)) {
        stop("bird_summary: incompatible raster geometry for ", selected_species)
      }
      missing <- terra::global(!is.na(reporting_mask) & is.na(prediction), "sum", na.rm = TRUE)[1, 1]
      if (missing > 0) {
        stop("bird_summary: simulation missing within fixed BAM reporting mask: ", row$rasterPath)
      }
      metrics <- birdRasterMetrics(terra::mask(prediction, reporting_mask), config$predictionDensityUnits)
      for (field in names(metrics)) data.table::set(updated, j, field, metrics[[field]])
    }
    raster <- terra::mask(raster, reporting_mask)
    rows[[i]] <- birdRasterMetrics(raster, "males/ha")
    rows[[i]][, `:=`(
      species = selected_species, year = 2020L, modelVersion = record$modelVersion,
      sourceURL = record$sourceURL, rasterPath = normalizePath(record$rasterPath, mustWork = TRUE),
      sourceDensityUnits = "males/ha", nBootstraps = record$nBootstraps,
      ensembleStatistic = "mean", alignment = alignment,
      reportingSupport = "bam_2020_valid", originalValidCells = original,
      excludedCells = original - retained, maskAlignment = "nearest_native_validity",
      scenario = config$scenario, climateScenario = config$climateScenario,
      studyAreaName = config$studyAreaName,
      processing = "official BAM public product; truncation and range/data/water masks"
    )]
  }
  updated[, reportingSupport := "bam_2020_valid"]
  list(pairs = updated, benchmark = data.table::rbindlist(rows), masks = masks)
}

#' Resolve and validate a species' current and reference values
#' @param pairs Validated pair metrics.
#' @param benchmark Optional BAM benchmark metrics.
#' @param bounds Per-species/metric NRV sample bounds.
#' @param species Species code.
#' @param metric Supported metric.
#' @param currentCondition BAM, SIM, or both.
#' @return Named groups of numeric values.
#' @keywords internal
birdCurrentNRVGroups <- function(pairs, benchmark, bounds, species, metric, currentCondition) {
  current_rows <- pairs$species == species & pairs$year == 2020L
  reference_rows <- pairs$species == species &
    pairs$year %in% c(2720L, 2820L, 2920L, 3020L)
  current <- pairs[current_rows][[metric]]
  reference <- pairs[reference_rows][[metric]]
  official <- if (currentCondition %in% c("BAM", "both")) {
    benchmark_rows <- benchmark$species == species
    benchmark[benchmark_rows][[metric]]
  } else numeric()
  bound_rows <- bounds$species == species & as.character(bounds$metric) == metric
  selected_bounds <- bounds[bound_rows]
  if ((currentCondition %in% c("SIM", "both") && !length(current)) ||
      (currentCondition %in% c("BAM", "both") && length(official) != 1L) ||
      !length(reference) || nrow(selected_bounds) != 1L ||
      any(!is.finite(c(current, reference, official))) ||
      !isTRUE(all.equal(range(reference), c(selected_bounds$lower, selected_bounds$upper)))) {
    stop("bird_summary: incomplete or inconsistent current/BAM/NRV plot data for ", species, ".")
  }
  groups <- list()
  if (currentCondition %in% c("SIM", "both")) groups[["Simulated 2020"]] <- current
  if (currentCondition %in% c("BAM", "both")) groups[["BAM 2020"]] <- official
  groups[["NRV"]] <- reference
  groups
}

#' Plot selected current condition and the NRV reference sample
#' @param pairs Validated simulation pair metrics.
#' @param benchmark Optional official benchmark metrics on the same reporting support.
#' @param bounds Validated per-species/metric NRV sample bounds.
#' @param destination Existing output directory.
#' @param metric Supported density or abundance metric.
#' @param currentCondition BAM, SIM, or both.
#' @return Saved per-species PNG paths.
#' @keywords internal
birdPlotCurrentBAMNRV <- function(pairs, benchmark, bounds, destination, metric,
                                  currentCondition = "both") {
  label <- birdMetricPlotLabel(metric)
  files <- character(length(unique(pairs$species)))
  species_codes <- sort(unique(pairs$species))
  for (i in seq_along(species_codes)) {
    selected_species <- species_codes[i]
    groups <- birdCurrentNRVGroups(pairs, benchmark, bounds, selected_species,
                                   metric, currentCondition)
    fileMode <- switch(currentCondition, BAM = "BAMOnly", SIM = "SIM", both = "BAM")
    file <- file.path(destination, paste0("birdCurrent", fileMode, "NRV",
                                         label$stem, "_", selected_species, ".png"))
    grDevices::png(file, width = 1200, height = 950, res = 120)
    device <- grDevices::dev.cur()
    on.exit(if (grDevices::dev.cur() == device) grDevices::dev.off(), add = TRUE)
    graphics::layout(matrix(1:2, ncol = 1), heights = c(1, 0.25))
    graphics::par(mar = c(4, 6.1, 4.1, 2.1))
    colours <- c("Simulated 2020" = "red", "BAM 2020" = "#18724B",
                 "NRV" = "#1F4E79")[names(groups)]
    ticks <- pretty(c(0, max(c(unlist(groups), 1e-10)) * 1.08), n = 5)
    graphics::plot.default(NA, xlim = c(0.5, length(groups) + 0.5), ylim = range(ticks), xaxt = "n", yaxt = "n",
                           yaxs = "i", xlab = "", ylab = label$axis,
                           main = paste(selected_species, "-", label$title))
    graphics::axis(1, at = seq_along(groups), labels = names(groups))
    graphics::axis(2, at = ticks, las = 1)
    graphics::abline(h = ticks, col = "grey90", lty = "dotted")
    for (j in seq_along(groups)) {
      group <- groups[[j]]
      if (names(groups)[j] != "BAM 2020") {
        limits <- if (names(groups)[j] == "NRV")
          stats::quantile(group, c(0.025, 0.975), names = FALSE) else range(group)
        graphics::segments(j, limits[1], j, limits[2], col = colours[j], lwd = 2)
        graphics::segments(j - 0.07, limits, j + 0.07, limits, col = colours[j], lwd = 2)
        graphics::points(j + seq(-0.12, 0.12, length.out = length(group)), group,
                         col = grDevices::adjustcolor(colours[j], alpha.f = 0.5), pch = 16, cex = 0.8)
      }
      graphics::points(j, stats::median(group), pch = 23, bg = colours[j], cex = 1.6)
    }
    graphics::mtext("SIM: replicate range; NRV: 2.5-97.5 percentiles; not confidence intervals", side = 3,
                    line = 0.2, cex = 0.85)
    graphics::par(mar = c(0, 0, 0, 0))
    graphics::plot.new()
    graphics::legend("center", bty = "n", cex = 0.9,
      legend = c(if (currentCondition %in% c("SIM", "both"))
                   "Simulated: median and range across 2020 replicates",
                 if (currentCondition %in% c("BAM", "both"))
                   "BAM: metric of official bootstrap-mean surface; no landscape interval",
                 "NRV: median and percentile interval; snapshots are correlated",
                 if (currentCondition %in% c("BAM", "both"))
                   "All groups use the same fixed BAM-supported reporting area"))
    grDevices::dev.off()
    files[i] <- normalizePath(file, mustWork = TRUE)
  }
  files
}

#' Plot current and NRV density medians and sample ranges across species
#' @param pairs Validated pair metrics in males/ha.
#' @param benchmark Optional BAM benchmark metrics on the same reporting support.
#' @param bounds Per-species NRV sample bounds.
#' @param destination Existing output directory.
#' @param currentCondition BAM, SIM, or both.
#' @return Paths of saved multi-species PNG pages.
#' @keywords internal
birdPlotCurrentNRVWhiskers <- function(pairs, benchmark, bounds, destination,
                                        currentCondition) {
  species_codes <- sort(unique(pairs$species))
  pages <- split(species_codes, ceiling(seq_along(species_codes) / 7))
  files <- character(length(pages))
  colours <- c("BAM 2020" = "#18724B", "Simulated 2020" = "#C23B3B", NRV = "#1F4E79")
  for (page in seq_along(pages)) {
    selected <- pages[[page]]
    groups <- lapply(selected, function(species) {
      values <- birdCurrentNRVGroups(pairs, benchmark, bounds, species,
                                     "meanDensity", currentCondition)
      values <- values[intersect(c("BAM 2020", "Simulated 2020", "NRV"), names(values))]
      lapply(values, function(x) x * 100)
    })
    file <- file.path(destination, sprintf("birdCurrentNRVWhiskers_%s_page%02d.png",
                                           currentCondition, page))
    grDevices::png(file, width = 2600, height = 1500, res = 160)
    device <- grDevices::dev.cur()
    on.exit(if (grDevices::dev.cur() == device) grDevices::dev.off(), add = TRUE)
    graphics::layout(matrix(1:8, nrow = 2, byrow = TRUE))
    graphics::par(mar = c(3.4, 5.2, 3.2, 1), oma = c(4, 0, 3, 0))
    for (i in seq_along(selected)) {
      values <- groups[[i]]
      ticks <- pretty(c(0, max(c(unlist(values), 1e-10)) * 1.15), n = 5)
      graphics::plot.default(NA, xlim = c(0.5, length(values) + 0.5), ylim = range(ticks),
                             xaxt = "n", xlab = "", ylab = "Mean density (males/km2)",
                             main = selected[i])
      graphics::abline(h = ticks, col = "grey90", lty = "dotted")
      graphics::axis(1, at = seq_along(values), labels = names(values))
      for (j in seq_along(values)) {
        group <- values[[j]]
        colour <- colours[[names(values)[j]]]
        if (names(values)[j] != "BAM 2020") {
          limits <- range(group)
          graphics::segments(j, limits[1], j, limits[2], col = colour, lwd = 2)
          graphics::segments(j - 0.06, limits, j + 0.06, limits, col = colour, lwd = 2)
          graphics::points(j + seq(-0.11, 0.11, length.out = length(group)), group,
                           col = grDevices::adjustcolor(colour, alpha.f = 0.38),
                           pch = 16, cex = 0.65)
        }
        graphics::points(j, stats::median(group), pch = if (names(values)[j] == "BAM 2020") 23 else 21,
                         bg = colour, col = colour, cex = 1.7)
      }
    }
    for (i in seq_len(7L - length(selected))) graphics::plot.new()
    graphics::plot.new()
    legend_groups <- intersect(c("BAM 2020", "Simulated 2020", "NRV"),
                               names(groups[[1]]))
    graphics::legend("center", bty = "n", cex = 1,
                     legend = c(if (currentCondition %in% c("BAM", "both"))
                                  "BAM: official 2020 mean product",
                                if (currentCondition %in% c("SIM", "both"))
                                  "SIM: median and min-max of 2020 replicates",
                                "NRV: median and min-max of reference values",
                                "Small points: individual modelled values"),
                     col = c(unname(colours[legend_groups]), "grey50"),
                     pt.bg = c(unname(colours[legend_groups]), "grey50"),
                     pch = c(if (currentCondition %in% c("BAM", "both")) 23,
                             if (currentCondition %in% c("SIM", "both")) 21, 21, 16))
    graphics::mtext(paste0(pairs$studyAreaName[1], " (", pairs$scenario[1],
                           "): current condition and sampled NRV"),
                    outer = TRUE, side = 3, line = 1, font = 2)
    graphics::mtext(paste0(if (currentCondition %in% c("BAM", "both"))
                             "Same BAM-supported cells within each species. ",
                           "NRV years: 2720, 2820, 2920, 3020; all selected replicates."),
                    outer = TRUE, side = 1, line = 1, cex = 0.85)
    graphics::mtext("Whiskers are sample min-max, not confidence intervals. BAM is one modelled estimate.",
                    outer = TRUE, side = 1, line = 2.4, cex = 0.85)
    grDevices::dev.off()
    files[page] <- normalizePath(file, mustWork = TRUE)
  }
  files
}
