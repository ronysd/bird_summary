test_that("reported per-hectare density is independent of grid resolution and raw units", {
  raster <- terra::rast(nrows = 2, ncols = 2, xmin = 0, xmax = 2000,
                        ymin = 0, ymax = 2000, crs = "EPSG:3857")
  terra::values(raster) <- c(1, 2, 3, NA)
  raw <- summary_env$birdRasterMetrics(raster, "males/km2")
  converted <- summary_env$birdRasterMetrics(raster, "males/ha")
  for (field in c("meanDensity", "medianDensity", "totalAbundance")) {
    expect_equal(converted[[field]], raw[[field]] * 100)
  }
  expect_equal(converted$validAreaKm2, raw$validAreaKm2)
  expect_identical(converted$densityUnits, "males/ha")
  expect_equal(as.vector(terra::values(raster)), c(1, 2, 3, NA))
  expect_error(summary_env$birdRasterMetrics(raster, "males/pixel"), "densityUnits")
})

test_that("BAM acquisition downloads missing species once and reuses local products", {
  root <- withr::local_tempdir()
  dir.create(file.path(root, "scratch"))
  summary_env$inputPath <- function(sim) sim$input
  summary_env$scratchPath <- function(sim) sim$scratch
  original <- summary_env$birdDownloadBAMBenchmark
  withr::defer(summary_env$birdDownloadBAMBenchmark <- original)
  calls <- character()
  summary_env$birdDownloadBAMBenchmark <- function(species, region, file, scratchRoot) {
    calls <<- c(calls, paste(species, region))
    raster <- terra::rast(nrows = 2, ncols = 2, crs = "EPSG:3978")
    terra::values(raster) <- 1
    names(raster) <- "mean"
    terra::writeRaster(raster, file)
  }
  sim <- list(input = root, scratch = file.path(root, "scratch"))
  config <- list(species = c("GCKI", "RCKI"), bamRegion = "Canada")
  first <- summary_env$birdPrepareBAMBenchmark(sim, NULL, config)
  expect_equal(length(calls), 2)
  expect_true(all(file.exists(first$rasterPath)))
  expect_equal(summary_env$birdPrepareBAMBenchmark(sim, NULL, config), first)
  expect_equal(length(calls), 2)
  file.remove(first$rasterPath[2])
  summary_env$birdPrepareBAMBenchmark(sim, NULL, config)
  expect_equal(length(calls), 3)
  expect_equal(tail(calls, 1), "RCKI Canada")
  supplied <- data.table::data.table(species = "custom")
  expect_identical(summary_env$birdPrepareBAMBenchmark(sim, supplied, config), supplied)
  expect_equal(length(calls), 3)
})

test_that("BAM downloads are staged and never overwrite an existing product", {
  root <- withr::local_tempdir()
  testthat::local_mocked_bindings(
    bam_get_layer = function(version, spList, destfolder, crop_ext, bcrNM, year) {
      raster <- terra::rast(nrows = 2, ncols = 2, crs = "EPSG:3978")
      terra::values(raster) <- 1
      names(raster) <- "mean"
      file <- file.path(destfolder, paste0(spList, "_", bcrNM, "_", year, ".tif"))
      terra::writeRaster(raster, file)
      stats::setNames(list(terra::rast(file)), spList)
    },
    .package = "BAMexploreR"
  )
  file <- file.path(root, "GCKI_Canada_2020.tif")
  summary_env$birdDownloadBAMBenchmark("GCKI", "Canada", file, root)
  expect_true(file.exists(file))
  expect_identical(list.files(root), basename(file))
  hash <- tools::md5sum(file)
  expect_error(summary_env$birdDownloadBAMBenchmark("GCKI", "Canada", file, root),
               "without overwriting")
  expect_identical(tools::md5sum(file), hash)
})

test_that("current condition and BAM region are validated", {
  p <- summary_parameters("multi")
  p$currentCondition <- "BAM"
  p$years <- c(2720L, 2820L, 2920L, 3020L)
  sim <- list(parameters = list(bird_summary = p))
  expect_identical(summary_env$birdSummaryConfig(sim)$currentCondition, "BAM")
  sim$parameters$bird_summary$currentCondition <- "both"
  expect_error(summary_env$birdSummaryConfig(sim), "requires 2020")
  sim$parameters$bird_summary$currentCondition <- "other"
  expect_error(summary_env$birdSummaryConfig(sim), "currentCondition")
  sim$parameters$bird_summary$currentCondition <- "BAM"
  sim$parameters$bird_summary$bamRegion <- "other"
  expect_error(summary_env$birdSummaryConfig(sim), "bamRegion")
})

bam_fixture <- function(root, years = c(2020L, 2720L, 2820L, 2920L, 3020L)) {
  summaries <- bird_saved_summaries(root, bird_predictions(31, years = years))
  rows <- list()
  for (species in names(summaries)) {
    for (year in years) {
      file <- summaries[[species]][[as.character(year)]]$meanPath
      row <- summary_env$birdRasterMetrics(terra::rast(file))
      row[, `:=`(species = species, year = year, repID = "rep1", bootstrapID = 31,
                  scenario = "reference", climateScenario = "forcing", studyAreaName = "fixture",
                  rasterPath = file, sourceDensityUnits = "males/ha")]
      rows[[length(rows) + 1L]] <- row
    }
  }
  list(pairs = data.table::rbindlist(rows), config = list(
    species = names(summaries), reps = "rep1", predictionDensityUnits = "males/ha", scenario = "reference",
    climateScenario = "forcing", studyAreaName = "fixture"))
}

test_that("current comparison plots support BAM without simulated 2020", {
  root <- withr::local_tempdir()
  fixture <- bam_fixture(root)
  bounds <- summary_env$birdNRVBounds(fixture$pairs, fixture$config)
  benchmark <- fixture$pairs[year == 2020L]
  for (mode in c("BAM", "SIM", "both")) {
    file <- summary_env$birdPlotCurrentBAMNRV(
      fixture$pairs, if (mode == "SIM") NULL else benchmark,
      bounds, root, "meanDensity", mode
    )
    expect_true(all(file.exists(file)))
  }
  without2020 <- fixture$pairs[year != 2020L]
  files <- summary_env$birdPlotCurrentBAMNRV(without2020, benchmark, bounds,
                                             root, "meanDensity", "BAM")
  expect_true(all(file.exists(files)))
  expect_error(summary_env$birdPlotCurrentBAMNRV(without2020, benchmark, bounds,
                                                  root, "meanDensity", "both"),
               "incomplete")
})

test_that("multi-species whiskers show converted medians and sample extrema", {
  root <- withr::local_tempdir()
  fixture <- bam_fixture(root)
  fixture$config$species <- "GCKI"
  pairs <- fixture$pairs[species == "GCKI"]
  bounds <- summary_env$birdNRVBounds(pairs, fixture$config)
  benchmark <- pairs[year == 2020L]
  captured <- list()
  original_segments <- graphics::segments
  original_points <- graphics::points
  testthat::local_mocked_bindings(
    segments = function(x0, y0, x1, y1, ...) {
      if (length(x0) == 1L && length(x1) == 1L && x0 == x1)
        captured[[paste0("range", x0)]] <<- c(y0, y1)
      original_segments(x0, y0, x1, y1, ...)
    },
    points = function(x, y, ...) {
      if (length(x) == 1L && length(y) == 1L)
        captured[[paste0("median", x)]] <<- y
      original_points(x, y, ...)
    },
    .package = "graphics"
  )
  files <- summary_env$birdPlotCurrentNRVWhiskers(pairs, benchmark, bounds, root, "both")
  expect_length(files, 1)
  expect_true(file.exists(files[1]))
  current <- pairs[year == 2020L]$meanDensity * 100
  reference <- pairs[year %in% c(2720L, 2820L, 2920L, 3020L)]$meanDensity * 100
  expect_equal(captured$median1, benchmark$meanDensity * 100)
  expect_equal(captured$median2, median(current))
  expect_equal(captured$median3, median(reference))
  expect_equal(captured$range2, range(current))
  expect_equal(captured$range3, range(reference))
})

test_that("whisker pages support more than seven species", {
  root <- withr::local_tempdir()
  fixture <- bam_fixture(root)
  pairs <- data.table::rbindlist(lapply(1:8, function(i) {
    x <- data.table::copy(fixture$pairs[species == "GCKI"])
    x[, species := paste0("SP", i)]
    x
  }))
  config <- fixture$config
  config$species <- paste0("SP", 1:8)
  bounds <- summary_env$birdNRVBounds(pairs, config)
  benchmark <- pairs[year == 2020L]
  files <- summary_env$birdPlotCurrentNRVWhiskers(pairs, benchmark, bounds, root, "BAM")
  expect_length(files, 2)
  expect_true(all(file.exists(files)))
})

test_that("whisker values remain specific to each species", {
  root <- withr::local_tempdir()
  fixture <- bam_fixture(root)
  pairs <- data.table::copy(fixture$pairs)
  pairs[species == "RCKI", meanDensity := meanDensity + 100]
  bounds <- summary_env$birdNRVBounds(pairs, fixture$config)
  benchmark <- pairs[year == 2020L]
  groups <- summary_env$birdCurrentNRVGroups(pairs, benchmark, bounds, "RCKI",
                                              "meanDensity", "both")
  expect_equal(groups[["BAM 2020"]], pairs[species == "RCKI" & year == 2020L]$meanDensity)
  expect_equal(groups[["NRV"]], pairs[species == "RCKI" & year != 2020L]$meanDensity)
})

test_that("legacy records regenerate without changing source files and are not scaled twice", {
  root <- withr::local_tempdir()
  fixture <- bam_fixture(root)
  pairs <- data.table::copy(fixture$pairs)
  pairs[, `:=`(sourceDensityUnits = NULL, densityUnits = "males/km2",
               totalAbundance = totalAbundance / 100)]
  original <- data.table::copy(pairs)
  hashes <- tools::md5sum(pairs$rasterPath)
  expect_message(
    corrected <- summary_env$birdValidatePairRasters(pairs, "males/ha"),
    "regenerating.*legacy"
  )
  expect_equal(corrected$meanDensity, pairs$meanDensity)
  expect_equal(corrected$totalAbundance, pairs$totalAbundance * 100)
  expect_true(all(corrected$sourceDensityUnits == "males/ha"))
  expect_equal(summary_env$birdValidatePairRasters(corrected, "males/ha"), corrected)
  expect_error(summary_env$birdValidatePairRasters(corrected, "males/km2"), "disagree")
  expect_identical(tools::md5sum(pairs$rasterPath), hashes)
  expect_equal(pairs, original)
  previously_normalized <- data.table::copy(fixture$pairs)
  previously_normalized[, `:=`(densityUnits = "males/km2", meanDensity = meanDensity * 100,
                                medianDensity = medianDensity * 100)]
  migrated <- summary_env$birdValidatePairRasters(previously_normalized, "males/ha")
  expect_equal(migrated$meanDensity, fixture$pairs$meanDensity)
  expect_equal(migrated$totalAbundance, fixture$pairs$totalAbundance)
  expect_true(all(migrated$densityUnits == "males/ha"))
})

test_that("official mean band uses a fixed BAM-supported mask only when enabled", {
  root <- withr::local_tempdir()
  fixture <- bam_fixture(root)
  manifest <- data.table::data.table(
    species = fixture$config$species, year = 2020L,
    rasterPath = file.path(root, paste0(fixture$config$species, "_2020.tif")),
    modelVersion = "v5", densityUnits = "males/ha",
    sourceURL = "https://example.org/verified-fixture", nBootstraps = 32L,
    ensembleStatistic = "mean"
  )
  coarse <- terra::rast(nrows = 2, ncols = 2, xmin = -1000, xmax = 3000,
                        ymin = -1000, ymax = 3000, crs = "EPSG:3857")
  terra::values(coarse) <- 0.1
  product <- c(coarse, coarse * 9, coarse * 99)
  names(product) <- c("mean", "standard_deviation", "detection_distance")
  for (file in manifest$rasterPath) terra::writeRaster(product, file)
  original <- data.table::copy(manifest)
  expect_message(
    comparison <- summary_env$birdBAMComparison(manifest, fixture$pairs, fixture$config),
    "bilinearly aligning"
  )
  benchmark <- comparison$benchmark
  expect_equal(benchmark$meanDensity, rep(0.1, 2), tolerance = 1e-6)
  expect_equal(benchmark$nValidCells, rep(3, 2))
  expect_equal(benchmark$totalAbundance, benchmark$validAreaKm2 * 10, tolerance = 1e-6)
  expect_true(all(benchmark$ensembleStatistic == "mean" & benchmark$nBootstraps == 32))
  expect_equal(manifest, original)
  disabled <- summary_env$birdBAMComparison(NULL, fixture$pairs, fixture$config)
  expect_identical(disabled$pairs, fixture$pairs)
  expect_null(disabled$benchmark)
  expect_length(disabled$masks, 0)
  expect_error(summary_env$birdBAMComparison(manifest[1], fixture$pairs, fixture$config),
               "one BAM benchmark")
  manifest$densityUnits <- "males/km2"
  expect_error(summary_env$birdBAMComparison(manifest, fixture$pairs, fixture$config),
               "males/ha")
  manifest$densityUnits <- "males/ha"
  original_pairs <- data.table::copy(fixture$pairs)
  missing <- terra::rast(fixture$pairs$rasterPath[1])
  terra::values(missing) <- c(0, 0.1, NA, 99)
  names(missing) <- "mean"
  terra::writeRaster(missing, manifest$rasterPath[1], overwrite = TRUE)
  masked <- summary_env$birdBAMComparison(manifest, fixture$pairs, fixture$config)
  expect_equal(masked$benchmark[species == "GCKI"]$nValidCells, 2)
  expect_equal(masked$benchmark[species == "GCKI"]$meanDensity, 0.05, tolerance = 1e-6)
  expect_true(all(masked$pairs[species == "GCKI"]$nValidCells == 2))
  expect_true(all(masked$pairs[species == "RCKI"]$nValidCells == 3))
  expect_equal(as.vector(terra::values(masked$masks$GCKI)), c(1, 1, NA, NA))
  expect_equal(fixture$pairs, original_pairs)
  expect_equal(masked$pairs[species == "GCKI" & year == 2020L]$meanDensity,
               1.5, tolerance = 1e-5)
  expect_lt(masked$pairs[species == "GCKI" & year == 2720L]$meanDensity,
            fixture$pairs[species == "GCKI" & year == 2720L]$meanDensity)
  terra::values(missing) <- NA_real_
  terra::writeRaster(missing, manifest$rasterPath[1], overwrite = TRUE)
  expect_error(summary_env$birdBAMComparison(manifest, fixture$pairs, fixture$config),
               "no BAM-supported reporting cells")
})

test_that("bilinear density interpolation cannot fill the native BAM validity footprint", {
  root <- withr::local_tempdir()
  fixture <- bam_fixture(root)
  fixture$config$species <- "GCKI"
  fixture$pairs <- fixture$pairs[species == "GCKI"]
  native <- terra::shift(terra::rast(fixture$pairs$rasterPath[1]), dx = -100, dy = -100)
  terra::values(native) <- c(0, 0.1, NA, 99)
  file <- file.path(root, "GCKI_2020.tif")
  terra::writeRaster(native, file)
  manifest <- data.table::data.table(
    species = "GCKI", year = 2020L, rasterPath = file, densityUnits = "males/ha",
    modelVersion = "v5", sourceURL = "https://example.org/verified-fixture",
    nBootstraps = 32L, ensembleStatistic = "mean"
  )
  comparison <- summary_env$birdBAMComparison(manifest, fixture$pairs, fixture$config)
  expect_equal(as.vector(terra::values(comparison$masks$GCKI)), c(1, 1, NA, NA))
  path <- fixture$pairs[year == 2720L]$rasterPath
  prediction <- terra::rast(path)
  terra::values(prediction) <- c(NA, 2, 3, NA)
  terra::writeRaster(prediction, path, overwrite = TRUE)
  expect_error(summary_env$birdBAMComparison(manifest, fixture$pairs, fixture$config),
               "simulation missing within fixed BAM reporting mask")
})

test_that("benchmark baselines affect axis limits and comparison plots retain three distinct groups", {
  root <- withr::local_tempdir()
  fixture <- bam_fixture(root)
  pairs <- fixture$pairs
  bounds <- summary_env$birdNRVBounds(pairs, fixture$config)
  benchmark <- pairs[year == 2020L]
  benchmark[, `:=`(meanDensity = 2000, medianDensity = 1500, totalAbundance = 8e6)]
  lines <- numeric()
  limits <- list()
  groups <- list()
  intervals <- numeric()
  original_abline <- graphics::abline
  original_plot <- graphics::plot.default
  original_axis <- graphics::axis
  original_segments <- graphics::segments
  testthat::local_mocked_bindings(
    abline = function(..., h = NULL, col) {
      if (col == "#18724B") lines <<- c(lines, h)
      original_abline(..., h = h, col = col)
    },
    plot.default = function(..., ylim) {
      limits[[length(limits) + 1L]] <<- ylim
      original_plot(..., ylim = ylim)
    },
    axis = function(side, ..., labels = TRUE) {
      if (side == 1) groups[[length(groups) + 1L]] <<- labels
      original_axis(side, ..., labels = labels)
    },
    segments = function(x0, ...) {
      intervals <<- c(intervals, x0)
      original_segments(x0, ...)
    },
    .package = "graphics"
  )
  files <- summary_env$birdPlotAcrossReplicateMetrics(pairs, root, "meanDensity", bounds, benchmark)
  expect_equal(lines, c(2000, 2000))
  for (limit in limits) expect_gt(limit[2], 2000)
  files <- c(files, summary_env$birdPlotCurrentBAMNRV(pairs, benchmark, bounds, root, "meanDensity"))
  expect_true(all(file.info(files)$size > 0))
  expect_identical(tail(groups, 1)[[1]], c("Simulated 2020", "BAM 2020", "NRV"))
  expect_false(any(intervals == 2))
})
