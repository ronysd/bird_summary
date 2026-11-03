test_that("area weighting and abundance use actual km2 cell areas", {
  raster <- terra::rast(nrows = 2, ncols = 1, xmin = 0, xmax = 1,
                        ymin = 0, ymax = 60, crs = "EPSG:4326")
  terra::values(raster) <- c(10, 2)
  area <- as.vector(terra::values(terra::cellSize(raster, unit = "km")))
  out <- summary_env$birdRasterMetrics(raster)
  expect_equal(out$totalAbundance, sum(c(10, 2) * area))
  expect_equal(out$meanDensity, sum(c(10, 2) * area) / sum(area))
  expect_false(isTRUE(all.equal(out$meanDensity, 6)))
  expect_equal(out$medianDensity, 6)
  expect_identical(out$densityUnits, "males/km2")
  expect_identical(out$abundanceUnits, "males")
  expect_equal(out$nValidCells, 2)
})

test_that("invalid and all-missing response rasters fail clearly", {
  raster <- bird_predictions()$GCKI[["2020"]]
  out <- summary_env$birdRasterMetrics(raster)
  expect_equal(out$nValidCells, 3)
  for (values in list(rep(NA_real_, 4), c(-1, 1, 2, NA), c(Inf, 1, 2, NA))) {
    terra::values(raster) <- values
    expect_error(summary_env$birdRasterMetrics(raster), "no valid|finite and nonnegative")
  }
})

test_that("single pairing and multi expected selection are validated", {
  p <- summary_parameters()
  sim <- list(parameters = list(bird_summary = p))
  config <- summary_env$birdSummaryConfig(sim)
  expect_identical(config$years, 2020L:2021L)
  expect_identical(config$repID, "rep1")
  predictions <- bird_predictions()
  expect_length(summary_env$birdPairRasters(predictions, config, 2020), 2)
  predictions$GCKI[["2020"]] <- c(predictions$GCKI[["2020"]], predictions$GCKI[["2020"]])
  expect_error(summary_env$birdPairRasters(predictions, config, 2020), "single-bootstrap")
  sim$parameters$bird_summary$bootstrapID <- c(1, 2)
  expect_error(summary_env$birdSummaryConfig(sim), "exactly one")
  sim$parameters$bird_summary <- summary_parameters("multi")
  sim$parameters$bird_summary$years <- NA_integer_
  expect_error(summary_env$birdSummaryConfig(sim), "years must be explicit")
  sim$parameters$bird_summary$years <- 2020L:2021L
  sim$parameters$bird_summary$species <- character(0)
  expect_error(summary_env$birdSummaryConfig(sim), "species must be explicit")
})

test_that("single defaults resolve the predictor calendar and declared loader species", {
  p <- summary_parameters()
  p$years <- NA_integer_
  p$species <- character(0)
  p$repID <- NA_character_
  p$bootstrapID <- NA_real_
  sim <- list(parameters = list(
    bird_summary = p,
    bird_modelPredict = list(repID = "rep1", bootstrapIDs = 31,
                             predictStartYear = 2020, predictEndYear = 2021, predictInterval = 1),
    bird_modelLoad = list(species = "GCKI")
  ), modelIndex = data.table::data.table(species = c("GCKI", "RCKI")))
  config <- summary_env$birdSummaryConfig(sim)
  expect_identical(config$species, "GCKI")
  expect_equal(config$bootstrapID, 31)
  expect_identical(config$years, 2020L:2021L)
  sim$parameters$bird_summary$bootstrapID <- 7
  expect_error(summary_env$birdSummaryConfig(sim), "must match the predictor")
})
