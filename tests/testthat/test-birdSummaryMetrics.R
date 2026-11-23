test_that("per-hectare density integrates with actual cell areas in hectares", {
  raster <- terra::rast(nrows = 2, ncols = 1, xmin = 0, xmax = 1,
                        ymin = 0, ymax = 60, crs = "EPSG:4326")
  terra::values(raster) <- c(10, 2)
  area <- as.vector(terra::values(terra::cellSize(raster, unit = "km")))
  out <- summary_env$birdRasterMetrics(raster)
  expect_equal(out$totalAbundance, sum(c(10, 2) * area * 100))
  expect_equal(out$meanDensity, sum(c(10, 2) * area) / sum(area))
  expect_false(isTRUE(all.equal(out$meanDensity, 6)))
  expect_equal(out$medianDensity, 6)
  expect_identical(out$densityUnits, "males/ha")
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
  root <- withr::local_tempdir()
  summaries <- bird_saved_summaries(root)
  expect_length(summary_env$birdPairRasters(summaries, config, 2020)$rasters, 2)
  summaries$GCKI[["2020"]]$meanPath <- NULL
  expect_error(summary_env$birdPairRasters(summaries, config, 2020), "saved meanPath")
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

test_that("only omitted single-mode species are deferred during initialization", {
  p <- summary_parameters()
  p$species <- character(0)
  sim <- list(parameters = list(bird_summary = p))
  expect_identical(summary_env$.inputObjects(sim), sim)
  expect_identical(summary_env$birdSummaryConfig(sim, resolveSpecies = FALSE)$species,
                   character(0))
  expect_error(summary_env$birdSummaryConfig(sim), "modelIndex")
  sim$parameters$bird_summary$species <- "invalid/code"
  expect_error(summary_env$.inputObjects(sim), "species must be explicit")
  sim$parameters$bird_summary <- summary_parameters("multi")
  sim$parameters$bird_summary$species <- character(0)
  expect_error(summary_env$.inputObjects(sim), "species must be explicit")
})

test_that("NRV bounds pool individual reference snapshots, not annual means or other years", {
  pairs <- data.table::CJ(species = c("GCKI", "RCKI"), year = c(2020L, 2720L, 2820L, 2920L, 3020L),
                          repID = paste0("rep", 1:5))
  pairs[, `:=`(scenario = "reference", climateScenario = "forcing", studyAreaName = "fixture")]
  pairs[, meanDensity := as.numeric(sub("rep", "", repID)) + (year - 2720) / 100 +
          ifelse(species == "RCKI", 100, 0)]
  pairs[year == 2020, meanDensity := 10000]
  pairs[, `:=`(medianDensity = meanDensity / 2, totalAbundance = meanDensity * 10)]
  original <- data.table::copy(pairs)
  config <- list(species = c("GCKI", "RCKI"), reps = paste0("rep", 1:5))
  bounds <- summary_env$birdNRVBounds(pairs, config)
  expect_equal(nrow(bounds), 6)
  expect_true(all(bounds$nObservations == 20 & bounds$nPairs == 5 & bounds$nYears == 4))
  expect_true(all(bounds$referenceYears == "2720;2820;2920;3020"))
  expect_true(all(bounds$method == "sample_min_max"))
  density <- bounds[metric == "meanDensity"]
  expect_equal(density$lower, c(1, 101))
  expect_equal(density$upper, c(8, 108))
  expect_true(all(density$units == "males/ha"))
  abundance <- bounds[metric == "totalAbundance"]
  expect_equal(abundance$lower, c(10, 1010))
  expect_equal(abundance$upper, c(80, 1080))
  expect_true(all(abundance$units == "males"))
  expect_equal(bounds[metric == "medianDensity"]$upper, c(4, 54))
  expect_equal(pairs, original)

  config$reps <- c("rep1", "rep2")
  preliminary <- pairs[repID %in% config$reps]
  two <- summary_env$birdNRVBounds(preliminary, config)
  expect_true(all(two$nObservations == 8 & two$nPairs == 2))
  expect_equal(two[metric == "meanDensity"]$upper, c(5, 105))
  expect_error(summary_env$birdNRVBounds(preliminary[-6], config), "every species/reference year")
  expect_error(summary_env$birdNRVBounds(data.table::rbindlist(list(preliminary, preliminary[6])), config),
               "every species/reference year")
})

test_that("NRV calculation is explicit, multi-only and requires all exact reference years", {
  p <- summary_parameters("multi")
  p$calculateNRV <- TRUE
  sim <- list(parameters = list(bird_summary = p))
  expect_error(summary_env$birdSummaryConfig(sim), "requires years")
  p$years <- c(2720L, 2820L, 2920L, 3020L)
  sim$parameters$bird_summary <- p
  expect_true(summary_env$birdSummaryConfig(sim)$calculateNRV)
  sim$parameters$bird_summary$mode <- "single"
  expect_error(summary_env$birdSummaryConfig(sim), "only in multi")
  for (value in list(NA, 1, c(TRUE, FALSE), NULL)) {
    sim$parameters$bird_summary["calculateNRV"] <- list(value)
    expect_error(summary_env$birdSummaryConfig(sim), "must be TRUE or FALSE")
  }
})
