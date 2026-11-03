test_that("file registry is authoritative and its caller-owned table is unchanged", {
  root <- withr::local_tempdir()
  config <- list(simOutputPath = root)
  file <- file.path(root, "birdPairSummary_year2020_rep1_boot-31.csv")
  writeLines("fixture", file)
  registry <- data.table::data.table(file = c(file, file))
  original <- data.table::copy(registry)
  expect_identical(summary_env$birdSummaryFiles(list(outputsDF = registry), config),
                   normalizePath(file))
  expect_equal(registry, original)
  unrelated <- data.table::data.table(file = file.path(root, "unrelated.csv"))
  expect_error(summary_env$birdSummaryFiles(list(outputsDF = unrelated), config),
               "no registered or on-disk pair")
  expect_error(summary_env$birdSummaryFiles(list(outputsDF = data.table::data.table(x = "bad")), config),
               "outputsDF must")
})

test_that("spatial support and grids must be comparable across pair-year rasters", {
  root <- withr::local_tempdir()
  predictions <- bird_predictions()
  rows <- list()
  for (year in 2020:2021) {
    raster <- predictions$GCKI[[as.character(year)]]
    file <- file.path(root, paste0(year, "_GCKI_mosaic_rep1_boot-31.tif"))
    terra::writeRaster(raster, file, overwrite = TRUE)
    row <- summary_env$birdRasterMetrics(terra::rast(file))
    row[, `:=`(species = "GCKI", year = year, repID = "rep1", bootstrapID = 31,
               rasterPath = file)]
    rows[[length(rows) + 1L]] <- row
  }
  pairs <- data.table::rbindlist(rows)
  expect_true(summary_env$birdValidatePairRasters(pairs))
  raster <- terra::rast(pairs$rasterPath[2])
  terra::values(raster) <- c(NA, 3, 4, 5)
  terra::writeRaster(raster, pairs$rasterPath[2], overwrite = TRUE)
  expect_error(summary_env$birdValidatePairRasters(pairs), "reporting support differs")
  raster <- terra::shift(predictions$GCKI[["2021"]], dx = 1000)
  terra::writeRaster(raster, pairs$rasterPath[2], overwrite = TRUE)
  expect_error(summary_env$birdValidatePairRasters(pairs), "incompatible raster geometry")
  pairs$rasterPath[2] <- file.path(root, "missing.tif")
  expect_error(summary_env$birdValidatePairRasters(pairs), "missing pair raster")
})

test_that("multi mode rejects inconsistent bootstrap mapping and corrupted metrics", {
  root <- withr::local_tempdir()
  config <- summary_env$birdSummaryConfig(list(
    parameters = list(bird_summary = summary_parameters("multi"))
  ))
  keys <- data.table::CJ(species = config$species, year = config$years, repID = config$reps)
  keys[, `:=`(bootstrapID = 31, scenario = "reference", climateScenario = "CanESM5_SSP370",
              studyAreaName = "fixture", meanDensity = 2, medianDensity = 2,
              totalAbundance = 6, nValidCells = 3, validAreaKm2 = 3,
              densityUnits = "males/km2", abundanceUnits = "males", rasterPath = "unused")]
  file <- file.path(root, "pairs.csv")
  data.table::fwrite(keys, file)
  expect_equal(nrow(summary_env$birdReadPairTables(file, config)), 8)
  keys$bootstrapID[1] <- 32
  data.table::fwrite(keys, file)
  expect_error(summary_env$birdReadPairTables(file, config), "exactly one bootstrapID")
  keys$bootstrapID[1] <- 31
  keys$densityUnits[1] <- "males/ha"
  data.table::fwrite(keys, file)
  expect_error(summary_env$birdReadPairTables(file, config), "response units")
  keys$densityUnits[1] <- "males/km2"
  keys$totalAbundance[1] <- NA_real_
  data.table::fwrite(keys, file)
  expect_error(summary_env$birdReadPairTables(file, config), "missing values")
})
