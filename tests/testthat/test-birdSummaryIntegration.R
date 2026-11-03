test_that("single and standalone multi modes communicate through saved registered files", {
  local_summary_runtime()
  root <- withr::local_tempdir()
  registry <- list()
  for (rep in 1:2) {
    p <- summary_parameters(rep = paste0("rep", rep), boot = 30 + rep)
    sim <- run_summary_sim(root, p, list(mosaicedPredictions = bird_predictions(30 + rep, rep - 1)))
    expect_s4_class(sim, "simList")
    expect_null(sim$birdPairSummary)
    expect_null(sim$birdAcrossPairSummary)
    files <- SpaDES.core::outputs(sim)$file
    expect_length(grep("\\.csv$", files), 2)
    expect_length(grep("\\.tif$", files), 4)
    expect_true(all(file.exists(files)))
    registry[[rep]] <- SpaDES.core::outputs(sim)
  }
  combined <- data.table::rbindlist(registry, fill = TRUE)
  p <- summary_parameters("multi")
  multi <- run_summary_sim(root, p, list(outputsDF = combined), label = "postprocess")
  expect_null(multi$mosaicedPredictions)
  output_files <- SpaDES.core::outputs(multi)$file
  annual_file <- output_files[basename(output_files) == "birdAcrossPairSummary.csv"]
  expect_length(annual_file, 1)
  annual <- data.table::fread(annual_file)
  expect_equal(nrow(annual), 12)
  expect_true(all(annual$nPairs == 2))
  expect_equal(annual[metric == "meanDensity" & year == 2020]$mean, c(2.5, 2.5), tolerance = 0.01)
  expect_equal(annual[metric == "meanDensity" & year == 2021]$median, c(3.5, 3.5), tolerance = 0.01)
  expect_false(any(grepl("NRV|Q05|SD", names(annual))))
  scan <- run_summary_sim(root, p, label = "scan")
  scanned <- data.table::fread(SpaDES.core::outputs(scan)$file[
    basename(SpaDES.core::outputs(scan)$file) == "birdAcrossPairSummary.csv"])
  expect_equal(as.data.frame(scanned), as.data.frame(annual))

  config <- summary_env$birdSummaryConfig(list(parameters = list(bird_summary = p)))
  tables <- combined$file[grepl("\\.csv$", combined$file)]
  pairs <- summary_env$birdReadPairTables(tables, config)
  expect_equal(nrow(pairs), 8)
  expect_error(summary_env$birdReadPairTables(tables[-1], config), "missing expected")
  duplicate <- data.table::fread(tables[1])
  duplicate_file <- file.path(root, "duplicate.csv")
  data.table::fwrite(duplicate, duplicate_file)
  expect_error(summary_env$birdReadPairTables(c(tables, duplicate_file), config), "duplicate")
  removed_registry <- data.table::data.table(file = file.path(root, "birdPairSummary_year2020_rep1_boot-31.csv"))
  expect_error(summary_env$birdSummaryFiles(list(outputsDF = removed_registry), config),
               "registered pair files are missing")

  raster <- terra::rast(pairs$rasterPath[1])
  terra::values(raster) <- c(1, 2, NA, NA)
  terra::writeRaster(raster, pairs$rasterPath[1], overwrite = TRUE)
  expect_error(summary_env$birdValidatePairRasters(pairs), "persisted metrics|reporting support")
})
