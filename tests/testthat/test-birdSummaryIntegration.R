test_that("inferred species work when an upstream module publishes modelIndex in init", {
  local_summary_runtime()
  root <- withr::local_tempdir()
  module_path <- file.path(root, "modules")
  dir.create(module_path)
  expect_true(file.symlink(normalizePath(testthat::test_path("..", "..")),
                           file.path(module_path, "bird_summary")))
  expect_true(file.copy(testthat::test_path("fixtures", "birdSummaryProvider"),
                         module_path, recursive = TRUE))
  p <- summary_parameters()
  p$species <- character(0)
  sim <- SpaDES.core::simInit(
    times = list(start = 2020, end = 2021),
    modules = c("birdSummaryProvider", "bird_summary"),
    params = list(bird_summary = p, .globals = list(verbose = FALSE)),
    objects = list(fixtureSummaries = bird_saved_summaries(root)),
    paths = list(modulePath = module_path, inputPath = file.path(root, "inputs"),
                 outputPath = file.path(root, "outputs"), cachePath = file.path(root, "cache"),
                 scratchPath = file.path(root, "scratch"))
  )
  expect_null(sim$modelIndex)
  sim <- SpaDES.core::spades(sim, .plots = "none", progress = FALSE)
  expect_setequal(sim$modelIndex$species, c("GCKI", "RCKI"))
  files <- SpaDES.core::outputs(sim)$file
  tables <- files[grepl("\\.csv$", files)]
  expect_length(tables, 2)
  for (file in tables) {
    expect_setequal(data.table::fread(file)$species, c("GCKI", "RCKI"))
  }
})

test_that("single and standalone multi modes communicate through saved registered files", {
  local_summary_runtime()
  root <- withr::local_tempdir()
  registry <- list()
  for (rep in 1:2) {
    p <- summary_parameters(rep = paste0("rep", rep), boot = 30 + rep)
    summaries <- bird_saved_summaries(root, bird_predictions(30 + rep, rep - 1),
                                      p$repID, 30 + rep)
    predictor_files <- list.files(file.path(root, "bird_modelPredict", p$repID),
                                  pattern = "\\.tif$", recursive = TRUE, full.names = TRUE)
    original_hashes <- tools::md5sum(predictor_files)
    sim <- run_summary_sim(root, p, list(summaries = summaries))
    expect_s4_class(sim, "simList")
    expect_null(sim$birdPairSummary)
    expect_null(sim$birdAcrossPairSummary)
    files <- SpaDES.core::outputs(sim)$file
    expect_length(grep("\\.csv$", files), 2)
    expect_length(grep("\\.tif$", files), 0)
    expect_true(all(file.exists(files)))
    expect_identical(tools::md5sum(predictor_files), original_hashes)
    expect_length(list.files(file.path(root, p$repID), pattern = "\\.tif$", recursive = TRUE), 0)
    for (file in files[grepl("\\.csv$", files)]) {
      table <- data.table::fread(file)
      expect_true(all(table$rasterPath %in% normalizePath(predictor_files)))
    }
    registry[[rep]] <- SpaDES.core::outputs(sim)
  }
  combined <- data.table::rbindlist(registry, fill = TRUE)
  p <- summary_parameters("multi")
  multi <- run_summary_sim(root, p, list(outputsDF = combined), label = "postprocess")
  expect_null(multi$mosaicedPredictions)
  expect_null(multi$summaries)
  output_files <- SpaDES.core::outputs(multi)$file
  annual_file <- output_files[basename(output_files) == "birdAcrossPairSummary.csv"]
  expect_length(annual_file, 1)
  plot_files <- output_files[grepl("\\.png$", output_files)]
  expect_length(plot_files, 6)
  expect_true(all(file.exists(plot_files)))
  expect_setequal(basename(plot_files),
                  paste0("birdAcrossReplicate",
                         rep(c("Density", "TotalAbundance", "MedianDensity"), each = 2),
                         "_", c("GCKI", "RCKI"), ".png"))
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
  p$plotMetrics <- "totalAbundance"
  selected <- run_summary_sim(root, p, list(outputsDF = combined), label = "abundance-only")
  selected_files <- SpaDES.core::outputs(selected)$file
  expect_setequal(basename(selected_files[grepl("\\.png$", selected_files)]),
                  paste0("birdAcrossReplicateTotalAbundance_", c("GCKI", "RCKI"), ".png"))
  p$plotMetrics <- character(0)
  no_plots <- run_summary_sim(root, p, list(outputsDF = combined), label = "no-plots")
  expect_false(any(grepl("\\.png$", SpaDES.core::outputs(no_plots)$file)))

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
