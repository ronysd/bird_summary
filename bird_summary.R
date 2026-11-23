defineModule(sim, list(
  name = "bird_summary",
  description = paste("Save pair/year metrics referencing existing single-bootstrap bird mosaics.",
                      "In standalone multi mode, summarize completed runs and plot replicate",
                      "metric trajectories from registered files or disk."),
  keywords = c("birds", "BAM", "summary"),
  authors = person("Sourav", "Das", role = c("aut", "cre")),
  childModules = character(0),
  version = list(bird_summary = "0.1.15"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.md", "bird_summary.Rmd"),
  loadOrder = list(after = "bird_modelPredict"),
  reqdPkgs = list("SpaDES.core (>= 3.2.1.9037)", "data.table", "terra",
                  "reproducible", "borealbirds/BAMexploreR@e500693af07104e096642ac31aaa731497718f19"),
  parameters = bindrows(
    defineParameter("mode", "character", "single", NA, NA,
                    "single saves the current replicate; multi postprocesses files from completed runs."),
    defineParameter("simOutputPath", "character", outputPath(sim), NA, NA,
                    "Root to scan recursively in multi mode when combined outputsDF is not supplied."),
    defineParameter("reps", "integer", 1L:10L, 1L, NA_integer_,
                    "Expected numeric replicate identifiers in multi mode, corresponding to rep1, rep2, etc."),
    defineParameter("years", "integer", NA_integer_, NA, NA,
                    paste("Exact simulation years to summarize, not a start/end range.",
                          "Required in multi mode. Single mode may use bird_modelPredict's calendar.")),
    defineParameter("species", "character", character(0), NA, NA,
                    paste("Expected species codes. Required in multi mode.",
                          "Single mode may use modelIndex filtered by bird_modelLoad's species selection.")),
    defineParameter("predictionDensityUnits", "character", "males/ha", NA, NA,
                    paste("Raw saved prediction units. Unmodified BAM V5 models return males/ha.",
                          "Use males/km2 only for predictions already converted to that scale.",
                          "All reported density metrics use males/ha, matching BAM.",
                          "Legacy pair records are checked and regenerated from existing rasters in multi mode.")),
    defineParameter("currentCondition", "character", "BAM", NA, NA,
                    paste("Multi-mode current-condition comparison: BAM, SIM, or both.",
                          "BAM uses the official 2020 product and SIM uses 2020 replicate predictions.",
                          "BAM and both download missing benchmark data to inputPath; SIM never downloads BAM.")),
    defineParameter("bamRegion", "character", "Canada", NA, NA,
                    paste("BAM V5 regional product to use: Canada, Alaska, or Lower48.",
                          "Select the region containing the reporting area; products are not combined.")),
    defineParameter("plotMetrics", "character",
                    c("meanDensity", "totalAbundance", "medianDensity"), NA, NA,
                    paste("Multi-mode metrics to plot per species: meanDensity (area-weighted males/ha),",
                          "totalAbundance (predicted males), medianDensity (unweighted pixel males/ha).",
                          "Default plots all three; character(0) disables plots. CSV metrics are unchanged.",
                          "Plots shade 2720-3020 as FireSense NRV reference years where visible.",
                          "SIM and both show a red dashed line for the 2020 replicate median when selected.",
                          "Y-axis ticks and gridlines use equal intervals within each plot.",
                          "NRV bound lines appear only when calculateNRV is TRUE.",
                          "With NRV enabled, meanDensity also produces multi-species median and",
                          "sample min-max whisker plots in males/km2.")),
    defineParameter("calculateNRV", "logical", FALSE, NA, NA,
                    paste("Multi mode only: calculate sample min-max NRV bounds for all three metrics",
                          "from exactly 2720, 2820, 2920, 3020 across every selected replicate.",
                          "All four years must be selected and complete. Saves bounds and adds plot lines.",
                          "FALSE retains descriptive summaries only; no confidence limits are calculated.")),
    defineParameter("repID", "character", NA_character_, NA, NA,
                    "Single-mode replicate label, e.g. rep1; NA uses bird_modelPredict's repID."),
    defineParameter("bootstrapID", "numeric", NA_real_, NA, NA,
                    "Single-mode BAM bootstrap ID; NA uses the predictor's sole bootstrapIDs value."),
    defineParameter("scenario", "character", NA_character_, NA, NA,
                    "Required scenario label, kept separate in all summaries and output paths."),
    defineParameter("climateScenario", "character", NA_character_, NA, NA,
                    "Required climate forcing label, e.g. CanESM5_SSP370."),
    defineParameter(".studyAreaName", "character", NA_character_, NA, NA,
                    "Required reporting-area label; use the same label and raster support across runs."),
    defineParameter(".useCache", "logical", FALSE, NA, NA,
                    "Caching is disabled by default so saved files are checked and registered on every run.")
  ),
  inputObjects = bindrows(
    expectsInput("summaries", "list",
                 paste("Single mode only: bird_modelPredict's [[species]][[year]]$meanPath",
                       "references an existing single-bootstrap density GeoTIFF in predictionDensityUnits.",
                       "Files must remain accessible for multi mode; no raster copies are written."), sourceURL = NA),
    expectsInput("modelIndex", "data.table",
                 "Optional single-mode species selection source from bird_modelLoad.", sourceURL = NA),
    expectsInput("outputsDF", "data.table",
                 "Optional multi-mode combined outputs(sim) records with a file column; never assembled automatically.", sourceURL = NA),
    expectsInput("bamBenchmarkFiles", "data.table",
                 paste("Optional multi-mode V5 2020 benchmark manifest: species, year, rasterPath,",
                       "densityUnits (males/ha), modelVersion (v5), sourceURL, nBootstraps,",
                       "ensembleStatistic (mean). Each raster must have a mean band.",
                       "In BAM or both mode, supplied files are used instead of downloads.",
                       "A fixed per-species BAM valid-data mask is applied",
                       "to every simulation year/replicate and the benchmark; valid zeros remain.",
                       "Density uses bilinear alignment; validity uses nearest neighbour.",
                       "Masked outputs and mask rasters are saved under multi/bam_supported.",
                       "Missing benchmarks are downloaded when no manifest is supplied.",
                       "Source rasters are not changed."), sourceURL = NA)
  ),
  outputObjects = data.frame()
))

doEvent.bird_summary <- function(sim, eventTime, eventType) {
  switch(eventType,
    init = {
      config <- birdSummaryConfig(sim, resolveSpecies = FALSE)
      if (config$mode == "single") {
        for (year in config$years) {
          sim <- scheduleEvent(sim, year, "bird_summary", "save_single", .last())
        }
      } else {
        sim <- InitMulti(sim)
      }
    },
    save_single = {
      sim <- SaveSingle(sim)
    },
    stop("bird_summary: unknown event type: ", eventType)
  )
  return(invisible(sim))
}

SaveSingle <- function(sim) {
  config <- birdSummaryConfig(sim)
  year <- as.integer(time(sim))
  if (!year %in% config$years) stop("bird_summary: unexpected save year: ", year)
  saved <- birdPairRasters(sim$summaries, config, year)
  rasters <- saved$rasters
  rows <- vector("list", length(rasters))
  for (i in seq_along(rasters)) {
    rows[[i]] <- birdRasterMetrics(rasters[[i]], config$predictionDensityUnits)
    rows[[i]][, `:=`(
      species = config$species[i], year = year, repID = config$repID,
      bootstrapID = config$bootstrapID, scenario = config$scenario,
      climateScenario = config$climateScenario, studyAreaName = config$studyAreaName,
      rasterPath = saved$paths[i], sourceDensityUnits = config$predictionDensityUnits,
      reportingSupport = "full_prediction"
    )]
  }
  destination <- birdSummaryDirectory(sim, config, config$repID, paste0("year_", year))
  table <- data.table::rbindlist(rows)
  file <- file.path(destination, paste0("birdPairSummary_year", year, "_",
                                       config$repID, "_boot-", config$bootstrapID, ".csv"))
  data.table::fwrite(table, file)
  sim <- registerOutputs(file, sim)
  message("bird_summary: saved ", nrow(table), " species for ", config$repID, " year ", year)
  return(invisible(sim))
}

InitMulti <- function(sim) {
  config <- birdSummaryConfig(sim)
  files <- birdSummaryFiles(sim, config)
  pairs <- birdReadPairTables(files, config)
  pairs <- birdValidatePairRasters(pairs, config$predictionDensityUnits)
  manifest <- if (config$currentCondition == "SIM") NULL else
    birdPrepareBAMBenchmark(sim, sim$bamBenchmarkFiles, config)
  comparison <- birdBAMComparison(manifest, pairs, config)
  pairs <- comparison$pairs
  benchmark <- comparison$benchmark
  annual <- birdAcrossPairMetrics(pairs)
  bounds <- if (config$calculateNRV) birdNRVBounds(pairs, config) else NULL
  destination <- if (is.null(benchmark)) birdSummaryDirectory(sim, config, "multi") else
    birdSummaryDirectory(sim, config, "multi", "bam_supported")
  for (species in names(comparison$masks)) {
    file <- file.path(destination, paste0("birdBAMReportingMask_", species, ".tif"))
    terra::writeRaster(comparison$masks[[species]], file, overwrite = TRUE)
    selected_pairs <- which(pairs$species == species)
    selected_benchmark <- which(benchmark$species == species)
    data.table::set(pairs, selected_pairs, "reportingMaskPath", file)
    data.table::set(benchmark, selected_benchmark, "reportingMaskPath", file)
    sim <- registerOutputs(file, sim)
  }
  if (!is.null(benchmark)) {
    support <- benchmark[, .(species, reportingSupport, reportingMaskPath, nValidCells, validAreaKm2)]
    annual <- merge(annual, support, by = "species", sort = FALSE)
    if (!is.null(bounds)) bounds <- merge(bounds, support, by = "species", sort = FALSE)
  }
  if (!is.null(bounds)) {
    file <- file.path(destination, "birdNRVBounds.csv")
    data.table::fwrite(bounds, file)
    sim <- registerOutputs(file, sim)
  }
  file <- file.path(destination, "birdAcrossPairSummary.csv")
  data.table::fwrite(annual, file)
  sim <- registerOutputs(file, sim)
  file <- file.path(destination, "birdPairSummary_all.csv")
  data.table::fwrite(pairs, file)
  sim <- registerOutputs(file, sim)
  if (!is.null(benchmark)) {
    file <- file.path(destination, "birdBAMBenchmark2020.csv")
    data.table::fwrite(benchmark, file)
    sim <- registerOutputs(file, sim)
  }
  for (metric in config$plotMetrics) {
    plot_files <- birdPlotAcrossReplicateMetrics(pairs, destination, metric, bounds, benchmark,
                                                 config$currentCondition)
    for (file in plot_files) sim <- registerOutputs(file, sim)
    if (!is.null(bounds) && (config$currentCondition == "BAM" || 2020L %in% config$years)) {
      comparison_files <- birdPlotCurrentBAMNRV(pairs, benchmark, bounds, destination, metric,
                                               config$currentCondition)
      for (file in comparison_files) sim <- registerOutputs(file, sim)
      if (metric == "meanDensity") {
        whisker_files <- birdPlotCurrentNRVWhiskers(pairs, benchmark, bounds, destination,
                                                   config$currentCondition)
        for (file in whisker_files) sim <- registerOutputs(file, sim)
      }
    }
  }
  if (!is.null(bounds) && config$currentCondition == "SIM" && !2020L %in% config$years) {
    message("bird_summary: 2020 is not selected; current-versus-NRV plots were omitted.")
  }
  message("bird_summary: combined ", length(config$reps), " replicates; ",
          if (config$calculateNRV) "sample min-max NRV bounds calculated from four reference years."
          else "NRV is not calculated.")
  return(invisible(sim))
}

.inputObjects <- function(sim) {
  birdSummaryConfig(sim, resolveSpecies = FALSE)
  return(invisible(sim))
}
