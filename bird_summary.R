defineModule(sim, list(
  name = "bird_summary",
  description = paste("Save pair/year metrics referencing existing single-bootstrap bird mosaics.",
                      "In standalone multi mode, summarize completed runs and plot replicate",
                      "metric trajectories from registered files or disk."),
  keywords = c("birds", "BAM", "summary"),
  authors = person("Sourav", "Das", role = c("aut", "cre")),
  childModules = character(0),
  version = list(bird_summary = "0.1.10"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.md", "bird_summary.Rmd"),
  loadOrder = list(after = "bird_modelPredict"),
  reqdPkgs = list("SpaDES.core (>= 3.2.1.9037)", "data.table", "terra"),
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
    defineParameter("plotMetrics", "character",
                    c("meanDensity", "totalAbundance", "medianDensity"), NA, NA,
                    paste("Multi-mode metrics to plot per species: meanDensity (area-weighted males/km2),",
                          "totalAbundance (predicted males), medianDensity (unweighted pixel males/km2).",
                          "Default plots all three; character(0) disables plots. CSV metrics are unchanged.",
                          "Plots shade 2720-3020 as FireSense NRV reference years where visible.",
                          "NRV bound lines appear only when calculateNRV is TRUE.")),
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
                       "references an existing single-bootstrap density GeoTIFF in males/km2.",
                       "Files must remain accessible for multi mode; no raster copies are written."), sourceURL = NA),
    expectsInput("modelIndex", "data.table",
                 "Optional single-mode species selection source from bird_modelLoad.", sourceURL = NA),
    expectsInput("outputsDF", "data.table",
                 "Optional multi-mode combined outputs(sim) records with a file column; never assembled automatically.", sourceURL = NA)
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
    rows[[i]] <- birdRasterMetrics(rasters[[i]])
    rows[[i]][, `:=`(
      species = config$species[i], year = year, repID = config$repID,
      bootstrapID = config$bootstrapID, scenario = config$scenario,
      climateScenario = config$climateScenario, studyAreaName = config$studyAreaName,
      rasterPath = saved$paths[i]
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
  birdValidatePairRasters(pairs)
  annual <- birdAcrossPairMetrics(pairs)
  bounds <- if (config$calculateNRV) birdNRVBounds(pairs, config) else NULL
  destination <- birdSummaryDirectory(sim, config, "multi")
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
  for (metric in config$plotMetrics) {
    plot_files <- birdPlotAcrossReplicateMetrics(pairs, destination, metric, bounds)
    for (file in plot_files) sim <- registerOutputs(file, sim)
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
