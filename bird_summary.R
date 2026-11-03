defineModule(sim, list(
  name = "bird_summary",
  description = paste("Save and register single-bootstrap bird mosaics and pair/year metrics.",
                      "In standalone multi mode, summarize completed runs from registered",
                      "files or disk without loading saved simulations."),
  keywords = c("birds", "BAM", "summary"),
  authors = person("Sourav", "Das", role = c("aut", "cre")),
  childModules = character(0),
  version = list(bird_summary = "0.1.0"),
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
    expectsInput("mosaicedPredictions", "list",
                 "Single mode only: [[species]][[year]] = one-layer SpatRaster in males/km2.", sourceURL = NA),
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
      config <- birdSummaryConfig(sim)
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
  rasters <- birdPairRasters(sim$mosaicedPredictions, config, year)
  rows <- vector("list", length(rasters))
  for (i in seq_along(rasters)) {
    rows[[i]] <- birdRasterMetrics(rasters[[i]])
    rows[[i]][, `:=`(
      species = config$species[i], year = year, repID = config$repID,
      bootstrapID = config$bootstrapID, scenario = config$scenario,
      climateScenario = config$climateScenario, studyAreaName = config$studyAreaName
    )]
  }
  destination <- birdSummaryDirectory(sim, config, config$repID, paste0("year_", year))
  for (i in seq_along(rasters)) {
    file <- file.path(destination, paste0(year, "_", config$species[i], "_mosaic_",
                                         config$repID, "_boot-", config$bootstrapID, ".tif"))
    terra::writeRaster(rasters[[i]], file, overwrite = TRUE)
    rows[[i]][, rasterPath := normalizePath(file, mustWork = TRUE)]
    sim <- registerOutputs(file, sim)
  }
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
  destination <- birdSummaryDirectory(sim, config, "multi")
  file <- file.path(destination, "birdAcrossPairSummary.csv")
  data.table::fwrite(annual, file)
  sim <- registerOutputs(file, sim)
  file <- file.path(destination, "birdPairSummary_all.csv")
  data.table::fwrite(pairs, file)
  sim <- registerOutputs(file, sim)
  message("bird_summary: combined ", length(config$reps), " replicates; NRV is not calculated.")
  return(invisible(sim))
}

.inputObjects <- function(sim) {
  birdSummaryConfig(sim)
  return(invisible(sim))
}
