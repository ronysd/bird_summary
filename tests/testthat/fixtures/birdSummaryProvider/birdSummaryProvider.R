defineModule(sim, list(
  name = "birdSummaryProvider",
  description = "Tiny local provider publishing its index only during init.",
  keywords = "test",
  authors = person("Sourav", "Das", role = "aut"),
  childModules = character(0),
  version = list(birdSummaryProvider = "0.0.1"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  reqdPkgs = list("data.table", "terra"),
  parameters = bindrows(defineParameter(".useCache", "logical", FALSE, NA, NA,
                                        "Disable caching in this fixture.")),
  inputObjects = bindrows(expectsInput("fixturePredictions", "list",
                                       "Small supplied prediction rasters.")),
  outputObjects = bindrows(
    createsOutput("modelIndex", "data.table", "Species codes published during init."),
    createsOutput("mosaicedPredictions", "list", "Small species/year rasters.")
  )
))

doEvent.birdSummaryProvider <- function(sim, eventTime, eventType) {
  if (eventType != "init") stop("Unexpected fixture event.")
  sim$modelIndex <- data.table::data.table(species = names(sim$fixturePredictions))
  sim$mosaicedPredictions <- sim$fixturePredictions
  invisible(sim)
}
