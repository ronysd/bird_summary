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
  inputObjects = bindrows(expectsInput("fixtureSummaries", "list",
                                       "Small supplied saved prediction summaries.")),
  outputObjects = bindrows(
    createsOutput("modelIndex", "data.table", "Species codes published during init."),
    createsOutput("summaries", "list", "Species/year entries with existing meanPath TIFFs.")
  )
))

doEvent.birdSummaryProvider <- function(sim, eventTime, eventType) {
  if (eventType != "init") stop("Unexpected fixture event.")
  sim$modelIndex <- data.table::data.table(species = names(sim$fixtureSummaries))
  sim$summaries <- sim$fixtureSummaries
  invisible(sim)
}
