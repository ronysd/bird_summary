summary_env <- new.env(parent = globalenv())
module_code <- parse(testthat::test_path("..", "..", "bird_summary.R"))
for (expr in module_code[-1L]) eval(expr, summary_env)
for (file in list.files(testthat::test_path("..", "..", "R"), full.names = TRUE)) {
  sys.source(file, summary_env)
}
summary_env$params <- function(sim) sim$parameters
summary_env$start <- function(sim) 2020
summary_env$end <- function(sim) 2021

summary_parameters <- function(mode = "single", rep = "rep1", boot = 31) {
  list(mode = mode, reps = 1L:2L, years = 2020L:2021L,
       species = c("GCKI", "RCKI"), repID = rep, bootstrapID = boot,
       scenario = "reference", climateScenario = "CanESM5_SSP370",
       .studyAreaName = "fixture")
}

bird_predictions <- function(boot = 31, offset = 0) {
  predictions <- list()
  for (species in c("GCKI", "RCKI")) {
    predictions[[species]] <- list()
    for (year in 2020:2021) {
      raster <- terra::rast(nrows = 2, ncols = 2, xmin = 0, xmax = 2000,
                            ymin = 0, ymax = 2000, crs = "EPSG:3857")
      terra::values(raster) <- c(1, 2, 3, NA) + offset + (year - 2020)
      names(raster) <- paste0("b", boot, "_", year)
      predictions[[species]][[as.character(year)]] <- raster
    }
  }
  predictions
}

local_summary_runtime <- function(envir = parent.frame()) {
  withr::local_options(repos = c(CRAN = "https://cloud.r-project.org"), .local_envir = envir)
  testthat::local_mocked_bindings(
    Require = function(...) {
      for (pkg in c("SpaDES.core", "data.table", "terra")) library(pkg, character.only = TRUE)
      invisible(TRUE)
    },
    .package = "Require", .env = envir
  )
}

run_summary_sim <- function(root, p, objects = list(), label = p$repID) {
  p$simOutputPath <- root
  sim <- SpaDES.core::simInit(
    times = list(start = 2020, end = 2021),
    params = list(bird_summary = p, .globals = list(verbose = FALSE)),
    modules = "bird_summary", objects = objects,
    paths = list(modulePath = testthat::test_path("..", "..", ".."),
                 inputPath = file.path(root, "inputs"),
                 outputPath = file.path(root, label),
                 cachePath = file.path(root, "cache"),
                 scratchPath = file.path(root, "scratch"))
  )
  SpaDES.core::spades(sim, .plots = "none", progress = FALSE)
}
