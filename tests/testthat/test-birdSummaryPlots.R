test_that("each species plot contains only its own replicate and annual mean values", {
  root <- withr::local_tempdir()
  pairs <- data.table::data.table(
    species = rep(c("GCKI", "RCKI"), each = 4),
    repID = rep(c("rep1", "rep1", "rep2", "rep2"), 2),
    year = rep(c(2020L, 2120L), 4),
    meanDensity = c(2, 4, 6, 8, 20, 40, 60, 80)
  )
  original <- data.table::copy(pairs)
  calls <- list()
  original_lines <- graphics::lines
  testthat::local_mocked_bindings(
    lines = function(x, y, ...) {
      calls[[length(calls) + 1L]] <<- list(x = x, y = y)
      original_lines(x, y, ...)
    },
    .package = "graphics"
  )
  files <- summary_env$birdPlotAcrossReplicateMetrics(pairs, root, "meanDensity")
  expect_identical(basename(files),
                   paste0("birdAcrossReplicateDensity_", c("GCKI", "RCKI"), ".png"))
  expect_true(all(file.info(files)$size > 0))
  expect_length(calls, 6)
  expected <- list(c(2, 4), c(6, 8), c(4, 6),
                   c(20, 40), c(60, 80), c(40, 60))
  for (i in seq_along(expected)) {
    expect_identical(calls[[i]]$x, c(2020L, 2120L))
    expect_equal(calls[[i]]$y, expected[[i]])
  }
  expect_equal(pairs, original)
})

test_that("all selectable plots use their metric values, species title and correct units", {
  root <- withr::local_tempdir()
  pairs <- data.table::data.table(
    species = "GCKI", repID = rep(c("rep1", "rep2"), each = 2),
    year = rep(c(2020L, 2120L), 2), meanDensity = c(2, 4, 6, 8),
    medianDensity = c(1, 2, 3, 4), totalAbundance = c(200, 400, 600, 800)
  )
  original_lines <- graphics::lines
  original_plot <- graphics::plot.default
  calls <- list()
  titles <- list()
  testthat::local_mocked_bindings(
    lines = function(x, y, ...) {
      calls[[length(calls) + 1L]] <<- y
      original_lines(x, y, ...)
    },
    plot.default = function(x, ..., main, ylab) {
      titles[[length(titles) + 1L]] <<- list(main = main, ylab = ylab)
      original_plot(x, ..., main = main, ylab = ylab)
    },
    .package = "graphics"
  )
  metrics <- c("meanDensity", "totalAbundance", "medianDensity")
  stems <- c("Density", "TotalAbundance", "MedianDensity")
  expected_means <- list(c(4, 6), c(400, 600), c(2, 3))
  for (i in seq_along(metrics)) {
    files <- summary_env$birdPlotAcrossReplicateMetrics(pairs, root, metrics[i])
    expect_identical(basename(files), paste0("birdAcrossReplicate", stems[i], "_GCKI.png"))
    expect_equal(calls[[3 * i]], expected_means[[i]])
    expect_match(titles[[i]]$main, "GCKI")
  }
  expect_identical(titles[[1]]$ylab,
                   expression("Area-weighted mean density (males/" * km^2 * ")"))
  expect_identical(titles[[2]]$ylab, "Total predicted abundance (males)")
  expect_identical(titles[[3]]$ylab,
                   expression("Median pixel density (males/" * km^2 * ")"))
  expect_error(summary_env$birdPlotAcrossReplicateMetrics(pairs, root, "unknown"),
               "unsupported plot metric")
})

test_that("plot metric choices fail clearly when unsupported or duplicated", {
  p <- summary_parameters("multi")
  for (metrics in list("SD", c("meanDensity", "meanDensity"), NA_character_, 1, NULL)) {
    p["plotMetrics"] <- list(metrics)
    expect_error(summary_env$birdSummaryConfig(list(parameters = list(bird_summary = p))),
                 "plotMetrics")
  }
  p$plotMetrics <- character(0)
  expect_identical(summary_env$birdSummaryConfig(
    list(parameters = list(bird_summary = p)))$plotMetrics, character(0))
})
