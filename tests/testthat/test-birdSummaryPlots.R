test_that("each species plot contains its replicate values and annual median", {
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
  expected_medians <- list(c(4, 6), c(400, 600), c(2, 3))
  for (i in seq_along(metrics)) {
    files <- summary_env$birdPlotAcrossReplicateMetrics(pairs, root, metrics[i])
    expect_identical(basename(files), paste0("birdAcrossReplicate", stems[i], "_GCKI.png"))
    expect_equal(calls[[3 * i]], expected_medians[[i]])
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

test_that("five-replicate plot trajectory is the median, not the mean", {
  root <- withr::local_tempdir()
  pairs <- data.table::data.table(
    species = "GCKI",
    repID = paste0("rep", 1:5),
    year = 2020L,
    meanDensity = c(1, 2, 3, 4, 100)
  )
  plotted <- list()
  labels <- character()
  original_lines <- graphics::lines
  original_legend <- graphics::legend
  testthat::local_mocked_bindings(
    lines = function(x, y, ...) {
      plotted[[length(plotted) + 1L]] <<- y
      original_lines(x, y, ...)
    },
    legend = function(x, legend, ...) {
      labels <<- legend
      original_legend(x, legend, ...)
    },
    .package = "graphics"
  )

  summary_env$birdPlotAcrossReplicateMetrics(pairs, root, "meanDensity")

  expect_length(plotted, 6L)
  expect_equal(plotted[[6L]], 3)
  expect_false(isTRUE(all.equal(plotted[[6L]], mean(pairs$meanDensity))))
  expect_true("Across-replicate median" %in% labels)
})

test_that("plot trajectories have no point markers and years have vertical guides", {
  root <- withr::local_tempdir()
  pairs <- data.table::data.table(
    species = "GCKI", repID = paste0("rep", 1:5),
    year = rep(c(2020L, 2120L), each = 5),
    meanDensity = seq_len(10)
  )
  point_calls <- 0L
  vertical_guides <- list()
  original_points <- graphics::points
  original_abline <- graphics::abline
  testthat::local_mocked_bindings(
    points = function(...) {
      point_calls <<- point_calls + 1L
      original_points(...)
    },
    legend = function(...) invisible(NULL),
    abline = function(..., h = NULL, v = NULL, col = NULL, lty = NULL, lwd = NULL) {
      if (!is.null(v)) {
        vertical_guides[[length(vertical_guides) + 1L]] <<- list(v = v, col = col, lty = lty)
      }
      args <- list(...)
      if (!is.null(h)) args$h <- h
      if (!is.null(v)) args$v <- v
      if (!is.null(col)) args$col <- col
      if (!is.null(lty)) args$lty <- lty
      if (!is.null(lwd)) args$lwd <- lwd
      do.call(original_abline, args)
    },
    .package = "graphics"
  )

  summary_env$birdPlotAcrossReplicateMetrics(pairs, root, "meanDensity")

  expect_identical(point_calls, 0L)
  expect_length(vertical_guides, 1L)
  expect_identical(vertical_guides[[1L]]$v, c(2020L, 2120L))
  expect_identical(vertical_guides[[1L]]$col, "grey90")
  expect_identical(vertical_guides[[1L]]$lty, "dotted")
})

test_that("FireSense reference shading uses exact visible years in every metric plot", {
  root <- withr::local_tempdir()
  pairs <- data.table::data.table(
    species = "GCKI", repID = "rep1", year = c(2020L, 2720L, 3020L, 3120L),
    meanDensity = c(2, 4, 6, 8), medianDensity = c(1, 2, 3, 4),
    totalAbundance = c(200, 400, 600, 800)
  )
  original <- data.table::copy(pairs)
  bands <- list()
  labels <- list()
  draw_order <- character()
  original_rect <- graphics::rect
  original_legend <- graphics::legend
  original_lines <- graphics::lines
  testthat::local_mocked_bindings(
    rect = function(xleft, ybottom, xright, ytop, ..., col) {
      if (identical(col, "#EEE9D5")) {
        bands[[length(bands) + 1L]] <<- c(xleft, xright, ybottom, ytop)
        draw_order <<- c(draw_order, "band")
      }
      original_rect(xleft, ybottom, xright, ytop, ..., col = col)
    },
    legend = function(x, legend, ...) {
      labels[[length(labels) + 1L]] <<- legend
      original_legend(x, legend, ...)
    },
    lines = function(...) {
      draw_order <<- c(draw_order, "line")
      original_lines(...)
    },
    .package = "graphics"
  )
  for (metric in c("meanDensity", "totalAbundance", "medianDensity")) {
    files <- summary_env$birdPlotAcrossReplicateMetrics(pairs, root, metric)
    expect_true(all(file.info(files)$size > 0))
  }
  expect_length(bands, 3)
  for (i in seq_along(bands)) {
    expect_equal(bands[[i]][1:2], c(2720, 3020))
    expect_equal(bands[[i]][3], 0)
    expect_identical(tail(labels[[i]], 1), "FireSense NRV reference years")
  }
  expect_equal(bands[[1]][4], 8 * 1.08)
  expect_equal(bands[[2]][4], 800 * 1.08)
  expect_equal(bands[[3]][4], 4 * 1.08)
  expect_identical(draw_order, rep(c("band", "line", "line"), 3))
  expect_equal(pairs, original)

  partial <- data.table::copy(pairs[1:2])
  partial$year <- c(2820L, 2920L)
  summary_env$birdPlotAcrossReplicateMetrics(partial, root, "meanDensity")
  expect_equal(bands[[4]][1:2], c(2820, 2920))
  partial$year <- c(2020L, 2120L)
  summary_env$birdPlotAcrossReplicateMetrics(partial, root, "meanDensity")
  expect_length(bands, 4)
  expect_false("FireSense NRV reference years" %in% labels[[5]])
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

test_that("NRV plot lines show species and metric specific sample bounds and counts", {
  root <- withr::local_tempdir()
  pairs <- data.table::CJ(species = c("GCKI", "RCKI"), repID = c("rep1", "rep2"),
                          year = c(2720L, 2820L, 2920L, 3020L))
  pairs[, `:=`(scenario = "reference", climateScenario = "forcing", studyAreaName = "fixture")]
  pairs[, meanDensity := (year - 2720) / 100 + ifelse(repID == "rep1", 1, 5) +
          ifelse(species == "RCKI", 100, 0)]
  pairs[, `:=`(medianDensity = meanDensity / 2, totalAbundance = meanDensity * 10)]
  bounds <- summary_env$birdNRVBounds(pairs, list(species = c("GCKI", "RCKI"),
                                                reps = c("rep1", "rep2")))
  lines <- list()
  subtitles <- character()
  original_abline <- graphics::abline
  original_mtext <- graphics::mtext
  testthat::local_mocked_bindings(
    abline = function(..., h, col) {
      if (col == "#9A4B25") lines[[length(lines) + 1L]] <<- h
      original_abline(..., h = h, col = col)
    },
    mtext = function(text, ...) {
      subtitles <<- c(subtitles, text)
      original_mtext(text, ...)
    },
    .package = "graphics"
  )
  for (metric in c("meanDensity", "totalAbundance", "medianDensity")) {
    summary_env$birdPlotAcrossReplicateMetrics(pairs, root, metric, bounds)
  }
  expect_equal(lines, list(c(1, 8), c(101, 108), c(10, 80), c(1010, 1080),
                          c(0.5, 4), c(50.5, 54)))
  expect_length(subtitles, 6)
  expect_true(all(grepl("2 replicates, 8 observations; not confidence limits", subtitles)))
  expect_error(summary_env$birdPlotAcrossReplicateMetrics(pairs, root, "meanDensity", bounds[species == "GCKI"]),
               "missing or invalid NRV plot bounds")
})

test_that("legends and their text remain in a separate panel outside the data", {
  root <- withr::local_tempdir()
  pairs <- data.table::CJ(species = "REVI", repID = c("rep1", "rep2"),
                          year = c(2720L, 2820L, 2920L, 3020L))
  pairs[, `:=`(scenario = "reference", climateScenario = "forcing", studyAreaName = "fixture",
               meanDensity = 0.4, medianDensity = 0.3, totalAbundance = 400)]
  bounds <- summary_env$birdNRVBounds(pairs, list(species = "REVI", reps = c("rep1", "rep2")))
  data_panels <- list()
  legend_panels <- list()
  original_plot <- graphics::plot.default
  original_legend <- graphics::legend
  testthat::local_mocked_bindings(
    plot.default = function(...) {
      result <- original_plot(...)
      data_panels[[length(data_panels) + 1L]] <<- graphics::par("fig")
      result
    },
    legend = function(...) {
      result <- original_legend(...)
      figure <- graphics::par("fig")
      user <- graphics::par("usr")
      expect_gte(result$rect$left, user[1])
      expect_lte(result$rect$left + result$rect$w, user[2])
      expect_lte(result$rect$top, user[4])
      expect_gte(result$rect$top - result$rect$h, user[3])
      legend_panels[[length(legend_panels) + 1L]] <<- figure
      result
    },
    .package = "graphics"
  )
  for (metric in c("meanDensity", "totalAbundance", "medianDensity")) {
    for (with_bounds in c(FALSE, TRUE)) {
      files <- summary_env$birdPlotAcrossReplicateMetrics(
        pairs, root, metric, if (with_bounds) bounds else NULL
      )
      expect_true(all(file.info(files)$size > 0))
    }
  }
  expect_length(data_panels, 6)
  expect_length(legend_panels, 6)
  for (i in seq_along(data_panels)) {
    expect_lte(legend_panels[[i]][4], data_panels[[i]][3])
  }
})
