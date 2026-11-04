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
  files <- summary_env$birdPlotAcrossReplicateDensity(pairs, root)
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
