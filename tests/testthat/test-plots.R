vdiffr_skip_stale()

context("Testing plots")

test_that("Volcano plot", {
  data("sample_data")
  
  result <- compute_amanida(sample_data)
  
  volc <- volcano_plot(result)
  
  expect_doppelganger("Volcano plot of adapted meta-analysis results", volc)
})

test_that("Vote plot", {
  data("sample_data")
  
  result <- compute_amanida(sample_data)
  
  votep <- vote_plot(result)
  
  expect_doppelganger("Vote plot", votep)
})

test_that("Explore plot", {
  data("sample_data")
  
  result <- compute_amanida(sample_data)
  
  expl <- explore_plot(sample_data, counts = 2)
  
  expect_doppelganger("Explore plot", expl)
})
