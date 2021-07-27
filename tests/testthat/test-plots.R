
context("Testing plots")

test_that("Volcano plot", {
  set.seed(123)
  data("sample_data")
  
  result <- compute_amanida(sample_data)
  
  volc <- volcano_plot(result, cutoff = c(0.05, 4))
  
  vdiffr::expect_doppelganger("Volcano plot of adapted meta-analysis results", volc)
})

test_that("Vote plot", {
  set.seed(123)
  data("sample_data")
  
  result <- compute_amanida(sample_data)
  
  votep <- vote_plot(result)
  
  vdiffr::expect_doppelganger("Vote plot", votep)
})

test_that("Explore plot", {
  set.seed(123)
  data("sample_data")
  
  result <- compute_amanida(sample_data)
  
  expl_a <- explore_plot(sample_data, type = "all", counts = 3)
  
  vdiffr::expect_doppelganger("Explore plot", expl_a)
  
  expl_m <- explore_plot(sample_data, type = "mix", counts = 2)
  
  vdiffr::expect_doppelganger("Explore plot", expl_m)
  
  expl_s <- explore_plot(sample_data, type = "sub", counts = 2)
  
  vdiffr::expect_doppelganger("Explore plot", expl_s)
})


