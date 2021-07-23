
context("Quantitative meta-analysis")

test_that("Quantitative analysis of data", {
  data("sample_data")
  
  result <- compute_amanida(sample_data)
  
  expect_equal(result@stat$fc[8], 0.6389)
  expect_equal(result@stat$N_total[6], 342)
  expect_equal(result@stat$pval[2], 0.0189)
  
  expect_equal(result@vote$votes[100], -1)
  expect_equal(result@vote$articles[48], 1)
  expect_equal(result@vote$vote_counting[62], -1)
})

