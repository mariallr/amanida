
context("Quantitative meta-analysis")

test_that("Quantitative analysis of data", {
 
  data("sample_data")
  
  result <- compute_amanida(sample_data)
  
  expect_equal(result@stat %>% filter(id == "2-Ethyl-5-methylfuran") %>% 
                 pull(fc), 0.6389, tolerance = 0.01)
  expect_equal(result@stat %>% filter(id == "Fumaric acid") %>% 
                 pull(N_total), 204)
  expect_equal(result@stat %>% filter(id == "Phenylacetic acid") %>% 
                 pull(pval), 4.636927e-05)
  
  expect_equal(result@vote %>% filter(id == "3-Phosphoglycerate") %>% 
                 pull(votes), -1)
  expect_equal(result@vote %>% filter(id == "L-Tryptophan") %>% 
                 pull(articles), 3)
  expect_equal(result@vote %>% filter(id == "Tiglic aldehyde") %>% 
                 pull(vote_counting), -1)
})

