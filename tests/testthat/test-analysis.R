
context("Quantitative meta-analysis")

test_that("Quantitative analysis of data", {
 
  data("sample_data")
  
  result <- compute_amanida(sample_data, comp.inf = F)
  
  expect_equal(result@stat %>% dplyr::filter(id == "2-Ethyl-5-methylfuran") %>% 
                 dplyr::pull(fc), 0.6389, tolerance = 0.01)
  expect_equal(result@stat %>% dplyr::filter(id == "Fumaric acid") %>% 
                 dplyr::pull(N_total), 204)
  expect_equal(result@stat %>% dplyr::filter(id == "Phenylacetic acid") %>% 
                 dplyr::pull(pval), 4.636927e-05)
  
  expect_equal(result@vote %>% dplyr::filter(id == "3-Phosphoglycerate") %>% 
                 dplyr::pull(votes), -1)
  expect_equal(result@vote %>% dplyr::filter(id == "L-Tryptophan") %>% 
                 dplyr::pull(articles), 3)
  expect_equal(result@vote %>% dplyr::filter(id == "Tiglic aldehyde") %>% 
                 dplyr::pull(vote_counting), -1)
})

