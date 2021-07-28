context("Qualitative meta-analysis")

test_that("Vote-counting function", {
  set.seed(123)
  coln = c("Compound Name", "Behaviour", "References")
  input_file <- system.file("extdata", "dataset2.csv", package = "amanida")
  data_votes <- amanida_read(input_file, mode = "qual", coln, separator = ";")
  
  vote_result <- amanida_vote(data_votes)
  
  expect_equal(vote_result@vote$votes[9], -1)
  expect_equal(vote_result@vote$vote_counting[4], -1)
})
