context("Report meta-analysis")

test_that("Report function quan", {
  local_edition(3)
  column_id = c("Compound Name", "P-value", "Fold-change", "N total", "References")
  input_file <- system.file("extdata", "dataset2.csv", package = "amanida")
  
  expect_snapshot(amanida_report(input_file, separator = ";", column_id, 
                                 analysis_type = "quan", pvalue_cutoff = 0.05, 
                                 fc_cutoff = 4, votecount_lim = 2))
  expect_snapshot(amanida_report(input_file, separator = ";", column_id, 
                                 analysis_type = "quan", pvalue_cutoff = 0.05, 
                                 fc_cutoff = 4, votecount_lim = 2))
})

test_that("Report function qual", {
  local_edition(3)
  column_id = c("Compound Name", "Behaviour", "References")
  input_file <- system.file("extdata", "dataset2.csv", package = "amanida")
  
  expect_snapshot(amanida_report(input_file, separator = ";", 
                                 column_id, analysis_type = "qual", 
                                 votecount_lim = 2))
  expect_snapshot(amanida_report(input_file, separator = ";", 
                                 column_id, analysis_type = "qual", 
                                 votecount_lim = 2))
})


