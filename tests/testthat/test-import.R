
context("Import data")

test_that("Read table function works", {
  coln = c("Compound Name", "P-value", "Fold-change", "N total", "References")
  input_file <- system.file("extdata", "dataset2.csv", package = "amanida")
  datafile <- amanida_read(input_file, mode = "quan", coln, separator=";")
  
  expect_equal(datafile$pvalue[5], 1e-3)
  
})

