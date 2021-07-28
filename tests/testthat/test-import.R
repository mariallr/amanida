
context("Import dataset")

test_that("Read table quantitative function works", {
  set.seed(123)
  coln = c("Compound Name", "P-value", "Fold-change", "N total", "References")
  input_file <- system.file("extdata", "dataset2.csv", package = "amanida")
  datafile <- amanida_read(input_file, mode = "quan", coln, separator=";")
  
  expect_equal(datafile$pvalue[5], 1e-3)
  expect_equal(datafile$foldchange[10], 0.34)
  expect_equal(datafile$N[55], 135)
  expect_equal(datafile$trend[15], 1)
})

test_that("Read table qualitative function works", {
  set.seed(123)
  coln = c("Compound Name", "N total", "References")
  input_file <- system.file("extdata", "dataset2.csv", package = "amanida")
  datafile <- amanida_read(input_file, mode = "qual", coln, separator=";")
  
  expect_equal(datafile$trend[15], 1)
  expect_equal(datafile$id[75], "3-Methylhistidine")
})


