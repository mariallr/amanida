
context("data type")

test_that("invalid data type stops when", {
  expect_error(data.read("bla.R", coln = c("id", "p-val")))
  })

