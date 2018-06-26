test_that('cyclic_loess works', {
  library(testthat)
  library(HiCcompare2)
  data("r1", "r2", "r3", "r4", "hicexp")
  groups <- c(1, 1, 2, 2)
  hicexp_new <- make_hicexp(r1, r2, r3, r4, groups = groups)
  hicexp_new <- cyclic_loess(hicexp_new, span = 0.5, verbose = FALSE, Plot = FALSE)
  expect_equal(hicexp_new@normalized, TRUE)
  
  # test for errors on wrong input
  expect_error(cyclic_loess(hicexp_new), "Data has already been normalized.")
  expect_error(cyclic_loess(hicexp, span = 0), "span must be set to NA or a value between 0 and 1")
  expect_error(cyclic_loess(hicexp, span = 2), "span must be set to NA or a value between 0 and 1")
})


test_that('fastlo works', {
  library(testthat)
  library(HiCcompare2)
  data("r1", "r2", "r3", "r4", "hicexp")
  groups <- c(1, 1, 2, 2)
  hicexp_new <- make_hicexp(r1, r2, r3, r4,groups = groups)
  hicexp_new <- fastlo(hicexp_new, verbose = FALSE, Plot = FALSE)
  expect_equal(hicexp_new@normalized, TRUE)
  
  # test for errors on wrong input
  expect_error(fastlo(hicexp_new), "Data has already been normalized.")
  expect_error(fastlo(hicexp, span = 0), "span must be set to NA or a value between 0 and 1")
  expect_error(fastlo(hicexp, span = 2), "span must be set to NA or a value between 0 and 1") 
})

