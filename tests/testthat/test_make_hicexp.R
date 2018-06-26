test_that('make_hicexp works', {
  library(testthat)
  library(HiCcompare2)
  data("sparse_mat")
  data("r1", "r2", "r3", "r4", "r5","r6", "r7", "hicexp")
  groups <- c(1, 1, 1, 2, 2, 2, 2)
  covariates <- data.frame(enzyme = c('mobi', 'mboi', 'mboi', 'dpnii', 'dpnii', 'dpnii', 'dpnii'), batch = c(1, 2, 1, 2, 1, 2, 2))
  hicexp_new <- make_hicexp(r1, r2, r3, r4, r5, r6, r7, groups = groups, covariates = covariates)
  expect_is(hicexp_new, "hicexp")
  
  # test for errors on wrong input
  expect_error(make_hicexp(r1, r2, r3, r4, r5, r6, groups = groups, covariates = covariates), "Length of groups must equal the number of Hi-C data objects entered")
  expect_error(make_hicexp(r1, r2, r3, r4, r5, r6, r7, groups = c(1,1,1,2,2,2), covariates = covariates), "Length of groups must equal the number of Hi-C data objects entered")
  
  # providing list or data give same results
  dat <- list(r1, r2, r3, r4, r5, r6, r7)
  hicexp1 <- make_hicexp(r1, r2, r3, r4, r5, r6, r7, groups = groups, covariates = covariates)
  hicexp2 <- make_hicexp(data_list = dat, groups = groups, covariates = covariates)
  expect_equal(hicexp1, hicexp2)
})
