test_that('make_hicexp works', {
  data("HCT116_r1", "HCT116_r2", "HCT116_r3", "HCT116_r4", "HCT116_r5","HCT116_r6", "HCT116_r7")
  groups <- c(1, 1, 1, 2, 2, 2, 2)
  covariates <- data.frame(enzyme = c('mobi', 'mboi', 'mboi', 'dpnii',
                                      'dpnii', 'dpnii', 'dpnii'), 
                           batch = c(1, 2, 1, 2, 1, 2, 2))
  hicexp_new <- make_hicexp(HCT116_r1, HCT116_r2, HCT116_r3, HCT116_r4, HCT116_r5, HCT116_r6, HCT116_r7, 
                            groups = groups, covariates = covariates)
  expect_is(hicexp_new, "Hicexp")
  
  # test for errors on wrong input
  expect_error(make_hicexp(HCT116_r1, HCT116_r2, HCT116_r3, HCT116_r4, HCT116_r5, HCT116_r6, 
                           groups = groups, covariates = covariates),
               "Length of groups must equal the number of Hi-C data objects entered")
  expect_error(make_hicexp(HCT116_r1, HCT116_r2, HCT116_r3, HCT116_r4, HCT116_r5, HCT116_r6, HCT116_r7, 
                           groups = c(1,1,1,2,2,2), 
                           covariates = covariates), 
               "Length of groups must equal the number of Hi-C data objects entered")
  
  # providing list or data give same results
  dat <- list(HCT116_r1, HCT116_r2, HCT116_r3, HCT116_r4, HCT116_r5, HCT116_r6, HCT116_r7)
  hicexp1 <- make_hicexp(HCT116_r1, HCT116_r2, HCT116_r3, HCT116_r4, HCT116_r5, HCT116_r6, HCT116_r7, 
                         groups = groups, covariates = covariates)
  hicexp2 <- make_hicexp(data_list = dat, groups = groups, 
                         covariates = covariates)
  expect_equal(hicexp1, hicexp2)
})
