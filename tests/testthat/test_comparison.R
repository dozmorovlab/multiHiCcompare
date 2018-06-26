test_that('exact test works', {
  library(testthat)
  library(HiCcompare2)
  data("r1", "r2", "r3", "r4", "r5", "r6", "hicexp")
  groups <- c(1, 1, 2, 2)
  hicexp_new <- make_hicexp(r1, r2, r3, r4,groups = groups)
  hicexp_new <- fastlo(hicexp_new, verbose = FALSE, Plot = FALSE)
  hicexp_new <- hic_exactTest(hicexp_new, Plot = FALSE)
  expect_gt(nrow(hicexp_new@comparison), 1)
  
  # test for errors on wrong input
  hicexp_new <- make_hicexp(r1, r2, r3, r4, r5, r6, groups = c(1,1,2,2,3,3))
  hicexp_new <- fastlo(hicexp_new, verbose = FALSE, Plot = FALSE)
  expect_error(hic_exactTest(hicexp_new))
  
})


test_that('glm works', {
  library(testthat)
  library(HiCcompare2)
  data("r1", "r2", "r3", "r4", "r5", "r6", "hicexp")
  groups <- c(1, 1, 2, 2)
  hicexp_new <- make_hicexp(r1, r2, r3, r4,groups = groups)
  hicexp_new <- fastlo(hicexp_new, verbose = FALSE, Plot = FALSE)
  d <- model.matrix(~factor(hicexp_new@metadata$group))
  hicexp_new <- hic_glm(hicexp_new, design = d, coef = 2, method = "QLFTest", Plot = FALSE)
  expect_gt(nrow(hicexp_new@comparison), 1)
  
  hicexp_new <- hic_glm(hicexp_new, design = d, coef = 2, method = "LRTest", Plot = FALSE)
  expect_gt(nrow(hicexp_new@comparison), 1)
  
  hicexp_new <- hic_glm(hicexp_new, design = d, coef = 2, method = "Treat", Plot = FALSE)
  expect_gt(nrow(hicexp_new@comparison), 1)
  
  # test for errors on wrong input
  expect_error(hic_glm(hicexp_new, design = d, coef = 2,contrast = 1, method = "QLFTest", Plot = FALSE), 
               "Please enter either a value for contrast or coef but NOT both.")
  expect_error(hic_glm(hicexp_new, design = d, coef = NA,contrast = NA, method = "QLFTest", Plot = FALSE), 
               "You must enter a value for contrast or a coef, but not both")
  
})
