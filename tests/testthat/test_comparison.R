test_that('exact test works', {
  data("hicexp_diff", "HCT116_r1", "HCT116_r2", "HCT116_r3", "HCT116_r4")
  hicexp_new <- hic_exactTest(hicexp_diff)
  expect_gt(nrow(hicexp_new@comparison), 1)
  
  # test for errors on wrong input
  hicexp_new <- suppressWarnings(make_hicexp(HCT116_r1, 
                                             HCT116_r2, HCT116_r3,
                                             HCT116_r4, groups = c(1,2,3,3)))
  slot(hicexp_new, "normalized") <- TRUE
  # hicexp_new <- fastlo(hicexp_new, verbose = FALSE)
  expect_error(hic_exactTest(hicexp_new))
  
})


test_that('glm works', {
   data("hicexp2")
  hicexp_new <- fastlo(hicexp2, verbose = FALSE)
  d <- model.matrix(~factor(hicexp_new@metadata$group))
  hicexp_new <- hic_glm(hicexp_new, design = d, coef = 2, 
                        method = "QLFTest")
  expect_gt(nrow(hicexp_new@comparison), 1)
  
  # hicexp_new <- hic_glm(hicexp_new, design = d, coef = 2, 
  #                       method = "LRTest")
  # expect_gt(nrow(hicexp_new@comparison), 1)
  # 
  # hicexp_new <- hic_glm(hicexp_new, design = d, coef = 2, 
  #                       method = "Treat")
  # expect_gt(nrow(hicexp_new@comparison), 1)
  
  # test for errors on wrong input
  expect_error(hic_glm(hicexp_new, design = d, coef = 2,contrast = 1, 
                       method = "QLFTest"), 
               "You must enter a value for contrast or a coef, but not both")
  expect_error(hic_glm(hicexp_new, design = d, coef = NA,contrast = NA,
                       method = "QLFTest"), 
               "You must enter a value for contrast or a coef, but not both")
  
})
