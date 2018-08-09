test_that('cyclic_loess works', {
  data("hicexp2")
  hicexp_new <- cyclic_loess(hicexp2, 
                             span = 0.5, verbose = FALSE)
  expect_equal(hicexp_new@normalized, TRUE)
  
  # test for errors on wrong input
  expect_error(cyclic_loess(hicexp_new), 
               "Data has already been normalized.")
  data('hicexp2')
  expect_error(cyclic_loess(hicexp2, span = 0),
               "span must be set to NA or a value between 0 and 1")
  expect_error(cyclic_loess(hicexp2, span = 2), 
               "span must be set to NA or a value between 0 and 1")
})


test_that('fastlo works', {
  data("hicexp2")
  hicexp_new <- fastlo(hicexp2, verbose = FALSE)
  expect_equal(hicexp_new@normalized, TRUE)
  
  # test for errors on wrong input
  expect_error(fastlo(hicexp_new),
               "Data has already been normalized.")
  data('hicexp2')
  expect_error(fastlo(hicexp2, span = 0), 
               "span must be set to NA or a value between 0 and 1")
  expect_error(fastlo(hicexp2, span = 2), 
               "span must be set to NA or a value between 0 and 1") 
})

