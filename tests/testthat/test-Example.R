test_that("This is just an example test which should not fail!", {
  set.seed(1242)
  x <- iris[1,1]
  expect_equal(x, 5.1, tolerance = 0.00001)
})
