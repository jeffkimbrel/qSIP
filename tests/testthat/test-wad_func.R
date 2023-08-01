test_that("WAD calculation works", {
  expect_equal(WAD_func(c(1,2,3), c(1,2,3)), 2.33333333)
  expect_equal(WAD_func(c(1,1,1), c(1,2,3)), 2)
})


test_that("vectors have different sizes", {
  expect_error(WAD_func(c(1, 2), c(1)))
})

test_that("non-numeric vectors trigger warning", {
  expect_error(WAD_func(c("A", "B"), c(1, 2)))
  expect_error(WAD_func(c(1, 2), c("A", "B")))
})
