test_that("get_nrule", {
  testthat::expect_true(get_nrule(5, 125.10))
  testthat::expect_true(get_nrule(5, 125.50))
  testthat::expect_true(get_nrule(5, 125.70))
  testthat::expect_false(get_nrule(6, 125.10))
  testthat::expect_false(get_nrule(6, 125.50))
  testthat::expect_false(get_nrule(6, 125.70))
})
