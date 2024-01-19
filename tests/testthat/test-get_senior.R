test_that("get_senior", {
  expect_true(Spec2Annot::get_senior("C6H12O3"))
  temp <- Spec2Annot::get_senior("C6H12O3", global = FALSE)
  expect_true(all(temp))
  expect_true(is.vector(temp))
  expect_true(is.logical(temp))
  expect_false(Spec2Annot::get_senior("H2O"))
})
