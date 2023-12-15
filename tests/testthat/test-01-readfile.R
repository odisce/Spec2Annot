test_that("Spec2Annot::test", {
  testthat::expect_true(TRUE)
})

test_that("Spec2Annot::test", {
  expect_true(is.data.table(Spec2Annot::get_iso_from_annot("C6H12_13C2")))
})
