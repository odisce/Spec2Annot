test_that("Spec2Annot::get_iso_from_annot", {
  temp <- Spec2Annot::get_iso_from_annot("C6H12_13C2")
  expect_true(data.table::is.data.table(temp))
  expect_true(temp[element == "13C", elmt_nb] == 2)
  
  temp <- Spec2Annot::get_iso_from_annot("C6H12_13C4_18O1")
  expect_true(data.table::is.data.table(temp))
  expect_true(temp[element == "13C", elmt_nb] == 4)
  expect_true(temp[element == "18O", elmt_nb] == 1)

  expect_false(Spec2Annot::get_iso_from_annot("C6H12"))
})

