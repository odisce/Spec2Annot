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

test_that("Spec2Annot::mz_calc_ion", {
  expect_true(Spec2Annot::mz_calc_ion(351.2564, form = "-H") %between% (350.2491 + c(-0.0001, +0.0001)))
  expect_true(Spec2Annot::mz_calc_ion(351.2564, form = "+H") %between% (352.2637 + c(-0.0001, +0.0001)))
  expect_true(Spec2Annot::mz_calc_ion(351.2564, form = "+NH4") %between% (369.2902 + c(-0.0001, +0.0001)))
})

