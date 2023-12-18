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

test_that("Spec2Annot::gen_formula_from_compo", {
  temp_form <- Spec2Annot::gen_formula_from_compo("C6H12O3")
  expect_true(is.character(temp_form))
  expect_true(temp_form == "(6*C+12*H+3*O)")
})

test_that("Spec2Annot::element_from_formula", {
  temp_dt <- Spec2Annot::element_from_formula(Spec2Annot::gen_formula_from_compo("C6H12O3"))
  expect_true(data.table::is.data.table(temp_dt))
  expect_true(temp_dt[element == "C", elmt_nb] == 6)
  expect_true(temp_dt[element == "H", elmt_nb] == 12)
  expect_true(temp_dt[element == "O", elmt_nb] == 3)
  temp_dt <- Spec2Annot::element_from_formula("C6H12O3")
  expect_true(data.table::is.data.table(temp_dt))
  expect_true(temp_dt[element == "C", elmt_nb] == 6)
  expect_true(temp_dt[element == "H", elmt_nb] == 12)
  expect_true(temp_dt[element == "O", elmt_nb] == 3)
  temp_dt <- Spec2Annot::element_from_formula(Spec2Annot::gen_formula_from_compo("C6H12O3_13C3"))
  expect_true(data.table::is.data.table(temp_dt))
  expect_true(temp_dt[element == "C", elmt_nb] == 3)
  expect_true(temp_dt[element == "H", elmt_nb] == 12)
  expect_true(temp_dt[element == "O", elmt_nb] == 3)
  expect_true(temp_dt[element == "13C", elmt_nb] == 3)
  temp_dt <- Spec2Annot::element_from_formula("C6H12O3_13C3")
  expect_true(data.table::is.data.table(temp_dt))
  expect_true(temp_dt[element == "C", elmt_nb] == 3)
  expect_true(temp_dt[element == "H", elmt_nb] == 12)
  expect_true(temp_dt[element == "O", elmt_nb] == 3)
  expect_true(temp_dt[element == "13C", elmt_nb] == 3)
})


test_that("Spec2Annot::mz_from_string", {
  require(data.table)
  expect_true(Spec2Annot::mz_from_string("C6H2O3") == 122.00039392317000874754739925265312194824218750)
  expect_true(Spec2Annot::mz_from_string("C6H3O3+") %between% (123.0077 + c(-0.0001, +0.0001)))
  expect_true(Spec2Annot::mz_from_string("C6H1O3+C6-O+H10") %between% (187.0759 + c(-0.0001, +0.0001)))
  expect_true(Spec2Annot::mz_from_string("C6H1O3+C6-O+H10+Ca2+K") %between% (305.9648 + c(-0.0001, +0.0001)))
  expect_true(Spec2Annot::mz_from_string("C6H1O3+C6-O+H10+Ca2+K+H+") %between% (306.9721 + c(-0.0001, +0.0001)))
  expect_true(Spec2Annot::mz_from_string("C6H1O3+C6-O+H10+Ca2+K-H-") %between% (304.9575 + c(-0.0001, +0.0001)))
  expect_true(Spec2Annot::mz_from_string("-(C6H2O)") == -Spec2Annot::mz_from_string("C6H2O"))
  expect_true(Spec2Annot::mz_from_string("-(C6H2O)-(H2O)") == -Spec2Annot::mz_from_string("C6H2O")-Spec2Annot::mz_from_string("H2O"))
  expect_true(Spec2Annot::mz_from_string("-(C6H2O)+(H2O)") == Spec2Annot::mz_from_string("-C6H2O+H2O"))
  expect_true(Spec2Annot::mz_from_string("-(C6H2O)+(H2O)") == Spec2Annot::mz_from_string("-(C6H2O)+H2O"))
  expect_true(Spec2Annot::mz_from_string("-(C6H2O)+(H2O)+") == Spec2Annot::mz_from_string("-(C6H2O)+H2O+"))
  expect_true(Spec2Annot::mz_from_string("-(C6H2O)+(H2O)-") == Spec2Annot::mz_from_string("-(C6H2O)+H2O-"))
  ## mz from Annotation
  expect_true(Spec2Annot::mz_from_string("[C6H2O+H-H2O]+") +
                Spec2Annot::Element[paste0(mass_nb,atomic_symb) == "13C", atomic_mass] -
                Spec2Annot::Element[paste0(mass_nb,atomic_symb) == "12C", atomic_mass] ==
                Spec2Annot::mz_from_string("[C6H2O+H-H2O]+_13C1"))
  expect_true(Spec2Annot::mz_from_string("[C6H2O+H-H2O]+") +
                2*(Spec2Annot::Element[paste0(mass_nb,atomic_symb) == "13C", atomic_mass] - Spec2Annot::Element[paste0(mass_nb,atomic_symb) == "12C", atomic_mass]) ==
                Spec2Annot::mz_from_string("[C6H2O+H-H2O]+_13C2"))
  expect_true(
    round(Spec2Annot::mz_from_string("[C6H2O+H-H2O]") + 2*Spec2Annot::Isotopes_db[isotope == "13C", mass_diff] + 1*Spec2Annot::Isotopes_db[isotope == "18O", mass_diff], 5) ==
      round(Spec2Annot::mz_from_string("[C6H2O+H-H2O]_13C2_18O"), 5)
  )
  expect_true(
    round(Spec2Annot::mz_from_string("[C6H2O+H-H2O]+") + 1*Spec2Annot::Isotopes_db[isotope == "18O", mass_diff], 5) ==
      round(Spec2Annot::mz_from_string("[C6H2O+H-H2O]+_18O"), 5)
  )
  expect_true(
    round(Spec2Annot::mz_from_string("[C6H2O+H-H2O]+") + 2*Spec2Annot::Isotopes_db[isotope == "13C", mass_diff] + 1*Spec2Annot::Isotopes_db[isotope == "18O", mass_diff], 5) ==
      round(Spec2Annot::mz_from_string("[C6H2O+H-H2O]+_13C2_18O"), 5)
  )
})

test_that("Spec2Annot::mz_ppm", {
  require(magrittr)
  expect_true(Spec2Annot::mz_ppm(135.2235, 135.2238) == c(135.2235, 135.2238) %>% {(((max(.) - min(.)) / mean(.)) * 10^6)})
  expect_true(Spec2Annot::mz_ppm(654.3256, 654.5256) == c(654.3256, 654.5256) %>% {(((max(.) - min(.)) / mean(.)) * 10^6)})
})

test_that("Spec2Annot::mz_range", {
  temp <- Spec2Annot::mz_range(100, 10)
  expect_true(is.vector(temp))
  expect_true(is.numeric(temp))
  expect_true(length(temp) == 2)
  Spec2Annot::mz_ppm(temp[1], temp[2]) %between% (10 + c(-0.0001, 0.0001))

  temp <- Spec2Annot::mz_range(652.2253, 10)
  expect_true(is.vector(temp))
  expect_true(is.numeric(temp))
  expect_true(length(temp) == 2)
  Spec2Annot::mz_ppm(temp[1], temp[2]) %between% (10 + c(-0.0001, 0.0001))

  temp <- Spec2Annot::mz_range(1256.24532, 10)
  expect_true(is.vector(temp))
  expect_true(is.numeric(temp))
  expect_true(length(temp) == 2)
  Spec2Annot::mz_ppm(temp[1], temp[2]) %between% (10 + c(-0.0001, 0.0001))
})

