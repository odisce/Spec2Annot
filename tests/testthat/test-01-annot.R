test_that("Spec2Annot::get_iso_from_annot", {
  temp <- Spec2Annot::get_iso_from_annot("C6H12_13C2")
  expect_true(data.table::is.data.table(temp))
  expect_true(temp[element == "13C", elmt_nb] == 2)
  expect_false(Spec2Annot::get_iso_from_annot("ZSDR"))

  temp <- Spec2Annot::get_iso_from_annot("C6H12_13C4_18O1")
  expect_true(data.table::is.data.table(temp))
  expect_true(temp[element == "13C", elmt_nb] == 4)
  expect_true(temp[element == "18O", elmt_nb] == 1)

  expect_false(Spec2Annot::get_iso_from_annot("C6H12"))
  expect_false(Spec2Annot::get_iso_from_annot("SDMLX"))
  expect_false(Spec2Annot::get_iso_from_annot("SD_MLX"))
})

test_that("Spec2Annot::get_charge_from_compo", {
  tp <- lapply(
    c(
      "C6H12_13C4_18O1",
      "C6H12_13C2",
      "C6H12O3",
      "C6H3O3",
      "C6H1O3+C6",
      "-(C6H2O)-(H2O)",
      "[C6H2O+H-H2O]"
    ),
    function(x) {
      temp <- get_charge_from_compo(x)
      # message(sprintf("compo: %s, charges: %i", x, temp))
      expect_true(is.numeric(temp))
      expect_true(temp == 0)
    }
  )

  for (i in c("+", "++", "+++", "-", "--", "---")) {
    tp <- lapply(
      c(
        "C6H12Z_13C4_18O1",
        "C6H12Z_13C2",
        "C6H12O3Z",
        "C6H3O3Z",
        "C6H1O3+C6Z",
        "-(C6H2O)-(H2O)Z",
        "[C6H2O+H-H2O]Z",
        "[C6H2O+H-H2O]Z"
      ),
      function(input) {
        x <- gsub("Z", i, input)
        temp <- get_charge_from_compo(x)
        expect_true(is.numeric(temp))
        expect_true(abs(temp) == nchar(i))
        if (grepl("\\-", i)) {
          expect_true(temp < 0)
        } else {
          expect_true(temp > 0)
        }
      }
    )
  }

  lapply(
    c(
      "C6H12_13C4_18O1+++",
      "[C6H12-H2O]_13C4_18O1+++",
      "C6H12_13C2---",
      "C6H++12--O3"
    ),
    function(x) {
      temp <- get_charge_from_compo(x)
      expect_true(is.numeric(temp))
      expect_true(temp == 0)
    }
  )

})

test_that("Spec2Annot::mz_calc_ion", {
  expect_true(
    Spec2Annot::mz_calc_ion(351.2564, form = "-H") %between%
      (350.2491 + c(-0.0001, +0.0001))
  )
  expect_true(
    Spec2Annot::mz_calc_ion(351.2564, form = "+H") %between%
      (352.2637 + c(-0.0001, +0.0001))
  )
  expect_true(
    Spec2Annot::mz_calc_ion(351.2564, form = "+NH4") %between%
      (369.2902 + c(-0.0001, +0.0001))
  )
  testthat::expect_warning(is.na(Spec2Annot::mz_calc_ion(351.2564, form = "+NHCMDL")))
})

test_that("Spec2Annot::gen_formula_from_compo", {
  temp_form <- Spec2Annot::gen_formula_from_compo("C6H12O3")
  expect_true(is.character(temp_form))
  expect_true(temp_form == "(6*C+12*H+3*O)")
})

test_that("Spec2Annot::element_from_formula", {
  temp_dt <- Spec2Annot::element_from_formula(
    Spec2Annot::gen_formula_from_compo("C6H12O3")
  )
  expect_true(data.table::is.data.table(temp_dt))
  expect_true(temp_dt[element == "C", elmt_nb] == 6)
  expect_true(temp_dt[element == "H", elmt_nb] == 12)
  expect_true(temp_dt[element == "O", elmt_nb] == 3)
  temp_dt <- Spec2Annot::element_from_formula("C6H12O3")
  expect_true(data.table::is.data.table(temp_dt))
  expect_true(temp_dt[element == "C", elmt_nb] == 6)
  expect_true(temp_dt[element == "H", elmt_nb] == 12)
  expect_true(temp_dt[element == "O", elmt_nb] == 3)
  temp_dt <- Spec2Annot::element_from_formula(
    Spec2Annot::gen_formula_from_compo("C6H12O3_13C3")
  )
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
  expect_error(Spec2Annot::element_from_formula("H12O3_13C1"))
})


test_that("Spec2Annot::mz_from_string", {
  require(data.table)
  expect_true(
    sprintf(
      "%0.6f",
      Spec2Annot::mz_from_string("-H-")
    ) == -1.007276
  )
  expect_true(
    sprintf(
      "%0.10f",
      Spec2Annot::mz_from_string("C6H2O3")
    ) == 122.0003939232
  )
  expect_true(
    Spec2Annot::mz_from_string("C6H3O3+") %between%
      (123.0077 + c(-0.0001, +0.0001))
  )
  expect_true(
    Spec2Annot::mz_from_string("C6H1O3+C6-O+H10") %between%
      (187.0759 + c(-0.0001, +0.0001))
  )
  expect_true(
    Spec2Annot::mz_from_string("C6H1O3+C6-O+H10+Ca2+K") %between%
      (305.9648 + c(-0.0001, +0.0001))
  )
  expect_true(
    Spec2Annot::mz_from_string("C6H1O3+C6-O+H10+Ca2+K+H+") %between%
      (306.9721 + c(-0.0001, +0.0001))
  )
  expect_true(
    Spec2Annot::mz_from_string("C6H1O3+C6-O+H10+Ca2+K-H-") %between%
      (304.9575 + c(-0.0001, +0.0001))
  )
  expect_true(
    Spec2Annot::mz_from_string("-(C6H2O)") ==
      - Spec2Annot::mz_from_string("C6H2O")
  )
  expect_true(
    Spec2Annot::mz_from_string("-(C6H2O)-(H2O)") ==
      - Spec2Annot::mz_from_string("C6H2O") -
        Spec2Annot::mz_from_string("H2O")
  )
  expect_true(
    Spec2Annot::mz_from_string("-(C6H2O)+(H2O)") ==
      Spec2Annot::mz_from_string("-C6H2O+H2O")
  )
  expect_true(
    Spec2Annot::mz_from_string("-(C6H2O)+(H2O)") ==
      Spec2Annot::mz_from_string("-(C6H2O)+H2O")
  )
  expect_true(
    Spec2Annot::mz_from_string("-(C6H2O)+(H2O)+") ==
      Spec2Annot::mz_from_string("-(C6H2O)+H2O+")
  )
  expect_true(
    Spec2Annot::mz_from_string("-(C6H2O)+(H2O)-") ==
      Spec2Annot::mz_from_string("-(C6H2O)+H2O-")
  )
  ## mz from Annotation
  expect_true(
    Spec2Annot::mz_from_string("[C6H2O+H-H2O]+") +
      Spec2Annot::Element[paste0(mass_nb, atomic_symb) == "13C", atomic_mass] -
      Spec2Annot::Element[paste0(mass_nb, atomic_symb) == "12C", atomic_mass] ==
      Spec2Annot::mz_from_string("[C6H2O+H-H2O]+_13C1")
  )
  expect_true(
    Spec2Annot::mz_from_string("[C6H2O+H-H2O]+") +
      2 * (Spec2Annot::Element[
        paste0(mass_nb, atomic_symb) == "13C", atomic_mass
      ] -
        Spec2Annot::Element[
          paste0(mass_nb, atomic_symb) == "12C", atomic_mass
        ]) ==
      Spec2Annot::mz_from_string("[C6H2O+H-H2O]+_13C2")
  )
  expect_true(
    round(
      Spec2Annot::mz_from_string("[C6H2O+H-H2O]") +
        2 * Spec2Annot::Isotopes_db[isotope == "13C", mass_diff] +
        1 * Spec2Annot::Isotopes_db[isotope == "18O", mass_diff], 5
    ) ==
      round(
        Spec2Annot::mz_from_string("[C6H2O+H-H2O]_13C2_18O"), 5
      )
  )
  expect_true(
    round(
      Spec2Annot::mz_from_string("[C6H2O+H-H2O]+") +
        1 * Spec2Annot::Isotopes_db[isotope == "18O", mass_diff], 5
    ) ==
      round(
        Spec2Annot::mz_from_string("[C6H2O+H-H2O]+_18O"), 5
      )
  )
  expect_true(
    round(
      Spec2Annot::mz_from_string("[C6H2O+H-H2O]+") +
        2 * Spec2Annot::Isotopes_db[isotope == "13C", mass_diff] +
        1 * Spec2Annot::Isotopes_db[isotope == "18O", mass_diff], 5
    ) ==
      round(
        Spec2Annot::mz_from_string("[C6H2O+H-H2O]+_13C2_18O"), 5
      )
  )
})

test_that("Spec2Annot::mz_ppm", {
  require(magrittr)
  expect_true(
    Spec2Annot::mz_ppm(135.2235, 135.2238) ==
      c(135.2235, 135.2238) %>%
        {
          (((max(.) - min(.)) / mean(.)) * 10^6)
        }
  )
  expect_true(
    Spec2Annot::mz_ppm(654.3256, 654.5256) ==
      c(654.3256, 654.5256) %>%
        {
          (((max(.) - min(.)) / mean(.)) * 10^6)
        }
  )
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

test_that("Spec2Annot::get_charge_from_compo", {
  tp <- mapply(
    function(x, y) {
      temp <- Spec2Annot::get_charge_from_compo(x)
      expect_true(temp == y)
    },
    c(
      "C6H2O3+++",
      "[M+H2+Na]+",
      "[M+H2+2(Na)]++",
      "[M+H2+2(HCOO)]--",
      "C6H12O3-"
    ),
    c(3, 1, 2, -2, -1)
  )
})

test_that("Spec2Annot::get_charge_from_compo", {
  tp <- mapply(
    function(x, y) {
      temp <- Spec2Annot::get_charge_from_compo(x)
      expect_true(temp == y)
    },
    c(
      "C6H2O3+3",
      "[M+H2+Na]+1",
      "[M+H2+2(Na)]+2",
      "[M+H2+2(HCOO)]-2",
      "C6H12O3-1"
    ),
    c(3, 1, 2, -2, -1)
  )
})