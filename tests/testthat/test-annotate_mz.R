testthat::test_that("annotat_mz()", {
  lapply(
    c(
      "C10H12N5O6P1",
      "C10H12O6P1"
    ),
    function(x) {
      spec_res <- annotate_mz(
        spectra_ms2,
        ppm = 3,
        polarity = 1,
        compo = x
      )
      testthat::expect_true(
        all(
          names(spec_res) %in%
            c("mz", "i", "irel", "formula", "ppm", "DBE", "nrule", "senior")
        )
      )
      testthat::expect_true(
        all(
          spec_res$mz %in% spec_res$mz
        )
      )

      testthat::expect_true(
        all(
          spec_res$i %in% spec_res$i
        )
      )

      testthat::expect_true(nrow(spec_res) == nrow(spectra_ms2))
    }
  )
})

testthat::test_that("annotat_mz(): no annot found", {
  spec_res <- annotate_mz(
    spectra_ms2[1:5,],
    ppm = 3,
    polarity = 1,
    compo = "C10H12N5O6P1"
  )
  expect_true(is.data.table(spec_res))
  expect_true(nrow(spec_res) == 5)
})

testthat::test_that("annotat_mz(): negative mode", {
  temp <- copy(Spec2Annot::spectra_ms2)
  ## Shift mz to negative mode
  mz_shift <- Spec2Annot::mz_calc_ion(0, form = "-H")*2
  temp[, mz+mz_shift]
  spec_res <- annotate_mz(
    temp,
    ppm = 5,
    polarity = 0,
    compo = "C10H12N5O6P1"
  )
  expect_true(is.data.table(spec_res))
  expect_true(nrow(spec_res) == nrow(temp))
  expect_true(all(spec_res$mz %in% temp$mz))
})

testthat::test_that("add_formula_to_annot(): no annot found", {
  expect_message(add_formula_to_annot(data.table()))
  expect_true(is.data.table(add_formula_to_annot(data.table())))
  expect_true(nrow(add_formula_to_annot(data.table())) == 0)
})
