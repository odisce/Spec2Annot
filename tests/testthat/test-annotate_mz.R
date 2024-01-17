testthat::test_that("annotat_mz()", {
  spec_res <- annotate_mz(
    spectra_ms2,
    ppm = 3,
    polarity = 1,
    compo = "C10H12N5O6P1"
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

})
