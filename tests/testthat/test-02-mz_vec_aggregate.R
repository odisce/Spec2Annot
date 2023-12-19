test_that("Spec2Annot::mz_vec_aggregate", {
  data(spectra_full)
  spectrab <- spectra_full$mz %>% {
    . + sample(
      seq(
        from = -0.01,
        to = 0.01,
        by = 0.00002
      ),
      size = length(.),
      replace = TRUE
    )
  }
  merged_spec <- Spec2Annot::mz_vec_aggregate(
    sort(
      c(
        spectra_full$mz,
        sort(spectrab)
      )
    ),
    0.015
  )
  expect_true(is.numeric(merged_spec))
  expect_true(is.vector(merged_spec))
  expect_true(
    length(merged_spec) ==
      length(spectrab) + length(spectra_full$mz)
  )
  expect_true(length(unique(merged_spec)) < length(spectrab))
})
