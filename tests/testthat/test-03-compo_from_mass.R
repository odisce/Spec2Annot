testthat::test_that(
  "Spec2Annot::find_compo_from_mass", {
    mapply(
      function(x, y) {
        ap_res <- find_compo_from_mass(mass_target = x, ppm = 1)
        query_dt <- Spec2Annot::element_from_formula(y) %>%
          {
            .[elmt_nb > 0, .(element, elmt_nb)]
          } %>%
          transpose(., make.names = 1)
        ## Does the real formula exist in the results ?
        testthat::expect_true(nrow(ap_res[query_dt, on = names(query_dt)]) == 1)
      },
      c(360.1936, 89.0476, 105.04259, 180.0633),
      c("C21H28O5", "C3H7NO2", "C3H7NO3", "C6H12O6")
    )
  }
)
