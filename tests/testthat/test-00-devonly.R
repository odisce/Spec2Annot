testthat::skip("For development only")

test_that("DEV_A", {
  ## annotate_mz
  # devtools::test()
  # devtools::document()
  # devtools::load_all()
  # renv::install(".", prompt = FALSE)
  # require(testthat)
  # covr::report()

  ## Generate list of ions to search
  ion_mass <- 125.2135
  pol_val <- "neg"
  121 - Spec2Annot::mz_calc_ion(121, form = "+H")
  
  ##  
})
