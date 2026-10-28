test_that("brute_force_const()", {
  in_formula <- "C2H4O1"
  mass <- mz_from_string(in_formula)

  mass_dt <- get_element_from_mass(mass) %>%
    gen_formula_from_compo() %>%
    element_from_formula()
  mass_dt[elmt_nb > 0, ][order(-mass), ]

  ## Test debug argument
  tp <- lapply(c(3, 5, 10, 50), function(x) {
    temp <- brute_force_const(
      mass = mass,
      ppm = x,
      mass_vc = mass_dt[elmt_nb > 0, ][order(-mass), mass],
      name_vc = mass_dt[elmt_nb > 0, ][order(-mass), element],
      maxiter_vc_ = mass_dt[elmt_nb > 0, ][order(-mass), elmt_nb],
      debugl = FALSE
    )
    expect_true(nrow(temp) >= 1)
    temp_form <- add_formula_to_annot(
      data.table::as.data.table(temp)
    )
    expect_true(temp_form$formula[1] == in_formula)
  })

  ## Test iter
  temp <- brute_force_const(
      mass = mass,
      ppm = 0.1,
      mass_vc = mass_dt[elmt_nb > 0, ][order(-mass), mass],
      name_vc = mass_dt[elmt_nb > 0, ][order(-mass), element],
      maxiter_vc_ = mass_dt[elmt_nb > 0, ][order(-mass), elmt_nb],
      debugl = 0,
      debugit = 1
  )
  expect_true(nrow(temp) == 0)
})
