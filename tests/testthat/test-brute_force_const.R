test_that("brute_force_const()", {
  formula <- "C2H4O1"
  mass <- mz_from_string(formula)

  mass_dt <- get_element_from_mass(mass) %>%
    gen_formula_from_compo() %>%
    element_from_formula()
  mass_dt[elmt_nb > 0, ][order(-mass), ]

  ## Test debug argument
  tp <- lapply(c(0,1,2,10), function(x) {
    temp <- brute_force_const(
      mass = mass,
      ppm = 0.1,
      mass_vc = mass_dt[elmt_nb > 0, ][order(-mass), mass],
      name_vc = mass_dt[elmt_nb > 0, ][order(-mass), element],
      maxiter_vc_ = mass_dt[elmt_nb > 0, ][order(-mass), elmt_nb],
      debugl = x
    )
    expect_true(nrow(temp) == 1)
    temp_form <- add_formula_to_annot(
      data.table::as.data.table(temp)
    )
    expect_true(temp_form[1, formula] == formula)
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
