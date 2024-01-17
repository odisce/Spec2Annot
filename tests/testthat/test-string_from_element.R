test_that(
  "Spec2Annot::string_from_element",
  {
    "C6H12O3_18O1_13C2" %>%
      gen_formula_from_compo() %>%
      element_from_formula() %>%
      string_from_element() %>%
      testthat::expect_equal(., "C6H12O3_13C2_18O1")

    "C20H15N2" %>%
      gen_formula_from_compo() %>%
      element_from_formula() %>%
      string_from_element() %>%
      testthat::expect_equal(., "C20H15N2")

    "C20H15N2_14N2" %>%
      gen_formula_from_compo() %>%
      element_from_formula() %>%
      string_from_element() %>%
      testthat::expect_equal(., "C20H15N2_14N2")
  }
)
