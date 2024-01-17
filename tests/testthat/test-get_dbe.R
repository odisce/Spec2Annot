testthat::test_that("Spec2Annot::get_dbe", {
  temp <- mapply(
    function(x, y) {
      testthat::expect_equal(
        get_dbe(x),
        y
      )
    },
    c(
      "C3H8",
      "C6H12O2",
      "C6H11NO",
      "C32H25N3O2S",
      "C8H20N2",
      "C7H9N",
      "C3H6BrN",
      "C9H11ClN4O2",
      "C17H32Cl2",
      "C17H32Cl2_13C7",
      "C17H32Cl2_13C2",
      "C17H32Cl2_13C2_37Cl1"
    ),
    c(0, 1, 2, 22, 0, 4, 1, 6, 1, 1, 1, 1)
  )
})
