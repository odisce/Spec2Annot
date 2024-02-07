test_that("gen_isotopes()", {
  tp <- mapply(
    function(x, y) {
      temp_db <- gen_isotopes(x, y)
      expect_true(is.data.table(temp_db))
      expect_true(all(c("ID", "isotope", "label", "mz_query") %in% names(temp_db)))
      expect_true(all(temp_db[, grepl(y, label)]))
    },
    c(100, 150, 250.1235, 10.1235),
    c("C6H2O3", "XXX", "Isomlqkf", "fdecd-lkj")
  )
})

test_that("gen_monocharge()", {
  tp <- mapply(
    function(x, y) {
      temp_db <- gen_monocharge(x, y)
      expect_true(is.data.table(temp_db))
      expect_true(all(c("ID", "label", "lossL", "mz_query") %in% names(temp_db)))
      expect_true(all(temp_db[, grepl(y, label)]))
    },
    c(100, 150, 250.1235, 10.1235),
    c("C6H2O3", "XXX", "Isomlqkf", "fdecd-lkj")
  )
  
  testthat::expect_error(gen_monocharge(100, "XX", "pos", "pos", "c"))
  testthat::expect_error(gen_monocharge(100, "XX", "qsdf"))
  testthat::expect_error(gen_monocharge(100, "XX", "pos", "mlkjqs"))
  testthat::expect_error(gen_monocharge("mlk"))
  testthat::expect_true(
    all(
      grepl(
        "X",
        gen_monocharge(100, NULL, "pos", "pos")$label
      )
    )
  )
})

test_that("gen_adduct()", {
  tp <- mapply(
    function(x, y) {
      temp_db <- gen_adduct(x, y)
      expect_true(is.data.table(temp_db))
      expect_true(all(c("ID", "label", "mz_query") %in% names(temp_db)))
      expect_true(all(temp_db[, grepl(y, label)]))
    },
    c(100, 150, 250.1235, 10.1235),
    c("C6H2O3", "XXX", "Isomlqkf", "fdecd-lkj")
  )
  testthat::expect_error(gen_adduct(125, adduct_db = "c"))
  testthat::expect_error(gen_adduct("x"))
  testthat::expect_true(
    all(
      grepl(
        "X",
        gen_adduct(100, NULL)$label
      )
    )
  )
})

test_that("gen_losses()", {
  tp <- mapply(
    function(x, y) {
      temp_db <- gen_losses(x, y)
      expect_true(is.data.table(temp_db))
      expect_true(all(c("ID", "label", "mz_query") %in% names(temp_db)))
      expect_true(all(temp_db[, grepl(y, label)]))
    },
    c(100, 150, 250.1235, 10.1235),
    c("C6H2O3", "XXX", "Isomlqkf", "fdecd-lkj")
  )
  testthat::expect_error(gen_losses("c"))
  testthat::expect_error(gen_losses(100, loss_db = "c"))
  testthat::expect_true(
    all(
      grepl(
        "X",
        gen_losses(100, NULL)$label
      )
    )
  )
})