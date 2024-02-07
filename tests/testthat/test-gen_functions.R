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
})