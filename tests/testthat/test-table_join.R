
test_that(
  "DEV_A", {
    require(magrittr)
    require(data.table)

    exp_dt <- Spec2Annot::spectra_full
    exp_dt[, id := ID]
    exp_dt[, rt := rep(1:10, length.out = .N)]
    db_dt <- exp_dt[c(1,1,2,3,4,5,6,7,8,9,10), ][, .(id = paste0('DB', seq_len(.N)), rt, mz)]
    
    output <- search_db(db_dt = db_dt, exp_dt = exp_dt)
    testthat::expect_true('data.table' %in% class(output))
    testthat::expect_true(output[, .N] >= 11)

    output <- search_db(db_dt = db_dt, exp_dt = exp_dt, rttol = NULL, ppmtol = 100)
    testthat::expect_true('data.table' %in% class(output))
    testthat::expect_true(output[, .N] == 17)
    testthat::expect_true(all(output[, range(mz)] == db_dt[, range(mz)]))

    output <- search_db(db_dt = db_dt, exp_dt = exp_dt, rttol = 10, ppmtol = NULL)
    testthat::expect_true('data.table' %in% class(output))
    testthat::expect_true(output[, .N] == 27324)
    testthat::expect_true(all(output[, range(rt)] == db_dt[, range(rt)]))
  }
)