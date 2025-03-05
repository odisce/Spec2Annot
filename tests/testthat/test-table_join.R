
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

test_that(
  "Rcpp function match_tables()", {
    require(magrittr)
    require(data.table)

    exp_dt <- Spec2Annot::spectra_full
    exp_dt[, id := ID]
    exp_dt[, rt := rep(1:10, length.out = .N) %>% as.numeric()]
    db_dt <- exp_dt[c(1,1,2,3,4,5,6,7,8,9,10), ][, .(id = paste0('DB', seq_len(.N)), rt, mz)]
    
    output <- search_db_cpp(db_dt, exp_dt, ppmtol = 10, rttol = 10)
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

test_that(
  "Rcpp function _()", {
    require(magrittr)
    require(data.table)
    exp_dt <- Spec2Annot::spectra_full
    exp_dt[, id := ID]
    exp_dt[, rt := rep(1:10, length.out = .N) %>% as.numeric()]
    exp_dt <- exp_dt[sample(seq_len(.N), size = 100)]
    db_dt <- exp_dt[sample(seq_len(.N), 10),
    ][
      , .(id = paste0("DB", seq_len(.N)), rt, mz)
    ]

    ## Unsort data to test
    for (db_size in c(5, 10, 100, 200)) {
      for (ppmi in c(1, 2, 3, 5)) {
        for (rttoli in c(1, 2, 3)) {
          db_dt_i <- db_dt[sample(seq_len(.N), size = db_size, replace = TRUE), ]
          exp_dt_i <- exp_dt[sample(seq_len(.N)), ]
          output <- search_db_cpp(
            in_db = db_dt_i,
            in_exp = exp_dt_i,
            ppmtol = ppmi,
            rttol = rttoli
          )
          testthat::expect_true("data.table" %in% class(output))
          testthat::expect_true(output[, .N] == db_dt_i[, .N])
          testthat::expect_true(
            all(
              output[, .(
                ppm_check = mz_ppm(db_mz, exp_mz) <= ppmi
              ), by = .(dbid, expid)][, ppm_check]
            )
          )
          testthat::expect_true(
            all(
              output[, .(
                rt_check = abs(db_rt - exp_rt) <= rttoli
              ), by = .(dbid, expid)][, rt_check]
            )
          )
        }
      }
    }
  }
)
