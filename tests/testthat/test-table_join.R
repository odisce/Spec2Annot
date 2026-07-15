
require(magrittr)
require(data.table)

test_that(
  "search_db()", {
    exp_dt <- spectra_full
    exp_dt[, id := ID]
    exp_dt[, rt := rep(1:10, length.out = .N)]
    db_dt <- exp_dt[c(1,1,2,3,4,5,6,7,8,9,10), ][, .(id = paste0('DB', seq_len(.N)), rt, mz)]
    for (rttol in c(0, 1.5, 10)) {
      for (mztol in c(0, 3, 10, 100)) {
        output <- search_db(db_dt = db_dt, exp_dt = exp_dt, rttol = rttol, ppmtol = mztol)
        testthat::expect_true('data.table' %in% class(output))
        testthat::expect_true(all(output[, expid] %in% exp_dt[, ID]))
        testthat::expect_true(all(output[, dbid] %in% db_dt[, id]))
        testthat::expect_true(all(output[, range(mz)] %between% exp_dt[, range(mz)]))
        output_complete <- merge(
          output,
          db_dt[, .(dbid = id, db_rt = rt, db_mz = mz)],
          by = "dbid"
        )
        output_complete[, ppm := mz_ppm(mz, db_mz), by = .(mz, db_mz)]
        output_complete[, rt_diff := abs(rt - db_rt), by = .(rt, db_rt)]
        testthat::expect_true(all(output_complete[, ppm <= mztol]))
        testthat::expect_true(all(output_complete[, rt_diff <= rttol]))
      }
    }
  }
)

test_that(
  "Rcpp function search_db_cpp()", {
    exp_dt <- spectra_full
    exp_dt[, id := ID]
    exp_dt[, rt := rep(1:10, length.out = .N)]
    db_dt <- exp_dt[c(1,1,2,3,4,5,6,7,8,9,10), ][, .(id = paste0('DB', seq_len(.N)), rt, mz)]
    for (rttol in c(0, 1.5, 10)) {
      for (mztol in c(0, 3, 10, 100)) {
        output <- search_db_cpp(in_db = db_dt, in_exp = exp_dt, rttol = rttol, ppmtol = mztol)
        testthat::expect_true('data.table' %in% class(output))
        testthat::expect_true(all(output[, expid] %in% exp_dt[, ID]))
        testthat::expect_true(all(output[, db_id] %in% db_dt[, id]))
        testthat::expect_true(all(output[, range(exp_mz)] %between% exp_dt[, range(mz)]))
        output[, ppm := mz_ppm(exp_mz, db_mz), by = .(exp_mz, db_mz)]
        output[, rt_diff := abs(exp_rt - db_rt), by = .(db_rt, db_rt)]
        testthat::expect_true(all(output[, ppm <= mztol]))
        testthat::expect_true(all(output[, rt_diff <= rttol]))
      }
    }
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
            in_db = db_dt_i[, .(exp_id = id, rt, mz)],
            in_exp = exp_dt_i[, .(db_id = id, rt, mz)],
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
