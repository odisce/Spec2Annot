#' Search DB egaint EXP
#'
#' @param db_dt a data.table with id, mz, rt columns
#' @param exp_dt a data.table with id, mz, rt columns
#' @param ppmtol numeric value to set ppm search in m/Z dimension.
#'               Can be set to NULL to use only rt search.
#' @param rttol numeric value to set rt tolerance in rt dimension.
#'              Can be set to NULL to use only mz search.
#' @return
#' Return the exp_dt data.table entries matching the criteria with
#' 2 new columns: dbid (ids in db_dt) and expid (ids in exp_dt).
#' 
#' @import magrittr data.table
#' @export
#'
search_db <- function(db_dt, exp_dt, ppmtol = 5, rttol = 5) {
  db_in <- copy(db_dt)
  exp_in <- copy(exp_dt)
  setnames(exp_in, 'id', 'expid')
  setnames(db_in, 'id', 'dbid')

  if (is.null(ppmtol)) {
    ## no ppm conditions = create NA column in db
    db_in[, mz := as.numeric(NA)]
  }
  if (is.null(rttol)) {
    ## no ppm conditions = create NA column in db
    db_in[, rt := as.numeric(NA)]
  }
  
  ## Matching
  x <- exp_in[, .(expid, rt, mz)] %>% unique()
  y <- db_in[, {
    if (!is.na(mz)) {
      mzrange <- Spec2Annot::mz_range(mz, ppmtol)
    } else {
      mzrange <- rep(as.numeric(NA), 2)
    }
    .(
      rtmin = rt - rttol,
      rtmax = rt + rttol,
      mzmin = mzrange[1],
      mzmax = mzrange[2])
  }, by = .(dbid)]
  data.table::setkey(y, "rtmin", "rtmax", "mzmin", "mzmax")
  output <- y[, {
    x[
      mz %between% c(mzmin, mzmax) &
        rt %between% c(rtmin, rtmax),
      .(expid)
    ]
  }, by = .(dbid)]
  if (output[, .N] <= 0) {
    return(FALSE)
  } else {
    return(
      merge(
        exp_in,
        output,
        by = "expid",
        allow.cartesian = TRUE
      )
    )
  }
}

#' Search DB egaint EXP using cpp function
#'
#' @param db_dt a data.table with id, mz, rt columns
#' @param in_exp a data.table with id, mz, rt columns
#' @param ppmtol numeric value to set ppm search in m/Z dimension.
#' @param rttol numeric value to set rt tolerance in rt dimension.
#' @return
#' Return the exp_dt data.table entries matching the criteria with
#' 2 new columns: dbid (ids in db_dt) and expid (ids in exp_dt).
#' 
#' @import magrittr data.table
#' @export
#'
search_db_cpp <- function(in_db, in_exp, ppmtol = 5, rttol = 5) {
  db_dt <- copy(in_db)
  exp_dt <- copy(in_exp)
  dt_n <- lapply(list(db_dt, exp_dt), names)
  if (!is.null(ppmtol) && !all(sapply(dt_n, function(x) {"mz" %in% x}))) {
    stop("missing mz column in A and/or B")
  }
  if (!is.null(rttol) && !all(sapply(dt_n, function(x) {"rt" %in% x}))) {
    stop("missing rt column in A and/or B")
  }
  
  ## Cpp optimized search function
  db_dt[, dbid := seq_len(.N)]
  setkey(db_dt, "dbid")
  exp_dt[, expid := seq_len(.N)]
  setkey(exp_dt, "expid")
  res <- match_tables(db_dt, exp_dt, ppmtol = ppmtol, rttol = rttol) %>%
    as.data.table()
  old_names <- names(db_dt[, -c("dbid")])
  new_names <- paste0("db_", old_names)
  setnames(db_dt, old_names, new_names)
  old_names <- names(exp_dt[, -c("expid")])
  new_names <- paste0("exp_", old_names)
  setnames(exp_dt, old_names, new_names)
  output <- merge(
    res,
    db_dt,
    by = "dbid"
  ) %>%
    merge(
      .,
      exp_dt,
      by = "expid"
    )
  setkey(output, expid, dbid)
  return(output[])
}
