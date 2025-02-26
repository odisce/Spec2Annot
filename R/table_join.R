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