# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Group m/Z from a vector based on tolerance
#'
#' This function returns a matrix with m/Z groups
#' based on a tolerance.
#'
#' @param xx A numeric and sorted vector of m/Z
#' @param tt Tolerance in absolute m/Z to group peaks
#' @export
#'
mz_vec_aggregate <- function(xx, tt) {
    .Call(`_Spec2Annot_mz_vec_aggregate`, xx, tt)
}
