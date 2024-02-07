#' Convert annotation to html string
#'
#' @param annotation String corresponding to the annotation
#' @param compo Elemental composition to replace to in the annotation string
#' @param compo_replace String to replace with `compo` in the annotation
#'
#' @return
#' An html string of the annotation.
#'
#' @export
#'
#' @examples
#' annot_to_html("[M+H]+_13C2", "C6H12O2", "M")
#' annot_to_html("[M+H-H2O]+_13C2", "C6H12O2", "M")
#' annot_to_html("[M+H-H2O]+_13C2_18O", "C6H12O2", "M")
#' annot_to_html("C6H13O2-H2O+_13C2_18O", "", "M")
#' annot_to_html("[M+H-(OH•)]+_13C2_18O", "C6H12O2", "M")
annot_to_html <- function(
  annotation,
  compo = "",
  compo_replace = "X"
) {
  ## Replace by compo
  form <- gsub(compo_replace, paste0("(", compo, ")"), annotation)
  ## Get element list
  temp <- gen_formula_from_compo(form)
  elem <- element_from_formula(temp)
  ## Get element for sorting
  elem[, elmt := gsub("[0-9]{1,2}([A-Z]{1}[a-z]{0,2})", "\\1", element)]
  ## Create html string without e-
  elem[, ID := seq_len(.N)]
  elem_sub <- elem[elmt != "e-" & elmt_nb > 0][order(elmt, mass)]
  out_string <- elem_sub[, {
    output <- ""
    ## Add isotope upper
    if (!is.na(isotope)) {
      output <- paste0("<sup>", isotope, "</sup>")
    }
    ## Add element
    output <- paste0(output, elmt)
    ## If element > 1 add nb sub
    if (elmt_nb > 1) {
      output <- paste0(output, "<sub>", elmt_nb, "</sub>")
    }
    data.table("text" = output)
  }, by = ID][, paste0(text, collapse = "")]

  ## Add sign if presents
  if (nrow(elem[element == "e-"]) > 0) {
    out_string <- ifelse(
      elem[element == "e-", elmt_nb < 0],
      paste0(out_string, "<sup>+</sup>"),
      paste0(out_string, "<sup>-</sup>")
    )
  }
  ## Add radical
  if (grepl("\u2022", annotation)) {
    out_string <- paste0(out_string, "<sup>&#x2022</sup>")
  }
  return(out_string)
}
