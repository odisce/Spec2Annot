#' Add formula to annotate_mz()
#'
#'
#' @param mzannot_dt A `data.table` containing the elemental composition of
#'                   ions as returned by `annotate_mz()`.
#' @import data.table magrittr
#' @return
#'   Return the mzannot_dt table with a new `formula` column
#' containing a string with the elemental composition.
#'
#' @export
add_formula_to_annot <- function(mzannot_dt) {
  annot_dt <- copy(mzannot_dt)

  if (nrow(annot_dt) <= 0) {
    message("0 annotation, returning empty table")
    mzannot_dt
    return(mzannot_dt)
  }
  annot_dt[, iter_annot := seq_len(.N)]

  ## Get column corresponding to elements
  names_vc <- gsub("^([0-9]{0,3})([A-Z]{1}[a-z]{0,2})$", "\\2", names(annot_dt))
  names(names_vc) <- names(annot_dt)
  col_elem <- names_vc[
    names_vc %in% Spec2Annot::Element[, unique(atomic_symb)]
  ] %>%
    names()

  annot_dt[, "formula" := {
    temp <- transpose(.SD, keep.names = "element")[V1 > 0, ]
    setnames(temp, "V1", "elmt_nb")
    temp[,
      isotope := ifelse(
        grepl("^[0-9]{1,2}[A-Z]{1}.*$", element),
        gsub("^([0-9]{1,2})([A-Z]{1}.*)$", "\\1", element),
        as.logical(NA)
      )
    ]
    string_from_element(temp)
  }, by = iter_annot, .SDcols = col_elem]
  annot_dt[, iter_annot := FALSE]
  return(annot_dt[])
}


#' Annotate a mass spectrum
#'
#' This function annotate a full mass spectrum using a full brute force
#' search optionnaly restricted by an elemental composition `compo`.
#' Since generating the space to research can take some times, the
#' parameter `search_space` can be pre-calculated in advance.
#'
#' @param input_spectrum A mass spectrum as a data.table with `mz` and
#'                       `i ` columns
#' @param ppm mass tolerance for the search in ppm.
#' @param polarity Polarity is used to convert the mass spectrum to neutral mass
#'                 before running the search.
#' @inheritParams find_compo_from_mass
#' @param compo (optional) Elemental composition as a string (`"C6H12O3"`) of
#'              the neutral form.
#' @import data.table magrittr
#' @return
#'   Return a data.table with the annotated spectrum.
#' @export
#' @examples
#' annotate_mz(spectra_ms2, ppm = 3, polarity = 1, compo = "C10H12N5O6P1")
annotate_mz <- function(
  input_spectrum,
  ppm = 5,
  polarity = NULL,
  compo = NULL,
  use_golden_ratio = TRUE
) {
  ## If ion mass, use polarity to add/remove electron and proton
  mz_dt <- copy(input_spectrum)
  mz_dt[, ion_iter := seq_len(.N)]
  mz_dt[, irel := i / max(i, na.rm = TRUE)]
  mz_dt[, mze := ifelse(
    is.null(polarity),
    mz,
    ifelse(
      polarity == 0,
      mz - electron_mass(),
      mz + electron_mass()
    )
  ), by = ion_iter]

  ## Add H to neutral formula for the search
  compo_ion <- compo
  compo_elem_dt <- gen_formula_from_compo(compo) %>%
    element_from_formula()
  if (!is.null(polarity) && polarity == 0) {
    compo_ion <- compo_elem_dt[element == "H", elmt_nb := elmt_nb - 1][] %>%
      string_from_element()
  } else if (!is.null(polarity) && polarity == 1) {
    compo_ion <- compo_elem_dt[element == "H", elmt_nb := elmt_nb + 1][] %>%
      string_from_element()
  }

  ## Annotate
  mz_dt_annot <- mz_dt[
    ,
    {
      find_compo_from_mass(
        mass_target = mze,
        ppm = ppm,
        use_golden_ratio = use_golden_ratio,
        elements_vc = compo_ion,
        debugl = 0
      )
    },
    by = .(ion_iter, mz, i, irel)
  ]

  if (nrow(mz_dt_annot) <= 0) {
    output <- mz_dt_annot
    output[, formula := character()]
    output[, DBE := numeric()]
    output[, nrule := logical()]
    output[, senior := logical()]
  } else {
    mz_dt_annot <- add_formula_to_annot(mz_dt_annot)
    ## Scores
    ##   DBE
    get_dbe_vec <- Vectorize(get_dbe, "element_dt")
    mz_dt_annot[, DBE := get_dbe_vec(formula)]

    ## Parity [ check with Annelaure]

    ## nrule
    if ("N" %in% names(mz_dt_annot)) {
      get_nrule_vec <- Vectorize(get_nrule, c("n", "mass"))
      mz_dt_annot[, nrule := get_nrule_vec(N, mass)]
    } else {
      mz_dt_annot[, nrule := as.logical(NA)]
    }

    ## Senior
    get_senior_vec <- Vectorize(get_senior, c("element_dt"))
    mz_dt_annot[, senior := get_senior_vec(element_dt = formula, global = TRUE)]

    ## Assemble data
    col_to_concatenate <- c("formula", "ppm", "DBE", "nrule", "senior")
    output <- mz_dt_annot[order(-nrule, -senior),
      lapply(.SD, function(x) {
        if (is.numeric(x)) {
          x <- sprintf("%.2f", x)
        }
        paste0(x, collapse = "/")
      }),
      by = .(ion_iter, mz, i, irel),
      .SDcols = col_to_concatenate
    ]
  }
  output_bind <- rbindlist(
    list(
      mz_dt[
        !ion_iter %in% output[, ion_iter],
        .SD,
        .SDcols = names(mz_dt) %>% {
          .[. %in% names(output)]
        }
      ],
      output
    ),
    fill = TRUE
  )
  output_bind[, ion_iter := NULL]
  return(
    output_bind[order(mz), ][]
  )
}
