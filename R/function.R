#' Get isotopes from annotation
#'
#' @param annotation Annotation in the form of `[M+H]+_13C1`
#'
#' @return
#' Return a data.table with the following isotopes informations:
#'   - `element`: element type (13C, 18O, ...)
#'   - `text`: label to use for the element
#'   - `elmt_nb`: Element count
#'   - `isotope`: Isotope number (13, 18, ...)
#'   - `mass`: Mass of the isotope
#'   - `ID`: Unique identifier (`integer`)
#'
#' @export
#'
#' @examples
#' get_iso_from_annot("[M+H]+_13C_18O")
#' get_iso_from_annot("[M+H]+_13C2")
get_iso_from_annot <- function(annotation) {
  if (grepl("_", annotation)) {
    # Check if isotopes are present by detecting "_"
    temp <- tstrsplit(annotation, "_") %>% unlist()
    annotation <- temp[1]

    ## For each iso, get Elmt and nb
    data.table(
      "element" = gsub(
        "([0-9]{1,2})([A-Z]{1}[a-z]{0,1})([0-9]{0,3})",
        "\\1\\2",
        temp[-1]
      ),
      "text" = temp[-1],
      "elmt_nb" = gsub(
        "([0-9]{1,2})([A-Z]{1}[a-z]{0,1})([0-9]{0,3})",
        "\\3",
        temp[-1]
      )
    ) %>%
      merge(
        .,
        Spec2Annot::Element[,
          .(
            element = paste0(mass_nb, atomic_symb),
            isotope = mass_nb,
            mass = atomic_mass
          )
        ],
        by = "element"
      ) %>%
      {
        .[, ID := seq_len(.N)][elmt_nb == "", elmt_nb := 1][]
      }
  } else {
    return(FALSE)
  }
}

#' Get element count from formula
#'
#' @param formula Formula as returned by `gen_formula_from_compo`
#'
#' @import data.table magrittr
#' @importFrom stringr str_extract_all
#'
#' @return
#' Return a data.table with element count
#'
#' @export
#'
#' @examples
#' formula <- gen_formula_from_compo("C6H2O3Ca2K1")
#' element_from_formula(formula)
#' formula <- gen_formula_from_compo("[C6H12O2+H]+_13C2")
#' element_from_formula(formula)
element_from_formula <- function(formula) {
  if (!grepl("\\*", formula)) {
    formula <- gen_formula_from_compo(formula)
  }
  ## Isotopes
  iso_res <- get_iso_from_annot(formula)
  if (grepl("_", formula)) {
    formula <- strsplit(formula, "_") %>%
      unlist() %>%
      {
        .[1]
      }
  }

  temp <- stringr::str_extract_all(formula, "([A-Z][a-z]{0,1})") %>%
    unlist() %>%
    unique() %>%
    {
      Spec2Annot::Element[
        atomic_symb %in% .
      ][,
        .SD[which.max(isotopic_compo)],
        by = atomic_symb
      ][,
        .(element = atomic_symb, mass = atomic_mass, isotope = NA)
      ]
    } %>%
    rbind(., list("e-", 0.0005489, NA))

  ## Add atom count
  temp[, elmt_nb := {
    gsub(paste0("(", element, ")([^a-z]|$)"), "1\\2", formula) %>%
      gsub("e\\-|([A-Z][a-z]{0,2})", "0", .) %>%
      str2expression() %>%
      eval()
  }, by = .(element, mass)]

  ## Remove isotopes from elements
  if (!isFALSE(iso_res)) {
    iso_res[, ID := seq_len(.N)]
    iso_res[, {
      elem <- gsub("[0-9]{1,2}([A-Z]{1}[a-z]{0,2})", "\\1", element)
      nb_to_rem <- as.numeric(elmt_nb)
      in_tp <- temp[element == elem]
      if (nrow(in_tp) == 0) {
        stop("Isotopes present when main atom is not ?")
      } else {
        temp[element == elem, elmt_nb := elmt_nb - nb_to_rem]
      }
    }, by = ID]
    temp <- rbindlist(
      list(
        temp,
        iso_res
      ),
      fill = TRUE
    )
  }

  return(temp[])
}

#' Generate formula from composition
#'
#' @param compo Elemental composition as a string (ex.: "C6H12O3NH2")
#'
#' @return Return a string with an arithmetic formula
#' used to calculate the final mass. This function can
#' be used to check if the string parser behave correctly.
#'
#' @export
#' @import magrittr
#'
#' @examples
#' gen_formula_from_compo("C6H2O3NH4")
#' gen_formula_from_compo("-(C6H2O)-(H2O)")
#' gen_formula_from_compo("-C6H2O-H2O+Ca2-")
#' gen_formula_from_compo("[C6H12O2+H-H2O]+_13C")
#' gen_formula_from_compo("[2(C6H2O3)+NH4]+_13C3")
gen_formula_from_compo <- function(compo) {
  ## Add isotopes
  ### Search for isotopes
  iso_res <- get_iso_from_annot(compo)

  if (grepl("_", compo)) {
    ## If iso, remove them before extracting formula
    compo <- strsplit(compo, "_") %>%
      unlist() %>%
      {
        .[1]
      }
  }

  ## Searh multimeres
  output <- gsub("\u2022", "", compo) %>% ## Remove radical character
    gsub("\\[", "(", .) %>% ## replace brackets by parenthesis
    gsub("\\]", ")", .) %>% ## replace brackets by parenthesis
    paste0("(", ., ")") %>%
    gsub("(\\+|\\-)", "\\)\\1\\(", .) %>%
    gsub("(\\+|\\-)\\(\\)$", "\\1", .) %>%
    gsub("^\\(\\)", "", .) %>%
    gsub("([A-Z]{1}[a-z]{0,1})([0-9]{1,4})", "\\2\\*\\1", .) %>%
    gsub("([0-9]{1,4})([A-Z]{1}[a-z]{0,1})", "\\1\\*\\2", .) %>%
    gsub("([A-Z]{1}[a-z]{0,1})", "\\1\\+", .) %>%
    gsub("([+|-])\\)", ")", .) %>%
    gsub("([A-Z][a-z]{0,1})(\\+$)", "\\1", .) %>%
    {
      ifelse(
        grepl("(\\+)$", .), gsub("(\\+)$", "-e-", .),
        ifelse(
          grepl("(\\-)$", .), gsub("(\\-)$", "+e-", .),
          .
        )
      )
    } %>%
    gsub("([0-9])\\(", "\\1*(", .)

  ## If iso, add them
  if (!isFALSE(iso_res)) {
    output <- paste0(output, iso_res[, paste0("_", text, collapse = "")])
  }
  return(output)
}

#' Calculate mz from a string with signs
#'
#' @param string Formula as a string of the form "C6H5O3+" or with isotopes
#'               `[C6H5O3+H]+_13C1`
#'
#' @return
#' Return a numeric value corresponding to the mass of the `string` input
#'
#' @export
#' @import data.table magrittr
#' @importFrom stringr str_extract_all
#'
#' @examples
#' mz_from_string("C6H5O3+")
#' mz_from_string("[C6H4O3+H]+_13C1")
#' mz_from_string("[C6H4O3+H]+_13C2")
mz_from_string <- function(string) {
  ## Format string for calcul
  temp <- gen_formula_from_compo(string)
  ## Get all unique element and their mz
  elem <- element_from_formula(temp)
  ## Calculate final m/z
  output <- elem[, mass * as.numeric(elmt_nb)] %>% sum()
  return(output)
}

#' Calculate ion mass using form
#'
#' @param mass numerical mass to use as base
#' @param form chemical form to add or substract (on
#'             of `Spec2Annot::db_monocharge[, unique(Formula)]`)
#'
#' @return
#' A numeric value corresponding to the mass form.
#'
#' @export
#' @import data.table
#'
#' @examples
#' mz_calc_ion(142.5236, "-H")
mz_calc_ion <- function(mass, form = "-H") {
  temp_dt <- Spec2Annot::db_monocharge[Formula == form, ]
  if (nrow(temp_dt) == 1) {
    return(mass + Spec2Annot::db_monocharge[Formula == form, mz_query])
  } else if (nrow(temp_dt) == 0) {
    warning("form not found")
    return(as.numeric(NA))
  } else if (nrow(temp_dt) > 1) {
    warning(
      paste0(
        "Multiple form available, choose one: ",
        paste0(temp_dt$Formula, collapse = ", ")
      )
    )
    return(as.numeric(NA))
  }
}

#' Calculate ppm deviation between two masses
#'
#' @param massa First mass as `numeric()`
#' @param massb Second mass as `numeric()`
#'
#' @return
#' A numeric value corresponding to the ppm deviation between massa and massb
#'
#' @export
#'
#' @examples
#' mz_ppm(142.5236, 142.5241)
mz_ppm <- function(massa, massb) {
  data <- c(massa, massb)
  mass_error <- max(data) - min(data)
  return((mass_error / mean(data)) * 10^6)
}

#' Calculate mass range using a ppm deviation
#'
#' @param mass mass to use as center value
#' @param ppm ppm deviation to get the mass range
#'
#' @return
#' A numeric vector with the mass range corresponding to mass +- ppm/2.
#' The range difference is equal to the asked ppm
#'
#' @export
#'
#' @examples
#' mz_range(142.5236, 5)
mz_range <- function(mass, ppm) {
  mass + (ppm * mass / 10^6 / 2) * c(-1, 1)
}
