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
    output <- data.table(
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
    )

    output <- merge(
      output,
      Spec2Annot::Element[,
        .(
          element = paste0(mass_nb, atomic_symb),
          isotope = mass_nb,
          mass = atomic_mass
        )
      ],
      by = "element"
    ) %>% {
      .[, ID := seq_len(.N)][
        elmt_nb == "", elmt_nb := 1
      ][,
        elmt_nb := as.numeric(elmt_nb)
      ][]
    }

    if (nrow(output) <= 0) {
      return(FALSE)
    } else {
      return(output)
    }
  } else {
    return(FALSE)
  }
}

#' Return a string from an element table
#'
#' @param element_dt A data.table as returned by
#'                   `Spec2Annot::element_from_formula()`
#' @return
#' Return a string from an element table
#'
#' @import data.table
#'
#' @export
#'
#' @examples
#' "C6H12O3" %>%
#'   gen_formula_from_compo() %>%
#'   element_from_formula() %>%
#'   string_from_element()
string_from_element <- function(element_dt) {
  ## Add isotope elements
  element_dt[, iso_element := ifelse(
    is.na(isotope),
    element,
    gsub("[0-9]{1,2}([A-Z]{1}[a-z]{0,2})", "\\1", element)
  )]
  compoa <- element_dt[
    elmt_nb > 0, .(elmt_cnt = sum(elmt_nb)), by = iso_element
  ][
    order(iso_element),
    paste0(iso_element, elmt_cnt, collapse = "")
  ]

  if (element_dt[!is.na(isotope), .N] > 0) {
    compob <- element_dt[
      !is.na(isotope),
    ][
      order(iso_element),
    ][,
      .(text = paste0(isotope, iso_element, elmt_nb, collapse = "")),
      by = "element"
    ][,
      paste0(text, collapse = "_")
    ]
    compoa <- paste0(compoa, "_", compob)
  }
  return(compoa)
}

#' Calculate DBE from an elemental composition
#'
#' Calculation is performed as explained in
#' http://ms-textbook.com/chapter-6/answer-6-3/
#'
#' @param element_dt Either a data.table as returned by
#'                   `Spec2Annot::element_from_formula()`
#'                   or a string as "C18H12O3P".
#' @return
#' Return a string from an element table
#'
#' @import data.table magrittr
#'
#' @export
#'
#' @examples
#' get_dbe("C32H25N3O2S")
get_dbe <- function(element_dt) {
  if (
    !is.data.table(element_dt) &&
      length(element_dt) == 1 &&
      is.character(element_dt)
  ) {
    element_dt <- gen_formula_from_compo(element_dt) %>%
      element_from_formula()
  }
  ## Count isotope as element
  element_dt[
    ,
    elmt_dbe := gsub(
      "^([0-9]{0,3})([A-Z]{1}[a-z]{0,2})$",
      "\\2",
      element
    )
  ]

  element_dt <- element_dt[elmt_nb > 0, ]
  ## Get C nb
  elem_grp <- list(
    "C" = "C",
    "Y" = c("H", "F", "Cl", "Br", "I", "At"),
    "Z" = c("N", "P")
  )

  elem_res <- lapply(elem_grp, function(x) {
    temp_dt <- element_dt[elmt_dbe %in% x, ]
    if (nrow(temp_dt) > 0) {
      out <- temp_dt[, sum(elmt_nb)]
    } else {
      out <- 0
    }
    return(out)
  })

  dbe <- elem_res$C - elem_res$Y / 2 + elem_res$Z / 2 + 1
  return(dbe)
}

#' Check Nitrogen rule
#'
#' Calculation is performed as explained in
#' https://en.wikipedia.org/wiki/Nitrogen_rule
#'
#' @param n Nitrogen atom number
#' @param mass Ion mass
#' @return
#' A logical corresponding to `TRUE` if the N
#' rule is respected or `FALSE` if not.
#'
#' @export
#' @examples
#' get_nrule(5, 125.1235)
get_nrule <- function(n, mass) {
  meven <- ifelse(floor(mass) %% 2 == 0,
    TRUE,
    FALSE
  )
  neven <- ifelse(n %% 2 == 0,
    TRUE,
    FALSE
  )
  return(meven == neven)
}

#' Check Senior theorems
#'
#' Check the Senior theorems as described in Morikawa and
#' Newbold (2003) with the following:
#' i) the sum of valencies is an even number, or the total
#' number of atoms having odd valencies is even.
#' ii) the sum of valencies is greater than or equal to twice
#' the maximum valency.
#' iii) the sum of valencies is greater than or equal to twice
#' the number of atoms minus 1.
#'
#' @inheritParams get_dbe
#' @param global Logical to return a unique value if all the
#'               theorem are valid (`TRUE`) or a vector with
#'               the result of each theorem (`FALSE`).

#' @return
#' If global is set to `TRUE`, return a unique logical value
#' corresponding to `TRUE` if all the theorems pass or `FALSE`
#' if any of them isn't.
#' If global is set to `FALSE`, return a named logical vector with
#' the result of each theorem.
#'
#' @references
#' \enumerate{
#' \item Senior JK (1951) Partitions and Their Representative Graphs. Am J Math 73:663. doi: 10.2307/2372318
#' \item Morikawa T, Newbold BT (2003) Analogous odd-even parities in mathematics and chemistry. Chemistry 12:445–450
#'}
#'
#' @export
#' @examples
#' get_senior("C6H12O3")
get_senior <- function(element_dt, global = TRUE) {
  if (
    !is.data.table(element_dt) &&
      length(element_dt) == 1 &&
      is.character(element_dt)
  ) {
    element_dt <- gen_formula_from_compo(element_dt) %>%
      element_from_formula()
  }
  # Sum of valence minus twice number of atome minus one
  elem_sum <- element_dt[, .(elmt_nb = sum(elmt_nb)), by = elmt]
  elem_val <- merge(
    elem_sum[, .(elmt, elmt_nb)],
    Spec2Annot::valence_db[, .(elmt = atomic_symb, val = valence)],
    by = "elmt"
  )

  ## Store results of the 3 theorems
  output <- c(FALSE, FALSE, FALSE)
  names(output) <- c("Senior1", "Senior2", "Senior3")

  ## Theorem 1
  if (
    (elem_val[, sum(val * elmt_nb)] %% 2 == 0) ||
      (elem_val[val %% 2, sum(elmt_nb)] %% 2 == 0)
  ) {
    output[1] <- TRUE
  }

  ## Theorem 2
  if (
    (elem_val[, sum(val * elmt_nb)]) >=
      (2 * elem_val[, max(val)])
  ) {
    output[2] <- TRUE
  }

  ## Theorem 3
  if (
    (elem_val[, sum(val * elmt_nb)]) >=
      (2 * elem_val[, sum(elmt_nb)] - 1)
  ) {
    output[3] <- TRUE
  }

  ## Output
  if (isTRUE(global)) {
    return(all(output))
  } else {
    return(output)
  }
}

#' Return the mass of an electorn
#'
#' @return
#' Return the mass of an electron as
#' a numeric value
#'
#' @export
#'
#' @examples
#' electron_mass()
electron_mass <- function() {
  return(as.numeric(0.0005489))
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
#' # Exempl A
#' formula <- gen_formula_from_compo("C6H2O3Ca2K1")
#' element_from_formula(formula)
#'
#' # Exempl B
#' formula <- gen_formula_from_compo("[C6H12O2+H]+_13C2")
#' element_from_formula(formula)
#'
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
    rbind(., list("e-", electron_mass(), NA))

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
  temp[, elmt := gsub("^([0-9]{0,2})([A-Z]{1}[a-z]{0,2})$", "\\2", element)]
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
    return(mass + unique(Spec2Annot::db_monocharge[Formula == form, mz_query]))
  } else if (nrow(temp_dt) == 0) {
    warning("form not found")
    return(as.numeric(NA))
  }
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
