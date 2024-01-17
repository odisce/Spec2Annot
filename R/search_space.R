#' Get Element from a mass
#'
#' @param mass input mass
#' @param use_golden_ratio Logical to return the maximum
#'                       number of each element by using
#'                       the ratios from the seven golden
#'                       rules (`TRUE`) or not (`FALSE`).
#' @import magrittr data.table
#' @return
#' Return a string with the number maximum number
#' of each element to expect from the mass.
#' get_element_from_mass(125.2535, FALSE)
#' get_element_from_mass(125.2535, TRUE)
#' @export
#'
get_element_from_mass <- function(mass, use_golden_ratio = FALSE) {
  ## Find the maximum number of Carbon, range of n:n - 10
  if (isFALSE(use_golden_ratio)) {
    element_vec <- c(
      "C",
      "H",
      "F",
      "Cl",
      "Br",
      "O",
      "P",
      "S",
      "Si",
      "Na",
      "K",
      "Cl",
      "N"
    )
    elmt_cnt <- lapply(
      element_vec,
      function(x) {
        elmt_nb <- ceiling(
          mass / Spec2Annot::Element[
            grepl(paste0("^", x, "$"), atomic_symb),
          ][
            which.max(isotopic_compo), atomic_mass
          ]
        )
        paste0(x, elmt_nb)
      }
    )
    output <- paste0(elmt_cnt, collapse = "")
    return(output)
  } else {
    c_nb <- (mass / 12) %>% ceiling()
    h_nb <- c_nb * 3.1 %>% ceiling()
    f_nb <- c_nb * 1.5 %>% ceiling()
    cl_nb <- c_nb * 0.8 %>% ceiling()
    br_nb <- c_nb * 0.8 %>% ceiling()
    n_nb <- c_nb * 1.3 %>% ceiling()
    o_nb <- c_nb * 1.2 %>% ceiling()
    p_nb <- c_nb * 0.3 %>% ceiling()
    s_nb <- c_nb * 0.8 %>% ceiling()
    si_nb <- c_nb * 0.5 %>% ceiling()
  }
  compo <- paste0(
    "C", c_nb,
    "H", h_nb,
    "F", f_nb,
    "Cl", cl_nb,
    "Br", br_nb,
    "N", n_nb,
    "O", o_nb,
    "P", p_nb,
    "S", s_nb,
    "Si", si_nb
  )
  return(compo)
}

#' Find composition from mass
#'
#' @param mass_target mass to decompose
#' @param ppm mass error to filter the results (in ppm)
#' @param elements_vc (optional) character vector containing
#'                    the element to include
#' @inheritParams get_element_from_mass
#' @inheritParams brute_force_const
#' @import magrittr data.table
#' @return
#'   A data.table with one proposition by line
#' with it's elemental composition, the theoretical
#' mass and the mass deviation with the query in
#' ppm.
#' @export
#' @examples
#'   find_compo_from_mass(125.215, ppm = 10)
#'   find_compo_from_mass(528.125, ppm = 3)
#'   find_compo_from_mass(89.0476, ppm = 3, elements_vc = "C3H7NO2")
find_compo_from_mass <- function(
  mass_target,
  ppm = 5,
  use_golden_ratio = TRUE,
  elements_vc = NULL,
  debugl = 0
) {
  element_dt_max <- get_element_from_mass(
    mass_target,
    use_golden_ratio = use_golden_ratio
  ) %>%
    gen_formula_from_compo() %>%
    element_from_formula() %>%
    {
      .[elmt_nb > 0, ]
    }

  if (!is.null(elements_vc)) {
    if (length(elements_vc) > 1) {
      ## Filter element but do not limit numbers
      element_dt_max <- element_dt_max[element %in% elements_vc, ]
    } else {
      ## Filter element and limit numbers
      element_dt_max <- elements_vc %>%
        gen_formula_from_compo() %>%
        element_from_formula() %>%
        {
          .[elmt_nb > 0, .(element, mass, elmt_nb)]
        }
    }
  }
  element_dt_max <- element_dt_max[order(-mass)]

  res <- brute_force_const(
    mass = mass_target,
    ppm = ppm,
    mass_vc = element_dt_max[, mass],
    name_vc = element_dt_max[, element],
    maxiter_vc_ = element_dt_max[, elmt_nb],
    debugl = debugl,
    debugit = 0
  ) %>%
    as.data.table(.)

  return(res[order(ppm)])
}
