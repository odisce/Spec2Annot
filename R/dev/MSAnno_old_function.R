#' Calculate the full database
#'
#' @param losses
#' @param iso
#' @param charges
#' @param neutral
#' @import RMassBank data.table magrittr
#'
#' @return
#' @export
#'
#' @examples
#' fun_get_ion_list()
fun_get_ion_list <- function(losses = MSAnno::Losses_db,
                             iso = MSAnno::Isotopes_db,
                             charges = MSAnno::Monocharge_db,
                             neutral = MSAnno::Adduct_db) {
  ## In data
  dt_losses <- losses
  dt_iso <- iso
  dt_charges <- charges
  dt_neutral <- neutral

  ## Extract loss compo
  dt_losses[, ID := 1:.N][, Compo := fun_get_atom_count(loss, "string"), by = ID]
  ## Add exact mass
  dt_losses[, ID := 1:.N][, exact_mass := RMassBank::findMz.formula(Compo, mode = "", ppm = 0)$mzCenter, by = ID]
  ## Remove duplicat
  dt_losses <- dt_losses[!duplicated(exact_mass)]
  ## Add charges mass
  dt_charges[, ID := 1:.N]
  elem_list <- dt_charges[, fun_get_atom_count(Formula, output = "data.table")]
  elem_list[, exact_mass := RMassBank::findMz.formula(element, mode = "", ppm = 0)$mzCenter, by = element]
  elem_list <- rbind(elem_list, data.table(element = "e-", count = "1", exact_mass = 5.489*10^-4))
  mono_charge <- elem_list[element == "H", exact_mass] - elem_list[element == "e-", exact_mass]
  ### Replace atom by masses and Add electron
  dt_charges[, exact_mass := gsub("(.)([0-9])", "\\2\\*\\1", Formula), by = ID]
  dt_charges[, exact_mass := gsub("([A-Z][a-z]{0,1})", "\\1\\+", exact_mass), by = ID]
  dt_charges[, exact_mass := gsub("\\+$", "", exact_mass)]
  dt_charges[, exact_mass := ifelse(charge < 0, paste0(exact_mass, "+", abs(charge), "*e-"), paste0(exact_mass, "-", abs(charge), "*e-"))]

  r <- lapply(1:nrow(elem_list), function(x) {
    char <- elem_list[x, element]
    repl <- elem_list[x, exact_mass]
    dt_charges[, exact_mass := gsub(char, repl, exact_mass), by = ID]
  })
  rm(r)
  dt_charges[, exact_mass := eval(parse(text = exact_mass)), by = ID]
  ### Neutral calculate mass
  dt_neutral[, ID := 1:.N]
  elem_list <- dt_neutral[, fun_get_atom_count(adduct, output = "data.table")]
  elem_list[, exact_mass := RMassBank::findMz.formula(element, mode = "", ppm = 0)$mzCenter, by = element]
  dt_neutral[, exact_mass := adduct, by = ID]
  dt_neutral[, exact_mass := gsub("([0-9]{1,2})", "\\1\\*", exact_mass)]

  r <- lapply(1:nrow(elem_list), function(x) {
    char <- elem_list[x, element]
    repl <- elem_list[x, exact_mass]
    dt_neutral[, exact_mass := gsub(char, repl, exact_mass), by = ID]
  })
  rm(r)
  dt_neutral[, exact_mass := eval(parse(text = exact_mass)), by = ID]

  ## Generate list of adduct
  add_list <- CJ1(dt_charges[, .(Formula, charge_mass = exact_mass, charge, lossL)], dt_neutral[, .(adduct, adduct_mass = exact_mass)])
  add_list[, ID := 1:.N]
  add_list[, Attribution := paste0("[(M", adduct, ")", Formula, "]", ifelse(charge > 0, "+", "-")), by = ID]
  add_list[, Attribution_mass := paste0("M+", adduct_mass, "+", charge_mass)]
  ## Remove adduct with +H (which are not possible)
  add_list <- add_list[!grepl("\\+H$", Formula)]
  add_list[, Type := "Adduct"]
  ## Add single charge forms
  dt_single <- dt_charges[, .(Formula, charge_mass = exact_mass, charge, lossL, Type = "Monocharge")]
  dt_single[, Attribution := paste0("[M", Formula, "]", ifelse(charge < 0, "-", "+"))]
  dt_single[, Attribution_mass := paste0("M+", charge_mass)]
  # openxlsx::write.xlsx(add_list, file = "./data/AnnoteV3/20210319_external_charge.xlsx", asTable = T, colWidths = "auto")
  add_list <- list(dt_single,
                   add_list) %>%
    rbindlist(fill = T)
  add_list[, ID := 1:.N]
  ### Filter adduct for combinaisons with the same mass
  #### Calculate relative mass of adduct
  add_list[, relative_mass := gsub("M\\+", "", Attribution_mass) %>% {eval(parse(text = .))}, by = ID]
  add_list[, relative_mass := sprintf("%.8f", relative_mass) %>% as.numeric(), by = ID]
  add_list[, ID := 1:.N]
  #### Check duplicated relative mass and keep the shortest adduct form
  dup_add <- add_list[relative_mass %in% add_list[duplicated(relative_mass), unique(relative_mass)]][order(relative_mass)]
  #### If multiple H loss, keep shortest expression
  keep_id <- dup_add[grepl("[2-9]", Attribution), .SD[which.min(nchar(adduct))], by = relative_mass][, ID]
  #### If single H loss, keep form with "K" as charge
  keep_id_2 <- dup_add[!grepl("[2-9]", Attribution), .SD[Formula == "+K"], by = relative_mass][, ID]
  rem_ID <- dup_add[!ID %in% c(keep_id, keep_id_2), ID]
  add_list <- add_list[!ID %in% rem_ID]

  ## Generate list of losses
  ### Neutral losses only for neutral form with T for lossL
  List_add_neut_loss <- CJ1(add_list[Type == "Monocharge" & lossL == "T", .(Formula, adduct, lossL, adduct_mass, charge, charge_mass)], dt_losses[, .(loss, loss_mass = exact_mass)])
  ## Annotation using convention :
  ### [(ion)-neutral loss]-
  List_add_neut_loss[, Attribution := paste0("[(M", ifelse(is.na(adduct), "", adduct), Formula, ")", loss, "]", ifelse(charge < 0, "-", "+"))]
  List_add_neut_loss[, ID := 1:.N]
  List_add_neut_loss[, exact_mass := paste0("M+", ifelse(is.na(adduct_mass), "", adduct_mass), "+", charge_mass, "-", loss_mass), by = ID]
  List_add_neut_loss[, Type := "Loss"]

  ## M+H and M-H
  # List_mono <- data.table(Type = "Mono", Attribution = c("[M+H]+", "[M-H]-"), exact_mass = c(paste0("M+", mono_charge), paste0("M-", mono_charge)), charge = c(+1, -1))
  ## Combine
  output <- rbind(add_list[, .(Attribution, exact_mass = Attribution_mass, charge, Type, lossL)],
                  List_add_neut_loss[, .(Attribution, exact_mass, charge, Type, lossL)])
  output[, ID := 1:.N]

  ## Check duplicate
  # output[, relative_mass := gsub("M\\+", "", exact_mass) %>% {eval(parse(text = .))}, by = ID]
  # output[relative_mass %in% List_add_neut_loss[duplicated(relative_mass), unique(relative_mass)]][order(relative_mass)] %>% View()

  return(data.table(output))
}

#' Return the atom count of a string
#'
#' @param string Chemical Formula as string
#' @param output either 'data.table' or 'string' for export format
#' @import magrittr data.table stringr
#'
#' @return
#' @export
#'
#' @examples
#' fun_get_atom_count("C6H5O2", output = "data.table")
fun_get_atom_count <- function(string, output = c("data.table", "string")[[1]]) {
  # string <- "-(H3PO4)-(CHNO)"
  elem <- stringr::str_extract_all(string, "[A-Z][a-z]?[0-9]*") %>%
    unlist()
  unique_elem <- gsub("[0-9]*", "", elem)
  count <- gsub("[A-z][a-z]*", "", elem) %>% as.numeric()
  list_elem <- data.table("element" = unique_elem, "count" = count)
  list_elem[is.na(count), count := 1]
  ## Sum all unique elements
  elem_calc <- list_elem[, .(count = sum(count)), by = element]
  if (output == "string") {
    elem_calc[, count_label := ifelse(count == 1, "", count)]
    return(elem_calc[order(element), paste0(element, count_label)] %>% paste(., collapse = ""))
  }
  if (output == "data.table") {
    return(elem_calc)
  }
}

#' Cross join of two data.table
#'
#' @param DT1 A data.table
#' @param DT2 A data.table
#' @import data.table
#'
#' @return
#' @export
#'
#' @examples
#' CJ1(data.table(A = 1:10, B = "A"), data.table(C = 50:100, D = rep(c("C", "D", "E", "F"), length.out = 51)))
CJ1 <- function(DT1, DT2){
  DT1[ , c(.SD, DT2), by = seq_len(nrow(DT1))]
}


#' Function to generate isotopes from an exact mass
#'
#' @param mz Exact mass
#'
#' @return
#' @export
#' @import data.table magrittr
#'
#' @examples
#' gen_isotopes(125.53658, "Isovalerine")
gen_isotopes <- function(mz, label = NULL, db_iso = MSAnno::Isotopes_db) {
  # mz <- 253.125465 ; label = NULL ; db_iso = Isotopes_db
  temp_db <- copy(db_iso)
  temp_db[, ID := 1:.N]
  if (is.null(label)) {
    label <- "X"
  }
  temp_db[, .(label = paste0("[", label, "]_", isotope), mz_query = mz + mass_diff), by = .(ID, isotope)]
}

#' Function to generate monocharges from an exact mass
#'
#' @param mz Exact mass of targeted compound
#' @param mz_type Type of mass for targeted compound: pos, neg or neutral
#' @param ion_mode Output mode wanted: pos, neg
#' @param label Label to add in output attribution
#'
#' @return A data.table with ID, label and mz_query
#' @export
#' @import data.table magrittr
#'
#' @examples
#' gen_monocharge(125.53658, "Isovalerine", "neutral", "neg")
#' gen_monocharge(116.0837, "C6H12O2", "neutral", "neg")
#' gen_monocharge(116.0837, "C6H12O2", "neutral", "pos")
gen_monocharge <- function(mz, label = NULL, mz_type = c("neutral", "pos", "neg")[[1]], ion_mode = c("pos", "neg")[[1]], mono_db = MSAnno::Monocharge_db) {
  # mz <- 125.53658 ; label = "Isovalerine" ; mz_type = "neutral" ; ion_mode <- "neg" ; mono_db = MSAnno::Monocharge_db
  if (!"data.table" %in% class(mono_db) | !all(c("Formula", "charge", "lossL", "mz_query") %in% names(mono_db))) {stop("mono_db must be a data.table with Formula, charge and lossL columns")}
  if (!is.numeric(mz)) {stop("mz must be numeric")}
  ## convert mz to neutral form
  mz <- switch(mz_type,
               pos = mz - 1.007276132,
               neg = mz + 1.007276132,
               neutral = mz,
               stop("Check mz_type argument"))
  ## Get monocharge list for selected polarity
  temp_db <- copy(mono_db)
  temp_db <- switch(ion_mode,
                    pos = temp_db[charge == 1],
                    neg = temp_db[charge == -1],
                    stop("ion_mode argument must be 'pos' or 'neg'"))
  ## Calculate m/Z from formula and add mz
  temp_db[, ID := 1:.N]
  ## Add charge
  temp_db[, Formula_calc := ifelse(charge == 1, paste0(Formula, "-", abs(charge)), ifelse(charge == -1, paste0(Formula, "+", abs(charge)), Formula))]
  temp_db[, mz_query := mz + mz_query, by = ID]
  ## Add label
  if (is.null(label)) {label <- "X"}
  ## Output
  temp_db[, label := paste0("[", label, Formula, "]", ifelse(charge < 0, "-", ifelse(charge > 0, "+", "")))]
  temp_db[, .(ID, label, lossL, mz_query)]
}

#' Calculate mz from a string with signs
#'
#' @param string
#'
#' @return
#' @export
#' @import data.table magrittr RMassBank
#' @importFrom stringr str_extract_all
#'
#' @examples
#' MSAnno_mz_from_string("C6H5O3+")
#' MSAnno_mz_from_string("[C6H4O3+H]+_13C1")
#' MSAnno_mz_from_string("[C6H4O3+H]+_13C2")
MSAnno_mz_from_string <- function(string) {
  ## Format string for calcul
  temp <- gen_formula_from_compo(string)
  ## Get all unique element and their mz
  elem <- element_from_formula(temp)
  ## replace each element in string, starting with longest
  ## If isotopes, replace _ by +
  # temp <- gsub("_", "\\+", temp)
  #
  # elem[order(-char_nb), {temp <<- gsub(element, mass, temp)}, by = element] %>%
  #   invisible()

  ## Calculate final m/z
  output <- elem[, mass * as.numeric(elmt_nb)] %>% sum()
  return(output)
}

#' Convert annotation to html string
#'
#' @param annotation
#' @param compo
#' @param compo_replace
#' @param string_format
#'
#' @return
#' @export
#'
#' @examples
#' MSAnno_html_from_annotation("[M+H]+_13C2", "C6H12O2", "M")
#' MSAnno_html_from_annotation("[M+H-H2O]+_13C2", "C6H12O2", "M")
#' MSAnno_html_from_annotation("[M+H-H2O]+_13C2_18O", "C6H12O2", "M")
#' MSAnno_html_from_annotation("C6H13O2-H2O+_13C2_18O", "", "M")
#' MSAnno_html_from_annotation("[M+H-(OH•)]+_13C2_18O", "C6H12O2", "M")
MSAnno_html_from_annotation <- function(annotation, compo = "", compo_replace = "X") {
  # annotation <- "[((X-2H+2Na)+H)-C3H7O2NS]+" ; compo = "C26H43NO6" ; compo_replace = "X"
  ## Replace by compo
  form <- gsub(compo_replace, paste0("(", compo, ")"), annotation)
  ## Get element list
  temp <- gen_formula_from_compo(form)
  elem <- element_from_formula(temp)
  ## Get element for sorting
  elem[, elmt := gsub("[0-9]{1,2}([A-Z]{1}[a-z]{0,2})", "\\1", element)]
  ## Create html string without e-
  elem[, ID := 1:.N]
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

#' Convert annotation or formula to latex format
#'
#' @param annotation
#' @param compo
#' @param compo_replace
#'
#' @return
#' @export
#'
#' @examples
#' MSAnno_latex_from_annotation("[M+H]+_13C", "C6H12O3", "M")
MSAnno_latex_from_annotation <- function(annotation, compo = "", compo_replace = "X") {
  # annotation <- "[2X+NH4]+_13C3"
  # compo <- "C27H32O14"
  # compo_replace = "X"
  ## Replace by compo
  form <- gsub(compo_replace, paste0("(", compo, ")"), annotation)
  ## Get element list
  temp <- gen_formula_from_compo(form)
  elem <- element_from_formula(temp)
  ## Get element for sorting
  elem[, elmt := gsub("[0-9]{1,2}([A-Z]{1}[a-z]{0,2})", "\\1", element)]
  ## Create html string without e-
  elem[, ID := 1:.N]
  elem_sub <- elem[elmt != "e-"][order(elmt, mass)]
  out_string <- elem_sub[, {
    output <- ""
    ## Add isotope upper
    if (!is.na(isotope)) {
      output <- paste0("^{", isotope, "}")
    }
    ## Add element
    output <- paste0(output, elmt)
    ## If element > 1 add nb sub
    if (elmt_nb > 1) {
      output <- paste0(output, "_{", elmt_nb, "}")
    }
    data.table("text" = output)
  }, by = ID][, paste0(text, collapse = "")]

  ## Add sign if presents
  if (nrow(elem[element == "e-"]) > 0) {
    out_string <- ifelse(
      elem[element == "e-", elmt_nb < 0],
      paste0(out_string, "{+}"),
      paste0(out_string, "^{-}")
    )
  }
  ## Add radical if present
  if (grepl("\u2022", annotation)) {
    out_string <- paste0(out_string, "^{U+2022}")
  }
  return(out_string)
}

#' Generate adduct list from mz
#'
#' @param mz
#' @param label
#' @param mz_type
#' @param adduct_db
#'
#' @return
#' @export
#'
#' @examples
#' gen_adduct(125.53658, "Isovalerine", "neutral")
#' gen_adduct(116.0837, "C6H12O2", "neutral")
#' gen_adduct(116.0837, "C6H12O2", "pos")
gen_adduct <- function(mz, label = NULL, mz_type = c("neutral", "pos", "neg")[[1]], adduct_db = MSAnno::Adduct_db) {
  # mz= 125.53658
  # label = "[X+H]+"
  # mz_type = "pos"
  # adduct_db = MSAnno::Adduct_db

  if (!"data.table" %in% class(adduct_db) | !all(c("adduct", "charge", "mz_query") %in% names(adduct_db))) {stop("adduct_db must be a data.table with Formula, charge and lossL columns")}
  if (!is.numeric(mz)) {stop("mz must be numeric")}
  ## convert mz to neutral form
  mz <- switch(mz_type,
               pos = mz - 1.007276132,
               neg = mz + 1.007276132,
               neutral = mz,
               stop("Check mz_type argument"))
  temp_db <- copy(adduct_db)
  ## Calculate m/Z from formula and add mz
  temp_db[, ID := 1:.N]
  ## Add charge
  temp_db[, mz_query := mz + mz_query]
  ## Add label
  if (is.null(label)) {label <- "X"}
  ## Output
  ### Perf optimization
  # temp_db[, label := paste0("[", label, adduct, "]", ifelse(charge < 0, "-", ifelse(charge > 0, "+", "")))]
  temp_db[, label := paste0("[", label, adduct, "]")]
  temp_db[charge < 0, label := paste0(label, "-")]
  temp_db[charge > 0, label := paste0(label, "+")]
  temp_db[, .(ID, label, mz_query)]
}

#' Generate losses list from mz
#'
#' @param mz
#' @param label
#' @param mz_type
#' @param loss_db
#'
#' @return
#' @export
#'
#' @examples
#' gen_losses(125.53658, "Isovalerine", "neutral")
#' gen_losses(116.0837, "C6H12O2", "neutral")
#' gen_losses(116.0837, "C6H12O2", "pos")
gen_losses <- function(mz, label = NULL, mz_type = c("neutral", "pos", "neg")[[1]], loss_db = MSAnno::Losses_db) {
  if (!"data.table" %in% class(loss_db) | !all(c("loss", "mz_query") %in% names(loss_db))) {stop("loss_db must be a data.table with Formula, charge and lossL columns")}
  if (!is.numeric(mz)) {stop("mz must be numeric")}
  ## convert mz to neutral form
  mz <- switch(mz_type,
               pos = mz - 1.007276132,
               neg = mz + 1.007276132,
               neutral = mz,
               stop("Check mz_type argument"))
  ## Calculate m/Z from formula and add mz
  temp_db <- copy(loss_db)
  temp_db[, ID := 1:.N]
  ## Add charge
  temp_db[, mz_query := mz + mz_query]
  ## Remove losse < 0
  temp_db <- temp_db[mz_query > 0]
  ## Add label
  if (is.null(label)) {label <- "X"}
  ## Output
  temp_db[, label := paste0("[", label, loss, "]")]
  temp_db[, .(ID, label, mz_query)]
}


#' Title
#'
#' @param mass
#' @param ppm
#'
#' @return
#' @export
#'
#' @examples
mz_range <- function(mass, ppm) {
  mass + (ppm * mass / 10^6 / 2) * c(-1, 1)
}

#' Title
#'
#' @param massa
#' @param massb
#'
#' @return
#' @export
#'
#' @examples
mz_ppm <- function(massa, massb) {
  data <- c(massa, massb)
  mass_error <- max(data) - min(data)
  return((mass_error / mean(data)) * 10^6)
}

#' Get sign from annotation
#'
#' @param string
#'
#' @return
#' @export
#'
#' @examples
#' get_sign(string = "[M+H]-")
get_sign <- function(string) {
  string %>% gsub(".*\\]([0-9]{0,2}(+|-))", "\\1", .)
}

#' Get all between brackets
#'
#' @param string
#'
#' @return
#' @export
#'
#' @examples
#' get_brackets("[M+H]+")
get_brackets <- function(string) {
  string %>% gsub("\\[(.*)(\\].*)", "\\1", .)
}

#' Get isotopes from annotation
#'
#' @param annotation
#'
#' @return
#' @export
#'
#' @examples
#' get_iso_from_annot("[M+H]+_13C_18O")
#' get_iso_from_annot("[M+H]+_13C2")
get_iso_from_annot <- function(annotation) {
  # annotation <- "[M+H]+_13C2_18O"
  if (grepl("_", annotation)) { # Check if isotopes are present by detecting "_"
    temp <- tstrsplit(annotation, "_") %>% unlist()
    annotation <- temp[1]
    iso <- temp[-1] %>%
      lapply(., function(x) {
        # x <- temp[-1][2]
        tp <- gsub("([0-9]{1,2})([A-Z]{1}[A-z]{0,2})([0-9]{0,2})", "\\1_\\2_\\3", x) %>%
          data.table::tstrsplit(., "_") %>%
          data.table::as.data.table()
        if (ncol(tp) == 2) {
          tp[, elmt_nb := 1]
        }
        tp[, text := x]
        setnames(tp, c("isotope", "element", "elmt_nb", "text"))
        return(tp)
      }) %>%
      rbindlist(fill = T)
    ## Replace NA by 1
    iso[, ID := 1:.N]
    output <- iso[, {
      query <- paste0(isotope, element)
      data.table(.SD, data.table(mass = MSAnno::Element[paste0(mass_nb, atomic_symb) == query, atomic_mass]))
    }, by = ID]

    output[, element := paste0(isotope, element)]
    return(output[, -'ID'])
  } else {
    return(F)
  }
}

#' Get isotopes from annotation
#'
#' @param annotation
#'
#' @return
#' @export
#'
#' @examples
#' get_iso_from_annot("[M+H]+_13C_18O")
#' get_iso_from_annot("[M+H]+_13C2")
get_iso_from_annotV2 <- function(annotation) {
  # annotation <- "[M+H]+_13C2_18O"
  # annotation <- "[M+H]+_18O3_13C2_44Ca4"
  if (grepl("_", annotation)) { # Check if isotopes are present by detecting "_"
    temp <- tstrsplit(annotation, "_") %>% unlist()
    annotation <- temp[1]

    ## For each iso, get Elmt and nb
    data.table("element" = gsub("([0-9]{1,2})([A-Z]{1}[a-z]{0,1})([0-9]{0,3})", "\\1\\2", temp[-1]),
               "text" = temp[-1],
               "elmt_nb" = gsub("([0-9]{1,2})([A-Z]{1}[a-z]{0,1})([0-9]{0,3})", "\\3", temp[-1])) %>%
      merge(., MSAnno::Element[, .(element = paste0(mass_nb, atomic_symb), isotope = mass_nb, mass = atomic_mass)], by = "element") %>%
      {.[, ID := 1:.N][elmt_nb == "", elmt_nb := 1][]}
  } else {
    return(F)
  }
}

#' Title
#'
#' @param compo
#'
#' @return
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
  # compo <- "[C6H12O2+H]+_13C2_18O"
  # compo <- "C6H8O10-2(H2O•)"
  # compo <- "[2X+NH4]+_13C3"
  iso_res <- get_iso_from_annotV2(compo)

  if (grepl("_", compo)) {
    ## If iso, remove them before extracting formula
    compo <- strsplit(compo, "_") %>% unlist() %>% {.[1]}
  }

  ## Searh multimeres
  output <- gsub("\u2022", "", compo) %>% ## Remove radical character
    gsub("\\[", "(", .) %>% ## replace brackets by parenthesis
    gsub("\\]", ")", .) %>% ## replace brackets by parenthesis
    paste0("(", ., ")") %>%
    gsub("(\\+|\\-)", "\\)\\1\\(", .) %>% ## Add parenthesis by compo groups
    gsub("(\\+|\\-)\\(\\)$", "\\1", .) %>% ## Remove empty ending parenthesis if charge
    gsub("^\\(\\)", "", .) %>% ## Remove empty starting parenthesis if starting sign
    gsub("([A-Z]{1}[a-z]{0,1})([0-9]{1,4})", "\\2\\*\\1", .) %>% ## transform operation ex: C2 to 2*C
    gsub("([0-9]{1,4})([A-Z]{1}[a-z]{0,1})", "\\1\\*\\2", .) %>% ## transform operation ex: 2C to 2*C
    gsub("([A-Z]{1}[a-z]{0,1})", "\\1\\+", .) %>% ## Extract 1 character elements
    gsub("([+|-])\\)", ")", .) %>% ## Resolve parenthesis signs
    gsub("([A-Z][a-z]{0,1})(\\+$)", "\\1", .) %>% {  ## Delete ending artefact + sign
      ifelse(
        grepl("(\\+)$", .), gsub("(\\+)$", "-e-", .), ## if string terminate by + = -e-
        ifelse(
          grepl("(\\-)$", .), gsub("(\\-)$", "+e-", .), ## if string terminate by - = +e-
          .)
      )
    } %>%
    gsub("([0-9])\\(", "\\1*(", .)

  ## If iso, add them
  if (!isFALSE(iso_res)) {
    output <- paste0(output, iso_res[, paste0("_", text, collapse = "")])
  }
  return(output)
}

#' Title
#'
#' @param formula
#'
#' @return
#' @export
#' @import data.table magrittr
#' @importFrom RMassBank findMz.formula
#' @importFrom stringr str_extract_all
#'
#' @examples
#' formula <- gen_formula_from_compo("C6H2O3Ca2K1")
#' element_from_formula(formula)
#' formula <- gen_formula_from_compo("[C6H12O2+H]+_13C2")
#' element_from_formula(formula)
element_from_formula <- function(formula) {
  # formula <- gen_formula_from_compo("[((C6H12-2H+2Na)+H)-C3H7O2NS]+_13C2")
  ## Isotopes
  iso_res <- get_iso_from_annotV2(formula)

  if (grepl("_", formula)) {
    formula <- strsplit(formula, "_") %>% unlist() %>% {.[1]}
  }

  temp <- stringr::str_extract_all(formula, "([A-Z][a-z]{0,1})") %>%
    unlist() %>%
    unique() %>%
    {MSAnno::Element[atomic_symb %in% .][, .SD[which.max(isotopic_compo)], by = atomic_symb][, .(element = atomic_symb, mass = atomic_mass, isotope = NA)]} %>%
    rbind(., list("e-", 0.0005489, NA))

  ## Add atom count
  temp[, elmt_nb := {
    # gsub(paste0(element, "(?![A-z])"), "1", formula, perl = T) %>% # Replace current element by 1
      gsub(paste0("(", element, ")([^a-z]|$)"), "1\\2", formula) %>%
      gsub("e\\-|([A-Z][a-z]{0,2})", "0", .) %>% # replace all other element by 0 %>%
      str2expression() %>%
      eval()
  }, by = .(element, mass)]

  ## Remove isotopes from elements
  if (!isFALSE(iso_res)) {
    iso_res[, ID := 1:.N]
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
      fill = T
    )
  }

  return(temp)
}

#' Title
#'
#' @param string
#' @param compo
#' @param compo_replace
#'
#' @return
#' @export
#'
#' @examples
fun_convert_annot_to_html <- function(string, compo, compo_replace) {
  string %>%
    strsplit(., "\\|") %>%
    unlist() %>%
    stringr::str_trim() %>%
    sapply(., function(x) {
      tryCatch({
        MSAnno_html_from_annotation(x, compo = compo, compo_replace = compo_replace)
      }, error = function(e) {message(e) ; return(" ")}
      )
    }) %>%
    paste0(., collapse = " | ")
}

#' Find elemental formula of ions in a spectra
#'
#' @param spectrum Spectrum matrix, data.table with the column: mz (named mz OR first column)
#' @param prec_compo Precursor composition in the form: [C5H12O3+H]+
#' @param ppm tolerance for the search in ppm
#' @import magrittr data.table
#' @importFrom Rdisop initializeElements
#' @importFrom MassTools calcMF
#'
#' @return
#' @export
#'
fun_spectra_elem_from_compo <- function(spectrum, prec_compo, ppm, polarity = c(0,1)[1]) {
  if ("data.table" %in% class(spectrum) && "mz" %in% names(spectrum)) {
    spectrum_dt <- spectrum
    spectrum_dt[, `:=`(ion_id, 1:.N)]
  } else {
    spectrum_dt <- as.data.table(spectrum)
    setnames(spectrum_dt, 1, "mz")
    spectrum_dt[, `:=`(ion_id, 1:.N)]
  }
  element_count <- prec_compo %>% MSAnno::gen_formula_from_compo() %>%
    MSAnno::element_from_formula() %>% {
      .[!grepl("^e", element), .(element, count = elmt_nb)]
    }
  precursor_elemental_compo <- paste(element_count[order(element), paste0(element, count)], collapse = "")
  temp_elem <- Rdisop::initializeElements(element_count[order(element), element])
  prec_mass <- MSAnno::MSAnno_mz_from_string(prec_compo)
  charge_val <- ifelse(polarity == 1, +1, -1)
  spectrum_dt[, c("calc_formula_nb", "calc_formula", "error_mz", "error_ppm", "RDBE") := {
    res <- .(as.integer(NA), as.character(NA), as.numeric(NA), as.numeric(NA), as.numeric(NA))
    temp_res <- MassTools::calcMF(
      mz = mz,
      ppm = 6,
      z = charge_val,
      elements = temp_elem,
      Filters = list(maxElements = precursor_elemental_compo)
    )
    if (!is.null(temp_res)) {
      res <- .(nrow(temp_res), temp_res[1, ]$MF, temp_res[1, ]$error, temp_res[1, ]$ppm, RMassBank::dbe(temp_res[1, ]$MF))
    }
    res
  }, by = ion_id]
  setnames(element_count, "count", "prec")
  spectrum_dt[!is.na(calc_formula), `:=`("loss_formula", {
    ion_atoms <- MSAnno::fun_get_atom_count(calc_formula)
    final_count <- merge(element_count, ion_atoms, by = "element",
                         all = T)
    final_count[is.na(count), `:=`(count, 0)]
    final_count[, `:=`(loss_count, prec - count)]
    final_loss <- final_count[!is.na(loss_count) & loss_count >
                                0][, paste0(element, loss_count)] %>% paste(., collapse = "")
    if (final_loss == "") {
      .(as.character(NA))
    }
    else {
      .(paste0("-", final_loss))
    }
  }), by = ion_id]
  spectrum_dt[, `:=`(loss_formula, gsub("([A-z])([1]{1}(?![0-9]))",
                                        "\\1", loss_formula, perl = T))]
  spectrum_dt[!is.na(calc_formula), `:=`(calc_formula, paste0(calc_formula,
                                                              ifelse(polarity == 0, "-", "+")))]
  return(spectrum_dt[])
}


