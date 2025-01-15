#' Calculate the ions database
#'
#' @param losses losses as data.table (see Spec2Annot::Losses_db)
#' @param charges charges as data.table (see Spec2Annot::db_monocharge)
#' @param adducts adducts as data.table (see Spec2Annot::Adduct_db)
#' @import data.table magrittr
#'
#' @return
#' A data.table containing the full list of possible ions with the formula to calculate them
#' 
#' @export
#'
#' @examples
#' fun_get_ion_list()
fun_get_ion_list <- function(
  losses = Spec2Annot::Losses_db,
  charges = Spec2Annot::db_monocharge,
  adducts = Spec2Annot::Adduct_db
) {
  ## In data
  if (is.character(losses) || (is.vector(losses) && all(is.character(losses)))) {
    dt_losses <- data.table(loss = gsub("(^[A-z]{1}).*$", "-\\1", losses)) %>%
      unique()
    dt_losses[, mz_query := mz_from_string(loss) %>% {ifelse(. < 0, ., . * -1)}, by = .(loss)]
    dt_losses[, encod := "unknown"]
  } else if (is.data.table(losses)) {
    dt_losses <- losses
  } else if (isFALSE(losses)) {
    dt_losses <- FALSE
  } else {
    stop("fun_get_ion_list(): losses argument with unexpected format")
  }

  dt_charges <- charges
  dt_neutral <- adducts
  
  dt_losses <- dt_losses[!duplicated(mz_query)]
  ## Generate list of adduct
  add_list <- CJ1(dt_charges[, .(Formula, charge_mass = mz_query, charge, lossL)], dt_neutral[, .(adduct, adduct_mass = mz_query)])
  add_list[, ID := 1:.N]
  add_list[, Attribution := paste0("[(M", adduct, ")", Formula, "]", ifelse(charge > 0, "+", "-")), by = ID]
  add_list[, Attribution_mass := paste0("M+", adduct_mass, "+", charge_mass)]
  ## Remove adduct with +H (which are not possible)
  add_list <- add_list[!grepl("\\+H$", Formula)]
  add_list[, Type := "Adduct"]
  ## Add single charge forms
  dt_single <- dt_charges[, .(Formula, charge_mass = mz_query, charge, lossL, Type = "Monocharge")]
  dt_single[, Attribution := paste0("[M", Formula, "]", ifelse(charge < 0, "-", "+"))]
  dt_single[, Attribution_mass := paste0("M+", charge_mass)]
  add_list <- list(dt_single,
                   add_list) %>%
    rbindlist(fill = TRUE)
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
  List_add_neut_loss <- CJ1(
    add_list[Type == "Monocharge" & lossL == TRUE, .(Formula, adduct, lossL, adduct_mass, charge, charge_mass)],
    dt_losses[, .(loss, loss_mass = mz_query)]
  )
  ## Annotation using convention :
  ### [(ion)-neutral loss]-
  List_add_neut_loss[, Attribution := paste0("[(M", ifelse(is.na(adduct), "", adduct), Formula, ")", loss, "]", ifelse(charge < 0, "-", "+"))]
  List_add_neut_loss[, ID := 1:.N]
  List_add_neut_loss[, exact_mass := paste0("M+", ifelse(is.na(adduct_mass), "", adduct_mass), "+", charge_mass, loss_mass), by = ID]
  List_add_neut_loss[, Type := "Loss"]

  ## M+H and M-H
  # List_mono <- data.table(Type = "Mono", Attribution = c("[M+H]+", "[M-H]-"), exact_mass = c(paste0("M+", mono_charge), paste0("M-", mono_charge)), charge = c(+1, -1))
  ## Combine
  output <- rbind(
    add_list[, .(Attribution, exact_mass = Attribution_mass, charge, Type, lossL)],
    List_add_neut_loss[, .(Attribution, exact_mass, charge, Type, lossL)]
  )
  output[, ID := 1:.N]

  ## Check duplicate
  # output[, relative_mass := gsub("M\\+", "", exact_mass) %>% {eval(parse(text = .))}, by = ID]
  # output[relative_mass %in% List_add_neut_loss[duplicated(relative_mass), unique(relative_mass)]][order(relative_mass)] %>% View()

  return(data.table(output))
}

#' Calculate the ions database
#'
#' @inheritParams fun_get_ion_list
#' @inherit mz_calc_ion
#' @import data.table magrittr
#'
#' @return
#' A data.table containing the full list of possible ions with the formula to calculate them
#' 
#' @export
#'
#' @examples
#' fun_get_ion_list()
fun_generate_ions_from_mz <- function(
  mass = 252.2534,
  losses = Spec2Annot::Losses_db,
  iso = Spec2Annot::Isotopes_db,
  charges = Spec2Annot::db_monocharge,
  adducts = Spec2Annot::Adduct_db
) {
  ions_full_list <- do.call(fun_get_ion_list, list(losses, charges, adducts))
  if (is.null(losses) || isFALSE(losses)) {
    ions_full_list <- ions_full_list[Type != "Loss", ]
  }
  if (is.null(adducts) || isFALSE(adducts)) {
    ions_full_list <- ions_full_list[Type != "Adduct", ]
  }
  ##Calculate masses
  ions_full_list[, ID := 1:.N]
  ions_full_list[, exact_mass := gsub("M", mass, exact_mass)]
  ions_full_list[
    ,
    mz_query := exact_mass %>%
      {eval(parse(text = .))},
    by = .(ID)
  ]
  ions_full_list[, mz_query := as.numeric(mz_query)]
  
  if (!(isFALSE(iso) || is.null(iso))) {
    isotopes_dt <- lapply(seq_len(nrow(ions_full_list)), function(x) {
      it_dt <- ions_full_list[x,]
      out <- gen_isotopes(
        mz = it_dt$mz_query,
        label = "X",
        db_iso = Spec2Annot::Isotopes_db
      )
      ## Create custom label
      new_label <- paste0(it_dt$Attribution, "_", out$isotope)
      out[, label := NULL]
      out[, Attribution := new_label]
      out[, Type := "Isotopes"]
      out[, lossL := FALSE]
      out[, isotope := NULL]
      out[, charge := it_dt$charge]
      data.table(out)
    }) %>%
      rbindlist()
    
    ions_full_list <- rbindlist(
      list(
        isotopes_dt,
        ions_full_list
      ),
      fill = TRUE
    )
  }

  return(ions_full_list[])
}
