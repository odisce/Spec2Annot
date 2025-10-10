#' Calculate the ions database
#'
#' @param losses losses as data.table (see Spec2Annot::Losses_db) or NULL
#' @param charges charges as data.table (see Spec2Annot::db_monocharge)
#' @param adducts adducts as data.table (see Spec2Annot::Adduct_db) or NULL
#' @param polarity Polarity to subset the forms (1 for positive, 0 for negative mode
#' or NULL for both)
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
  adducts = Spec2Annot::Adduct_db,
  polarity = NULL
) {
  ## First check charges (mandatory input)
  if (is.null(charges) || isFALSE(charges)) {
    stop("Cannot calculate ion form without any charges, please provide the charges argument.")
  } else {
    dt_charges <- charges[, .(Formula, charge_mass = mz_query, charge, lossL)]
    if (is.null(polarity)) {
      charge_sel <- c(-1, +1)
    } else if (polarity == 0) {
      charge_sel <- -1
    } else if (polarity == 1) {
      charge_sel <- +1
    } else {
      stop("polarity should be 1 (positive), 0 (negative) or NULL (default)")
    }
    dt_charges <- dt_charges[charge %in% charge_sel, ]
    ## Add single charges
    dt_single <- dt_charges[, .(Formula, charge_mass, charge, lossL, Type = "Monocharge")]
    dt_single[, Attribution := paste0("[M", Formula, "]", ifelse(charge < 0, "-", "+"))]
    dt_single[, Attribution_mass := paste0("M+", charge_mass)]
    add_list <- dt_single
  }

  if (is.null(adducts) || isFALSE(adducts)) {
    add_list[, adduct_mass := as.numeric(NA)]
  } else {
    ## If adduct present, add combination of adduct and charges
    dt_neutral <- adducts[, .(adduct, adduct_mass = mz_query)]
    ## Remove adduct with H (which are not possible)
    adduct_list <- CJ1(dt_charges[!Formula %in% c("-H", "+H"), .(Formula, charge_mass, charge, lossL, Type = "Adduct")], dt_neutral)
    adduct_list[, Attribution := paste0("[(M", adduct, ")", Formula, "]", ifelse(charge > 0, "+", "-"))]
    adduct_list[, Attribution_mass := paste0("M+", adduct_mass, "+", charge_mass)]
    ## Add to the list
    add_list <- rbindlist(
      list(add_list, adduct_list),
      fill = TRUE,
      use.names = TRUE
    )
  }
  if (is.null(losses) || isFALSE(losses)) {
    add_list[, loss_mass := as.numeric(NA)]
  } else {
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
    dt_losses <- dt_losses[!duplicated(mz_query), ]
    ## Generate losses ions
    ## Generate list of losses
    ### Neutral losses only for neutral form with T for lossL and monocharges
    charge_sub <- add_list[Type == "Monocharge" & lossL == TRUE, ]
    if (charge_sub[, .N] <= 0) {
      warning("No charges elligible to calculate losses")
    } else {
      List_add_neut_loss <- CJ1(
        charge_sub[, -c("Type")],
        dt_losses[, .(loss, loss_mass = mz_query, Type = "Loss")]
      )
      ## Annotation using convention :
      ### [(ion)-neutral loss]-
      # List_add_neut_loss[, Attribution := paste0("[(M", ifelse(is.na(adduct), "", adduct), Formula, ")", loss, "]", ifelse(charge < 0, "-", "+"))]
      # List_add_neut_loss[, ID := 1:.N]
      ## Add to the list
      add_list <- rbindlist(
        list(add_list, List_add_neut_loss),
        fill = TRUE,
        use.names = TRUE
      )
    }
  }

  if (add_list[, .N] > 0) {
    add_list[, ID := seq_len(.N)]
    ### Filter adduct for combinaisons with the same mass
    #### Calculate relative mass of adduct
    add_list[, relative_mass := gsub("M\\+", "", Attribution_mass) %>% {eval(parse(text = .))}, by = ID]
    if (add_list[Type == "Adduct", .N] > 0) {
      # add_list[, relative_mass := sprintf("%.8f", relative_mass) %>% as.numeric(), by = ID]
      #### Check duplicated relative mass and keep the shortest adduct form
      dup_add <- add_list[Type == "Adduct" & relative_mass %in% add_list[duplicated(relative_mass), unique(relative_mass)]][order(relative_mass)]
      if (dup_add[, .N] > 0) {
        #### If multiple H loss, keep shortest expression
        keep_id <- dup_add[grepl("[2-9]", Attribution), .SD[which.min(nchar(adduct))], by = relative_mass][, ID]
        #### If single H loss, keep form with "K" as charge
        keep_id_2 <- dup_add[!grepl("[2-9]", Attribution), .SD[Formula == "+K"], by = relative_mass][, ID]
        if (length(c(keep_id, keep_id_2)) > 0) {
          rem_ID <- dup_add[!ID %in% c(keep_id, keep_id_2), ID]
          add_list <- add_list[!ID %in% rem_ID]
        }
      }
    }
  } else {
    return(NULL)
  }
  
  add_list[, exact_mass := paste0(
    "M",
    ifelse(is.na(adduct_mass), "", paste0("+", adduct_mass)),
    ifelse(is.na(charge_mass), "", paste0("+", charge_mass)),
    ifelse(is.na(loss_mass), "", paste0("+", loss_mass))
  ), by = ID]

  return(data.table(add_list))
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
  adducts = Spec2Annot::Adduct_db,
  polarity = NULL
) {
  ions_full_list <- do.call(fun_get_ion_list, list(losses, charges, adducts, polarity))
  if (ions_full_list[, .N] <= 0) {
    stop("No ion form found")
  }
  ##Calculate masses
  ions_full_list[, ID := seq_len(.N)]
  ions_full_list[,
    relative_mass := gsub("M", 0, exact_mass) %>%
      {eval(parse(text = .))},
    by = .(ID)
  ]
  ions_full_list[, exact_mass := gsub("M", mass, exact_mass)]
  ions_full_list[,
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
