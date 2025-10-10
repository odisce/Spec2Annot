#' Calculate list of targetd ions from neutral mz
#'
#' @param neutral_mz Neutral Mass
#' @param polarity Polarity
#' @param iso Logical
#' @param multi Multimers
#' @param losses Logical
#' @param adducts Logical
#' @param mono_db Database
#' @param loss_db Database
#' @param adduct_db Database
#' @param db_iso Database
#'
#' @return a data.table
#' @export
#'
#' @examples
#' gen_ionlist(153.5125, 1, TRUE, 2, TRUE, TRUE)
gen_ionlist <- function(
  neutral_mz = 153.5125,
  polarity = c(0,1)[1],
  iso = c(TRUE, FALSE)[1],
  multi = 0,
  losses = c(TRUE, FALSE)[1],
  adducts = c(TRUE, FALSE)[1],
  mono_db = Spec2Annot::db_monocharge,
  loss_db = Spec2Annot::Losses_db,
  adduct_db = Spec2Annot::Adduct_db,
  db_iso = Spec2Annot::Isotopes_db
) {
  ## Check
  if (!is.numeric(c(neutral_mz, multi))) {stop("neutral_mz and multi must be numeric")}
  if (!is.logical(c(iso, losses, adducts))) {stop("iso, multi, losses and adducts must be logicals")}
  if (!polarity %in% c(0,1)) {stop("polarity must be 0 (NEG) or 1 (POS)")}

  ## Calculate monocharges
  mono_dt <- gen_monocharge(mz = neutral_mz, mz_type = "neutral", ion_mode = ifelse(polarity == 1, "pos", "neg"), mono_db = mono_db)
  mono_dt[, type := "Charge"]
  ion_list <- mono_dt
  ## Calculate multimeres
  if (multi > 0) {
    multi_dt <- lapply((1:multi)+1, function(x) {
      out <- gen_monocharge(mz = neutral_mz * x, label = paste0(x, "X"), mz_type = "neutral", ion_mode = ifelse(polarity == 1, "pos", "neg"), mono_db = mono_db)
      out[, lossL := FALSE]
      data.table(out)
    }) %>%
      rbindlist()
    multi_dt[, type := "multimer"]
    ion_list <- rbind(ion_list, multi_dt)
  }

  ## Calculate losses
  if (isTRUE(losses)) {
    ion_list_sub <- ion_list[lossL == TRUE,]
    loss_dt <- lapply(1:nrow(ion_list_sub), function(x) {
      # x <- 1
      dt_iter <- ion_list_sub[x,]
      out <- gen_losses(
        mz = dt_iter$mz_query,
        label = "X",
        mz_type = "neutral",
        loss_db = loss_db
      )
      bet_bracket <- dt_iter$label %>% gsub("\\[(.*)(\\].*)", "\\1", .)
      final_sign <- dt_iter$label %>% gsub(".*\\](.*)", "\\1", .)
      out[, label := gsub("X", bet_bracket, label) %>% paste0(., final_sign)]
      out[, type := "Loss"]
      out[, lossL := FALSE]
      data.table(out)
    }) %>%
      rbindlist()
    ion_list <- rbind(ion_list, loss_dt)
  }
  ## Calculate adducts only for monocharge
  if (isTRUE(adducts)) {
    add_dt <- lapply(seq_len(mono_dt[, .N]), function(x) {
      it_dt <- mono_dt[x, ]
      out <- gen_adduct(
        mz = it_dt$mz_query,
        label = "X",
        mz_type = "neutral",
        adduct_db = adduct_db
      )
      ## Create custom label
      bet_bracket <- out$label %>%
        gsub("\\[(.*)(\\].*)", "\\1", .) %>%
        paste0("(", ., ")")
      new_label <- sapply(bet_bracket, function(x) {
        gsub("X", x, it_dt$label)
      })
      out[, label := new_label]
      out[, type := "Adduct"]
      out[, lossL := FALSE]
    }) %>%
      rbindlist()
    ion_list <- rbind(ion_list, add_dt)
  }
  ## Calculate isotopes
  if (isTRUE(db_iso)) {
    isotopes_dt <- lapply(
      seq_len(ion_list[, .N]),
      function(x) {
        it_dt <- ion_list[x,]
        out <- gen_isotopes(mz = it_dt$mz_query, label = "X", db_iso = db_iso)
        ## Create custom label
        new_label <- paste0(it_dt$label, "_", out$isotope)
        out[, label := new_label]
        out[, type := "Isotopes"]
        out[, lossL := FALSE]
        out[, isotope := NULL]
        data.table(out)
      }
    ) %>%
      rbindlist()
    ion_list <- rbind(ion_list, isotopes_dt)
  }
  ## OUTPUT
  return(data.table(ion_list))
}
