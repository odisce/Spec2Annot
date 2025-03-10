% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gen_ionlist.R
\name{gen_ionlist}
\alias{gen_ionlist}
\title{Calculate list of targetd ions from neutral mz}
\usage{
gen_ionlist(
  neutral_mz = 153.5125,
  polarity = c(0, 1)[1],
  iso = c(TRUE, FALSE)[1],
  multi = 0,
  losses = c(TRUE, FALSE)[1],
  adducts = c(TRUE, FALSE)[1],
  mono_db = Spec2Annot::db_monocharge,
  loss_db = Spec2Annot::Losses_db,
  adduct_db = Spec2Annot::Adduct_db,
  db_iso = Spec2Annot::Isotopes_db
)
}
\arguments{
\item{adducts}{}
}
\description{
Calculate list of targetd ions from neutral mz
}
\examples{
gen_ionlist(153.5125, 1, T, 2, T, T)
}
