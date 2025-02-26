% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/function.R
\name{gen_formula_from_compo}
\alias{gen_formula_from_compo}
\title{Generate formula from composition}
\usage{
gen_formula_from_compo(compo)
}
\arguments{
\item{compo}{Elemental composition as a string (ex.: "C6H12O3NH2")}
}
\value{
Return a string with an arithmetic formula
used to calculate the final mass. This function can
be used to check if the string parser behave correctly.
}
\description{
Generate formula from composition
}
\examples{
gen_formula_from_compo("C6H2O3NH4")
gen_formula_from_compo("-(C6H2O)-(H2O)")
gen_formula_from_compo("-C6H2O-H2O+Ca2-")
gen_formula_from_compo("[C6H12O2+H-H2O]+_13C")
gen_formula_from_compo("[2(C6H2O3)+NH4]+_13C3")
gen_formula_from_compo("[C6H12O2+H-H2O]++")
gen_formula_from_compo("[C6H12O2+H-H2O]--_13C1")
}
