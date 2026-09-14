#' Periodic table as a data.table
#'
#' @references \url{https://physics.nist.gov/cgi-bin/Compositions/
#' stand_alone.pl?ele=&ascii=ascii2&isotype=some}
#' @format A `data.table` with the Periodic Table Elements fields
#' \describe{
#'   \item{atomic_nb}{Atomic number}
#'   \item{atomic_symb}{Atomic symbol}
#'   \item{mass_nb}{Mass number}
#'   \item{atomic_mass}{Exact atomic mass}
#'   \item{isotopic_compo}{Relative isotopic abundance}
#' }
"Element"

#' Mono charged ions
#' 
#' @format A `data.table` with the common charges in LC-ESI-MS
#' \describe{
#'   \item{Formula}{Formula}
#'   \item{charge}{charge} 
#'   \item{lossL}{Does the resulting ion can produce in-source losses ?}
#'   \item{mz_query}{Exact mass value}
#' }
"db_monocharge"

#' adduct data
#'
#' @format A `data.table` with the common adducts in LC-ESI-MS
#' \describe{
#'   \item{adduct}{Formula}
#'   \item{charge}{charge} 
#'   \item{mz_query}{Exact mass value}
#' }
"Adduct_db"

#' isotopes data
#'
#' @format A `data.table` with the common isotopes found in LC-ESI-MS
#' \describe{
#'   \item{isotope}{Isotope formula}
#'   \item{mass}{Exact mass} 
#'   \item{abundance}{Natural abundance (in percent)}
#'   \item{mass_diff}{Exact mass difference from non-isotopic element}
#' }
"Isotopes_db"

#' losses data
#' 
#' @format A `data.table` with the common losses found in LC-ESI-MS
#' \describe{
#'   \item{loss}{loss formula}
#'   \item{mz_query}{Exact mass} 
#'   \item{encod}{Character encoding}
#' }
"Losses_db"

#' ions data
#'
#' @format A `data.table` with the common ions found in LC-ESI-MS
#' \describe{
#'   \item{Attribution}{Encoded formula}
#'   \item{exact_mass}{Exact mass formula as a character string} 
#'   \item{charge}{charge}
#'   \item{Type}{Molecular form (Monocharge, Loss, Adduct, ...)}
#' }
"ions_database"

#' Spectra example
#'
#' @format A `data.table` with a centroided MS spectra used as example
#' \describe{
#'   \item{mz}{m/Z value}
#'   \item{intensity}{Measured intensity} 
#'   \item{ID}{Unique identifier}
#' }
"spectra_full"

#' Spectra example
#'
#' MS2 spectrum of C10H12N5O6P1
#' Precursor mass: 330.0597
#' Mode: ESI +
#' Collision mode: HCD
#' Energy: 30
#'
#' @format A `data.table` with a centroided MS2 spectra used as example
#' \describe{
#'   \item{mz}{m/Z value}
#'   \item{i}{Measured intensity} 
#' }
"spectra_ms2"

#' Valence data
#'
#' @format A `data.table` with element valence
#' \describe{
#'   \item{atomic_symb}{Atomic symbol}
#'   \item{valence}{valence value} 
#' }
"valence_db"
