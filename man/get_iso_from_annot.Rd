% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/function.R
\name{get_iso_from_annot}
\alias{get_iso_from_annot}
\title{Get isotopes from annotation}
\usage{
get_iso_from_annot(annotation)
}
\arguments{
\item{annotation}{Annotation in the form of \verb{[M+H]+_13C1}}
}
\value{
Return a data.table with the following isotopes informations:
\itemize{
\item \code{element}: element type (13C, 18O, ...)
\item \code{text}: label to use for the element
\item \code{elmt_nb}: Element count
\item \code{isotope}: Isotope number (13, 18, ...)
\item \code{mass}: Mass of the isotope
\item \code{ID}: Unique identifier (\code{integer})
}
}
\description{
Get isotopes from annotation
}
\examples{
get_iso_from_annot("[M+H]+_13C_18O")
get_iso_from_annot("[M+H]+_13C2")
}
