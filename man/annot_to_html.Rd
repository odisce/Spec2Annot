% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annot_to_html.R
\name{annot_to_html}
\alias{annot_to_html}
\title{Convert annotation to html string}
\usage{
annot_to_html(annotation, compo = "", compo_replace = "X")
}
\arguments{
\item{annotation}{String corresponding to the annotation}

\item{compo}{Elemental composition to replace to in the annotation string}

\item{compo_replace}{String to replace with \code{compo} in the annotation}
}
\value{
An html string of the annotation.
}
\description{
Convert annotation to html string
}
\examples{
annot_to_html("[M+H]+_13C2", "C6H12O2", "M")
annot_to_html("[M+H-H2O]+_13C2", "C6H12O2", "M")
annot_to_html("[M+H-H2O]+_13C2_18O", "C6H12O2", "M")
annot_to_html("C6H13O2-H2O+_13C2_18O", "", "M")
annot_to_html("[M+H-(OH•)]+_13C2_18O", "C6H12O2", "M")
}
