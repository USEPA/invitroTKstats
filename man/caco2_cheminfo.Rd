% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{caco2_cheminfo}
\alias{caco2_cheminfo}
\title{Caco-2 Chemical Information Example Data set}
\format{
A chemical info data.frame with 554 rows and 7 variables: \describe{
\item{\code{DTXSID}}{DSSTox Substance Identifier (CompTox Chemicals Dashboard)}
\item{\code{PREFERRED_NAME}}{Preferred compound name from the CompTox Chemicals Dashboard (CCD)}
\item{\code{CASRN}}{CAS Registry Number of the test compound}
\item{\code{MOLECULAR_FORMULA}}{Molecular formula of the test compound}
\item{\code{AVERAGE_MASS}}{Molecular weight of the compound in daltons}
\item{\code{QSAR_READY_SMILES}}{SMILES (Simplified molecular-input line-entry system) chemical structure description.}
\item{\code{test_article}}{Compound ID used in the laboratory}
}
}
\usage{
caco2_cheminfo
}
\description{
The chemical ID mapping information from tandem mass spectrometry (MS/MS) measurements 
of Caco-2 assay-specific data \insertCite{honda2025impact}{invitroTKstats} . 
This data set contains 520 unique compounds/chemicals.
}
\references{
\insertRef{honda2025impact}{invitroTKstats}
}
\keyword{datasets}
