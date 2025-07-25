% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_method_table.R
\name{create_method_table}
\alias{create_method_table}
\title{Creates a Standardized Data Table for Chemical Analysis Methods}
\usage{
create_method_table(
  input.table,
  dtxsid.col = "DTXSID",
  compound.col = "Compound.Name",
  istd.name.col = "ISTD.Name",
  analysis.method.col = "Analysis.Method",
  analysis.instrument.col = "Analysis.Instrument",
  analysis.parameters.col = "Analysis.Parameters"
)
}
\arguments{
\item{input.table}{(Data Frame) A level-1 or level-2 data frame containing mass-spectrometry peak areas,
indication of chemical identity, and analytical chemistry methods.
It should contain columns with names specified by the following arguments:}

\item{dtxsid.col}{(Character) Column name of \code{input.table} containing EPA's DSSTox Structure 
ID (\url{http://comptox.epa.gov/dashboard}). (Defaults to "DTXSID".)}

\item{compound.col}{(Character) Column name of \code{input.table} containing the test compound.
(Defaults to "Compound.Name".)}

\item{istd.name.col}{(Character) Column name of \code{input.table} containing identity of the 
internal standard. (Defaults to "ISTD.Name".)}

\item{analysis.method.col}{(Character) Column name of \code{input.table} containing the analytical
chemistry analysis method, typically "LCMS" or "GCMS", liquid or gas chromatography 
mass spectrometry, respectively. (Defaults to "Analysis.Method".)}

\item{analysis.instrument.col}{(Character) Column name of \code{input.table} containing the 
instrument used for chemical analysis. For example, "Agilent 6890 GC with model 5973 MS". 
(Defaults to "Analysis.Instrument".)}

\item{analysis.parameters.col}{(Character) Column name of \code{input.table} containing the 
parameters used to identify the compound on the chemical analysis instrument.
For example, "Negative Mode, 221.6/161.6, -DPb=26, FPc=-200, EPd=-10, CEe=-20, CXPf=-25.0".
(Defaults to "Analysis.Parameters".)}
}
\value{
A data frame with one row per chemical-method pair containing 
information on analysis parameters, instruments, internal standards, 
and compound identifiers used for each pair.
}
\description{
This function extracts the chemical analysis methods from a set of MS data and 
returns a data frame with each row representing a unique chemical-method
pair. (Unique chemical identified by DTXSID.) Each row contains all compound names, analysis parameters, analysis instruments, 
and internal standards used for each chemical-method pair.
}
\examples{
library(invitroTKstats)
# Smeltz et al. (2020) data:
##  Clint ##
create_method_table(
  input.table = invitroTKstats::clint_L1,
  dtxsid.col = "DTXSID",
  compound.col = "Compound.Name",
  istd.name.col = "ISTD.Name",
  analysis.method.col = "Analysis.Method",
  analysis.instrument.col = "Analysis.Instrument",
  analysis.parameters.col = "Analysis.Parameters"
  )
## Fup RED ##
create_method_table(
  input.table = invitroTKstats::fup_red_L1,
  dtxsid.col = "DTXSID",
  compound.col = "Compound.Name",
  istd.name.col = "ISTD.Name",
  analysis.method.col = "Analysis.Method",
  analysis.instrument.col = "Analysis.Instrument",
  analysis.parameters.col = "Analysis.Parameters"
  )
## Fup UC ##
create_method_table(
  input.table = invitroTKstats::fup_uc_L1,
  dtxsid.col = "DTXSID",
  compound.col = "Compound.Name",
  istd.name.col = "ISTD.Name",
  analysis.method.col = "Analysis.Method",
  analysis.instrument.col = "Analysis.Instrument",
  analysis.parameters.col = "Analysis.Parameters"
  )
# Honda et al. () data:
## Caco2 ##
create_method_table(
  input.table = invitroTKstats::caco2_L1,
  dtxsid.col = "DTXSID",
  compound.col = "Compound.Name",
  istd.name.col = "ISTD.Name",
  analysis.method.col = "Analysis.Method",
  analysis.instrument.col = "Analysis.Instrument",
  analysis.parameters.col = "Analysis.Parameters"
  )

}
\author{
John Wambaugh
}
