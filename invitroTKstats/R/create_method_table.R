#' Creates a standardized data table for chemical analysis methods
#'
#' This function extracts the chemical analysis methods from a set of MS data.
#'
#' @param input.table A data frame containing mass-spectrometry peak areas,
#' indication of chemical identiy, and analytical chemistry methods.
#' It should contain columns with names specified by the following arguments:
#' 
#' @param dtxsid.col Which column of input.table indicates EPA's DSSTox Structure 
#' ID (\url{http://comptox.epa.gov/dashboard}) (Defaults to "DTXSID")
#' 
#' @param compound.col Which column of input.table indicates the test compound
#' (Defaults to "Compound.Name")
#' 
#' @param istd.name.col Which column of PPB.data indicates identity of the 
#' internal standard (Defaults to "ISTD.Name")
#' 
#' @param analysis.method.col Which column of PPB.data indicates the analytical
#' chemistry analysis method, typically "LCMS" or "GCMS" (Defaults to 
#' "Analysis.Method")
#'
#' @param analysis.instrument.col Which column of PPB.data indicates the 
#' instrument used for chemical analysis, for example 
#' "Agilent 6890 GC with model 5973 MS" (Defaults to 
#' "Analysis.Instrument")
#'
#' @param analysis.parameters.col Which column of PPB.data indicates the 
#' parameters used to identify the compound on the chemical analysis instrument,
#' for example 
#' "Negative Mode, 221.6/161.6, -DPb=26, FPc=-200, EPd=-10, CEe=-20, CXPf=-25.0"
#' (Defaulys to "Analysis.Paramaters"). 
#'
#' @return A data.frame with one row per method 
#'
#' @author John Wambaugh
#'
#' @export create_method_table
create_method_table <- function(input.table,
  dtxsid.col="DTXSID",
  compound.col="Compound.Name",
  istd.name.col="ISTD.Name",
  analysis.method.col="Analysis.Method",
  analysis.instrument.col="Analysis.Instrument",
  analysis.parameters.col="Analysis.Parameters"
  )
{
  N.chems <- length(unique(input.table[,dtxsid.col]))
  N.methods <- 0
  
  out.table <- NULL
  for (this.chem in sort(unique(input.table[,dtxsid.col])))
  {
    this.subset <- subset(input.table,input.table[,dtxsid.col]==this.chem)
    for (this.method in analysis.method.col)
    {
      this.method.subset <- subset(this.subset,
        this.subset[,analysis.method.col]==this.method)
      this.row <- data.frame(
        Compound.Name=paste(unique(this.method.subset[,compound.col]),
          collapse=", "),
        DTXSID=paste(unique(this.method.subset[,dtxisd.col]),
          collapse=", "),
        Analysis.Method=paste(unique(this.method.subset[,analysis.method.col]),
          collapse=", "),
        Analysis.Instrument=paste(unique(
          this.method.subset[,analysis.instrument.col]),
          collapse=", "),
        Analysis.Parameters=paste(unique(
          this.method.subset[,analysis.parameters.col]),
          collapse=", "),
        ISTD.Name=paste(unique(
          this.method.subset[,istd.name.col]),
          collapse=", "),
          )
        out.table <- rbind(out.table,this.row)
        N.methods <- N.methods + 1
      }
    }

  cat(paste(N.methods,"analytical methods for",N.chems,"chemicals.\n"))
  
  return(out.table)
}