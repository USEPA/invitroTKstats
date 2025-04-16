#' Creates a Standardized Data Table for Chemical Analysis Methods
#'
#' This function extracts the chemical analysis methods from a set of MS data and 
#' returns a data frame with each row representing a unique chemical-method
#' pair. (Unique chemical identified by DTXSID's.) Each row contains all compound names, analysis parameters, analysis instruments, 
#' and internal standards used for each chemical-method pair. 
#'
#' @param input.table (Data Frame) A Level-1 or Level-2 data frame containing mass-spectrometry peak areas,
#' indication of chemical identity, and analytical chemistry methods.
#' It should contain columns with names specified by the following arguments:
#' 
#' @param dtxsid.col (Character) Column name of \code{input.table} containing EPA's DSSTox Structure 
#' ID (\url{http://comptox.epa.gov/dashboard}). (Defaults to "DTXSID".)
#' 
#' @param compound.col (Character) Column name of \code{input.table} containing the test compound.
#' (Defaults to "Compound.Name".)
#' 
#' @param istd.name.col (Character) Column name of \code{input.table} containing identity of the 
#' internal standard. (Defaults to "ISTD.Name".)
#' 
#' @param analysis.method.col (Character) Column name of \code{input.table} containing the analytical
#' chemistry analysis method, typically "LCMS" or "GCMS", liquid or gas chromatography 
#' mass spectrometry, respectively. (Defaults to "Analysis.Method".)
#'
#' @param analysis.instrument.col (Character) Column name of \code{input.table} containing the 
#' instrument used for chemical analysis. For example, "Agilent 6890 GC with model 5973 MS". 
#' (Defaults to "Analysis.Instrument".)
#'
#' @param analysis.parameters.col (Character) Column name of \code{input.table} containing the 
#' parameters used to identify the compound on the chemical analysis instrument.
#' For example, "Negative Mode, 221.6/161.6, -DPb=26, FPc=-200, EPd=-10, CEe=-20, CXPf=-25.0".
#' (Defaults to "Analysis.Parameters".) 
#'
#' @return A data frame with one row per chemical-method pair containing 
#' information on analysis parameters, instruments, internal standards, 
#' and compound identifiers used for each pair.
#'
#' @author John Wambaugh
#'
#' @examples
#' library(invitroTKstats)
#' # Smeltz et al. (2020) data:
#' ##  Clint ##
#' create_method_table(
#'   input.table = invitroTKstats::clint_L1,
#'   dtxsid.col = "DTXSID",
#'   compound.col = "Compound.Name",
#'   istd.name.col = "ISTD.Name",
#'   analysis.method.col = "Analysis.Method",
#'   analysis.instrument.col = "Analysis.Instrument",
#'   analysis.parameters.col = "Analysis.Parameters"
#'   )
#' ## Fup RED ##
#' create_method_table(
#'   input.table = invitroTKstats::fup_red_L1,
#'   dtxsid.col = "DTXSID",
#'   compound.col = "Compound.Name",
#'   istd.name.col = "ISTD.Name",
#'   analysis.method.col = "Analysis.Method",
#'   analysis.instrument.col = "Analysis.Instrument",
#'   analysis.parameters.col = "Analysis.Parameters"
#'   )
#' ## Fup UC ##
#' create_method_table(
#'   input.table = invitroTKstats::fup_uc_L1,
#'   dtxsid.col = "DTXSID",
#'   compound.col = "Compound.Name",
#'   istd.name.col = "ISTD.Name",
#'   analysis.method.col = "Analysis.Method",
#'   analysis.instrument.col = "Analysis.Instrument",
#'   analysis.parameters.col = "Analysis.Parameters"
#'   )
#' # Honda et al. () data:
#' ## Caco2 ##
#' create_method_table(
#'   input.table = invitroTKstats::caco2_L1,
#'   dtxsid.col = "DTXSID",
#'   compound.col = "Compound.Name",
#'   istd.name.col = "ISTD.Name",
#'   analysis.method.col = "Analysis.Method",
#'   analysis.instrument.col = "Analysis.Instrument",
#'   analysis.parameters.col = "Analysis.Parameters"
#'   )
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
# We need all these columns in input.table
  cols <-c(
    compound.col,
    dtxsid.col,
    istd.name.col,
    analysis.method.col,
    analysis.instrument.col,
    analysis.parameters.col
    )
  
  if (!(all(cols %in% colnames(input.table))))
  {
    stop(paste("Missing columns named:",
      paste(cols[!(cols%in%colnames(input.table))],collapse=", ")))
  }
  
# Get rid of blank methods:
  input.table <- subset(input.table, !is.na(input.table[,analysis.method.col]))
  input.table <- subset(input.table, input.table[,analysis.method.col] != "")
  
  N.chems <- length(unique(input.table[,dtxsid.col]))
  N.methods <- 0
  
  out.table <- NULL
  for (this.chem in sort(unique(input.table[,dtxsid.col])))
  {
    this.subset <- subset(input.table,input.table[,dtxsid.col]==this.chem)
    for (this.method in sort(unique(this.subset[,analysis.method.col])))
    {
      this.method.subset <- subset(this.subset,
        this.subset[,analysis.method.col]==this.method)
      this.row <- data.frame(
        Compound.Name=paste(unique(this.method.subset[,compound.col]),
          collapse=", "),
        DTXSID=paste(unique(this.method.subset[,dtxsid.col]),
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
          collapse=", ")
          )
        out.table <- rbind(out.table,this.row)
        N.methods <- N.methods + 1
      }
    }

  cat(paste(N.methods,"analytical methods for",N.chems,"chemicals.\n"))
  
  return(out.table)
}