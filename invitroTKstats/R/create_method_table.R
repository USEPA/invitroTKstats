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
#' @examples
#' library(invitroTKstats)
#'
#' create_method_table(smeltz2020)
#'
#' create_method_table(kreutz2020,compound.col="Name")
#'
#' # Create the method columns invitroTKstats is expecting from the
#' # Wambaugh et al. (2019) supplemental table 10:
#' wambaugh2019.methods$Method <- ""
#' # Need to convert NA's to something to allow logic to work on whole column:
#' wambaugh2019.methods[is.na(wambaugh2019.methods$LC),"LC"] <- ""
#' # Describe the chemcials with LC as such:
#' wambaugh2019.methods[wambaugh2019.methods$LC=="Y","Method"] <- "LC"
#' # Need to convert NA's to something to allow logic to work on whole column:
#' wambaugh2019.methods[is.na(wambaugh2019.methods$GC),"GC"] <- ""
#' # Describe the chemcials with GC as such:
#' wambaugh2019.methods[wambaugh2019.methods$GC=="Y","Method"] <- "GC"
#' # Remove the non GC/LC-able chemicals:
#' 
#' # Set instruments and notes:
#' # Agilent QQQ:
#' # Need to convert NA's to something to allow logic to work on whole column:
#' wambaugh2019.methods[is.na(wambaugh2019.methods$Agilent.QQQ),"Agilent.QQQ"] <-
#'   "Failed"
#' wambaugh2019.methods[wambaugh2019.methods$Agilent.QQQ!="Failed","Instrument"] <-
#'   "Agilent QQQ"
#' wambaugh2019.methods[wambaugh2019.methods$Agilent.QQQ!="Failed",
#'   "Analysis.Notes"] <- 
#'   wambaugh2019.methods[wambaugh2019.methods$Agilent.QQQ!="Failed","Agilent.QQQ"]
#' 
#' # Water's Xevo:
#' # Need to convert NA's to something to allow logic to work on whole column:
#' wambaugh2019.methods[is.na(wambaugh2019.methods$Water.s.Xevo),"Water.s.Xevo"] <-
#'   "Failed"
#' wambaugh2019.methods[wambaugh2019.methods$Water.s.Xevo!="Failed",
#'   "Instrument"] <- "Waters Xevo"
#' wambaugh2019.methods[wambaugh2019.methods$Water.s.Xevo!="Failed",
#'   "Analysis.Notes"] <- 
#'   wambaugh2019.methods[wambaugh2019.methods$Water.s.Xevo!="Failed",
#'   "Water.s.Xevo"]
#'   
#' # AB.Sciex.Qtrap:
#' # Need to convert NA's to something to allow logic to work on whole column:
#' wambaugh2019.methods[is.na(wambaugh2019.methods$AB.Sciex.Qtrap),
#'   "AB.Sciex.Qtrap"] <-  "Failed"
#' wambaugh2019.methods[wambaugh2019.methods$AB.Sciex.Qtrap!="Failed",
#'   "Instrument"] <- "AB Sciex QTRAP"
#' wambaugh2019.methods[wambaugh2019.methods$AB.Sciex.Qtrap!="Failed",
#'   "Analysis.Notes"] <- 
#'   wambaugh2019.methods[wambaugh2019.methods$AB.Sciex.Qtrap!="Failed",
#'   "AB.Sciex.Qtrap"]  
#' 
#' # Agilent.GCMS:
#' # Need to convert NA's to something to allow logic to work on whole column:
#' wambaugh2019.methods[is.na(wambaugh2019.methods$Agilent.GCMS),
#'   "Agilent.GCMS"] <-  "Failed"
#' wambaugh2019.methods[wambaugh2019.methods$Agilent.GCMS!="Failed",
#'   "Instrument"] <- "Agilent GCMS"
#' wambaugh2019.methods[wambaugh2019.methods$Agilent.GCMS!="Failed",
#'   "Analysis.Notes"] <- 
#'   wambaugh2019.methods[wambaugh2019.methods$Agilent.GCMS!="Failed",
#'   "Agilent.GCMS"]  
#' 
#' # GCTOF:
#' # Need to convert NA's to something to allow logic to work on whole column:
#' wambaugh2019.methods[is.na(wambaugh2019.methods$GCTOF),
#'   "GCTOF"] <-  "Failed"
#' wambaugh2019.methods[wambaugh2019.methods$GCTOF!="Failed",
#'   "Instrument"] <- "GCTOF"
#' wambaugh2019.methods[wambaugh2019.methods$GCTOF!="Failed",
#'   "Analysis.Notes"] <- 
#'   wambaugh2019.methods[wambaugh2019.methods$GCTOF!="Failed",
#'   "GCTOF"]  
#'   
#' # Add other needed columns:
#' wambaugh2019.methods$ISTD.Name <- "Bucetin and Diclofenac"
#' 
#' create_method_table(wambaugh2019.methods,
#'   compound.col="PREFERRED_NAME",
#'   analysis.method.col="Method",
#'   analysis.instrument.col="Instrument", 
#'   analysis.parameters.col="Analysis.Notes")
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