#' Creates a standardized data table reporting Caco2 data
#'
#' This function formats data describing mass spectrometry (MS) peak areas
#' from samples collected as part of in vitro measurement of membrane 
#' permeability using Caco2 cells (). 
#' unbound in plasma using ultracentrifugation (Redgrave 1975?).
#' An input dataframe is organized into a standard set of columns and is written
#' to a tab-separated text file. 
#'
#' The data frame of observations should be annotated according to direction 
#' (either apical to basal -- "AtoB" -- or basal to apical -- "BtoA") and type
#' of concentrtion measured:
#' \tabular{rr}{
#'   Blank with no chemical added \tab Blank \tab \cr
#'   Dosing vehicle (C0) at target concentration \tab Dosing\cr
#'   Donor compartment at end of experiment \tab Donor\cr
#'   Receiver compartment at end of experiment\tab Receiver\cr
#' }
#' Chemical concentration is calculated qualitatively as a response:
#'
#' Response <- AREA / ISTD.AREA * ISTD.CONC
#'
#' @param FILENAME A string used to identify outputs of the function call.
#' (defaults to "MYDATA")
#' 
#' @param input.data A data frame containing mass-spectrometry peak areas,
#' indication of chemical identiy, and measurment type. The data frame should
#' contain columns with names specified by the following arguments:
#' 
#' @param sample.col Which column of input.data indicates the unique mass 
#' spectrometry (MS) sample name used by the laboratory. (Defaults to 
#' "Lab.Sample.Name")
#' 
#' @param lab.compound.col Which column of input.data indicates The test compound 
#' name used by the laboratory (Defaults to "Lab.Compound.Name")
#' 
#' @param dtxsid.col Which column of input.data indicates EPA's DSSTox Structure 
#' ID (\url{http://comptox.epa.gov/dashboard}) (Defaults to "DTXSID")
#' 
#' @param date.col Which column of input.data indicates the laboratory measurment
#' date (Defaults to "Date")
#' 
#' @param compound.col Which column of input.data indicates the test compound
#' (Defaults to "Compound.Name")
#' 
#' @param area.col Which column of input.data indicates the target analyte (that 
#' is, the test compound) MS peak area (Defaults to "Area")
#' 
#' @param type.col Which column of input.data indicates the sample type (see table
#' above)(Defaults to "Type")
#' 
#' @param type.col Which column of input.data indicates the direction of the 
#' measurements (either "AtoB" for apical to basolateral or "BtoA" for vice 
#' versa) (Defaults to "Direction")
#' 
#' @param cal.col Which column of input.data indicates the MS calibration -- for
#' instance different machines on the same day or different days with the same
#' MS analyzer (Defaults to "Cal")
#' 
#' @param cal If this argument is used (defaults to NULL) every observation in
#' the table is assigned the value of the argument and the corresponding
#' column in input.table (if present) is ignored.
#' 
#' #param compound.conc.col Which column indictes the intended concentration 
#' of the test chemical for calibration curves (Defaults to "Standard.Conc")
#'
#' @param dilution.col Which column of input.data indicates how many times the
#' sample was diluted before MS analysis (Defaults to "Dilution.Factor")
#' 
#' @param dilution If this argument is used (defaults to NULL) every 
#' observation in the table is assigned the value of the argument and the 
#' corresponding column in input.table (if present) is ignored.
#' 
#' @param istd.col Which column of input.data indicates the MS peak area for the
#' internal standard (Defaults to "ISTD.Area")
#' 
#' @param istd.name.col Which column of input.data indicates identity of the 
#' internal standard (Defaults to "ISTD.Name")
#' 
#' @param istd.name If this argument is used (defaults to NULL) every 
#' observation in the table is assigned the value of the argument and the 
#' corresponding column in input.table (if present) is ignored.
#' 
#' @param istd.conc.col Which column of input.data indicates the concentration of
#' the internal standard (Defaults to "ISTD.Conc")
#' 
#' @param istd.conc If this argument is used (defaults to NULL) every 
#' observation in the table is assigned the value of the argument and the 
#' corresponding column in input.table (if present) is ignored.
#' 
#' @param nominal.test.conc.col Which column of input.data indicates the intended
#' test chemical concentration at time zero (Defaults to "Test.Target.Conc") 
#' 
#' @param nominal.test.conc If this argument is used (defaults to NULL) every 
#' observation in the table is assigned the value of the argument and the 
#' corresponding column in input.table (if present) is ignored.
#'
#' @param analysis.method.col Which column of input.data indicates the analytical
#' chemistry analysis method, typically "LCMS" or "GCMS" (Defaults to 
#' "Analysis.Method")
#' 
#' @param analysis.method If this argument is used (defaults to NULL) every 
#' observation in the table is assigned the value of the argument and the 
#' corresponding column in input.table (if present) is ignored.
#'
#' @param analysis.instrument.col Which column of input.data indicates the 
#' instrument used for chemical analysis, for example 
#' "Agilent 6890 GC with model 5973 MS" (Defaults to 
#' "Analysis.Instrument")
#' 
#' @param analysis.instrument If this argument is used (defaults to NULL) every 
#' observation in the table is assigned the value of the argument and the 
#' corresponding column in input.table (if present) is ignored.
#'
#' @param analysis.parameters.col Which column of input.data indicates the 
#' parameters used to identify the compound on the chemical analysis instrument,
#' for example 
#' "Negative Mode, 221.6/161.6, -DPb=26, FPc=-200, EPd=-10, CEe=-20, CXPf=-25.0"
#' (Defaulys to "Analysis.Paramaters"). 
#' 
#' @param analysis.parameters If this argument is used (defaults to NULL) every 
#' observation in the table is assigned the value of the argument and the 
#' corresponding column in input.table (if present) is ignored.
#'
#' @return data.frame A data.frame in standardized "level1" format 
#'
#' @author John Wambaugh
#' 
#' @examples
#' library(invitroTKstats)
#' level0 <- kreutz2020
#' level1 <- format_fup_uc(level0,
#'   FILENAME="Kreutz2020",
#'   compound.col="Name",
#'   compound.conc.col="Standard.Conc",
#'   area.col="Chem.Area"
#'   )
#' 
#' @references
#' Hubatsch, Ina, Eva GE Ragnarsson, and Per Artursson. 
#' "Determination of drug permeability and prediction of drug absorption in 
#' Caco-2 monolayers." Nature protocols 2.9 (2007): 2111.
#' 
#' @export format_caco2
format_caco2 <- function(input.data,
  FILENAME = "MYDATA",
  sample.col="Lab.Sample.Name",
  lab.compound.col="Lab.Compound.Name",
  dtxsid.col="DTXSID",
  date.col="Date",
  compound.col="Compound.Name",
  area.col="Area",
  istd.col="ISTD.Area",
  type.col="Type",
  direction.col="Direction",
  compound.conc=NULL,
  compound.conc.col="Nominal.Conc",
  cal=NULL,
  cal.col="Cal",
  dilution=NULL,
  dilution.col="Dilution.Factor",
  istd.name=NULL,
  istd.name.col="ISTD.Name",
  istd.conc=NULL,
  istd.conc.col="ISTD.Conc",
  nominal.test.conc=NULL,
  nominal.test.conc.col="Test.Target.Conc",
  analysis.method=NULL,
  analysis.method.col="Analysis.Method",
  analysis.instrument=NULL,
  analysis.instrument.col="Analysis.Instrument",
  analysis.parameters=NULL,
  analysis.parameters.col="Analysis.Parameters" 
  )
{  
  output.data <- as.data.frame(input.data)

# These arguments allow the user to specify a single value for every obseration 
# in the table:  
  if (!is.null(cal)) output.data$Cal <- cal
  if (!is.null(dilution)) output.data$Dilution.Factor <- dilution
  if (!is.null(istd.name)) output.data$ISTD.Name <- istd.name
  if (!is.null(istd.conc)) output.data$ISTD.Conc <- istd.conc
  if (!is.null(nominal.test.conc)) output.data$Test.Target.Conc <- 
    nominal.test.conc
  if (!is.null(analysis.method)) output.data$Analysis.Method <- analysis.method
  if (!is.null(analysis.instrument)) output.data$Analysis.Instrument <- 
    analysis.instrument
  if (!is.null(analysis.parameters)) output.data$Analysis.Parameters <- 
    analysis.parameters
  
  
  
  
  
# We need all these columns in input.data
  cols <-c(
    sample.col,
    date.col,
    compound.col,
    dtxsid.col,
    lab.compound.col,
    type.col,
    dilution.col,
    cal.col,
    compound.conc.col,
    nominal.test.conc.col,
    istd.name.col,
    istd.conc.col,
    istd.col,
    series.col,
    area.col,
    analysis.method.col,
    analysis.instrument.col,
    analysis.parameters.col
    )
  
  if (!(all(cols %in% colnames(input.data))))
  {
    stop(paste("Missing columns named:",
      paste(cols[!(cols%in%colnames(input.data))],collapse=", ")))
  }

  # Only include the data types used:
  input.data <- subset(input.data,input.data[,type.col] %in% c("CC","T1","T5","AF"))
  
  # Organize the columns:
  input.data <- input.data[,cols]
    
  # Standardize the column names:
    sample.col <- "Lab.Sample.Name"
    date.col <- "Date"
    compound.col <- "Compound.Name"
    dtxsid.col <- "DTXSID"
    lab.compound.col <- "Lab.Compound.Name"
    type.col <- "Sample.Type"
    dilution.col <- "Dilution.Factor"
    cal.col <- "Calibration"
    compound.conc.col <- "Standard.Conc"
    nominal.test.conc.col <- "Test.Target.Conc"
    istd.name.col <- "ISTD.Name"
    istd.conc.col <- "ISTD.Conc"
    istd.col <- "ISTD.Area"
    series.col <- "Series"
    area.col <- "Area"
    analysis.method.col <- "Analysis.Method"
    analysis.instrument.col <- "Analysis.Instrument"
    analysis.parameters.col <- "Analysis.Parameters" 
    
  colnames(input.data) <- c(
    sample.col,
    date.col,
    compound.col,
    dtxsid.col,
    lab.compound.col,
    type.col,
    dilution.col,
    cal.col,
    compound.conc.col,
    nominal.test.conc.col,
    istd.name.col,
    istd.conc.col,
    istd.col,
    series.col,
    area.col,
    analysis.method.col,
    analysis.instrument.col,
    analysis.parameters.col
    )
  
  # calculate the reponse:
  input.data[,"Response"] <- input.data[,area.col] /
     input.data[,istd.col] *  input.data[,istd.conc.col]
  
# Write out a "level 1" file (data organized into a standard format):  
  write.table(input.data, 
    file=paste(FILENAME,"-PPB-UC-Level1.tsv",sep=""),
    sep="\t",
    row.names=F,
    quote=F)
    
  summarize_table(input.data,
    req.types=c("CC","T1","T5","AF"))

  return(input.data)  
}


