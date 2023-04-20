#' Creates a standardized data table reporting RED PPB data
#'
#' This function formats data describing mass spectrometry (MS) peak areas
#' from samples collected as part of in vitro measurement of chemical fraction
#' unbound in plasma using rapid equilibrium dialysis (Waters, et al, 2008).
#' An input dataframe is organized into a standard set of columns and is written
#' to a tab-separated text file. 
#'
#' The data frame of observations should be annotated according to
#' of these types:
#' \tabular{rrrrr}{
#'   Blank (ignored) \tab Blank\cr
#'   Plasma well concentration \tab Plasma\cr
#'   Phosphate-buffered well concentration\tab PBS\cr
#'   Time zero plasma concentration \tab T0\cr
#'   Plasma stability sample \tab Stability\cr
#'   Equilibrium Control Well 1 \tab EC1\cr
#'   Equilibrium Control Well 2 \tab EC2\cr
#'   Calibration Curve \tab CC\cr
#' }
#' Chemical concentration is calculated qualitatively as a response:
#'
#' Response <- AREA / ISTD.AREA * ISTD.CONC
#'
#' @param FILENAME A string used to identify outputs of the function call.
#' (defaults to "MYDATA")
#' 
#' @param data.in A data frame containing mass-spectrometry peak areas,
#' indication of chemical identiy, and measurment type. The data frame should
#' contain columns with names specified by the following arguments:
#' 
#' @param sample.col Which column of data.in indicates the unique mass 
#' spectrometry (MS) sample name used by the laboratory. (Defaults to 
#' "Lab.Sample.Name")
#' 
#' @param lab.compound.col Which column of data.in indicates The test compound 
#' name used by the laboratory (Defaults to "Lab.Compound.Name")
#' 
#' @param dtxsid.col Which column of data.in indicates EPA's DSSTox Structure 
#' ID (\url{http://comptox.epa.gov/dashboard}) (Defaults to "DTXSID")
#' 
#' @param date.col Which column of data.in indicates the laboratory measurment
#' date (Defaults to "Date")
#' 
#' @param compound.col Which column of data.in indicates the test compound
#' (Defaults to "Compound.Name")
#' 
#' @param area.col Which column of data.in indicates the target analyte (that 
#' is, the test compound) MS peak area (Defaults to "Area")
#' 
#' @param series.col Which column of data.in indicates the "series", that is
#' a simultaneous replicate (Defaults to "Series")
#' 
#' @param type.col Which column of data.in indicates the sample type (see table
#' above)(Defaults to "Sample.Type")
#' 
#' @param cal.col Which column of data.in indicates the MS calibration -- for
#' instance different machines on the same day or different days with the same
#' MS analyzer (Defaults to "Cal")
#' 
#' @param dilution.col Which column of data.in indicates how many times the
#' sample was diluted before MS analysis (Defaults to "Dilution.Factor")
#' 
#' @param istd.col Which column of data.in indicates the MS peak area for the
#' internal standard (Defaults to "ISTD.Area")
#' 
#' @param istd.name.col Which column of data.in indicates identity of the 
#' internal standard (Defaults to "ISTD.Name")
#' 
#' @param istd.conc.col Which column of data.in indicates the concentration of
#' the internal standard (Defaults to "ISTD.Conc")
#' 
#' @param nominal.test.conc.col Which column of data.in indicates the intended
#' test chemical concentration at time zero (Defaults to "Test.Target.Conc") 
#'
#' @param analysis.method.col Which column of data.in indicates the analytical
#' chemistry analysis method, typically "LCMS" or "GCMS" (Defaults to 
#' "Analysis.Method")
#'
#' @param analysis.instrument.col Which column of data.in indicates the 
#' instrument used for chemical analysis, for example 
#' "Agilent 6890 GC with model 5973 MS" (Defaults to 
#' "Analysis.Instrument")
#'
#' @param analysis.parameters.col Which column of data.in indicates the 
#' parameters used to identify the compound on the chemical analysis instrument,
#' for example 
#' "Negative Mode, 221.6/161.6, -DPb=26, FPc=-200, EPd=-10, CEe=-20, CXPf=-25.0"
#' (Defaulys to "Analysis.Paramaters"). 
#'
#' @return \item{data.frame}{A data.frame in standardized "level1" format} 
#'
#' @author John Wambaugh
#' 
#' @examples
#' library(invitroTKstats)
#'red <- wambaugh2019.red
#'red$Date <- "2019"
#'red$Sample.Type <- "Blank"
#'red <- subset(red,!is.na(SampleName))
#'red[regexpr("PBS",red$SampleName)!=-1,"Sample.Type"] <- "PBS"
#'red[regexpr("Plasma",red$SampleName)!=-1,"Sample.Type"] <- "Plasma"
#'red$Dilution.Factor <- NA
#'red$Dilution.Factor <- as.numeric(red$Dilution.Factor)
#'red[red$Sample.Type=="PBS","Dilution.Factor"] <- 2
#'red[red$Sample.Type=="Plasma","Dilution.Factor"] <- 5
#'red[regexpr("T0",red$SampleName)!=-1,"Sample.Type"] <- "T0"
#'
#'red$Test.Target.Conc <- 5
#'red$ISTD.Name <- "Bucetin and Diclofenac"
#'red$ISTD.Conc <- 1
#'red$Series <- 1
#'
#'level1 <- format_fup_red(red,
#'  FILENAME="Wambaugh2019",
#'  sample.col="SampleName",
#'  compound.col="Preferred.Name",
#'  lab.compound.col="CompoundName",
#'  cal.col="RawDataSet")
#'
#' @references
#' Waters, Nigel J., et al. "Validation of a rapid equilibrium dialysis 
#' approach for the measurement of plasma protein binding." Journal of 
#' Pharmaceutical Sciences 97.10 (2008): 4586-4595.
#'
#' @export format_fup_red
format_fup_red <- function(data.in,
  FILENAME = "MYDATA",
  sample.col="Lab.Sample.Name",
  date=NULL,
  date.col="Date",
  compound.col="Compound.Name",
  dtxsid.col="DTXSID",
  lab.compound.col="Lab.Compound.Name",
  type.col="Sample.Type",
  cal=NULL,
  cal.col="Cal",
  replicate=NULL,
  replicate.col="Replicate",
  dilution=NULL,
  dilution.col="Dilution.Factor",
  time.col="Time",
  time = 4,
  istd.col="ISTD.Area",
  istd.name=NULL,
  istd.name.col="ISTD.Name",
  istd.conc=NULL,
  istd.conc.col="ISTD.Conc",
  nominal.test.conc=NULL,
  nominal.test.conc.col="Test.Target.Conc",
  plasma.percent=NULL,
  plasma.percent.col="Plasma.Percent",
  std.conc=NULL,
  std.conc.col="Standard.Conc",
  area.col="Area",
  analysis.method=NULL,
  analysis.method.col="Analysis.Method",
  analysis.instrument=NULL,
  analysis.instrument.col="Analysis.Instrument",
  analysis.parameters=NULL,
  analysis.parameters.col="Analysis.Parameters",
  note.col="Note",
  level0.file.col="Level0.File",
  level0.file=NULL,
  level0.sheet.col="Level0.Sheet",
  level0.sheet=NULL
  )
{
  data.in <- as.data.frame(data.in)

  # Write out a "level 0" file (data as the function received it):  
  write.table(data.in, 
    file=paste(FILENAME,"-PPB-RED-Level0.tsv",sep=""),
    sep="\t",
    row.names=F,
    quote=F)
    
  if (is.null(note.col)) 
  {
    data.in[,"Note"] <- ""
    note.col <- "Note"
  }

# These arguments allow the user to specify a single value for every obseration 
# in the table:  
  if (!is.null(cal)) data.in[,cal.col] <- cal
  if (!is.null(replicate)) data.in[,replicate.col] <- replicate
  if (!is.null(dilution)) data.in[,dilution.col] <- dilution
  if (!is.null(istd.name)) data.in[,istd.name.col] <- istd.name
  if (!is.null(istd.conc)) data.in[,istd.conc.col] <- istd.conc
  if (!is.null(std.conc)) data.in[,std.conc.col] <- 
    std.conc
  if (!is.null(nominal.test.conc)) data.in[,nominal.test.conc.col] <- 
    nominal.test.conc
  if (!is.null(plasma.percent)) data.in[,plasma.percent.col] <- 
    plasma.percent
  if (!is.null(analysis.method)) data.in[,analysis.method.col]<- analysis.method
  if (!is.null(analysis.instrument)) data.in[,analysis.instrument.col] <- 
    analysis.instrument
  if (!is.null(analysis.parameters)) data.in[,analysis.parameters.col] <- 
    analysis.parameters
  if (!is.null(level0.file)) data.in[,level0.file.col] <- level0.file
  if (!is.null(level0.sheet)) data.in[,level0.sheet.col] <- level0.sheet
  
# We need all these columns in data.in
  cols <-c(
    sample.col,
    date.col,
    compound.col,
    dtxsid.col,
    lab.compound.col,
    type.col,
    dilution.col,
    replicate.col,
    cal.col,
    istd.name.col,
    istd.conc.col,
    istd.col,
    std.conc.col,
    nominal.test.conc.col,
    plasma.percent.col,
    time.col,
    area.col,
    analysis.method.col,
    analysis.instrument.col,
    analysis.parameters.col,
    note.col,
    level0.file.col,
    level0.sheet.col
    )
  
  if (!(all(cols %in% colnames(data.in))))
  {
    stop(paste("Missing columns named:",
      paste(cols[!(cols%in%colnames(data.in))],collapse=", ")))
  }

  # Only include the data types used:
  data.out <- subset(data.in,data.in[,type.col] %in% c(
    "Plasma","PBS","T0","Plasma.Blank","NoPlasma.Blank","CC","Stability","EQ1","EQ2"))
  # Force code to throw error if data.in accessed after this point:
  rm(data.in)
  
  # Organize the columns:
  data.out <- data.out[,cols]
    
  # Standardize the column names:
  sample.col <- "Lab.Sample.Name"
  date.col <- "Date"
  compound.col <- "Compound.Name"
  dtxsid.col <- "DTXSID"
  lab.compound.col <- "Lab.Compound.Name"
  type.col <- "Sample.Type"
  dilution.col <- "Dilution.Factor"
  replicate.col <- "Replicate"
  cal.col <- "Calibration"
  istd.name.col <- "ISTD.Name"
  istd.conc.col <- "ISTD.Conc"
  istd.col <- "ISTD.Area"
  std.conc.col <- "Std.Conc"
  nominal.test.conc.col <- "Nominal.Test.Conc"
  plasma.percent.col <- "Percent.Physiologic.Plasma"
  time.col <- "Time"
  area.col <- "Area"
  analysis.method.col <- "Analysis.Method"
  analysis.instrument.col <- "Analysis.Instrument"
  analysis.parameters.col <- "Analysis.Parameters" 
  note.col <- "Note"
  level0.file.col <- "Level0.File"
  level0.sheet.col <- "Level0.Sheet"

  colnames(data.out) <- c(
    sample.col,
    date.col,
    compound.col,
    dtxsid.col,
    lab.compound.col,
    type.col,
    dilution.col,
    replicate.col,
    cal.col,
    istd.name.col,
    istd.conc.col,
    istd.col,
    std.conc.col,
    nominal.test.conc.col,
    plasma.percent.col,
    time.col,
    area.col,
    analysis.method.col,
    analysis.instrument.col,
    analysis.parameters.col,
    note.col,
    level0.file.col,
    level0.sheet.col
    )

  # Set reasonable sig figs:
  for (this.col in c("Area", "ISTD.Area"))
    data.out[,this.col] <- signif(data.out[,this.col], 5)
  
  # calculate the reponse:
  data.out[,"Response"] <- signif(data.out[,area.col] /
     data.out[,istd.col] * data.out[,istd.conc.col], 4)
  
# Write out a "level 1" file (data organized into a standard format):  
  write.table(data.out, 
    file=paste(FILENAME,"-PPB-RED-Level1.tsv",sep=""),
    sep="\t",
    row.names=F,
    quote=F)

  summarize_table(data.out,
    req.types=c("Plasma","PBS","T0","Blank"))

  return(data.out)  
}


