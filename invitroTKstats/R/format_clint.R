#' Creates a standardized data table reporting hepatocyte clearance data
#'
#' This function formats data describing mass spectrometry (MS) peak areas
#' from samples collected as part of in vitro measurement of chemical stability
#' when incubated with suspended hepatocytes (Shibata, 2000). Disappearance of 
#' the chemical over time is assumed to be due to metabolism by the hepatocytes.
#' An input data frame is organized into a standard set of columns and is written
#' to a tab-separated text file. 
#'
#' The data frame of observations should be annotated according to
#' of these types:
#' \tabular{rrrrr}{
#'   Blank \tab Blank\cr
#'   Hepatocyte incubation concentration \tab Cvst\cr
#'   Inactivated Hepatocytes \tab Inactive\cr
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
#' indication of chemical identity, and measurement type. The data frame should
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
#' @param date.col Which column of data.in indicates the laboratory measurement
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
#' @param density.col Which column of data.in indicates the density (units of
#' millions of hepatocytes per mL) hepatocytes in the in vitro incubation 
#' (Defaults to "Hep.Density" )
#' 
#' @param istd.col Which column of data.in indicates the MS peak area for the
#' internal standard (Defaults to "ISTD.Area")
#' 
#' @param istd.name.col Which column of data.in indicates identity of the 
#' internal standard (Defaults to "ISTD.Name")
#' 
#' @param istd.conc.col Which column of data.in indicates the concentration 
#' (units if uM) of
#' the internal standard (Defaults to "ISTD.Conc")
#' 
#' @param conc.col Which column of data.in indicates the intended
#' test chemical concentration 
#' (units if uM) of
#' at time zero (Defaults to "Conc") 
#'
#' @param time.col Which column of data.in indicates the intended
#' time of the measurment (in minutes) since the test chemical was introduced
#' into the hepatocyte incubation (Defaults to "Time") 
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
#' (Defaults to "Analysis.Parameters"). 
#'
#' @param FILENAME A string used to identify outputs of the function call.
#' (defaults to "MYDATA")
#' 
#' @param input.data A data frame containing mass-spectrometry peak areas,
#' indication of chemical identity, and measurement type. The data frame should
#' contain columns with names specified by the following arguments:
#' 
#' @param sample.col Which column of input.data indicates the unique mass 
#' spectrometry (MS) sample name used by the laboratory. (Defaults to 
#' "Lab.Sample.Name")
#' 
#' @param density.col Which column of input.data indicates the density of 
#' hepatocytes in suspension (10^6 hepatocytes / mL) (Defaults to "Hep.Density")
#' 
#' @param density.col A single value to be used for all samples indicating
#' the density of hepatocytes in suspension (10^6 hepatocytes / mL) 
#' (Defaults to NULL)
#' 
#' @param lab.compound.col Which column of input.data indicates The test compound 
#' name used by the laboratory (Defaults to "Lab.Compound.Name")
#' 
#' @param dtxsid.col Which column of input.data indicates EPA's DSSTox Structure 
#' ID (\url{http://comptox.epa.gov/dashboard}) (Defaults to "DTXSID")
#' 
#' @param date.col Which column of input.data indicates the laboratory measurement
#' date (Defaults to "Date")
#' 
#' @param series.col Which column of PPB.data indicates the "series", that is
#' a simultaneous replicate with the same analytical chemistry 
#' (Defaults to "Series")
#' 
#' @param series If this argument is used (defaults to NULL) every observation 
#' in the table is assigned the value of the argument and the corresponding
#' column in input.table (if present) is ignored.
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
#' #param compound.conc.col Which column indicates the intended concentration 
#' of the test chemical for calibration curves (Defaults to "Standard.Conc")
#'
#' @param dilution If this argument is used (defaults to NULL) every 
#' observation in the table is assigned the value of the argument and the 
#' corresponding column in input.table (if present) is ignored.
#' 
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
#' test chemical concentration at time zero in the dosing solution (added to the
#' donor side of the Caco-2 test well) (Defaults to "Test.Target.Conc") 
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
#' @param level0.file.col Which column of PPB.data indicates the file from
#' which the data were obtained (for example "MyWorkbook.xlsx").
#'  
#' @param level0.file If this argument is used (defaults to NULL) every
#' observation in the table is assigned the value of the argument and the
#' corresponding column in input.table (if present) is ignored.
#' 
#' @param level0.sheet.col Which column of PPB.data indicates the specific 
#' sheet containing the data if the file is an Excel workbook
#'  
#' @param level0.sheet If this argument is used (defaults to NULL) every
#' observation in the table is assigned the value of the argument and the
#' corresponding column in input.table (if present) is ignored.
#' 
#' @return \item{data.frame}{A data.frame in standardized "level1" format} 
#'
#' @author John Wambaugh
#' 
#' @examples
#' 
#' library(invitroTKstats)
#' 
#' clint <- wambaugh2019.clint
#' clint$Date <- "2019"
#' clint$Sample.Type <- "Blank"
#' clint$Time..mins. <- as.numeric(clint$Time..mins.)
#' clint[!is.na(clint$Time..mins.),"Sample.Type"] <- "Cvst"
#' clint$ISTD.Name <- "Bucetin, Propranolol, and Diclofenac"
#' clint$ISTD.Conc <- 1
#' clint$Dilution.Factor <- 1
#' clint[is.na(clint$FileName),"FileName"]<-"Wambaugh2019"
#' clint$Hep.Density <- 0.5
#' 
#' level1 <- format_clint(clint,
#'   FILENAME="Wambaugh2019",
#'   sample.col="Sample.Name",
#'   compound.col="Preferred.Name",
#'   lab.compound.col="Name",
#'   time.col="Time..mins.",
#'   cal.col="FileName")
#'
#' @references
#' \insertRef{shibata2002prediction}{invitroTKstats}
#'
#' @export format_clint
format_clint <- function(data.in,
  FILENAME = "MYDATA",
  sample.col="Lab.Sample.Name",
  date=NULL,
  date.col="Date",
  compound.col="Compound.Name",
  dtxsid.col="DTXSID",
  lab.compound.col="Lab.Compound.Name",
  type.col="Sample.Type",
  density=NULL,
  density.col="Hep.Density",
  compound.conc=NULL,
  compound.conc.col="Nominal.Conc",
  cal=NULL,
  cal.col="Cal",
  dilution=NULL,
  dilution.col="Dilution.Factor",
  time.col="Time",
  time = 2,
  istd.col="ISTD.Area",
  istd.name=NULL,
  istd.name.col="ISTD.Name",
  istd.conc=NULL,
  istd.conc.col="ISTD.Conc",
  std.conc=NULL,
  std.conc.col="Standard.Conc",
  clint.assay.conc=NULL,
  clint.assay.conc.col="Clint.Assay.Conc",
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
    file=paste(FILENAME,"-Clint-Level0.tsv",sep=""),
    sep="\t",
    row.names=F,
    quote=F)
    
  if (is.null(note.col)) 
  {
    data.in[,"Note"] <- ""
    note.col <- "Note"
  }
  
  if (!(std.conc.col %in% colnames(data.in)))
  {
    if (is.null(std.conc))
    {
      data.in[,std.conc.col] <- NA
    } else {
      data.in[,std.conc.col] <- std.conc
    }
  }
  
# These arguments allow the user to specify a single value for every obseration 
# in the table:  
  if (!is.null(cal)) data.in[,cal.col] <- cal
  if (!is.null(dilution)) data.in[,dilution.col] <- dilution
  if (!is.null(density)) data.in[,density.col] <- density
  if (!is.null(istd.name)) data.in[,istd.name.col] <- istd.name
  if (!is.null(istd.conc)) data.in[,istd.conc.col] <- istd.conc
  if (!is.null(std.conc)) data.in[,std.conc.col] <- 
    std.conc
  if (!is.null(clint.assay.conc)) data.in[,clint.assay.conc.col] <- 
    clint.assay.conc
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
    cal.col,
    istd.name.col,
    istd.conc.col,
    istd.col,
    density.col,
    std.conc.col,
    clint.assay.conc.col,
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
    "Blank","Cvst","CC","Inactive"))
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
  cal.col <- "Calibration"
  istd.name.col <- "ISTD.Name"
  istd.conc.col <- "ISTD.Conc"
  istd.col <- "ISTD.Area"
  density.col <- "Hep.Density"
  std.conc.col <- "Std.Conc"
  clint.assay.conc.col <- "Clint.Assay.Conc"
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
    cal.col,
    istd.name.col,
    istd.conc.col,
    istd.col,
    density.col,
    std.conc.col,
    clint.assay.conc.col,
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
    file=paste(FILENAME,"-Clint-Level1.tsv",sep=""),
    sep="\t",
    row.names=F,
    quote=F)

  summarize_table(data.out,
    req.types=c("Blank","Cvst"))

  return(data.out)  
}


