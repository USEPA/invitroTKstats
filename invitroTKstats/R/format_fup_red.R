#' Creates a Standardized Data Frame with Rapid Equilibrium Dialysis (RED)
#' Plasma Protein Binding (PPB) Data (Level-1)
#'
#' This function formats data describing mass spectrometry (MS) peak areas
#' from samples collected as part of in vitro measurement of chemical fraction
#' unbound in plasma using rapid equilibrium dialysis
#' \insertCite{waters2008validation}{invitroTKstats}.
#' The input data frame is organized into a standard set of columns and written
#' to a tab-separated text file.
#'
#' The data frame of observations should be annotated according to these types:
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
#' Chemical concentration is calculated qualitatively as a response and 
#' returned as a column in the output data frame:
#'
#' Response <- AREA / ISTD.AREA * ISTD.CONC
#'
#' @param FILENAME (Character) A string used to identify the output Level-1 file.
#' "<FILENAME>-fup-RED-Level1.tsv". (Defaults to "MYDATA".)
#'
#' @param data.in (Data Frame) A Level-0 data frame containing mass-spectrometry peak areas,
#' indication of chemical identity, and measurement type. The data frame should
#' contain columns with names specified by the following arguments:
#'
#' @param sample.col (Character) Column name of \code{data.in} containing the unique mass
#' spectrometry (MS) sample name used by the laboratory. (Defaults to
#' "Lab.Sample.Name".)
#' 
#param date 
#' 
#' @param date.col (Character) Column name of \code{data.in} containing the laboratory measurement
#' date. (Defaults to "Date".)
#' 
#' @param compound.col (Character) Column name of \code{data.in} containing the test compound.
#' (Defaults to "Compound.Name".)
#' 
#' @param dtxsid.col (Character) Column name of \code{data.in} containing EPA's DSSTox Structure
#' ID (\url{http://comptox.epa.gov/dashboard}). (Defaults to "DTXSID".)
#' 
#' @param lab.compound.col (Character) Column name of \code{data.in} containing the test compound
#' name used by the laboratory. (Defaults to "Lab.Compound.Name".)
#' 
#' @param type.col (Character) Column name of \code{data.in} containing the sample type (see table
#' under Details). (Defaults to "Sample.Type".)
#' 
#' @param cal (Character) MS calibration the samples were based on. Typically, this uses 
#' indices or dates to represent if the analyses were done on different machines on 
#' the same day or on different days with the same MS analyzer. (Defaults to \code{NULL}.) 
#' (Note: Single entry only, 
#' use only if all data were collected based on the same calibration.)
#' 
#' @param cal.col (Character) Column name containing \code{cal} 
#' information. (Defaults to "Cal".) (Note: \code{data.in} does not
#' necessarily have this field. If this field is missing, it can be auto-filled with the value 
#' specified in \code{cal}.)
#'
#' @param replicate (Numeric) Simultaneous replicates with the same analytical chemistry. 
#' (Defaults to \code{NULL}.) (Note: Single entry only, use only if all tested compounds 
#' use the same number of replicates.)
#' 
#' @param replicate.col (Character) Column name containing \code{replicate} 
#' information. (Defaults to "Replicate".) (Note: \code{data.in} does not
#' necessarily have this field. If this field is missing, it can be auto-filled with the value 
#' specified in \code{replicate}.)
#' 
#' @param dilution (Numeric) Number of times the sample was diluted before MS 
#' analysis. (Defaults to \code{NULL}.) (Note: Single entry only, use only if all 
#' samples underwent the same number of dilutions.)
#' 
#' @param dilution.col (Character) Column name containing \code{dilution} 
#' information. (Defaults to "Dilution.Factor".) (Note: \code{data.in} does not
#' necessarily have this field. If this field is missing, it can be auto-filled with the value 
#' specified in \code{dilution}.)
#'
#' @param time.col (Character) Column name of \code{data.in} containing the time (in hours) 
#' from the start of incubation to when the measurements were taken.
#' (Defaults to "Time".)
#' 
#param time 
#'
#' @param istd.col (Character) Column name of \code{data.in} containing the
#' MS peak area for the internal standard. (Defaults to "ISTD.Area".)
#' 
#' @param istd.name (Character) The identity of the internal standard. (Defaults to \code{NULL}.) 
#' (Note: Single entry only, use only if all tested compounds use the same internal standard.) 
#'
#' @param istd.name.col (Character) Column name containing \code{istd.name} information. 
#' (Defaults to "ISTD.Name".) (Note: \code{data.in} does not
#' necessarily have this field. If this field is missing, it can be auto-filled with the value 
#' specified in \code{istd.name}.)
#' 
#' @param istd.conc (Numeric) The concentration for the internal standard. (Defaults to \code{NULL}.) 
#' (Note: Single entry only, use only if all tested compounds have the same 
#' internal standard concentration.) 
#'
#' @param istd.conc.col (Character) Column name containing \code{istd.conc} information. 
#' (Defaults to "ISTD.Conc".) (Note: \code{data.in} does not
#' necessarily have this field. If this field is missing, it can be auto-filled with the value 
#' specified in \code{istd.conc}.)
#' 
#' @param test.nominal.conc (Numeric) The intended test chemical concentration 
#' at time zero. (Defaults to \code{NULL}.) (Note: Single entry only, use only 
#' if all tested compounds used the same concentration at time zero.)
#'
#' @param test.nominal.conc.col (Character) Column name containing \code{test.nominal.conc} 
#' information. (Defaults to "Test.Target.Conc".) (Note: \code{data.in} does not
#' necessarily have this field. If this field is missing, it can be auto-filled with the value 
#' specified in \code{test.nominal.conc}.)
#'
#' @param plasma.percent (Numeric) The percent of the physiological plasma concentration 
#' used in RED assay. (Defaults to \code{NULL}.) (Note: Single entry only, use only 
#' if all compounds were tested with the same plasma percent.)
#' 
#' @param plasma.percent.col (Character) Column name containing \code{plasma.percent} 
#' information. (Defaults to "Plasma.Percent".) (Note: \code{data.in} does not
#' necessarily have this field. If this field is missing, it can be auto-filled with the value 
#' specified in \code{plasma.percent}.)
#'
#' @param std.conc (Numeric) The standard test chemical concentration for 
#' the intrinsic clearance assay. (Defaults to \code{NULL}.) (Note: Single entry only, 
#' use only if the same standard concentration was used for all tested compounds.)
#' 
#' @param std.conc.col (Character) Column name containing \code{std.conc} 
#' information. (Defaults to "Standard.Conc".) (Note: \code{data.in} does not
#' necessarily have this field. If this field is missing, it can be auto-filled with the value 
#' specified in \code{std.conc}.)
#'
#' @param area.col (Character) Column name of \code{data.in} containing the target analyte (that
#' is, the test compound) MS peak area. (Defaults to "Area".)
#' 
#' @param analysis.method (Character) The analytical chemistry analysis method, 
#' typically "LCMS" or "GCMS", liquid chromatography or gas chromatography–mass spectrometry, respectively. 
#' (Defaults to \code{NULL}.) (Note: Single entry only, 
#' use only if the same method was used for all tested compounds.)
#'
#' @param analysis.method.col (Character) Column name containing \code{analysis.method} 
#' information. (Defaults to "Analysis.Method".) (Note: \code{data.in} does not
#' necessarily have this field. If this field is missing, it can be auto-filled with the value 
#' specified in \code{analysis.method}.)
#' 
#' @param analysis.instrument (Character) The instrument used for chemical analysis, 
#' for example "Waters ACQUITY I-Class UHPLC - Xevo TQ-S uTQMS". (Defaults to \code{NULL}.) 
#' (Note: Single entry only, use only if the same instrument was used for all tested compounds.)
#'
#' @param analysis.instrument.col (Character) Column name containing \code{analysis.instrument} 
#' information. (Defaults to "Analysis.Instrument".) (Note: \code{data.in} does not
#' necessarily have this field. If this field is missing, it can be auto-filled with the value 
#' specified in \code{analysis.instrument}.)
#' 
#' @param analysis.parameters (Character) The parameters used to identify the 
#' compound on the chemical analysis instrument. (Defaults to \code{NULL}.) 
#' (Note: Single entry only, use only if the same parameters were used for all tested compounds.) 
#'
#' @param analysis.parameters.col (Character) Column name containing \code{analysis.parameters} 
#' information. (Defaults to "Analysis.Parameters".) (Note: \code{data.in} does not
#' necessarily have this field. If this field is missing, it can be auto-filled with the value 
#' specified in \code{analysis.parameters}.)
#' 
#' @param note.col (Character) Column name of \code{data.in} containing additional notes on 
#' test compounds. (Defaults to "Note".)
#'
#' @param level0.file.col (Character) Column name containing \code{level0.file} information. 
#' (Defaults to "Level0.File".) (Note: \code{data.in} does not
#' necessarily have this field. If this field is missing, it can be auto-filled with the value 
#' specified in \code{level0.file}.)
#'
#' @param level0.file (Character) The Level-0 file from which the \code{data.in} were obtained.
#' (Defaults to \code{NULL}.) (Note: Single entry only, use only if all rows in \code{data.in}
#' were obtained from the same Level-0 file.) 
#' 
#' @param level0.sheet.col (Character) Column name containing \code{level0.sheet} information.
#' (Defaults to "Level0.Sheet".) (Note: \code{data.in} does not
#' necessarily have this field. If this field is missing, it can be auto-filled with the value 
#' specified in \code{level0.sheet}.)
#'
#' @param level0.sheet (Character) The specific sheet name of Level-0 file from which the 
#' \code{data.in} is obtained from, if the Level-0 file is an Excel workbook. 
#' (Defaults to \code{NULL}.) (Note: Single entry only, use only if all rows in \code{data.in}
#' were obtained from the same sheet in the same Level-0 file.)
#' 
#' @return A Level-1 data frame with a standardized format containing a  
#' standardized set of columns and column names with plasma protein
#' binding (PPB) data from an rapid equilibrium dialysis (RED) assay. 
#'
#' @author John Wambaugh
#'
#' @examples
#' library(invitroTKstats)
#' red <- wambaugh2019.red
#' red$Date <- "2019"
#' red$Sample.Type <- "Blank"
#' red <- subset(red,!is.na(SampleName))
#' red[regexpr("PBS",red$SampleName)!=-1,"Sample.Type"] <- "PBS"
#' red[regexpr("Plasma",red$SampleName)!=-1,"Sample.Type"] <- "Plasma"
#' red$Dilution.Factor <- NA
#' red$Dilution.Factor <- as.numeric(red$Dilution.Factor)
#' red[red$Sample.Type=="PBS","Dilution.Factor"] <- 2
#' red[red$Sample.Type=="Plasma","Dilution.Factor"] <- 5
#' red[regexpr("T0",red$SampleName)!=-1,"Sample.Type"] <- "T0"
#'
#' red$Test.Target.Conc <- 5
#' red$ISTD.Name <- "Bucetin and Diclofenac"
#' red$ISTD.Conc <- 1
#' red$Series <- 1
#'
#' level1 <- format_fup_red(red,
#'  FILENAME="Wambaugh2019",
#'  sample.col="SampleName",
#'  compound.col="Preferred.Name",
#'  lab.compound.col="CompoundName",
#'  cal.col="RawDataSet")
#'
#' @references
#' \insertRef{waters2008validation}{invitroTKstats}
#'
#' @export format_fup_red
format_fup_red <- function(
  FILENAME = "MYDATA",
  data.in,
  sample.col="Lab.Sample.Name",
  #date=NULL,
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
  #time = 4,
  istd.col="ISTD.Area",
  istd.name=NULL,
  istd.name.col="ISTD.Name",
  istd.conc=NULL,
  istd.conc.col="ISTD.Conc",
  test.nominal.conc=NULL,
  test.nominal.conc.col="Test.Target.Conc",
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
    file=paste(FILENAME,"-fup-RED-Level0.tsv",sep=""),
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
  if (!is.null(test.nominal.conc)) data.in[,test.nominal.conc.col] <-
    test.nominal.conc
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
    test.nominal.conc.col,
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
    "Plasma",
    "PBS",
    "T0",
    "Plasma.Blank",
    "NoPlasma.Blank",
    "CC",
    "Stability",
    "EQ1",
    "EQ2"))
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
  test.nominal.conc.col <- "Test.Nominal.Conc"
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
    test.nominal.conc.col,
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
  
  # calculate the response:
  data.out[,"Response"] <- signif(data.out[,area.col] /
     data.out[,istd.col] * data.out[,istd.conc.col], 4)

# Write out a "level 1" file (data organized into a standard format):
  write.table(data.out,
    file=paste(FILENAME,"-fup-RED-Level1.tsv",sep=""),
    sep="\t",
    row.names=F,
    quote=F)

  summarize_table(data.out,
    req.types=c("Plasma","PBS","Plasma.Blank","NoPlasma.Blank"))

  return(data.out)
}


