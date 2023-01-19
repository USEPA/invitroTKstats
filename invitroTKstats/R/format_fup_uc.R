#' Creates a standardized data table reporting UC PPB data
#'
#' This function formats data describing mass spectrometry (MS) peak areas
#' from samples collected as part of in vitro measurement of chemical fraction
#' unbound in plasma using ultracentrifugation (Redgrave 1975?).
#' An input dataframe is organized into a standard set of columns and is written
#' to a tab-separated text file. 
#'
#' The data frame of observations should be annotated according to
#' of these types:
#' \tabular{rrrrr}{
#'   Calibration Curve \tab CC\cr
#'   Ultra-centrifugation Aqueous Fraction \tab UC\cr
#'   Whole Plasma T1h Sample  \tab T1\cr
#'   Whole Plasma T5h Sample \tab T5\cr
#' }
#' Chemical concentration is calculated qualitatively as a response:
#'
#' Response <- AREA / ISTD.AREA * ISTD.CONC
#'
#' @param FILENAME A string used to identify outputs of the function call.
#' (defaults to "MYDATA")
#' 
#' @param PPB.data A data frame containing mass-spectrometry peak areas,
#' indication of chemical identiy, and measurment type. The data frame should
#' contain columns with names specified by the following arguments:
#' 
#' @param sample.col Which column of PPB.data indicates the unique mass 
#' spectrometry (MS) sample name used by the laboratory. (Defaults to 
#' "Lab.Sample.Name")
#' 
#' @param lab.compound.col Which column of PPB.data indicates The test compound 
#' name used by the laboratory (Defaults to "Lab.Compound.Name")
#' 
#' @param dtxsid.col Which column of PPB.data indicates EPA's DSSTox Structure 
#' ID (\url{http://comptox.epa.gov/dashboard}) (Defaults to "DTXSID")
#' 
#' @param date.col Which column of PPB.data indicates the laboratory measurment
#' date (Defaults to "Date")
#' 
#' @param compound.col Which column of PPB.data indicates the test compound
#' (Defaults to "Compound.Name")
#' 
#' @param area.col Which column of PPB.data indicates the target analyte (that 
#' is, the test compound) MS peak area (Defaults to "Area")
#' 
#' @param series.col Which column of PPB.data indicates the "series", that is
#' a simultaneous replicate with the same analytical chemistry 
#' (Defaults to "Series")
#' 
#' @param type.col Which column of PPB.data indicates the sample type (see table
#' above)(Defaults to "Sample.Type")
#' 
#' @param cal If this argument is used (defaults to NULL) every 
#' observation in the table is assigned the value of the argument and the 
#' corresponding column in input.table (if present) is ignored.
#' 
#' @param cal.col Which column of PPB.data indicates the MS calibration -- for
#' instance different machines on the same day or different days with the same
#' MS analyzer (Defaults to "Cal")
#' 
#' @param std.conc.col Which column indictes the intended concentration 
#' of the test chemical for calibration curves in uM (Defaults to "Standard.Conc")
#' 
#' @param dilution If this argument is used (defaults to NULL) every 
#' observation in the table is assigned the value of the argument and the 
#' corresponding column in input.table (if present) is ignored.
#'
#' @param dilution.col Which column of PPB.data indicates how many times the
#' sample was diluted before MS analysis (Defaults to "Dilution.Factor")
#' 
#' @param istd.col Which column of PPB.data indicates the MS peak area for the
#' internal standard (Defaults to "ISTD.Area")
#' 
#' @param istd.name If this argument is used (defaults to NULL) every 
#' observation in the table is assigned the value of the argument and the 
#' corresponding column in input.table (if present) is ignored.
#' 
#' @param istd.name.col Which column of PPB.data indicates identity of the 
#' internal standard (Defaults to "ISTD.Name")
#' 
#' @param istd.conc If this argument is used (defaults to NULL) every 
#' observation in the table is assigned the value of the argument and the 
#' corresponding column in input.table (if present) is ignored.
#' 
#' @param istd.conc.col Which column of PPB.data indicates the concentration of
#' the internal standard in uM (Defaults to "ISTD.Conc")
#' 
#' @param uc.assay.conc If this argument is used (defaults to NULL) every 
#' observation in the table is assigned the value of the argument and the 
#' corresponding column in input.table (if present) is ignored.
#'
#' @param uc.assay.conc.col Which column indicates the intended initial
#' test chemical concentration in the UC assay in uM (Defaults to "Test.Target.Conc") 
#'
#' @param analysis.method If this argument is used (defaults to NULL) every 
#' observation in the table is assigned the value of the argument and the 
#' corresponding column in input.table (if present) is ignored.
#'
#' @param analysis.method.col Which column of PPB.data indicates the analytical
#' chemistry analysis method, typically "LCMS" or "GCMS" (Defaults to 
#' "Analysis.Method")
#'
#' @param analysis.instrument If this argument is used (defaults to NULL) every 
#' observation in the table is assigned the value of the argument and the 
#' corresponding column in input.table (if present) is ignored.
#'
#' @param analysis.instrument.col Which column of PPB.data indicates the 
#' instrument used for chemical analysis, for example 
#' "Agilent 6890 GC with model 5973 MS" (Defaults to 
#' "Analysis.Instrument")
#'
#' @param analysis.parameters If this argument is used (defaults to NULL) every 
#' observation in the table is assigned the value of the argument and the 
#' corresponding column in input.table (if present) is ignored.
#'
#' @param analysis.parameters.col Which column of PPB.data indicates the 
#' parameters used to identify the compound on the chemical analysis instrument,
#' for example 
#' "Negative Mode, 221.6/161.6, -DPb=26, FPc=-200, EPd=-10, CEe=-20, CXPf=-25.0"
#' (Defaulys to "Analysis.Paramaters"). 
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
#'   std.conc.col="Standard.Conc",
#'   area.col="Chem.Area"
#'   )
#' 
#' @references
#' Redgrave, T. G., D. C. K. Roberts, and C. E. West. "Separation of plasma 
#' lipoproteins by density-gradient ultracentrifugation." Analytical 
#' Biochemistry 65.1-2 (1975): 42-49.#' 
#' 
#' @export format_fup_uc
format_fup_uc <- function(PPB.data,
  FILENAME = "MYDATA",
  sample.col="Lab.Sample.Name",
  lab.compound.col="Lab.Compound.Name",
  dtxsid.col="DTXSID",
  date.col="Date",
  compound.col="Compound.Name",
  area.col="Area",
  series.col="Series",
  type.col="Sample.Type",
  std.conc=NULL,
  std.conc.col="Standard.Conc",
  cal=NULL,
  cal.col="Cal",
  dilution=NULL,
  dilution.col="Dilution.Factor",
  istd.col="ISTD.Area",
  istd.name=NULL,
  istd.name.col="ISTD.Name",
  istd.conc=NULL,
  istd.conc.col="ISTD.Conc",
  uc.assay.conc=NULL,
  uc.assay.conc.col="UC.Assay.Conc",
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
  PPB.data <- as.data.frame(PPB.data)

# Write out a "level 0" file (data the function received it):  
  write.table(PPB.data, 
    file=paste(FILENAME,"-PPB-UC-Level0.tsv",sep=""),
    sep="\t",
    row.names=F,
    quote=F)
     
  if (is.null(note.col)) PPB.data[,"Note"] <- ""

# These arguments allow the user to specify a single value for every obseration 
# in the table:  
  if (!is.null(cal)) PPB.data[,cal.col] <- cal
  if (!is.null(dilution)) PPB.data[,dilution.col] <- dilution
  if (!is.null(istd.name)) PPB.data[,istd.name.col] <- istd.name
  if (!is.null(istd.conc)) PPB.data[,istd.conc.col] <- istd.conc
  if (!is.null(std.conc)) PPB.data[,std.conc.col] <- std.conc
  if (!is.null(uc.assay.conc)) PPB.data[,uc.assay.conc.col] <- uc.assay.conc
  if (!is.null(analysis.method)) PPB.data[,analysis.method.col]<- analysis.method
  if (!is.null(analysis.instrument)) PPB.data[,analysis.instrument.col] <- 
    analysis.instrument
  if (!is.null(analysis.parameters)) PPB.data[,analysis.parameters.col] <- 
    analysis.parameters
  if (!is.null(level0.file)) PPB.data[,level0.file.col] <- level0.file
  if (!is.null(level0.sheet)) PPB.data[,level0.sheet.col] <- level0.sheet

# We need all these columns in PPB.data
  cols <-c(
    sample.col,
    date.col,
    compound.col,
    dtxsid.col,
    lab.compound.col,
    type.col,
    dilution.col,
    cal.col,
    std.conc.col,
    uc.assay.conc.col,
    istd.name.col,
    istd.conc.col,
    istd.col,
    series.col,
    area.col,
    analysis.method.col,
    analysis.instrument.col,
    analysis.parameters.col,
    note.col,
    level0.file.col,
    level0.sheet.col
    )
  
  if (!(all(cols %in% colnames(PPB.data))))
  {
    stop(paste("Missing columns named:",
      paste(cols[!(cols%in%colnames(PPB.data))],collapse=", ")))
  }

  # Only include the data types used:
  PPB.data <- subset(PPB.data,PPB.data[,type.col] %in% c("CC","T1","T5","AF"))
  
  # Organize the columns:
  PPB.data <- PPB.data[,cols]
    
  # Standardize the column names:
    sample.col <- "Lab.Sample.Name"
    date.col <- "Date"
    compound.col <- "Compound.Name"
    dtxsid.col <- "DTXSID"
    lab.compound.col <- "Lab.Compound.Name"
    type.col <- "Sample.Type"
    dilution.col <- "Dilution.Factor"
    cal.col <- "Calibration"
    std.conc.col <- "Standard.Conc"
    uc.assay.conc.col <- "UC.Assay.T1.Conc"
    istd.name.col <- "ISTD.Name"
    istd.conc.col <- "ISTD.Conc"
    istd.col <- "ISTD.Area"
    series.col <- "Series"
    area.col <- "Area"
    analysis.method.col <- "Analysis.Method"
    analysis.instrument.col <- "Analysis.Instrument"
    analysis.parameters.col <- "Analysis.Parameters" 
    note.col <- "Note"
    level0.file.col <- "Level0.File"
    level0.sheet.col <- "Level0.Sheet"
    
  colnames(PPB.data) <- c(
    sample.col,
    date.col,
    compound.col,
    dtxsid.col,
    lab.compound.col,
    type.col,
    dilution.col,
    cal.col,
    std.conc.col,
    uc.assay.conc.col,
    istd.name.col,
    istd.conc.col,
    istd.col,
    series.col,
    area.col,
    analysis.method.col,
    analysis.instrument.col,
    analysis.parameters.col,
    note.col,
    level0.file.col,
    level0.sheet.col
    )

  # Blanks don't always have internal standard -- add average ISTD.Area
  # First identify the blanks (have to deal with NA standard.concs:
  blanks <- PPB.data[,"Standard.Conc"]
  blanks[is.na(blanks)] <- -999
  blanks <- blanks == 0
  for (this.chem in unique(PPB.data[,"DTXSID"]))
  {
    this.subset <- subset(PPB.data, DTXSID==this.chem)
    for (this.cal in unique(PPB.data[,"Calibration"]))
    {
      this.cal.subset <- subset(this.subset, Calibration==this.cal)
      if (any(is.na(this.cal.subset[,"ISTD.Area"])))
      {
        this.mean.ISTD <- signif(mean(this.cal.subset$ISTD.Area,na.rm=TRUE))
        which.indices <- PPB.data[,"DTXSID"] == this.chem &
          PPB.data[,"Calibration"] == this.cal &
          is.na(PPB.data[,"ISTD.Area"]) & 
          blanks
        PPB.data[which.indices, 
                 "ISTD.Area"] <- this.mean.ISTD
        PPB.data[which.indices, 
                 "Area"] <- 0
      }
    }
  }

  # calculate the response:
  PPB.data[,"Response"] <- signif(as.numeric(PPB.data[,area.col]) /
     as.numeric(PPB.data[,istd.col]) * as.numeric(PPB.data[,istd.conc.col]),4)
  
# Write out a "level 1" file (data organized into a standard format):  
  write.table(PPB.data, 
    file=paste(FILENAME,"-PPB-UC-Level1.tsv",sep=""),
    sep="\t",
    row.names=F,
    quote=F)
    
  summarize_table(PPB.data,
    req.types=c("CC","T1","T5","AF"))

  return(PPB.data)  
}


