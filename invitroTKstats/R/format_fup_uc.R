#' Creates a Standardized Data Frame Reporting Ultracentrifugation (UC) PPB Data (Level-1)
#'
#' This function formats data describing mass spectrometry (MS) peak areas
#' from samples collected as part of in vitro measurement of chemical fraction
#' unbound in plasma using ultracentrifugation
#' \insertCite{redgrave1975separation}{invitroTKstats}.
#' The input data frame is organized into a standard set of columns and is written
#' to a tab-separated text file.
#'
#' The data frame of observations should be annotated according to
#' of these types:
#' \tabular{rrrrr}{
#'   Calibration Curve \tab CC\cr
#'   Ultracentrifugation Aqueous Fraction \tab AF\cr
#'   Whole Plasma T1h Sample  \tab T1\cr
#'   Whole Plasma T5h Sample \tab T5\cr
#' }
#' Chemical concentration is calculated qualitatively as a response and 
#' returned as a column in the output data frame:
#'
#' Response <- AREA / ISTD.AREA * ISTD.CONC
#'
#' @param FILENAME (Character) A string used to identify the output Level-1 file.
#' "<FILENAME>-fup-UC-Level1.tsv". (Defaults to "MYDATA").
#'
#' @param data.in A data frame containing mass-spectrometry peak areas,
#' indication of chemical identity, and measurement type. The data frame should
#' contain columns with names specified by the following arguments:
#'
#' @param sample.col (Character) Column name of data.in containing the unique mass
#' spectrometry (MS) sample name used by the laboratory. (Defaults to
#' "Lab.Sample.Name".)
#'
#' @param lab.compound.col (Character) Column name of data.in containing the test compound
#' name used by the laboratory. (Defaults to "Lab.Compound.Name".)
#'
#' @param dtxsid.col (Character) Column name of data.in containing EPA's DSSTox Structure
#' ID (\url{http://comptox.epa.gov/dashboard}). (Defaults to "DTXSID".)
#'
#' @param date.col (Character) Column name of data.in containing the laboratory measurement
#' date. (Defaults to "Date".)
#'
#' @param compound.col (Character) Column name of data.in containing the test compound.
#' (Defaults to "Compound.Name".)
#'
#' @param area.col (Character) Column name of data.in containing the target analyte (that
#' is, the test compound) MS peak area. (Defaults to "Area".)
#'
#' @param series.col (Character) Column name of data.in containing the number of 
#' simultaneous replicate with the same analytical chemistry. 
#' (Defaults to "Series".)
#'
#' @param type.col (Character) Column name of data.in containing the sample type (see table
#' under Details). (Defaults to "Sample.Type".)
#' 
#' @param std.conc.col (Character) Column name containing \code{std.conc} 
#' information. (Defaults to "Standard.Conc".)
#'
#' @param std.conc (Numeric) The standard test chemical concentration for 
#' the intrinsic clearance assay. (Defaults to \code{NULL}.) (Note: Single entry only, 
#' use only if the same standard concentration was used for all tested compounds.)
#'
#' @param cal (Character) MS calibration the samples were based on, typically uses 
#' indices or dates to represent if the analyses were done on different machines on 
#' the same day or on different days with the same MS analyzer. (Defaults to \code{NULL}.) 
#' (Note: Single entry only, 
#' use only if all data were collected based on the same calibration.)
#'
#' @param cal.col (Character) Column name containing \code{cal} 
#' information. (Defaults to "Cal".)
#' 
#' @param dilution  (Numeric) Number of times the sample was diluted before MS 
#' analysis. (Defaults to \code{NULL}.) (Note: Single entry only, use only if all 
#' samples underwent the same number of dilutions.)
#'
#' @param dilution.col (Character) Column name containing \code{dilution} 
#' information. (Defaults to "Dilution.Factor".)
#'
#' @param istd.col (Character) Column name of data.in containing the
#' MS peak area for the internal standard. (Defaults to "ISTD.Area".)
#'
#' @param istd.name (Character) The identity of the internal standard. (Defaults to \code{NULL}.) 
#' (Note: Single entry only, use only if all tested compounds use the same internal standard.) 
#'
#' @param istd.name.col (Character) Column name containing \code{istd.name} information. 
#' (Defaults to "ISTD.Name".)
#'
#' @param istd.conc (Numeric) The concentration for the internal standard. (Defaults to \code{NULL}.) 
#' (Note: Single entry only, use only if all tested compounds have the same 
#' internal standard concentration.) 
#'
#' @param istd.conc.col (Character) Column name containing \code{istd.conc} information. 
#' (Defaults to "ISTD.Conc".)
#'
#' @param uc.assay.conc Which column indicates the intended initial
#' test chemical concentration in the UC assay in uM. (Defaults to \code{NULL}.)
#' (Note: Single entry only, 
#' use only if the same initial concentration was used for all tested compounds.)
#'
#' @param uc.assay.conc.col (Character) Column name containing \code{uc.assay.conc} 
#' information. (Defaults to "UC.Assay.Conc".)
#'
#' @param analysis.method (Character) The analytical chemistry analysis method, 
#' typically "LCMS" or "GCMS", liquid chromatography or gas chromatography–mass spectrometry, respectively. 
#' (Defaults to \code{NULL}.) (Note: Single entry only, 
#' use only if the same method was used for all tested compounds.)
#'
#' @param analysis.method.col (Character) Column name containing \code{analysis.method} 
#' information. (Defaults to "Analysis.Method".)
#' 
#' @param analysis.instrument (Character) The instrument used for chemical analysis, 
#' for example "Waters Xevo TQ-S micro (QEB0036)". (Defaults to \code{NULL}.) 
#' (Note: Single entry only, use only if the same instrument was used for all tested compounds.) 
#'
#' @param analysis.instrument.col (Character) Column name containing \code{analysis.instrument} 
#' information. (Defaults to "Analysis.Instrument".)
#'
#' @param analysis.parameters (Character) The parameters used to identify the 
#' compound on the chemical analysis instrument. (Defaults to \code{NULL}.) 
#' (Note: Single entry only, use only if the same parameters were used for all tested compounds.) 
#'
#' @param analysis.parameters.col (Character) Column name containing \code{analysis.parameters} 
#' information. (Defaults to "Analysis.Parameters".)
#' 
#' @param note.col (Character) Column name of data.in containing additional notes on 
#' test compounds. (Defaults to "Note").
#'
#' @param level0.file.col (Character) Column name containing \code{level0.file} information. 
#' (Defaults to "Level0.File".)
#'
#' @param level0.file (Character) The Level-0 file from which the data.in were obtained.
#' (Defaults to \code{NULL}.) (Note: Single entry only, use only if all rows in data.in
#' were obtained from the same Level-0 file.) 
#'
#' @param level0.sheet.col (Character) Column name containing \code{level0.sheet} information.
#' (Defaults to "Level0.Sheet".)
#'
#' @param level0.sheet (Character) The specific sheet name of Level-0 file from which the 
#' data.in is obtained from, if the level-0 file is an Excel workbook. 
#' (Defaults to \code{NULL}.) (Note: Single entry only, use only if all rows in data.in
#' were obtained from the same sheet in the same Level-0 file.) 
#'
#' @return A data frame in standardized Level-1 format containing a  
#' standardized set of columns with standardized column names.  
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
#' \insertRef{redgrave1975separation}{invitroTKstats}
#'
#' @export format_fup_uc
format_fup_uc <- function(data.in,
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
  data.in <- as.data.frame(data.in)

# Write out a "level 0" file (data the function received it):
  write.table(data.in,
    file=paste(FILENAME,"-fup-UC-Level0.tsv",sep=""),
    sep="\t",
    row.names=F,
    quote=F)

  if (is.null(note.col)) data.in[,"Note"] <- ""

# These arguments allow the user to specify a single value for every obseration
# in the table:
  if (!is.null(cal)) data.in[,cal.col] <- cal
  if (!is.null(dilution)) data.in[,dilution.col] <- dilution
  if (!is.null(istd.name)) data.in[,istd.name.col] <- istd.name
  if (!is.null(istd.conc)) data.in[,istd.conc.col] <- istd.conc
  if (!is.null(std.conc)) data.in[,std.conc.col] <- std.conc
  if (!is.null(uc.assay.conc)) data.in[,uc.assay.conc.col] <- uc.assay.conc
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

  if (!(all(cols %in% colnames(data.in))))
  {
    stop(paste("Missing columns named:",
      paste(cols[!(cols%in%colnames(data.in))],collapse=", ")))
  }

  # Check for sample types we don't know what to do with:
  data.in.badtype <- subset(data.in,!(data.in[,type.col] %in%
                             c("CC","T1","T5","AF")))
  # Write out a "level 0" file identifying those observations we threw out:
  write.table(data.in.badtype,
    file=paste(FILENAME,"-fup-UC-Level0-badtype.tsv",sep=""),
    sep="\t",
    row.names=F,
    quote=F)
  # Only include the data types used:
  data.out <- subset(data.in,data.in[,type.col] %in% c("CC","T1","T5","AF"))
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

  colnames(data.out) <- c(
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
  blanks <- data.out[,"Standard.Conc"]
  blanks[is.na(blanks)] <- -999
  blanks <- blanks == 0
  for (this.chem in unique(data.out[,"DTXSID"]))
  {
    this.subset <- subset(data.out, DTXSID==this.chem)
    for (this.cal in unique(data.out[,"Calibration"]))
    {
      this.cal.subset <- subset(this.subset, Calibration==this.cal)
      if (any(is.na(this.cal.subset[,"ISTD.Area"])))
      {
        this.mean.ISTD <- signif(mean(this.cal.subset$ISTD.Area,na.rm=TRUE))
        which.indices <- data.out[,"DTXSID"] == this.chem &
          data.out[,"Calibration"] == this.cal &
          is.na(data.out[,"ISTD.Area"]) &
          blanks
        data.out[which.indices,
                 "ISTD.Area"] <- this.mean.ISTD
        data.out[which.indices,
                 "Area"] <- 0
      }
    }
  }

  # Set reasonable sig figs:
  for (this.col in c("Area", "ISTD.Area"))
    data.out[,this.col] <- signif(data.out[,this.col], 5)

  # calculate the response:
  data.out[,"Response"] <- signif(as.numeric(data.out[,area.col]) /
     as.numeric(data.out[,istd.col]) * as.numeric(data.out[,istd.conc.col]),4)

# Write out a "level 1" file (data organized into a standard format):
  write.table(data.out,
    file=paste(FILENAME,"-fup-UC-Level1.tsv",sep=""),
    sep="\t",
    row.names=F,
    quote=F)

  summarize_table(data.out,
    req.types=c("CC","T1","T5","AF"))

  return(data.out)
}


