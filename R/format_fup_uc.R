#' Creates a Standardized Data Frame with Ultracentrifugation (UC)
#' Plasma Protein Binding (PPB) Data (Level-1)
#'
#' This function formats data describing mass spectrometry (MS) peak areas
#' from samples collected as part of \emph{in vitro} measurements of chemical fraction
#' unbound in plasma using ultracentrifugation
#' \insertCite{redgrave1975separation}{invitroTKstats}.
#' The input data frame is organized into a standard set of columns and written
#' to a tab-separated text file.
#'
#' The data frame of observations should be annotated according to
#' these types:
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
#' If the output level-1 result table is chosen to be exported and an output 
#' directory is not specified, it will be exported to the user's R session
#' temporary directory. This temporary directory is a per-session directory 
#' whose path can be found with the following code: \code{tempdir()}. For more 
#' details, see \url{https://www.collinberke.com/til/posts/2023-10-24-temp-directories/}.
#' 
#' As a best practice, \code{INPUT.DIR} and/or \code{OUTPUT.DIR} should be 
#' specified to simplify the process of importing and exporting files. This 
#' practice ensures that the exported files can easily be found and will not be 
#' exported to a temporary directory.
#'
#' @param FILENAME (Character) A string used to identify the output level-1 file,
#' "<FILENAME>-fup-UC-Level1.tsv", and/or used to identify the input level-0 file,
#' "<FILENAME>-fup-UC-Level0.tsv" if importing from a .tsv file. (Defaults to "MYDATA".)
#'
#' @param data.in (Data Frame) A level-0 data frame containing mass-spectrometry
#' peak areas, indication of chemical identity, and measurement type.
#' The data frame should contain columns with names specified by the following arguments:
#'
#' @param sample.col (Character) Column name from \code{data.in} containing the unique mass
#' spectrometry (MS) sample name used by the laboratory. (Defaults to
#' "Lab.Sample.Name".)
#'
#' @param lab.compound.col (Character) Column name from \code{data.in} containing the test compound
#' name used by the laboratory. (Defaults to "Lab.Compound.Name".)
#'
#' @param dtxsid.col (Character) Column name from \code{data.in} containing EPA's DSSTox Structure
#' ID (\url{http://comptox.epa.gov/dashboard}). (Defaults to "DTXSID".)
#'
#' @param date (Character) The laboratory measurement date, format "MMDDYY" where 
#' "MM" = 2 digit month, "DD" = 2 digit day, and "YY" = 2 digit year. (Defaults to \code{NULL}.) 
#' (Note: Single entry only, use only if all data were collected on the same date.)
#'
#' @param date.col (Character) Column name containing \code{date} information. (Defaults to "Date".) (Note: \code{data.in} does not
#' necessarily have this field. If this field is missing, it can be auto-filled with the value 
#' specified in \code{date}.)
#'
#' @param compound.col (Character) Column name from \code{data.in} containing the test compound.
#' (Defaults to "Compound.Name".)
#'
#' @param area.col (Character) Column name from \code{data.in} containing the target analyte (that
#' is, the test compound) MS peak area. (Defaults to "Area".)
#' 
#' @param type.col (Character) Column name from \code{data.in} containing the sample type (see table
#' under Details). (Defaults to "Sample.Type".)
#' 
#' @param test.conc (Numeric) The standard test chemical concentration for 
#' the fup UC assay. (Defaults to \code{NULL}.) (Note: Single entry only, 
#' use only if the same standard concentration was used for all tested compounds.)
#'
#' @param test.conc.col (Character) Column name containing \code{test.conc} 
#' information. (Defaults to Test.Compound.Conc".) (Note: \code{data.in} does not
#' necessarily have this field. If this field is missing, it can be
#' auto-filled with the value specified in \code{test.conc}.)
#' 
#' @param cal (Character) MS calibration the samples were based on. Typically, this uses 
#' indices or dates to represent if the analyses were done on different machines on 
#' the same day or on different days with the same MS analyzer. (Defaults to \code{NULL}.) 
#' (Note: Single entry only, use only if all data were collected based on the
#' same calibration.)
#'
#' @param cal.col (Character) Column name containing \code{cal} 
#' information. (Defaults to "Cal".) (Note: \code{data.in} does not
#' necessarily have this field. If this field is missing, it can be
#' auto-filled with the value specified in \code{cal}.)
#' 
#' @param dilution (Numeric) Number of times the sample was diluted before MS 
#' analysis. (Defaults to \code{NULL}.) (Note: Single entry only, use only if all 
#' samples underwent the same number of dilutions.)
#'
#' @param dilution.col (Character) Column name containing \code{dilution} 
#' information. (Defaults to "Dilution.Factor".) (Note: \code{data.in} does not
#' necessarily have this field. If this field is missing, it can be
#' auto-filled with the value specified in \code{dilution}.)
#'
#' @param istd.col (Character) Column name of \code{data.in} containing the
#' MS peak area for the internal standard. (Defaults to "ISTD.Area".)
#'
#' @param istd.name (Character) The identity of the internal standard.
#' (Defaults to \code{NULL}.) (Note: Single entry only, use only if all
#' tested compounds use the same internal standard.) 
#'
#' @param istd.name.col (Character) Column name containing \code{istd.name} information. 
#' (Defaults to "ISTD.Name".) (Note: \code{data.in} does not necessarily have
#' this field. If this field is missing, it can be auto-filled with the value 
#' specified in \code{istd.name}.)
#'
#' @param istd.conc (Numeric) The concentration for the internal standard.
#' (Defaults to \code{NULL}.) (Note: Single entry only, use only if all
#' tested compounds have the same internal standard concentration.) 
#'
#' @param istd.conc.col (Character) Column name containing \code{istd.conc}
#' information.  (Defaults to "ISTD.Conc".) (Note: \code{data.in} does not
#' necessarily have this field. If this field is missing, it can be
#' auto-filled with the value specified in \code{istd.conc}.)
#' 
#' @param test.nominal.conc (Numeric) The nominal concentration added to the UC assay
#' at time 0. (Defaults to \code{NULL}.)
#' (Note: Single entry only,  use only if all tested compounds used the same 
#' concentration at time 0.)
#'
#' @param test.nominal.conc.col (Character) Column name containing \code{test.nominal.conc} 
#' information. (Defaults to "Test.Target.Conc".) (Note: \code{data.in} does not
#' necessarily have this field. If this field is missing, it can be auto-filled
#' with the value specified in \code{test.nominal.conc}.)
#' 
#' @param biological.replicates (Character) Replicates with the same analyte. Typically, this uses 
#' numbers or letters to index. (Defaults to \code{NULL}.) (Note: Single entry only, 
#' use only if none of the test compounds have replicates.)
#' 
#' @param biological.replicates.col (Character) Column name of \code{data.in} containing the number or 
#' the indices of replicates with the same analyte. (Defaults to "Biological.Replicates".) 
#' (Note: \code{data.in} does not necessarily have this field. If this field is missing, it can be auto-filled
#' with the value specified in \code{biological.replicates}.)
#' 
#' @param technical.replicates (Character) Repeated measurements from one sample. Typically, this uses 
#' numbers or letters to index. (Defaults to \code{NULL}.) (Note: Single entry only, 
#' use only if none of the test compounds have replicates.)
#' 
#' @param technical.replicates.col (Character) Column name of \code{data.in} containing the number or 
#' the indices of replicates taken from the one sample. (Defaults to "Technical.Replicates".) 
#' (Note: \code{data.in} does not necessarily have this field. If this field is missing, it can be auto-filled
#' with the value specified in \code{technical.replicates}.)
#'
#' @param analysis.method (Character) The analytical chemistry analysis method, 
#' typically "LCMS" or "GCMS", liquid chromatography or gas chromatography–mass
#' spectrometry, respectively. (Defaults to \code{NULL}.) (Note: Single entry only, 
#' use only if the same method was used for all tested compounds.)
#'
#' @param analysis.method.col (Character) Column name containing \code{analysis.method} 
#' information. (Defaults to "Analysis.Method".) (Note: \code{data.in} does not
#' necessarily have this field. If this field is missing, it can be auto-filled
#' with the value specified in \code{analysis.method}.)
#' 
#' @param analysis.instrument (Character) The instrument used for chemical analysis, 
#' for example "Waters Xevo TQ-S micro (QEB0036)". (Defaults to \code{NULL}.) 
#' (Note: Single entry only, use only if the same instrument was used for all
#' tested compounds.) 
#'
#' @param analysis.instrument.col (Character) Column name containing
#' \code{analysis.instrument} information. (Defaults to "Analysis.Instrument".)
#' (Note: \code{data.in} does not necessarily have this field. If this field
#' is missing, it can be auto-filled with the value specified in
#' \code{analysis.instrument}.)
#'
#' @param analysis.parameters (Character) The parameters used to identify the 
#' compound on the chemical analysis instrument. (Defaults to \code{NULL}.) 
#' (Note: Single entry only, use only if the same parameters were used for all
#' tested compounds.) 
#'
#' @param analysis.parameters.col (Character) Column name containing
#' \code{analysis.parameters} information. (Defaults to "Analysis.Parameters".)
#' (Note: \code{data.in} does not necessarily have this field. If this field
#' is missing, it can be auto-filled with the value specified in
#' \code{analysis.parameters}.)
#' 
#' @param note.col (Character) Column name of \code{data.in} containing
#' additional notes on the test compounds. (Defaults to "Note").
#'
#' @param level0.file (Character) The level-0 file from which the \code{data.in}
#' were obtained. (Defaults to \code{NULL}.) (Note: Single entry only, use only
#' if all rows in \code{data.in} were obtained from the same level-0 file.) 
#' 
#' @param level0.file.col (Character) Column name containing \code{level0.file}
#' information. (Defaults to "Level0.File".) (Note: \code{data.in} does not
#' necessarily have this field. If this field is missing, it can be auto-filled
#' with the value specified in \code{level0.file}.)
#'
#' @param level0.sheet (Character) The specific sheet name of the level-0 file
#' where \code{data.in} is obtained from, if the level-0 file is an Excel workbook. 
#' (Defaults to \code{NULL}.) (Note: Single entry only, use only if all rows in
#' \code{data.in} were obtained from the same sheet in the same level-0 file.) 
#'
#' @param level0.sheet.col (Character) Column name containing \code{level0.sheet}
#' information. (Defaults to "Level0.Sheet".) (Note: \code{data.in} does not
#' necessarily have this field. If this field is missing, it can be auto-filled
#' with the value specified in \code{level0.sheet}.)
#' 
#' @param output.res (Logical) When set to \code{TRUE}, the result 
#' table (level-1) will be exported to the user's per-session temporary directory
#' or \code{OUTPUT.DIR} (if specified) as a .tsv file. 
#' (Defaults to \code{FALSE}.)
#' 
#' @param save.bad.types (Logical) When set to \code{TRUE}, export data removed 
#' due to inappropriate sample types. See the Detail section for the required sample types. 
#' (Defaults to \code{FALSE}.)
#' 
#' @param sig.figs (Numeric) The number of significant figures to round the exported result table (level-1). 
#' (Defaults to \code{5}.)
#' 
#' @param INPUT.DIR (Character) Path to the directory where the input level-0 file exists. 
#' If \code{NULL}, looking for the input level-0 file in the current working
#' directory. (Defaults to \code{NULL}.)
#' 
#' @param OUTPUT.DIR (Character) Path to the directory to save the output file. 
#' If \code{NULL}, the output file will be saved to the user's per-session temporary
#' directory or \code{INPUT.DIR} if specified. (Defaults to \code{NULL}.)
#'
#' @return A level-1 data frame with a standardized format containing a  
#' standardized set of columns and column names with plasma protein binding
#' (PPB) data from an ultracentrifugation (UC) assay.
#'
#' @author John Wambaugh
#'
#' @examples
#' 
#' ## Load the example level-0 data
#' level0 <- invitroTKstats::fup_uc_L0
#' 
#' ## Run it through level-1 processing function
#' ## This example shows the use of data.in argument which allows users to pass
#' ## in a data frame from the R session.
#' ## If the input level-0 data exists in an external file such as a .tsv file,
#' ## users may import it using INPUT.DIR to specify the path and FILENAME
#' ## to specify the file name. See documentation for details.
#' level1 <- format_fup_uc(data.in = level0,
#'                         sample.col="Sample",
#'                         compound.col="Compound",
#'                         test.conc.col ="Compound.Conc", 
#'                         lab.compound.col="Lab.Compound.ID", 
#'                         type.col="Sample.Type", 
#'                         istd.col="ISTD.Peak.Area",
#'                         cal.col = "Date",
#'                         area.col = "Peak.Area",
#'                         istd.conc = 1,
#'                         note.col = NULL,
#'                         test.nominal.conc = 10,
#'                         analysis.method = "UPLC-MS/MS",
#'                         analysis.instrument = "Waters Xevo TQ-S micro (QEB0036)",
#'                         analysis.parameters.col = "Analysis.Params",
#'                         technical.replicates.col = "Replicate",
#'                         biological.replicates = 1,
#'                         output.res = FALSE
#'                         )
#'
#' @references
#' \insertRef{redgrave1975separation}{invitroTKstats}
#'
#' @import Rdpack
#'
#' @export format_fup_uc
format_fup_uc <- function(
  FILENAME = "MYDATA",  
  data.in,
  sample.col="Lab.Sample.Name",
  lab.compound.col="Lab.Compound.Name",
  dtxsid.col="DTXSID",
  date=NULL,
  date.col="Date",
  compound.col="Compound.Name",
  area.col="Area",
  type.col="Sample.Type",
  test.conc=NULL,
  test.conc.col="Test.Compound.Conc",
  cal=NULL,
  cal.col="Cal",
  dilution=NULL,
  dilution.col="Dilution.Factor",
  istd.col="ISTD.Area",
  istd.name=NULL,
  istd.name.col="ISTD.Name",
  istd.conc=NULL,
  istd.conc.col="ISTD.Conc",
  test.nominal.conc=NULL,
  test.nominal.conc.col="Test.Target.Conc",
  biological.replicates = NULL,
  biological.replicates.col = "Biological.Replicates",
  technical.replicates = NULL,
  technical.replicates.col = "Technical.Replicates",
  analysis.method=NULL,
  analysis.method.col="Analysis.Method",
  analysis.instrument=NULL,
  analysis.instrument.col="Analysis.Instrument",
  analysis.parameters=NULL,
  analysis.parameters.col="Analysis.Parameters",
  note.col="Note",
  level0.file=NULL,
  level0.file.col="Level0.File",
  level0.sheet=NULL,
  level0.sheet.col="Level0.Sheet",
  output.res = FALSE,
  save.bad.types = FALSE,
  sig.figs = 5, 
  INPUT.DIR = NULL,
  OUTPUT.DIR = NULL
  )
{

  #assigning global variables
  DTXSID <- Calibration <- NULL
  
  if (!missing(data.in)) {
    data.in <- as.data.frame(data.in)
  } else if (!is.null(INPUT.DIR)) {
    data.in <- read.csv(file=paste0(INPUT.DIR, "/", FILENAME,"-fup-UC-Level0.tsv"),
                        sep="\t",header=T)
    } else {
    data.in <- read.csv(file=paste0(FILENAME,"-fup-UC-Level0.tsv"),
                        sep="\t",header=T)
    }
  
  if (is.null(note.col)) {
    data.in[,"Note"] <- ""
    note.col <- "Note"
  }
  
  # determine the path for output files 
  if (!is.null(OUTPUT.DIR)) {
    file.path <- OUTPUT.DIR
  } else if (!is.null(INPUT.DIR)) {
    file.path <- INPUT.DIR
  } else {
    file.path <- tempdir()
  }

# These arguments allow the user to specify a single value for every observation
# in the table:
  if (!is.null(date)){
    # if numeric, convert to string and ensuring leading zero is kept for single digit months
    if (is.numeric(date)) date <- base::sprintf("%06d", date)
    data.in[,date.col] <- date
  }
  if (!is.null(cal)) data.in[,cal.col] <- cal
  if (!is.null(dilution)) data.in[,dilution.col] <- dilution
  if (!is.null(istd.name)) data.in[,istd.name.col] <- istd.name
  if (!is.null(istd.conc)) data.in[,istd.conc.col] <- istd.conc
  if (!is.null(test.conc)) data.in[,test.conc.col] <- test.conc
  if (!is.null(test.nominal.conc)) data.in[,test.nominal.conc.col] <- test.nominal.conc
  if (!is.null(analysis.method)) data.in[,analysis.method.col]<- analysis.method
  if (!is.null(analysis.instrument)) data.in[,analysis.instrument.col] <-
    analysis.instrument
  if (!is.null(analysis.parameters)) data.in[,analysis.parameters.col] <-
    analysis.parameters
  if (!is.null(level0.file)) data.in[,level0.file.col] <- level0.file
  if (!is.null(level0.sheet)) data.in[,level0.sheet.col] <- level0.sheet
  if (!is.null(biological.replicates)) data.in[,biological.replicates.col]<- biological.replicates
  if (!is.null(technical.replicates)) data.in[,technical.replicates.col]<- technical.replicates

  # We need all these columns in data.in
  fup.uc.cols <- c(L1.common.cols,
                   test.conc.col = "Test.Compound.Conc",
                   test.nominal.conc.col = "Test.Nominal.Conc"
  )
  
  ## allow either one of the two, or both replicate columns in the data
  if (biological.replicates.col %in% colnames(data.in))
    fup.uc.cols <- c(fup.uc.cols, 
                    biological.replicates.col = "Biological.Replicates")
  if (technical.replicates.col %in% colnames(data.in))
    fup.uc.cols <- c(fup.uc.cols, 
                    technical.replicates.col = "Technical.Replicates")
  if (!any(c(biological.replicates.col, technical.replicates.col) %in% colnames(data.in)))
    stop(paste("Missing columns, need to specify/auto-fill least one replicate columns:", 
               paste(c(biological.replicates.col, technical.replicates.col),collapse = ", ")))
  
  cols <- unlist(mget(names(fup.uc.cols)))
  if (!(all(cols %in% colnames(data.in))))
  {
    stop(paste("Missing columns named:",
      paste(cols[!(cols%in%colnames(data.in))],collapse=", ")))
  }

  # Check for sample types we don't know what to do with:
  req.types=c("CC","T1","T5","AF")
  # Only include the data types used:
  data.out <- subset(data.in,data.in[,type.col] %in% req.types)
  data.in.badtype <- subset(data.in,!(data.in[,type.col] %in% req.types))
  
  # Force code to throw error if data.in accessed after this point:
  rm(data.in)
  
  # Option to export data with bad types
  if (nrow(data.in.badtype) != 0) {
    if (save.bad.types) {
      write.table(data.in.badtype,
                  file=paste0(file.path, "/", FILENAME,"-fup-UC-Level0-badtype.tsv"),
                  sep="\t",
                  row.names=F,
                  quote=F)
      cat(paste0("Data with inappropriate sample types were removed. Removed samples were exported to ",
                 FILENAME,"-fup-UC-Level0-badtype.tsv", " in the following directory: ", file.path), "\n")
    } else {
      warning("Data with inappropriate sample types were removed.")
    }
  }

  # Organize the columns:
  data.out <- data.out[,cols]

  colnames(data.out) <- fup.uc.cols

  # Blanks don't always have internal standard -- add average ISTD.Area
  # First identify the blanks (have to deal with NA Test.Compound.Conc:
  blanks <- data.out[,"Test.Compound.Conc"]
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
        this.mean.ISTD <- mean(this.cal.subset$ISTD.Area,na.rm=TRUE)
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

  # calculate the response:
  data.out[,"Response"] <- as.numeric(data.out[,"Area"]) /
     as.numeric(data.out[,"ISTD.Area"]) * as.numeric(data.out[,"ISTD.Conc"])

  if (output.res) {
    
    rounded.data.out <- data.out 
    
    # Round results to desired number of sig figs 
    if (!is.null(sig.figs)){
      rounded.data.out[,"Area"] <- signif(rounded.data.out[,"Area"], sig.figs)
      rounded.data.out[,"ISTD.Area"] <- signif(rounded.data.out[,"ISTD.Area"], sig.figs)
      rounded.data.out[,"Response"] <- signif(rounded.data.out[,"Response"], sig.figs)
      cat(paste0("\nData to export has been rounded to ", sig.figs, " significant figures.\n"))
    }
    
    # Write out a "level-1" file (data organized into a standard format):
    write.table(rounded.data.out,
                file=paste0(file.path, "/", FILENAME,"-fup-UC-Level1.tsv"),
                sep="\t",
                row.names=F,
                quote=F)
    cat(paste0("A level-1 file named ",FILENAME,"-fup-UC-Level1.tsv", 
                " has been exported to the following directory: ", file.path), "\n")
  }


  summarize_table(data.out,
    req.types=c("CC","T1","T5","AF"))

  return(data.out)
}


