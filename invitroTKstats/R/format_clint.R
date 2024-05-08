#' Creates a Standardized Data Table with Hepatocyte Clearance Data (Level-1)
#'
#' This function formats data describing mass spectrometry (MS) peak areas
#' from samples collected as part of in vitro measurement of chemical stability
#' when incubated with suspended hepatocytes \insertCite{shibata2002prediction}{invitroTKstats}.
#' Disappearance of the chemical over time is assumed to be due to metabolism
#' by the hepatocytes.
#' The input data frame is organized into a standard set of columns and is written
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
#' 
#' Chemical concentration is calculated qualitatively as a response and 
#' returned as a column in the output data frame:
#'
#' Response <- AREA / ISTD.AREA * ISTD.CONC
#'
#' @param FILENAME (Character) A string used to identify the output Level-1 file.
#' "<FILENAME>-Clint-Level1.tsv", and/or used to identify the input Level-0 file,
#' "<FILENAME>-Clint-Level0.tsv" if importing from a .tsv file. (Defaults to "MYDATA").
#'
#' @param data.in (Data Frame) A Level-0 data frame or a matrix containing mass-spectrometry peak areas,
#' indication of chemical identity, and measurement type. The data frame should
#' contain columns with names specified by the following arguments:
#'
#' @param sample.col (Character) Column name of \code{data.in} containing the unique mass
#' spectrometry (MS) sample name used by the laboratory. (Defaults to
#' "Lab.Sample.Name".)
#' 
#' @param date (Numeric) The laboratory measurement date. (Defaults to \code{NULL}.) 
#' (Note: Single entry only, 
#' use only if all data were collected on the same date.)
#'
#' @param date.col (Character) Column name containing \code{date} information. (Defaults to "Date".)
#' (Note: \code{data.in} does not
#' necessarily have this field. If this field is missing, it can be auto-filled with the value 
#' specified in \code{date}.)
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
#' @param density (Numeric) The density (units of
#' millions of hepatocytes per mL) hepatocytes in the in vitro incubation. 
#' (Defaults to \code{NULL}.) (Note: Single entry only, 
#' use only if all tested compounds have the same density.)
#' 
#' @param density.col (Character) Column name containing \code{density} 
#' information. (Defaults to "Hep.Density".) (Note: \code{data.in} does not
#' necessarily have this field. If this field is missing, it can be auto-filled with the value 
#' specified in \code{density}.)
#' 
#' @param compound.conc (Numeric) The concentration
#' of the test chemical for calibration curves. (Defaults to \code{NULL}.) (Note: Single entry only,
#' use only if all tested compounds have the same concentration for calibration curves.)
#' 
#' @param compound.conc.col (Character) Column name containing \code{compound.conc}
#' information. (Defaults to "Nominal.Conc".) (Note: \code{data.in} does not
#' necessarily have this field. If this field is missing, it can be auto-filled with the value
#' specified in \code{compound.conc}.)
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
#' @param dilution (Numeric) Number of times the sample was diluted before MS 
#' analysis. (Defaults to \code{NULL}.) (Note: Single entry only, use only if all 
#' samples underwent the same number of dilutions.)
#' 
#' @param dilution.col (Character) Column name containing \code{dilution} 
#' information. (Defaults to "Dilution.Factor".) (Note: \code{data.in} does not
#' necessarily have this field. If this field is missing, it can be auto-filled with the value 
#' specified in \code{dilution}.)
#' 
#' @param time (Numeric) Time of the measurement (in minutes) since the test 
#' chemicals was introduced into the hepatocyte incubation. (Defaults to \code{NULL}.)
#' (Note: Single entry only, use only if all measurements were taken after the same amount of time.)
#' 
#' @param time.col (Character) Column name containing \code{time} 
#' information. (Defaults to "Time".) (Note: \code{data.in} does not
#' necessarily have this field. If this field is missing, it can be auto-filled with the value 
#' specified in \code{time}.)
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
#' @param test.conc (Numeric) The standard test chemical concentration for 
#' the intrinsic clearance assay. (Defaults to \code{NULL}.) (Note: Single entry only, 
#' use only if the same standard concentration was used for all tested compounds.)
#' 
#' @param test.conc.col (Character) Column name containing \code{test.conc} 
#' information. (Defaults to "Test.Compound.Conc".) (Note: \code{data.in} does not
#' necessarily have this field. If this field is missing, it can be auto-filled with the value 
#' specified in \code{test.conc}.)
#' 
#' @param clint.assay.conc (Numeric) The initial test chemical concentration for 
#' the intrinsic clearance assay. (Defaults to \code{NULL}.) (Note: Single entry only, 
#' use only if the same initial concentration was used for all tested compounds.)
#' 
#' @param clint.assay.conc.col (Character) Column name containing \code{clint.assay.conc} 
#' information. (Defaults to "Clint.Assay.Conc".) (Note: \code{data.in} does not
#' necessarily have this field. If this field is missing, it can be auto-filled with the value 
#' specified in \code{clint.assay.conc}.)
#' 
#' @param area.col (Character) Column name of \code{data.in} containing the target analyte (that
#' is, the test compound) MS peak area. (Defaults to "Area".)
#' 
#' @param biological.replicates (Character) Replicates with the same analyte. Typically, this uses 
#' numbers or letters to index. (Defaults to \code{NULL}.) (Note: Single entry only, 
#' use only if none of the test compounds have replicates.)
#' 
#' @param biological.replicates.col (Character) Column name of \code{data.in} containing the number or 
#' the indices of replicates with the same analyte. (Defaults to "Biological.Replicates".)
#' (Note: \code{data.in} does not necessarily have this field. If this field is missing, it can be auto-filled with the value 
#' specified in \code{biological.replicates}.)
#' 
#' @param technical.replicates (Character) Repeated measurements from one sample. Typically, this uses 
#' numbers or letters to index. (Defaults to \code{NULL}.) (Note: Single entry only, 
#' use only if none of the test compounds have replicates.)
#' 
#' @param technical.replicates.col (Character) Column name of \code{data.in} containing the number or 
#' the indices of replicates taken from the one sample. (Defaults to "Technical.Replicates".) 
#' (Note: \code{data.in} does not necessarily have this field. If this field is missing, it can be auto-filled with the value 
#' specified in \code{technical.replicates}.)
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
#' for example "Waters Xevo TQ-S micro (QEB0036)". (Defaults to \code{NULL}.) 
#' (Note: Single entry only, use only if the same instrument was used for all tested compounds.) 
#' 
#' @param analysis.instrument.col (Character) Column name containing \code{analysis.instrument} 
#' information. (Defaults to "Analysis.Instrument".) (Note: \code{data.in} does not
#' necessarily have this field. If this field is missing, it can be auto-filled with the value 
#' specified in \code{analysis.instrument}.)
#' 
#' @param analysis.parameters (Numeric) The parameters used to identify the 
#' compound on the chemical analysis instrument. (Defaults to \code{NULL}.) 
#' (Note: Single entry only, use only if the same parameters were used for all tested compounds.)
#'
#' @param analysis.parameters.col (Character) Column name containing \code{analysis.parameters} 
#' information. (Defaults to "Analysis.Parameters".) (Note: \code{data.in} does not
#' necessarily have this field. If this field is missing, it can be auto-filled with the value 
#' specified in \code{analysis.parameters}.)
#' 
#' @param note.col (Character) Column name of \code{data.in} containing additional notes on 
#' test compounds. (Defaults to "Note").
#' 
#' @param level0.file (Character) The Level-0 file from which the \code{data.in} were obtained.
#' (Defaults to \code{NULL}.) (Note: Single entry only, use only if all rows in data.in
#' were obtained from the same Level-0 file.) 
#' 
#' @param level0.file.col (Character) Column name containing \code{level0.file} information. 
#' (Defaults to "Level0.File".) (Note: \code{data.in} does not
#' necessarily have this field. If this field is missing, it can be auto-filled with the value 
#' specified in \code{level0.file}.)
#' 
#' @param level0.sheet (Character) The specific sheet name of Level-0 file from which the 
#' \code{data.in} is obtained from, if the Level-0 file is an Excel workbook. 
#' (Defaults to \code{NULL}.) (Note: Single entry only, use only if all rows in \code{data.in}
#' were obtained from the same sheet in the same Level-0 file.) 
#' 
#' @param level0.sheet.col (Character) Column name containing \code{level0.sheet} information.
#' (Defaults to "Level0.Sheet".) (Note: \code{data.in} does not
#' necessarily have this field. If this field is missing, it can be auto-filled with the value 
#' specified in \code{level0.sheet}.)
#' 
#' @param output.res (Logical) When set to \code{TRUE}, the result 
#' table (Level-1) will be exported the current directory as a .tsv file. 
#' (Defaults to \code{TRUE}.)
#' 
#' @param save.bad.types (Logical) When set to \code{TRUE}, export data removed 
#' due to inappropriate sample types. See the Detail section for the required sample types. 
#' (Defaults to \code{FALSE}.)
#' 
#' @param INPUT.DIR (Character) Path to the directory where the input level-0 file exists. 
#' If \code{NULL}, looking for the input level-0 file in the current working
#' directory. (Defaults to \code{NULL}.)
#' 
#' @param OUTPUT.DIR (Character) Path to the directory to save the output file. 
#' If \code{NULL}, the output file will be saved to the current working
#' directory or \code{INPUT.DIR} if specified. (Defaults to \code{NULL}.)
#'
#' @return A Level-1 data frame with a standardized format containing a  
#' standardized set of columns and column names with hepatic
#' clearance data for a variety of chemicals. 
#'
#' @author John Wambaugh
#'
#' @examples
#' ## Load the example level-0 data
#' level0 <- invitroTKstats::clint_L0
#' 
#' ## Run it through level-1 processing function
#' ## This example shows the use of the data.in argument which allows users to pass
#' ## in a data frame from the R session.
#' ## If the input level-0 data exists in an external file such as a .tsv file,
#' ## users may import it using INPUT.DIR to specify the path and FILENAME
#' ## to specify the file name. See documentation for details.
#' level1 <- format_clint(data.in = level0,
#'                        sample.col ="Sample",
#'                        date.col="Date",
#'                        compound.col="Compound",
#'                        lab.compound.col="Lab.Compound.ID",
#'                        type.col="Type",
#'                        dilution.col="Dilution.Factor",
#'                        cal=1,
#'                        istd.conc = 10/1000,
#'                        istd.col= "ISTD.Peak.Area",
#'                        area.col = "Peak.Area",
#'                        density = 0.5,
#'                        clint.assay.conc = 1,
#'                        biological.replicates = 1,
#'                        test.conc.col="Compound.Conc",
#'                        time.col = "Time",
#'                        analysis.method = "LCMS",
#'                        analysis.instrument = "Unknown",
#'                        analysis.parameters.col = "Analysis.Params",
#'                        note="Sample Text",
#'                        output.res = FALSE
#'                        )
#'
#' @references
#' \insertRef{shibata2002prediction}{invitroTKstats}
#'
#' @import Rdpack
#'
#' @export format_clint
format_clint <- function(
  FILENAME = "MYDATA",
  data.in,
  sample.col="Lab.Sample.Name",
  date=NULL,
  date.col="Date",
  compound.col="Compound.Name",
  dtxsid.col="DTXSID",
  lab.compound.col="Lab.Compound.Name",
  type.col="Sample.Type",
  density=NULL,
  density.col="Hep.Density",
  # compound.conc=NULL,
  # compound.conc.col="Nominal.Conc",
  cal=NULL,
  cal.col="Cal",
  dilution=NULL,
  dilution.col="Dilution.Factor",
  time = NULL,
  time.col="Time",
  istd.col="ISTD.Area",
  istd.name=NULL,
  istd.name.col="ISTD.Name",
  istd.conc=NULL,
  istd.conc.col="ISTD.Conc",
  test.conc=NULL,
  test.conc.col="Test.Compound.Conc",
  clint.assay.conc=NULL,
  clint.assay.conc.col="Clint.Assay.Conc",
  area.col="Area",
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
  output.res = TRUE,
  save.bad.types = FALSE,
  INPUT.DIR = NULL,
  OUTPUT.DIR = NULL
  )
{
  
  if (!missing(data.in)) {
    data.in <- as.data.frame(data.in)
  } else if (!is.null(INPUT.DIR)) {
    data.in <- read.csv(file=paste0(INPUT.DIR, "/", FILENAME,"-Clint-Level0.tsv"),
                        sep="\t",header=T)
    } else {
      data.in <- read.csv(file=paste0(FILENAME,"-Clint-Level0.tsv"),
                           sep="\t",header=T)
    }
  
  if (is.null(note.col))
  {
    data.in[,"Note"] <- ""
    note.col <- "Note"
  }

  if (!(test.conc.col %in% colnames(data.in)))
  {
    if (is.null(test.conc))
    {
      data.in[, test.conc.col] <- NA
    } else {
      data.in[, test.conc.col] <- test.conc
    }
  }

  # determine the path for output files
  if (!is.null(OUTPUT.DIR)) {
    file.path <- OUTPUT.DIR
  } else if (!is.null(INPUT.DIR)) {
    file.path <- INPUT.DIR
  } else {
    file.path <- getwd()
  }

# These arguments allow the user to specify a single value for every obseration
# in the table:
  if (!is.null(date)) data.in[,date.col] <- date
  #if (!is.null(compound.conc)) data.in[,compound.conc.col] <- compound.conc
  if (!is.null(time)) data.in[,time.col] <- time
  if (!is.null(cal)) data.in[,cal.col] <- cal
  if (!is.null(dilution)) data.in[,dilution.col] <- dilution
  if (!is.null(density)) data.in[,density.col] <- density
  if (!is.null(istd.name)) data.in[,istd.name.col] <- istd.name
  if (!is.null(istd.conc)) data.in[,istd.conc.col] <- istd.conc
  if (!is.null(test.conc)) data.in[,test.conc.col] <- test.conc
  if (!is.null(clint.assay.conc)) data.in[,clint.assay.conc.col] <-
    clint.assay.conc
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
  clint.cols <- c(L1.common.cols,
                  time.col = "Time",
                  test.conc.col = "Test.Compound.Conc",
                  clint.assay.conc.col = "Clint.Assay.Conc",
                  density.col = "Hep.Density"
                  )
  
  ## allow either one of the two, or both replicate columns in the data
  if (biological.replicates.col %in% colnames(data.in))
    clint.cols <- c(clint.cols, 
                    biological.replicates.col = "Biological.Replicates")
  if (technical.replicates.col %in% colnames(data.in))
    clint.cols <- c(clint.cols, 
                    technical.replicates.col = "Technical.Replicates")
  if (!any(c(biological.replicates.col, technical.replicates.col) %in% colnames(data.in)))
    stop(paste("Missing columns, need to specify/auto-fill least one replicate columns:", 
               paste(c(biological.replicates.col, technical.replicates.col),collapse = ", ")))
  
  cols <- unlist(mget(names(clint.cols)))
  if (!(all(cols %in% colnames(data.in))))
  {
    stop(paste("Missing columns named:",
      paste(cols[!(cols%in%colnames(data.in))],collapse=", ")))
  }

  # Only include the data types used:
  req.types=c("Blank","Cvst","CC","Inactive")
  data.out <- subset(data.in,data.in[,type.col] %in% req.types)
  data.in.badtype <- subset(data.in,!(data.in[,type.col] %in% req.types))

  # Force code to throw error if data.in accessed after this point:
  rm(data.in)

  # Option to export data with bad types
  if (nrow(data.in.badtype) != 0) {
    if (save.bad.types) {
      write.table(data.in.badtype,
                  file=paste0(file.path, "/", FILENAME,"-Clint-Level0-badtype.tsv"),
                  sep="\t",
                  row.names=F,
                  quote=F)
      cat(paste0("Data with inappropriate sample types were removed. Removed samples were exported to ",
                 FILENAME,"-Clint-Level0-badtype.tsv", " in the following directory: ", file.path), "\n")
    } else {
      warning("Data with inappropriate sample types were removed.")
    }
  }

  # Organize the columns:
  data.out <- data.out[,cols]

  colnames(data.out) <- clint.cols

  # Set reasonable significant figures:
  for (this.col in c("Area", "ISTD.Area"))
    data.out[,this.col] <- signif(data.out[,this.col], 5)

  # calculate the response:
  data.out[,"Response"] <- signif(data.out[,"Area"] /
     data.out[,"ISTD.Area"] * data.out[,"ISTD.Conc"], 4)

  if (output.res) {
    # Write out a "level 1" file (data organized into a standard format):
    write.table(data.out,
                file=paste0(file.path, "/", FILENAME,"-Clint-Level1.tsv"),
                sep="\t",
                row.names=F,
                quote=F)
    cat(paste0("A Level-1 file named ",FILENAME,"-Clint-Level1.tsv",
                " has been exported to the following directory: ", file.path), "\n")
  }

  summarize_table(data.out,
    req.types=c("Blank","Cvst"))

  return(data.out)
}


