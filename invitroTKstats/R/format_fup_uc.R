#' Creates a Standardized Data Frame with Ultracentrifugation (UC)
#' Plasma Protein Binding Data (Level-1)
#'
#' This function formats data describing mass spectrometry (MS) peak areas
#' from samples collected as part of in vitro measurement of chemical fraction
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
#' @param FILENAME (Character) A string used to identify the output Level-1 file.
#' "<FILENAME>-fup-UC-Level1.tsv". (Defaults to "MYDATA".)
#'
#' @param data.in (Data Frame) A Level-0 data frame containing mass-spectrometry
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
#' @param date.col (Character) Column name from \code{data.in} containing the laboratory measurement
#' date. (Defaults to "Date".)
#'
#' @param compound.col (Character) Column name from \code{data.in} containing the test compound.
#' (Defaults to "Compound.Name".)
#'
#' @param area.col (Character) Column name from \code{data.in} containing the target analyte (that
#' is, the test compound) MS peak area. (Defaults to "Area".)
#'
#' @param series.col (Character) Column name from \code{data.in} containing the number of 
#' simultaneous replicates with the same analytical chemistry. 
#' (Defaults to "Series".)
#'
#' @param type.col (Character) Column name from \code{data.in} containing the sample type (see table
#' under Details). (Defaults to "Sample.Type".)
#' 
#' @param std.conc (Numeric) The standard test chemical concentration for 
#' the intrinsic clearance assay. (Defaults to \code{NULL}.) (Note: Single entry only, 
#' use only if the same standard concentration was used for all tested compounds.)
#'
#' @param std.conc.col (Character) Column name containing \code{std.conc} 
#' information. (Defaults to "Standard.Conc".) (Note: \code{data.in} does not
#' necessarily have this field. If this field is missing, it can be
#' auto-filled with the value specified in \code{std.conc}.)
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
#' @param uc.assay.conc (Numeric) The intended initial test chemical
#' concentration in the UC assay in uM. (Defaults to \code{NULL}.)
#' (Note: Single entry only,  use only if the intended initial concentration
#' was the same for all tested compounds.)
#'
#' @param uc.assay.conc.col (Character) Column name containing \code{uc.assay.conc} 
#' information. (Defaults to "UC.Assay.Conc".) (Note: \code{data.in} does not
#' necessarily have this field. If this field is missing, it can be auto-filled
#' with the value specified in \code{uc.assay.conc}.)
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
#' @param level0.file (Character) The Level-0 file from which the \code{data.in}
#' were obtained. (Defaults to \code{NULL}.) (Note: Single entry only, use only
#' if all rows in \code{data.in} were obtained from the same Level-0 file.) 
#' 
#' @param level0.file.col (Character) Column name containing \code{level0.file}
#' information. (Defaults to "Level0.File".) (Note: \code{data.in} does not
#' necessarily have this field. If this field is missing, it can be auto-filled
#' with the value specified in \code{level0.file}.)
#'
#' @param level0.sheet (Character) The specific sheet name of the Level-0 file
#' where \code{data.in} is obtained from, if the Level-0 file is an Excel workbook. 
#' (Defaults to \code{NULL}.) (Note: Single entry only, use only if all rows in
#' \code{data.in} were obtained from the same sheet in the same Level-0 file.) 
#'
#' @param level0.sheet.col (Character) Column name containing \code{level0.sheet}
#' information. (Defaults to "Level0.Sheet".) (Note: \code{data.in} does not
#' necessarily have this field. If this field is missing, it can be auto-filled
#' with the value specified in \code{level0.sheet}.)
#' 
#' @param output.res (Logical) When set to \code{TRUE}, the result 
#' table (Level-1) will be exported the current directory as a .tsv file. 
#' (Defaults to \code{TRUE}.)
#' 
#' @param save.bad.types (Logical) When set to \code{TRUE}, export
#' any data being removed due to having inappropriate sample types. 
#' See the Detail section for the required sample types. 
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
#' standardized set of columns and column names with plasma protein binding
#' (PPB) data from an ultracentrifugation (UC) assay.
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
#' @import Rdpack
#'
#' @export format_fup_uc
format_fup_uc <- function(
  FILENAME = "MYDATA",  
  data.in,
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
    data.in <- read.csv(file=paste0(INPUT.DIR, "/", FILENAME,"-fup-UC-Level0.tsv"),
                        sep="\t",header=T)
    } else {
    data.in <- read.csv(file=paste0(FILENAME,"-fup-UC-Level0.tsv"),
                        sep="\t",header=T)
    }
  
  if (is.null(note.col)) data.in[,"Note"] <- ""
  
  # determine the path for output files 
  if (!is.null(OUTPUT.DIR)) {
    file.path <- OUTPUT.DIR
  } else if (!is.null(INPUT.DIR)) {
    file.path <- INPUT.DIR
  } else {
    file.path <- getwd()
  }

# These arguments allow the user to specify a single value for every observation
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
  req.types=c("CC","T1","T5","AF")
  # Only include the data types used:
  data.out <- subset(data.in,data.in[,type.col] %in% req.types))
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
      cat(paste0("Data with inappropriate sample types are being removed and exported as ",
                 FILENAME,"-fup-UC-Level0-badtype.tsv", " to the following directory: ", file.path), "\n")
    } else {
      warning("Some data with inappropriate sample types are being removed.\n")
    }
  }

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

  if (output.res) {
    # Write out a "level 1" file (data organized into a standard format):
    write.table(data.out,
                file=paste0(file.path, "/", FILENAME,"-fup-UC-Level1.tsv"),
                sep="\t",
                row.names=F,
                quote=F)
    cat(paste0("A Level-1 file named ",FILENAME,"-fup-UC-Level1.tsv", 
                " has been exported to the following directory: ", file.path), "\n")
  }


  summarize_table(data.out,
    req.types=c("CC","T1","T5","AF"))

  return(data.out)
}


