#' Creates a Standardized Data Table with Caco-2 Data (Level-1)
#'
#' This function formats data describing mass spectrometry (MS) peak areas
#' from samples collected as part of in vitro measurement of membrane
#' permeability using Caco-2 cells \insertCite{hubatsch2007determination}{invitroTKstats}.
#' The input data frame is organized into a standard set of columns and is
#' written to a tab-separated text file.
#'
#' In this experiment an
#' in vitro well is separated into two by a membrane composed of a monolayer of
#' Caco-2 cells. A test chemical is added to either the apical or basal side of
#' of the monolayer at time 0, and after a set time samples are taken from both
#' the "donor" (side where the test chemical was added) and the "receiver side.
#' Depending on the direction of the test the donor side can be either apical or
#' basal.
#'
#' The data frame of observations should be annotated according to direction
#' (either apical to basal -- "AtoB" -- or basal to apical -- "BtoA") and type
#' of concentration measured:
#' \tabular{rr}{
#'   Blank with no chemical added \tab Blank \cr
#'   Dosing vehicle (C0) at target concentration \tab D0\cr
#'   Donor compartment at end of experiment \tab D2\cr
#'   Receiver compartment at end of experiment\tab R2\cr
#' }
#'
#' Chemical concentration is calculated qualitatively as a response and 
#' returned as a column in the output data frame:
#'
#' Response <- AREA / ISTD.AREA * ISTD.CONC
#'
#' @param FILENAME (Character) A string used to identify the output Level-1 file.
#' "<FILENAME>-Caco-2-Level1.tsv". (Defaults to "MYDATA".) 
#'
#' @param data.in (Data Frame) A Level-0 data frame containing
#' mass-spectrometry peak areas, indication of chemical identity,
#' and measurement type. The data frame should
#' contain columns with names specified by the following arguments:
#'
#' @param sample.col (Character) Column name of \code{data.in} containing the unique mass
#' spectrometry (MS) sample name used by the laboratory. (Defaults to
#' "Lab.Sample.Name".)
#'
#' @param lab.compound.col (Character) Column name of \code{data.in} containing the test compound
#' name used by the laboratory. (Defaults to "Lab.Compound.Name".)
#'
#' @param dtxsid.col (Character) Column name of \code{data.in} containing EPA's DSSTox Structure
#' ID (\url{http://comptox.epa.gov/dashboard}). (Defaults to "DTXSID".)
#'
#' @param date.col (Character) Column name of \code{data.in} containing the laboratory measurement
#' date. (Defaults to "Date".)
#' 
#' @param series (Numeric) Index of simultaneous replicates with the same analytical chemistry. 
#' (Defaults to \code{NULL}.) (Note: Single entry only, use only if all tested compounds 
#' use the same number of replicates.)
#' 
#' @param series.col (Character) Column name containing \code{series} information.
#' (Defaults to "Series".) 
#' (Note: \code{data.in} does not necessarily have this field. 
#' If this field is missing, it can be auto-filled with the value 
#' specified in \code{series}.)
#'
#' @param compound.col (Character) Column name of \code{data.in} containing the test compound.
#' (Defaults to "Compound.Name".)
#'
#' @param area.col (Character) Column name of \code{data.in} containing the target analyte (that
#' is, the test compound) MS peak area. (Defaults to "Area".)
#' 
#' @param istd.col (Character) Column name of \code{data.in} containing the
#' MS peak area for the internal standard. (Defaults to "ISTD.Area".)
#'
#' @param type.col (Character) Column name of \code{data.in} containing the sample type (see table
#' under Details). (Defaults to "Type".)
#'
#' @param direction.col (Character) Column name of \code{data.in} containing the direction of
#' the Caco-2 permeability experiment: either apical donor to basal receiver (AtoB), or 
#' basal donor to apical receiver (BtoA). (Defaults to "Direction".)
#' 
#' @param membrane.area (Numeric) The area of the Caco-2 monolayer (in cm^2). 
#' (Defaults to \code{NULL}.) (Note: Single entry only, use only if all tested compounds 
#' have the same area for the Caco-2 monolayer.)
#' 
#' @param membrane.area.col (Character) Column name containing \code{membrane.area} 
#' information. (Defaults to "Membrane.Area".) (Note: \code{data.in} does not
#' necessarily have this field. If this field is missing, it can be auto-filled with the value 
#' specified in \code{membrane.area}.)
#'
#' @param receiver.vol.col (Character) Column name of \code{data.in} containing the volume
#' (in cm^3) of the receiver portion of the Caco-2 experimental well. 
#' (Defaults to "Vol.Receiver".)
#'
#' @param donor.vol.col (Character) Column name of \code{data.in} containing the volume
#' (in cm^3) of the donor portion of the Caco-2 experimental well where the
#' test chemical is added. (Defaults to "Vol.Donor".)
#' 
#' @param compound.conc (Numeric) The concentration
#' of the test chemical for calibration curves. (Defaults to \code{NULL}.) (Note: Single entry only, 
#' use only if all test compounds have the same concentration for calibration curves.) 
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
#' @param time (Numeric) The amount of time (in hours) before the receiver and donor 
#' compartments are measured. (Defaults to \code{NULL}.)
#' 
#' @param time.col (Character) Column name containing \code{meas.time} 
#' information. (Defaults to "Time".) (Note: \code{data.in} does not
#' necessarily have this field. If this field is missing, it can be auto-filled with the value 
#' specified in \code{meas.time}.)
#' 
#' @param istd.name (Character) The identity of the internal standard. (Defaults to \code{NULL}.) 
#' (Note: Single entry only, use only if all tested compounds use the same internal standard.)
#' 
#' @param istd.name.col (Character) Column name containing \code{istd.name} information. (Defaults to "ISTD.Name".) 
#' (Note: \code{data.in} does not necessarily have this field. If this field is missing, 
#' it can be auto-filled with the value specified in \code{istd.name}.)
#'
#' @param istd.conc (Numeric) The concentration for the internal standard. (Defaults to \code{NULL}.) 
#' (Note: Single entry only, use only if all tested compounds have the same 
#' internal standard concentration.) 
#'
#' @param istd.conc.col (Character) Column name containing \code{istd.conc} information. (Defaults to "ISTD.Conc".) 
#' (Note: \code{data.in} does not necessarily have this field. If this field is missing, 
#' it can be auto-filled with the value specified in \code{istd.conc}.)
#'
#' @param nominal.test.conc (Numeric) The test chemical concentration in the dosing solution
#' that is added to the donor side at time zero. (Defaults to \code{NULL}.)
#' (Note: Single entry only, use only if the same initial 
#' concentration was used for all tested compounds.)
#'
#' @param nominal.test.conc.col (Character) Column name containing \code{nominal.test.conc} 
#' information. (Defaults to "Test.Target.Conc".) (Note: \code{data.in} does not
#' necessarily have this field. If this field is missing, it can be auto-filled with the value 
#' specified in \code{nominal.test.conc}.)
#'
#' @param analysis.method (Character) The analytical chemistry analysis method, 
#' typically "LCMS" or "GCMS", liquid chromatography or gas chromatography–mass
#' spectrometry, respectively. (Defaults to \code{NULL}.)
#' (Note: Single entry only, use only if the same method was used for all tested compounds.)
#'
#' @param analysis.method.col (Character) Column name containing \code{analysis.method} 
#' information. (Defaults to "Analysis.Method".) (Note: \code{data.in} does not
#' necessarily have this field. If this field is missing, it can be auto-filled with the value 
#' specified in \code{analysis.method}.)
#'
#' @param analysis.instrument (Character) The instrument used for chemical analysis, 
#' for example "Agilent 6890 GC with model 5973 MS". (Defaults to \code{NULL}.) 
#' (Note: Single entry only, use only if the same instrument was used for all tested compounds.) 
#'
#' @param analysis.instrument.col (Character) Column name containing \code{analysis.instrument} 
#' information. (Defaults to "Analysis.Instrument".) (Note: \code{data.in} does not
#' necessarily have this field. If this field is missing, it can be auto-filled with the value 
#' specified in \code{analysis.instrument}.)
#' 
#' @param analysis.parameters (Character) The parameters used to identify the 
#' compound on the chemical analysis instrument, for example
#' "Negative Mode, 221.6/161.6, -DPb=26, FPc=-200, EPd=-10, CEe=-20, CXPf=-25.0". (Defaults to \code{NULL}.) 
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
#' standardized set of columns and column names with membrane permeability data
#' from a Caco-2 assay.
#'
#' @author John Wambaugh
#'
#' @examples
#' library(invitroTKstats)
#' level0 <- TO1caco2
#' level1 <- format_caco2(level0,
#'                        FILENAME="EPACyprotex2021",
#'                        sample.col="SampleName",
#'                        dtxsid.col="CompoundName",
#'                        lab.compound.col="CompoundName",
#'                        cal=1,
#'                        istd.conc.col="ISTD.Conc",
#'                        compound.col="CompoundName",
#'                        compound.conc.col="Test.Target.Conc",
#'                        membrane.area=0.11,
#'                        series=1,
#'                        analysis.parameters="Feature",
#'                        analysis.instrument="GC or LC",
#'                        analysis.method="Mass Spec"
#'                       )
#'
#' @references
#' \insertRef{hubatsch2007determination}{invitroTKstats}
#'
#' @import Rdpack
#'
#' @export format_caco2
format_caco2 <- function(
  FILENAME = "MYDATA",
  data.in,
  sample.col="Lab.Sample.Name",
  lab.compound.col="Lab.Compound.Name",
  dtxsid.col="DTXSID",
  date.col="Date",
  series=NULL,
  series.col="Series",
  compound.col="Compound.Name",
  area.col="Area",
  istd.col="ISTD.Area",
  type.col="Type",
  direction.col="Direction",
  membrane.area=NULL,
  membrane.area.col="Membrane.Area",
  receiver.vol.col="Vol.Receiver",
  donor.vol.col="Vol.Donor",
  compound.conc=NULL,
  compound.conc.col="Nominal.Conc",
  cal=NULL,
  cal.col="Cal",
  dilution=NULL,
  dilution.col="Dilution.Factor",
  time = NULL,
  time.col="Time",
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
  # These are the required data types as indicated by type.col.
  # In order to calculate the parameter a chemical must have peak areas for each
  # of these measurements:
  
  if (!missing(data.in)) {
    data.out <- as.data.frame(data.in)
    # Force code to throw error if data.in accessed after this point:
    rm(data.in)
    } else if (!is.null(INPUT.DIR)) {
    data.out <- read.csv(file=paste0(INPUT.DIR, "/", FILENAME,"-Caco-2-Level0.tsv"),
                         sep="\t",header=T)
    } else {
    data.out <- read.csv(file=paste0(FILENAME,"-Caco-2-Level0.tsv"),
                         sep="\t",header=T)
    }
  
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
  if (!is.null(cal)) data.out[,cal.col] <- cal
  if (!is.null(dilution)) data.out[,dilution.factor.col] <- dilution
  if (!is.null(istd.name)) data.out[,istd.name.col] <- istd.name
  if (!is.null(istd.conc)) data.out[,istd.conc.col] <- istd.conc
  if (!is.null(nominal.test.conc)) data.out[,nominal.test.conc.col] <-
    nominal.test.conc
  if (!is.null(analysis.method)) data.out[,analysis.method.col]<- analysis.method
  if (!is.null(analysis.instrument)) data.out[,analysis.instrument.col] <-
    analysis.instrument
  if (!is.null(analysis.parameters)) data.out[,analysis.parameters.col] <-
    analysis.parameters
  if (!is.null(membrane.area)) data.out[,membrane.area.col] <- membrane.area
  if (!is.null(series)) data.out[,series.col] <- series
  if (!is.null(time)) data.out[,time.col] <- time
  if (!is.null(compound.conc)) data.out[,compound.conc.col] <- compound.conc
  
  caco2.cols <- c(L1.common.cols, 
                  series.col="Series",
                  time.col = "Time",
                  direction.col="Direction",
                  compound.conc.col="Nominal.Conc",
                  nominal.test.conc.col="Test.Target.Conc",
                  membrane.area.col="Membrane.Area",
                  receiver.vol.col="Vol.Receiver",
                  donor.vol.col="Vol.Donor"
  )


  cols <- unlist(mget(names(caco2.cols)))
  if (!(all(cols %in% colnames(data.out))))
  {
    stop(paste("Missing columns named:",
      paste(cols[!(cols%in%colnames(data.out))],collapse=", ")))
  }

  # Only include the data types used:
  req.types=c("Blank","D0","D2","R2")
  data.out <- subset(data.out,data.out[,type.col] %in% req.types)
  data.in.badtype <- subset(data.out,!(data.out[,type.col] %in% req.types))
  
  # Option to export data with bad types
  if (nrow(data.in.badtype) != 0) {
    if (save.bad.types) {
      write.table(data.in.badtype,
                file=paste0(file.path, "/", FILENAME,"-Caco-2-Level0-badtype.tsv"),
                sep="\t",
                row.names=F,
                quote=F)
      cat(paste0("Data with inappropriate sample types were removed. Removed samples were exported to ",
                 FILENAME,"-Caco-2-Level0-badtype.tsv", " in the following directory: ", file.path), "\n")
    } else {
      warning("Data with inappropriate sample types were removed.")
    }
  }
  
                  
  # Organize the columns:
  data.out <- data.out[,cols]

  colnames(data.out) <- caco2.cols

  # calculate the response:
  data.out[,"Area"] <- signif(as.numeric(data.out[,"Area"]), 5)
  data.out[,"ISTD.Area"] <- signif(as.numeric(data.out[,"ISTD.Area"]), 5)
  data.out[,"ISTD.Conc"] <- as.numeric(data.out[,"ISTD.Conc"])
  data.out[,"Response"] <- signif(data.out[,"Area"] /
                                    data.out[,"ISTD.Area"] *  data.out[,"ISTD.Conc"], 4)
  
  if (output.res) {
    # Write out a "level 1" file (data organized into a standard format):
    write.table(data.out,
                file=paste0(file.path, "/", FILENAME,"-Caco-2-Level1.tsv"),
                sep="\t",
                row.names=F,
                quote=F)
    cat(paste0("A Level-1 file named ",FILENAME,"-Caco-2-Level1.tsv", 
               " has been exported to the following directory: ", file.path), "\n")
  }

  summarize_table(data.out,
    req.types=req.types)

  return(data.out)
}


