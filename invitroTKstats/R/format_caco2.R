#' Creates a Standardized Data Table Reporting Caco-2 Data (Level-1)
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
#' @param data.in (Data Frame) A data frame containing mass-spectrometry peak areas,
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
#' @param series.col (Character) Column name containing `series` information. (Defaults to "Series".)
#'
#' @param series (Numeric) Index of simultaneous replicates with the same analytical chemistry. 
#' (Defaults to \code{NULL}.) (Note: Single entry only, use 1 only if all test compounds have one replicate.) 
#'
#' @param compound.col (Character) Column name of data.in containing the test compound.
#' (Defaults to "Compound.Name".)
#'
#' @param area.col (Character) Column name of data.in containing the target analyte (that
#' is, the test compound) MS peak area. (Defaults to "Area".)
#'
#' @param type.col (Character) Column name of data.in containing the sample type (see table
#' under Details). (Defaults to "Type".)
#'
#' @param direction.col (Character) Column name of data.in containing the direction of
#' the Caco-2 permeability experiment: either apical donor to basal receiver (AtoB), or 
#' basal donor to apical receiver (BtoA). (Defaults to "Direction".)
#'
#' @param cal.col (Character) Column name containing `cal` 
#' information. (Defaults to "Cal".)
#'
#' @param cal (Character) MS calibration the samples were based on, typically uses 
#' indices or dates to represent if the analyses were done on different machines on 
#' the same day or on different days with the same MS analyzer. (Defaults to \code{NULL}.) 
#' (Note: Single entry only, 
#' use only if all data were collected based on the same calibration.)
#' 
#' @param compound.conc.col (Character) Column name containing `compound.conc` 
#' information. (Defaults to "Nominal.Conc".)
#' 
#param compound.conc (Numeric) The intended concentration
#of the test chemical for calibration curves. (Note: Single entry only,
#use only if all test compounds have the same intended concentration.)
#(Defaults to \code{NULL}.)
#'
#' @param dilution.col (Character) Column name containing `dilution` 
#' information. (Defaults to "Dilution.Factor".)
#'
#' @param dilution (Numeric) Number of times the sample was diluted before MS 
#' analysis. (Defaults to \code{NULL}.) (Note: Single entry only, use only if all 
#' samples underwent the same number of dilutions.)
#'
#' @param membrane.area.col (Character) Column name containing `membrane.area` 
#' information. (Defaults to "Membrane.Area".)
#'
#' @param membrane.area (Numeric) The area of the Caco-2 monolayer (in cm^2). 
#' (Defaults to \code{NULL}.) (Note: Single entry only, use only if all tested compounds 
#' have the same area for the Caco-2 monolayer.)
#'
#' @param receiver.vol.col (Character) Column name of data.in containing the volume
#' (in cm^3) of the receiver portion of the Caco-2 experimental well. 
#' (Defaults to "Vol.Receiver".)
#'
#' @param donor.vol.col (Character) Column name of data.in containing the volume
#' (in cm^3) of the donor portion of the Caco-2 experimental well where the
#' test chemical is added. (Defaults to "Vol.Donor".)
#'
#' @param meas.time.col (Character) Column name containing `meas.time` 
#' information. (Defaults to "Time".)
#'
#' @param meas.time (Numeric) The amount of time (in hours) before the receiver and donor 
#' compartments are measured. (Defaults to 2.)
#'
#' @param istd.col (Character) Column name of data.in containing the
#' MS peak area for the internal standard. (Defaults to "ISTD.Area".)
#'
#' @param istd.name.col (Character) Column name containing `istd.name` information. (Defaults to "ISTD.Name".)
#'
#' @param istd.name (Character) The identity of the internal standard. (Defaults to \code{NULL}.) 
#' (Note: Single entry only, use only if all tested compounds use the same internal standard.) 
#'
#' @param istd.conc.col (Character) Column name containing `istd.conc` information. (Defaults to "ISTD.Conc".)
#'
#' @param istd.conc (Numeric) The concentration for the internal standard. (Defaults to \code{NULL}.) 
#' (Note: Single entry only, use only if all tested compounds have the same 
#' concentration of the internal standard.) 
#'
#' @param nominal.test.conc.col (Character) Column name containing `nominal.test.conc` 
#' information. (Defaults to "Test.Target.Conc".)
#'
#' @param nominal.test.conc (Numeric) The test chemical concentration in the dosing solution
#' that added to the donor side at time zero. (Defaults to \code{NULL}.) (Note: Single entry only, use only if the same initial 
#' concentration was used for all tested compounds.)
#'
#' @param analysis.method.col (Character) Column name containing `analysis.method` 
#' information. (Defaults to "Analysis.Method".)
#'
#' @param analysis.method (Character) The analytical chemistry analysis method, 
#' typically "LCMS" or "GCMS", liquid chromatography or gas chromatography–mass spectrometry. 
#' (Defaults to \code{NULL}.) (Note: Single entry only, 
#' use only if the same method was used for all tested compounds.)
#'
#' @param analysis.instrument.col (Character) Column name containing `analysis.instrument` 
#' information. (Defaults to "Analysis.Instrument".)
#'
#' @param analysis.instrument (Character) The instrument used for chemical analysis, 
#' for example "Agilent 6890 GC with model 5973 MS". (Defaults to \code{NULL}.) 
#' (Note: Single entry only, use only if the same instrument was used for all tested compounds.) 
#'
#' @param analysis.parameters.col (Character) Column name containing `analysis.parameters` 
#' information. (Defaults to "Analysis.Paramaters".)
#'
#' @param analysis.parameters (Character) The parameters used to identify the 
#' compound on the chemical analysis instrument, for example
#' "Negative Mode, 221.6/161.6, -DPb=26, FPc=-200, EPd=-10, CEe=-20, CXPf=-25.0". (Defaults to \code{NULL}.) 
#' (Note: Single entry only, use only if the same parameters were used for all tested compounds.) 
#'
#' @return data.frame A data.frame in standardized "level1" format containing a  
#' standardized set of columns with standardized column names. 
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
#' @export format_caco2
format_caco2 <- function(data.in,
  FILENAME = "MYDATA",
  sample.col="Lab.Sample.Name",
  lab.compound.col="Lab.Compound.Name",
  dtxsid.col="DTXSID",
  date.col="Date",
  series.col="Series",
  series=NULL,
  compound.col="Compound.Name",
  area.col="Area",
  istd.col="ISTD.Area",
  type.col="Type",
  direction.col="Direction",
  membrane.area.col="Membrane.Area",
  membrane.area=NULL,
  receiver.vol.col="Vol.Receiver",
  donor.vol.col="Vol.Donor",
  #compound.conc=NULL,
  compound.conc.col="Nominal.Conc",
  cal=NULL,
  cal.col="Cal",
  dilution=NULL,
  dilution.col="Dilution.Factor",
  meas.time.col="Time",
  meas.time = 2,
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
  # These are the required data types as indicated by type.col.
  # In order to calculate the parameter a chemical must have peak areas for each
  # of these measurements:
  req.types=c("Blank","D0","D2","R2")

  data.out <- as.data.frame(data.in)
  # Force code to throw error if data.in accessed after this point:
  rm(data.in)

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
  if (!is.null(meas.time)) data.out[,meas.time.col] <- meas.time

# We need all these columns in data.out
  cols <-c(
    sample.col,
    date.col,
    compound.col,
    dtxsid.col,
    lab.compound.col,
    type.col,
    direction.col,
    dilution.col,
    cal.col,
    series.col,
    compound.conc.col,
    nominal.test.conc.col,
    meas.time.col,
    istd.name.col,
    istd.conc.col,
    istd.col,
    area.col,
    membrane.area.col,
    donor.vol.col,
    receiver.vol.col,
    analysis.method.col,
    analysis.instrument.col,
    analysis.parameters.col
    )

  if (!(all(cols %in% colnames(data.out))))
  {
    stop(paste("Missing columns named:",
      paste(cols[!(cols%in%colnames(data.out))],collapse=", ")))
  }

  # Only include the data types used:
  data.out <- subset(data.out,data.out[,type.col] %in% req.types)

  # Organize the columns:
  data.out <- data.out[,cols]

  # Standardize the column names:
    sample.col <- "Lab.Sample.Name"
    date.col <- "Date"
    compound.col <- "Compound.Name"
    dtxsid.col <- "DTXSID"
    lab.compound.col <- "Lab.Compound.Name"
    type.col <- "Sample.Type"
    direction.col <- "Direction"
    dilution.col <- "Dilution.Factor"
    cal.col <- "Calibration"
    series.col <- "Series"
    compound.conc.col <- "Standard.Conc"
    nominal.test.conc.col <- "Test.Target.Conc"
    meas.time.col <- "Time"
    istd.name.col <- "ISTD.Name"
    istd.conc.col <- "ISTD.Conc"
    istd.col <- "ISTD.Area"
    area.col <- "Area"
    membrane.area.col <- "Membrane.Area"
    donor.vol.col <- "Vol.Donor"
    recevier.vol.col <- "Vol.Receiver"
    analysis.method.col <- "Analysis.Method"
    analysis.instrument.col <- "Analysis.Instrument"
    analysis.parameters.col <- "Analysis.Parameters"

  colnames(data.out) <- c(
    sample.col,
    date.col,
    compound.col,
    dtxsid.col,
    lab.compound.col,
    type.col,
    direction.col,
    dilution.col,
    cal.col,
    series.col,
    compound.conc.col,
    nominal.test.conc.col,
    meas.time.col,
    istd.name.col,
    istd.conc.col,
    istd.col,
    area.col,
    membrane.area.col,
    donor.vol.col,
    receiver.vol.col,
    analysis.method.col,
    analysis.instrument.col,
    analysis.parameters.col
    )

  # calculate the response:
  data.out[,area.col] <- signif(as.numeric(data.out[,area.col]), 5)
  data.out[,istd.col] <- signif(as.numeric(data.out[,istd.col]), 5)
  data.out[,istd.conc.col] <- as.numeric(data.out[,istd.conc.col])
  data.out[,"Response"] <- signif(data.out[,area.col] /
     data.out[,istd.col] *  data.out[,istd.conc.col], 4)

# Write out a "level 1" file (data organized into a standard format):
  write.table(data.out,
    file=paste(FILENAME,"-Caco-2-Level1.tsv",sep=""),
    sep="\t",
    row.names=F,
    quote=F)

  summarize_table(data.out,
    req.types=req.types)

  return(data.out)
}


