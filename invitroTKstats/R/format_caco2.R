#' Creates a standardized data table reporting Caco-2 data
#'
#' This function formats data describing mass spectrometry (MS) peak areas
#' from samples collected as part of in vitro measurement of membrane
#' permeability using Caco-2 cells \insertCite{hubatsch2007determination}{invitroTKstats}. The input dataframe is
#' organized into a standard set of columns and is written to a tab-separated
#' text file.
#'
#' In this experiment an
#' in vitro well is separated into two by a membrane composed of a monolayer of
#' Caco-2 cells. A test chemical is added to eithe the apical or basal side of
#' of the monolayer at time 0, and after a set time samples are taken from both
#' the "donor" (side where the test chemical was added) and the "receiver side.
#' Depending on the direction of the test the donor side can be either apical or
#' basal.
#'
#' The data frame of observations should be annotated according to direction
#' (either apical to basal -- "AtoB" -- or basal to apical -- "BtoA") and type
#' of concentrtion measured:
#' \tabular{rr}{
#'   Blank with no chemical added \tab Blank \cr
#'   Dosing vehicle (C0) at target concentration \tab D0\cr
#'   Donor compartment at end of experiment \tab D2\cr
#'   Receiver compartment at end of experiment\tab R2\cr
#' }
#'
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
#' @param series.col Which column of PPB.data indicates the "series", that is
#' a simultaneous replicate with the same analytical chemistry
#' (Defaults to "Series")
#'
#' @param series If this argument is used (defaults to NULL) every observation
#' in the table is assigned the value of the argument and the corresponding
#' column in input.table (if present) is ignored.
#'
#' @param compound.col Which column of data.in indicates the test compound
#' (Defaults to "Compound.Name")
#'
#' @param area.col Which column of data.in indicates the target analyte (that
#' is, the test compound) MS peak area (Defaults to "Area")
#'
#' @param type.col Which column of data.in indicates the sample type (see table
#' above)(Defaults to "Type")
#'
#' @param type.col Which column of data.in indicates the direction of the
#' measurements (either "AtoB" for apical to basolateral or "BtoA" for vice
#' versa) (Defaults to "Direction")
#'
#' @param cal.col Which column of data.in indicates the MS calibration -- for
#' instance different machines on the same day or different days with the same
#' MS analyzer (Defaults to "Cal")
#'
#' @param cal If this argument is used (defaults to NULL) every observation in
#' the table is assigned the value of the argument and the corresponding
#' column in input.table (if present) is ignored.
#'
#' #param compound.conc.col Which column indictes the intended concentration
#' of the test chemical for calibration curves (Defaults to "Standard.Conc")
#'
#' @param dilution.col Which column of data.in indicates how many times the
#' sample was diluted before MS analysis (Defaults to "Dilution.Factor")
#'
#' @param dilution If this argument is used (defaults to NULL) every
#' observation in the table is assigned the value of the argument and the
#' corresponding column in input.table (if present) is ignored.
#'
#' @param membrane.area.col Which column of data.in indicates the area of the
#' Caco-2 monolayer (in cm^2) (Defaults to "Membrane.Area")
#'
#' @param membrane.area If this argument is used (defaults to NULL) every
#' observation in the table is assigned the value of the argument and the
#' corresponding column in input.table (if present) is ignored.
#'
#' @param receiver.vol.col Which column of data.in indicates the volume
#' (in cm^3) of the receiver portion of the Caco-2 experimental well
#' (Defaults to "Vol.Receiver")
#'
#' @param donor.vol.col Which column of data.in indicates the volume
#' (in cm^3) of the donor portion of the Caco-2 experimental well where the
#' test chemical is added
#' (Defaults to "Vol.Donor")
#'
#' @param direction.col Which column of data.in indicates the direction of
#' the Caco-2 permeability experiment, either apical to basal (AtoB) or basal
#' to aprical (BtoA). (Defaults to "Direction")
#'
#' @param meas.time.col Which column of data.in indicates the amount of time
#' before the receiver and donor compartments are measured (Defaults to "Time")
#'
#' @param meas.time If this argument is used (defaults to 2 h) every
#' observation in the table is assigned the value of the argument and the
#' corresponding column in input.table (if present) is ignored.
#'
#' @param istd.col Which column of data.in indicates the MS peak area for the
#' internal standard (Defaults to "ISTD.Area")
#'
#' @param istd.name.col Which column of data.in indicates identity of the
#' internal standard (Defaults to "ISTD.Name")
#'
#' @param istd.name If this argument is used (defaults to NULL) every
#' observation in the table is assigned the value of the argument and the
#' corresponding column in input.table (if present) is ignored.
#'
#' @param istd.conc.col Which column of data.in indicates the concentration of
#' the internal standard (Defaults to "ISTD.Conc")
#'
#' @param istd.conc If this argument is used (defaults to NULL) every
#' observation in the table is assigned the value of the argument and the
#' corresponding column in input.table (if present) is ignored.
#'
#' @param nominal.test.conc.col Which column of data.in indicates the intended
#' test chemical concentration at time zero in the dosing solution (added to the
#' donor side of the Caco-2 test well) (Defaults to "Test.Target.Conc")
#'
#' @param nominal.test.conc If this argument is used (defaults to NULL) every
#' observation in the table is assigned the value of the argument and the
#' corresponding column in input.table (if present) is ignored.
#'
#' @param analysis.method.col Which column of data.in indicates the analytical
#' chemistry analysis method, typically "LCMS" or "GCMS" (Defaults to
#' "Analysis.Method")
#'
#' @param analysis.method If this argument is used (defaults to NULL) every
#' observation in the table is assigned the value of the argument and the
#' corresponding column in input.table (if present) is ignored.
#'
#' @param analysis.instrument.col Which column of data.in indicates the
#' instrument used for chemical analysis, for example
#' "Agilent 6890 GC with model 5973 MS" (Defaults to
#' "Analysis.Instrument")
#'
#' @param analysis.instrument If this argument is used (defaults to NULL) every
#' observation in the table is assigned the value of the argument and the
#' corresponding column in input.table (if present) is ignored.
#'
#' @param analysis.parameters.col Which column of data.in indicates the
#' parameters used to identify the compound on the chemical analysis instrument,
#' for example
#' "Negative Mode, 221.6/161.6, -DPb=26, FPc=-200, EPd=-10, CEe=-20, CXPf=-25.0"
#' (Defaulys to "Analysis.Paramaters").
#'
#' @param analysis.parameters If this argument is used (defaults to NULL) every
#' observation in the table is assigned the value of the argument and the
#' corresponding column in input.table (if present) is ignored.
#'
#' @return data.frame A data.frame in standardized "level1" format
#'
#' @author John Wambaugh
#'
#' @examples
#' library(invitroTKstats)
#' level0 <- TO1caco2
#' level1 <- format_caco2(level0,
#'   FILENAME="EPACyprotex2021",
#'   compound.col="CompoundName",
#'   compound.conc.col="Standard.Conc",
#'   membrane.area.col=0.11
#'   )
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
  compound.conc=NULL,
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

# These arguments allow the user to specify a single value for every obseration
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

  # calculate the reponse:
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


