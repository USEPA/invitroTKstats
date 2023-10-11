#' Calculate a point estimate of apparent membrane permeability from Caco-2 data
#'
#' This function uses mass spectrometry (MS) peak areas
#' from samples collected as part of in vitro measurement of membrane
#' permeability using Caco-2 cells \insertCite{hubatsch2007determination}{invitroTKstats}.
#' Data are read from a "Level2" text file that should have been formatted and
#' created
#' by \code{\link{format_caco2}} (this is the "Level1" file). The Level1 file
#' should have been curated by adding a column with the value "Y" indicating
#' that each row is verified as usable for analysis (that is, the Level2 file).
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
#' Apparent membrane permeability (\eqn{P_{app}}) is calculated from MS responses as:
#'
#'
#' \eqn{P_{app} = \frac{dQ/dt}{c_0/A}}
#'
#' The rate of permeation, \eqn{\frac{dQ}{dt}}\eqn{\left(\frac{\text{peak area}}{\text{time (s)}} \right)} is calculated as:
#'
#' \eqn{\frac{dQ}{dt} = \max\left(0, \frac{\sum_{i=1}^n(r_P * c_{DF})}{n_P} - \frac{\sum_{i=1}^n(r_B * c_{DF})}{n_B}\right)}
#'
#' where \eqn{r_P} is PBS Response, \eqn{c_{DF}} is Dilution Factor, \eqn{r_B} is Blank Response,
#' \eqn{n_P} is the number of PBS Responses, and \eqn{n_B} is the number of Blank Responses.
#'
#' @param FILENAME A string used to identify the input file, whatever the
#' argument given, "-Caco-2-Level2.tsv" is appended (defaults to "MYDATA")
#'
#' @param good.col Name of a column indicating which rows have been verified for
#' analysis, indicated by a "Y" (Defaults to "Verified")
#'
#' @return \item{data.frame}{A data.frame in standardized format}
#' \tabular{rrr}{
#'   C0_A2B \tab Time zero donor concentration \tab Mass Spec Response Ratio (RR) \cr
#'   dQdt_A2B \tab Estimated rate of mass movement through membrane \tab RR*cm^3/s \cr
#'   Papp_A2B \tab Apparent membrane permeability \tab 10^-6 cm/s\cr
#'   C0_B2A \tab Time zero donor concentration \tab Mass Spec Response Ratio (RR) \cr
#'   dQdt_B2A \tab Estimated rate of mass movement through membrane \tab RR*cm^3/s \cr
#'   Papp_B2A \tab Apparent membrane permeability \tab 10^-6 cm/s\cr
#'   Refflux \tab Efflux ratio \tab unitless\cr
#' }
#'
#' @author John Wambaugh
#'
#' @examples
#' library(invitroTKstats)
#'
#' level0 <- TO1caco2
#' level1 <- format_caco2(level0,
#'    FILENAME="EPACyprotex2021",
#'    sample.col="SampleName",
#'    dtxsid.col="CompoundName",
#'    lab.compound.col="CompoundName",
#'    cal=1,
#'    istd.conc.col="ISTD.Conc",
#'    compound.col="CompoundName",
#'    compound.conc.col="Test.Target.Conc",
#'    membrane.area=0.11,
#'    series=1,
#'    analysis.parameters="Feature",
#'    analysis.instrument="GC or LC",
#'    analysis.method="Mass Spec"
#'    )
#'
#' level2 <- level1
#' level2$Verified <- "Y"
#'
#' write.table(level2,
#'   file="EPACyprotex2021-Caco-2-Level2.tsv",
#'   sep="\t",
#'   row.names=F,
#'   quote=F)
#'
#' level3 <- calc_caco2_point(FILENAME="EPACyprotex2021")
#'
#' write.table(level3,
#'   file="EPACyprotex2021-Caco-2-Level3.tsv",
#'   sep="\t",
#'   row.names=F,
#'   quote=F)
#'
#' @references
#' \insertRef{hubatsch2007determination}{invitroTKstats}
#'
#' @export calc_caco2_point
calc_caco2_point <- function(FILENAME, good.col="Verified")
{
  # These are the required data types as indicated by type.col.
  # In order to calculate the parameter a chemical must have peak areas for each
  # of these measurements:
  req.types=c("Blank","D0","D2","R2")

  input.table <- read.csv(file=paste(FILENAME,"-Caco-2-Level2.tsv",sep=""),
    sep="\t",header=T)
  input.table <- subset(input.table,!is.na(Compound.Name))
  input.table <- subset(input.table,!is.na(Response))

  # Standardize the column names:
    sample.col <- "Lab.Sample.Name"
    date.col <- "Date"
    compound.col <- "Compound.Name"
    dtxsid.col <- "DTXSID"
    lab.compound.col <- "Lab.Compound.Name"
    type.col <- "Sample.Type"
    dilution.col <- "Dilution.Factor"
    cal.col <- "Calibration"
    series.col <- "Series"
    compound.conc.col <- "Standard.Conc"
    nominal.test.conc.col <- "Test.Target.Conc"
    meas.time.col="Time"
    istd.name.col <- "ISTD.Name"
    istd.conc.col <- "ISTD.Conc"
    istd.col <- "ISTD.Area"
    series.col <- "Series"
    area.col <- "Area"
    membrane.area.col <- "Membrane.Area"
    donor.vol.col <- "Vol.Donor"
    receiver.vol.col <- "Vol.Receiver"
    analysis.method.col <- "Analysis.Method"
    analysis.instrument.col <- "Analysis.Instrument"
    analysis.parameters.col <- "Analysis.Parameters"

# For a properly formatted level 2 file we should have all these columns:
  cols <-c(
    sample.col,
    date.col,
    compound.col,
    dtxsid.col,
    lab.compound.col,
    type.col,
    dilution.col,
    cal.col,
    series.col,
    compound.conc.col,
    nominal.test.conc.col,
    meas.time.col,
    istd.name.col,
    istd.conc.col,
    istd.col,
    series.col,
    area.col,
    membrane.area.col,
    donor.vol.col,
    receiver.vol.col,
    analysis.method.col,
    analysis.instrument.col,
    analysis.parameters.col
    )

  if (!(all(cols %in% colnames(input.table))))
  {
    warning("Run format_fup_red first (level 1) then curate to (level 2).")
    stop(paste("Missing columns named:",
      paste(cols[!(cols%in%colnames(input.table))],collapse=", ")))
  }

  # Only include the data types used:
  input.table <- subset(input.table,input.table[,type.col] %in% req.types)

  # Only used verfied data:
  input.table <- subset(input.table, input.table[,good.col] == "Y")

  out.table <-NULL
  num.a2b <- 0
  num.b2a <- 0
  num.efflux <- 0
  for (this.chem in unique(input.table[,compound.col]))
  {
    this.subset <- subset(input.table, input.table[,compound.col]==this.chem)
    this.dtxsid <- this.subset$dtxsid[1]
    this.row <- cbind(this.subset[1,
      c(compound.col, dtxsid.col, meas.time.col, membrane.area.col)],
      data.frame(Calibration="All Data",
        C0_A2B = NaN, dQdt_A2B=NaN, Papp_A2B=NaN,
        C0_B2A = NaN, dQdt_B2A=NaN, Papp_B2A=NaN, Refflux=NaN))
    for (this.direction in c("AtoB","BtoA"))
    {
      this.blank <- subset(this.subset, Sample.Type=="Blank" &
        Direction==this.direction)
      this.dosing <- subset(this.subset, Sample.Type=="D0" &
        Direction==this.direction)
      this.donor <- subset(this.subset,Sample.Type=="D2" &
        Direction==this.direction)
      this.receiver <- subset(this.subset,Sample.Type=="R2" &
        Direction==this.direction)

    # Check to make sure there are data for PBS and plasma:
      if (dim(this.blank)[1]> 0 &
        dim(this.dosing)[1] > 0 &
        dim(this.donor)[1] > 0 &
        dim(this.receiver)[1] > 0)
      {
        if (this.direction == "AtoB")
        {
          dir.string <- "A2B"
          num.a2b <- num.a2b + 1
        } else {
          dir.string <- "B2A"
          num.b2a <- num.b2a+1
        }

        this.row[paste("C0",dir.string,sep="_")] <- max(0,
          mean(this.dosing$Response*this.dosing$Dilution.Factor) -
          mean(this.blank$Response*this.blank$Dilution.Factor))

        if (length(unique(this.dosing$Vol.Receiver))>1 |
          length(unique(this.dosing$Time))>1) browser()
        this.row[paste("dQdt",dir.string,sep="_")] <- max(0,
          mean(this.receiver$Response*this.receiver$Dilution.Factor) -
          mean(this.blank$Response*this.blank$Dilution.Factor)*
          this.dosing$Vol.Receiver[1] / this.dosing$Time[1] / 3600)

        this.row[paste("Papp",dir.string,sep="_")] <- max(0,
          as.numeric(this.row[paste("dQdt",dir.string,sep="_")]) /
          as.numeric(this.row[paste("C0",dir.string,sep="_")]) /
          as.numeric(this.row["Membrane.Area"]) * 1e6)
        }
      }

      if (!is.nan(unlist(this.row["Papp_A2B"])) &
        !is.nan(unlist(this.row["Papp_B2A"])))
      {
        num.efflux <- num.efflux + 1
        this.row["Refflux"] <- as.numeric(this.row["Papp_B2A"]) /
          as.numeric(this.row["Papp_A2B"])
      }
      out.table <- rbind(out.table, this.row)
      print(paste(this.row$Compound.Name,"Refflux =",
        signif(this.row$Refflux,3)))
  }

  rownames(out.table) <- make.names(out.table$Compound.Name, unique=TRUE)
  out.table[,"C0_A2B"] <- signif(as.numeric(out.table[,"C0_A2B"]),3)
  out.table[,"C0_B2A"] <- signif(as.numeric(out.table[,"C0_B2A"]),3)
  out.table[,"dQdt_A2B"] <- signif(as.numeric(out.table[,"dQdt_A2B"]),3)
  out.table[,"dQdt_B2A"] <- signif(as.numeric(out.table[,"dQdt_B2A"]),3)
  out.table[,"Papp_A2B"] <- signif(as.numeric(out.table[,"Papp_A2B"]),3)
  out.table[,"Papp_B2A"] <- signif(as.numeric(out.table[,"Papp_B2A"]),3)
  out.table[,"Refflux"] <- signif(as.numeric(out.table[,"Refflux"]),3)
  out.table <- as.data.frame(out.table)

# Write out a "level 3" file (data organized into a standard format):
  write.table(out.table,
    file=paste(FILENAME,"-Caco-2-Level3.tsv",sep=""),
    sep="\t",
    row.names=F,
    quote=F)

  print(paste("Apical to basal permeability calculated for",num.a2b,"chemicals."))
  print(paste("Basal to apical permeability calculated for",num.b2a,"chemicals."))
  print(paste("Efflux ratio calculated for",num.efflux,"chemicals."))

  return(out.table)
}


