#' Plot Mass Spec. Response for Measurement of Fraction Unbound in Plasma (UC)
#'
#' This function use describing mass spectrometry (MS) peak areas
#' from samples collected as part of in vitro measurement of chemical fraction
#' unbound in plasma using ultracentrifugation \insertCite{redgrave1975separation}{invitroTKstats}.
#' Data are read from a "Level2" text file that should have been formatted and created
#' by \code{\link{format_fup_red}} (this is the "Level1" file). The Level1 file
#' should have been curated and had a column added with the value "Y" indicating
#' that each row is verified as usable for analysis (that is, the Level2 file).
#'
#' The should be annotated according to
#' of these types:
#' \tabular{rrrrr}{
#'   Blank (ignored) \tab Blank\cr
#'   Plasma well concentration \tab Plasma\cr
#'   Phosphate-buffered well concentration\tab PBS\cr
#'   Time zero plasma concentration \tab T0\cr
#' }
#' @param level2 A data.frame containing level2 data for fraction unbound in
#' plasma measured by ultracentrifugation.
#' 
#' @param dtxsid Which chemical to be plotted.
#'
#' @return \item{ggplot2}{A figure of mass spec. response for different sample types}
#'
#' @author John Wambaugh
#'
#' @export plot_fup_uc
#' @import ggplot2
plot_fup_uc <- function(level2,dtxsid, good.col="Verified")
{
# We need all these columns in uc data
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

# For a properly formatted level 2 file we should have all these columns:
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
    "Response",
    good.col)
  if (!(all(cols %in% colnames(level2))))
  {
    warning("Is this UC fup data? Run format_fup_uc first (level 1) then curate to level 2.")
    stop(paste("Missing columns named:",
      paste(cols[!(cols%in%colnames(level2))],collapse=", ")))
  }

  level2 <- subset(level2, DTXSID==dtxsid & Verified=="Y")
  frac <- subset(level2, Sample.Type=="AF")
  for (this.cal in unique(frac$Calibration))
  {
    this.t5 <- subset(level2, Calibration==this.cal & Sample.Type=="T5")
    mean.t5 <- mean(this.t5$Response*this.t5$Dilution.Factor,na.rm=TRUE)
    frac[frac$Calibration == this.cal,"Response"] <-
      frac[frac$Calibration == this.cal,"Response"] *
      frac[frac$Calibration == this.cal,"Dilution.Factor"] /
      mean.t5
  }
  frac$Sample.Type = "Rough Fup"
  level2 <- rbind(level2,frac)

  out <- ggplot(level2, aes(x=factor(Sample.Type), y=Response*Dilution.Factor)) +
    geom_point(mapping = aes(
      fill = factor(Calibration),
      shape = factor(Calibration),
      color=factor(Calibration)), size = 5)+
    scale_y_log10() +
    ylab("Mass Spec. Intensity / Fraction Unbound") +
    xlab("Sample Type") +
    ggtitle(paste(level2[1,"Compound.Name"]," (",level2[1,"DTXSID"],")",sep="")) +
    guides(
      fill=guide_legend(title="Calibration"),
      shape=guide_legend(title="Calibration"),
      color=guide_legend(title="Calibration"))
  print(out)

  return(out)
}
