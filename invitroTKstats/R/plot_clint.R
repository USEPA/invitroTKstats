#' Plot Mass Spec. Response for Measurment of Fraction Unbound in Plasma (UC)
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
#' plasma measured by ultracentrifucation.
#'
#' @param dtxsid Which chemical to be plotted.
#'
#' @return \item{ggplot2}{A figure of mass spec. response for different sample types}
#'
#' @author John Wambaugh
#'
#' @export plot_clint
#' @import ggplot2
plot_clint <- function(level2,dtxsid)
{
# We need all these columns in clint.data
# Standardize the column names:
  sample.col <- "Lab.Sample.Name"
  date.col <- "Date"
  compound.col <- "Compound.Name"
  dtxsid.col <- "DTXSID"
  lab.compound.col <- "Lab.Compound.Name"
  type.col <- "Sample.Type"
  dilution.col <- "Dilution.Factor"
  cal.col <- "Calibration"
  istd.name.col <- "ISTD.Name"
  istd.conc.col <- "ISTD.Conc"
  istd.col <- "ISTD.Area"
  density.col <- "Hep.Density"
  conc.col <- "Conc"
  time.col <- "Time"
  area.col <- "Area"
  analysis.method.col <- "Analysis.Method"
  analysis.instrument.col <- "Analysis.Instrument"
  analysis.parameters.col <- "Analysis.Parameters"


  cols <-c(
    sample.col,
    date.col,
    compound.col,
    dtxsid.col,
    lab.compound.col,
    type.col,
    dilution.col,
    cal.col,
    istd.name.col,
    istd.conc.col,
    istd.col,
    density.col,
    conc.col,
    time.col,
    area.col)

  if (!(all(cols %in% colnames(level2))))
  {
    warning("Is this Clint data? Run format_clint first (level 1) then curate to (level 2).")
    stop(paste("Missing columns named:",
      paste(cols[!(cols%in%colnames(level2))],collapse=", ")))
  }

  level2 <- subset(level2, DTXSID==dtxsid)

  out <- ggplot(level2, aes(x=Time, y=Response)) +
    geom_point(mapping = aes(
      fill = factor(Sample.Type),
      shape = factor(Sample.Type),
      color=factor(Calibration)), size = 5)
  print(out)
  return(out)
}
