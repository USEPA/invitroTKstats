#' Plot Mass Spectrometry Responses from Measurements of Intrinsic
#' Hepatic Clearance
#'
#' This function generates a response-versus-time plot of mass spectrometry (MS) 
#' responses collected from measurements of intrinsic hepatic clearance for a chemical.
#' Responses from different measurements/calibrations are labeled with different colors, 
#' and responses from various sample types are labeled with different shapes.  
#'
#' The function requires "Level-2" data for plotting. Level-2 data is Level-1,
#' data formatted with the \code{\link{format_clint}} function, and curated
#' with a verification column. "Y" in the verification column indicates the
#' data row is valid for plotting.  
#' 
#' @param level2 (Data Frame) A data frame containing Level-2 data with a measure
#' of chemical clearance over time when incubated with suspended hepatocytes.
#' 
#' @param dtxsid (Character) EPA's DSSTox Structure ID for the chemical to be plotted.
#'
#' @return \item{ggplot2}{A figure of mass spectrometry responses over time for
#' various sample types.}
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
