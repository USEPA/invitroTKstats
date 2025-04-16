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
#' @param color.palette (Character) A character string indicating which 
#' \code{viridis} R package color map option to use. (Defaults to "viridis".) 
#'
#' @return \item{ggplot2}{A figure of mass spectrometry responses over time for
#' various sample types.}
#'
#' @author John Wambaugh
#' 
#' @examples
#' ## Load example level-2 data 
#' level2 <- invitroTKstats::clint_L2
#' plot_clint(level2, dtxsid = "DTXSID1021116")
#' 
#' @export plot_clint
#' @import ggplot2
plot_clint <- function(level2,dtxsid,color.palette = "viridis")
{
  #assigning global variables
  DTXSID <- Time <- Response <- Sample.Type <- Calibration <- NULL
  
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
  std.conc.col <- "Test.Compound.Conc"
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
    std.conc.col,
    time.col,
    area.col)

  if (!(all(cols %in% colnames(level2))))
  {
    warning("Is this Clint data? Run format_clint first (level 1) then curate to (level 2).")
    stop(paste("Missing columns named:",
      paste(cols[!(cols%in%colnames(level2))],collapse=", ")))
  }

  level2 <- subset(level2, DTXSID==dtxsid)
  # Remove data with missing Time (i.e. Blank and CC samples)
  level2 <- subset(level2, !is.na(Time))

  out <- ggplot(level2, aes(x=Time, y=Response)) +
    labs(title = level2[1,"DTXSID"], 
         caption = level2[1,"Compound.Name"]) +
    geom_point(mapping = aes(
      #fill = factor(Calibration),
      shape = factor(Sample.Type),
      color=factor(Calibration)), size = 3, alpha = 0.6) +
      #scale_color_brewer(palette = color.palette) +
    guides(color=guide_legend(title="Calibrations"),
           shape=guide_legend(title="Sample Types")) + 
    scale_colour_viridis_d(option = color.palette, end = 0.75) + 
    theme(plot.caption = element_text(hjust = 0.5))

  return(out)
}
