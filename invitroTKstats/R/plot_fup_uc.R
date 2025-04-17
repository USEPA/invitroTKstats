#' Plot Mass Spectrometry Responses for Fraction Unbound in Plasma Data from
#' Ultracentrifugation (UC)
#'
#' This function generates a scatter plot of mass spectrometry (MS) responses 
#' for one chemical collected from measurement of fraction unbound in plasma (Fup)
#' using ultracentrifugation (UC). The scatter plot displays the MS 
#' responses (y-axis) by sample types (x-axis). Responses from different
#' measurements/calibrations are labeled with different shapes and colors. 
#' 
#' This function requires "Level-2" data for plotting. Level-2 data is Level-1,
#' data formatted with the \code{\link{format_fup_uc}} function, and curated
#' with a verification column. "Y" in the verification column indicates the
#' data row is valid for plotting.  
#'
#' @param level2 (Data Frame) A data.frame containing Level-2 data for fraction
#' unbound in plasma (Fup) measured by ultracentrifugation (UC).
#' 
#' @param dtxsid (Character) EPA's DSSTox Structure ID for the chemical to be plotted.
#' 
#' @param compare (Character) A string indicating the plot is for 
#' comparing the responses across sample types ("type") or across calibrations ("cal").
#' (Defaults to "type".) 
#' 
#' @param good.col (Character) Column name containg verification information,
#' data rows valid for plotting are indicated with a "Y". (Defaults to "Verified".)
#' 
#' @param color.palette (Character) A character string indicating which 
#' \code{viridis} R package color map option to use. (Defaults to "viridis".) 
#'
#' @return \item{ggplot2}{A figure of mass spectrometry responses for
#' various sample types.}
#'
#' @author John Wambaugh
#' 
#' @examples
#' ## Load example level-2 data 
#' level2 <- invitroTKstats::fup_uc_L2
#' plot_fup_uc(level2, dtxsid = "DTXSID0059829")
#' 
#' @import ggplot2
#' 
#' @export plot_fup_uc
plot_fup_uc <- function(level2,dtxsid, compare = "type",good.col="Verified", color.palette = "viridis")
{
  #assigning global variables
  DTXSID <- Verified <- Response <- Sample.Type <- Calibration <- Dilution.Factor <- NULL
  
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
    std.conc.col <- "Test.Compound.Conc"
    uc.assay.conc.col <- "Test.Nominal.Conc"
    istd.name.col <- "ISTD.Name"
    istd.conc.col <- "ISTD.Conc"
    istd.col <- "ISTD.Area"
    series.col <- "Biological.Replicates"
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
  
  if (compare == "type"){
  out <- ggplot(level2, aes(x=factor(Sample.Type), y=Response*Dilution.Factor)) +
    geom_boxplot(mapping = aes(
      color=factor(Calibration)))+
    geom_point(aes(x = factor(Sample.Type), colour = factor(Calibration)), alpha = 0.6, position = position_jitterdodge()) +
    guides(
      color=guide_legend(title="Calibrations"),
      x =  guide_axis(angle = 45)) 
  } else if (compare == "cal") {
    out <- ggplot(level2, aes(x=factor(Calibration), y=Response*Dilution.Factor)) +
      geom_boxplot(mapping = aes(
        color=factor(Sample.Type))) +
      geom_point(aes(x = factor(Calibration), colour = factor(Sample.Type)), alpha = 0.6, position = position_jitterdodge()) +
      guides(
        color=guide_legend(title="Sample Types"),
        x =  guide_axis(angle = 45))
  }
  
  out <- out + 
    labs(title = level2[1,"DTXSID"], 
                    caption = level2[1,"Compound.Name"]) +
    scale_y_continuous(trans=scales::pseudo_log_trans(base = 10)) +
    ylab("Mass Spec. Intensity / Fraction Unbound") +
    xlab("Sample Type") +
    scale_colour_viridis_d(option = color.palette, end = 0.75) +
    theme(plot.caption = element_text(hjust = 0.5))
    
  return(out)
}
