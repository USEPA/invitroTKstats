#' Creates a standardized data table reporting UC PPB data
#'
#' This function formats data describing mass spectrometry (MS) peak areas
#' from samples collected as part of in vitro measurement of chemical fraction
#' unbound in plasma using ultracentrifugation (Redgrave 1975?).
#' An input dataframe is organized into a standard set of columns and is written
#' to a tab-separated text file. 
#'
#' The data frame of observations should be annotated according to
#' of these types:
#' \tabular{rrrrr}{
#'   Calibration Curve \tab CC\cr
#'   Ultra-centrifugation Aqueous Fraction \tab UC\cr
#'   Whole Plasma T1h Sample  \tab T1\cr
#'   Whole Plasma T5h Sample \tab T5\cr
#' }
#' Chemical concentration is calculated qualitatively as a response:
#'
#' Response <- AREA / ISTD.AREA * ISTD.CONC
#'
#' @param FILENAME A string used to identify outputs of the function call.
#' (defaults to "MYDATA")
#' 
#' @param PPB.data A data frame containing mass-spectrometry peak areas,
#' indication of chemical identiy, and measurment type. The data frame should
#' contain columns with names specified by the following arguments:
#' 
#' @param sample.col Which column of PPB.data indicates the unique mass 
#' spectrometry (MS) sample name used by the laboratory. (Defaults to 
#' "Lab.Sample.Name")
#' 
#' @param lab.compound.col Which column of PPB.data indicates The test compound 
#' name used by the laboratory (Defaults to "Lab.Compound.Name")
#' 
#' @param dtxsid.col Which column of PPB.data indicates EPA's DSSTox Structure 
#' ID (\url{http://comptox.epa.gov/dashboard}) (Defaults to "DTXSID")
#' 
#' @param date.col Which column of PPB.data indicates the laboratory measurment
#' date (Defaults to "Date")
#' 
#' @param compound.col Which column of PPB.data indicates the test compound
#' (Defaults to "Compound.Name")
#' 
#' @param area.col Which column of PPB.data indicates the target analyte (that 
#' is, the test compound) MS peak area (Defaults to "Area")
#' 
#' @param series.col Which column of PPB.data indicates the "series", that is
#' a simultaneous replicate (Defaults to "Series")
#' 
#' @param type.col Which column of PPB.data indicates the sample type (see table
#' above)(Defaults to "Sample.Type")
#' 
#' @param cal.col Which column of PPB.data indicates the MS calibration -- for
#' instance different machines on the same day or different days with the same
#' MS analyzer (Defaults to "Cal")
#' 
#' @param dilution.col Which column of PPB.data indicates how many times the
#' sample was diluted before MS analysis (Defaults to "Dilution.Factor")
#' 
#' @param istd.col Which column of PPB.data indicates the MS peak area for the
#' internal standard (Defaults to "ISTD.Area")
#' 
#' @param istd.name.col Which column of PPB.data indicates identity of the 
#' internal standard (Defaults to "ISTD.Name")
#' 
#' @param istd.conc.col Which column of PPB.data indicates the concentration of
#' the internal standard (Defaults to "ISTD.Conc")
#' 
#' @param nominal.test.conc.col Which column of PPB.data indicates the intended
#' test chemical concentration at time zero (Defaults to "Test.Target.Conc") 
#'
#' @return \item{data.frame}{A data.frame in standardized format} 
#'
#' @author John Wambaugh
#' 
#' @examples
#' library(invitroTKstats)
#' level0 <- kreutz2020
#' level1 <- format_fup_uc(level0,
#'   FILENAME="Kreutz2020",
#'   compound.col="Name",
#'   compound.conc.col="Standard.Conc",
#'   area.col="Chem.Area"
#'   )
#' 
#' @references
#' Redgrave, T. G., D. C. K. Roberts, and C. E. West. "Separation of plasma 
#' lipoproteins by density-gradient ultracentrifugation." Analytical 
#' Biochemistry 65.1-2 (1975): 42-49.#' 
#' 
#' @export format_fup_uc
format_fup_uc <- function(PPB.data,
  FILENAME = "MYDATA",
  sample.col="Lab.Sample.Name",
  lab.compound.col="Lab.Compound.Name",
  dtxsid.col="DTXSID",
  date.col="Date",
  compound.col="Compound.Name",
  area.col="Area",
  series.col="Series",
  type.col="Sample.Type",
  compound.conc.col="Nominal.Conc",
  cal.col="Cal",
  dilution.col="Dilution.Factor",
  istd.col="ISTD.Area",
  istd.name.col="ISTD.Name",
  istd.conc.col="ISTD.Conc",
  nominal.test.conc.col="Test.Target.Conc" 
  )
{  
  PPB.data <- as.data.frame(PPB.data)

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
    compound.conc.col,
    nominal.test.conc.col,
    istd.name.col,
    istd.conc.col,
    istd.col,
    series.col,
    area.col)
  
  if (!(all(cols %in% colnames(PPB.data))))
  {
    stop(paste("Missing columns named:",
      paste(cols[!(cols%in%colnames(PPB.data))],collapse=", ")))
  }

  # Only include the data types used:
  PPB.data <- subset(PPB.data,PPB.data[,type.col] %in% c("CC","T1","T5","AF"))
  
  # Organize the columns:
  PPB.data <- PPB.data[,cols]
    
  # Standardize the column names:
    sample.col <- "Lab.Sample.Name"
    date.col <- "Date"
    compound.col <- "Compound.Name"
    dtxsid.col <- "DTXSID"
    lab.compound.col <- "Lab.Compound.Name"
    type.col <- "Sample.Type"
    dilution.col <- "Dilution.Factor"
    cal.col <- "Calibration"
    compound.conc.col <- "Standard.Conc"
    nominal.test.conc.col <- "Test.Target.Conc"
    istd.name.col <- "ISTD.Name"
    istd.conc.col <- "ISTD.Conc"
    istd.col <- "ISTD.Area"
    series.col <- "Series"
    area.col <- "Area"
    
  colnames(PPB.data) <- c(
    sample.col,
    date.col,
    compound.col,
    dtxsid.col,
    lab.compound.col,
    type.col,
    dilution.col,
    cal.col,
    compound.conc.col,
    nominal.test.conc.col,
    istd.name.col,
    istd.conc.col,
    istd.col,
    series.col,
    area.col)
  
  # calculate the reponse:
  PPB.data[,"Response"] <- PPB.data[,area.col] /
     PPB.data[,istd.col] *  PPB.data[,istd.conc.col]
  
# Write out a "level 1" file (data organized into a standard format):  
  write.table(PPB.data, 
    file=paste(FILENAME,"-PPB-UC-Level1.tsv",sep=""),
    sep="\t",
    row.names=F,
    quote=F)
    
  summarize_table(PPB.data,
    req.types=c("CC","T1","T5","AF"))

  return(PPB.data)  
}


