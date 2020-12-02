#' Creates a standardized data table reporting hepatocye clearance data
#'
#' This function formats data describing mass spectrometry (MS) peak areas
#' from samples collected as part of in vitro measurement of chemical stability
#' when incubated with suspended hepatocytes (Shibata, 2000). Disappareance of 
#' the chemical over time is assumed to be due to metabolism by the hepatocytes.
#' An input dataframe is organized into a standard set of columns and is written
#' to a tab-separated text file. 
#'
#' The data frame of observations should be annotated according to
#' of these types:
#' \tabular{rrrrr}{
#'   Blank \tab Blank\cr
#'   Hepatocyte inciubation concentration \tab Cvst\cr
#' }
#' Chemical concentration is calculated qualitatively as a response:
#'
#' Response <- AREA / ISTD.AREA * ISTD.CONC
#'
#' @param FILENAME A string used to identify outputs of the function call.
#' (defaults to "MYDATA")
#' 
#' @param clint.data A data frame containing mass-spectrometry peak areas,
#' indication of chemical identiy, and measurment type. The data frame should
#' contain columns with names specified by the following arguments:
#' 
#' @param sample.col Which column of clint.data indicates the unique mass 
#' spectrometry (MS) sample name used by the laboratory. (Defaults to 
#' "Lab.Sample.Name")
#' 
#' @param lab.compound.col Which column of clint.data indicates The test compound 
#' name used by the laboratory (Defaults to "Lab.Compound.Name")
#' 
#' @param dtxsid.col Which column of clint.data indicates EPA's DSSTox Structure 
#' ID (\url{http://comptox.epa.gov/dashboard}) (Defaults to "DTXSID")
#' 
#' @param date.col Which column of clint.data indicates the laboratory measurment
#' date (Defaults to "Date")
#' 
#' @param compound.col Which column of clint.data indicates the test compound
#' (Defaults to "Compound.Name")
#' 
#' @param area.col Which column of clint.data indicates the target analyte (that 
#' is, the test compound) MS peak area (Defaults to "Area")
#' 
#' @param series.col Which column of clint.data indicates the "series", that is
#' a simultaneous replicate (Defaults to "Series")
#' 
#' @param type.col Which column of clint.data indicates the sample type (see table
#' above)(Defaults to "Sample.Type")
#' 
#' @param cal.col Which column of clint.data indicates the MS calibration -- for
#' instance different machines on the same day or different days with the same
#' MS analyzer (Defaults to "Cal")
#' 
#' @param dilution.col Which column of clint.data indicates how many times the
#' sample was diluted before MS analysis (Defaults to "Dilution.Factor")
#'
#' @param density.col Which column of clint.data indicates the density (units of
#' millions of hepatocytes per mL) hepatocytes in the in vitro incubation 
#' (Defaults to "Hep.Density" )
#' 
#' @param istd.col Which column of clint.data indicates the MS peak area for the
#' internal standard (Defaults to "ISTD.Area")
#' 
#' @param istd.name.col Which column of clint.data indicates identity of the 
#' internal standard (Defaults to "ISTD.Name")
#' 
#' @param istd.conc.col Which column of clint.data indicates the concentration 
#' (units if uM) of
#' the internal standard (Defaults to "ISTD.Conc")
#' 
#' @param conc.col Which column of clint.data indicates the intended
#' test chemical concentration 
#' (units if uM) of
#' at time zero (Defaults to "Conc") 
#'
#' @param time.col Which column of clint.data indicates the intended
#' time of the measurment (in minutes) since the test chemical was introduced
#' into the hepatocyte incubation (Defaults to "Time") 
#'
#' @return \item{data.frame}{A data.frame in standardized format} 
#'
#' @author John Wambaugh
#' 
#' @examples
#' library(invitroTKstats)
#'red <- wambaugh2019.red
#'red$Date <- "2019"
#'red$Sample.Type <- "Blank"
#'red <- subset(red,!is.na(SampleName))
#'red[regexpr("PBS",red$SampleName)!=-1,"Sample.Type"] <- "PBS"
#'red[regexpr("Plasma",red$SampleName)!=-1,"Sample.Type"] <- "Plasma"
#'red$Dilution.Factor <- NA
#'red$Dilution.Factor <- as.numeric(red$Dilution.Factor)
#'red[red$Sample.Type=="PBS","Dilution.Factor"] <- 2
#'red[red$Sample.Type=="Plasma","Dilution.Factor"] <- 5
#'red[regexpr("T0",red$SampleName)!=-1,"Sample.Type"] <- "T0"
#'
#'red$Test.Target.Conc <- 5
#'red$ISTD.Name <- "Bucetin and Diclofenac"
#'red$ISTD.Conc <- 1
#'red$Series <- 1
#'
#'level1 <- format_fup_red(red,
#'  FILENAME="Wambaugh2019",
#'  sample.col="SampleName",
#'  compound.col="Preferred.Name",
#'  lab.compound.col="CompoundName",
#'  cal.col="RawDataSet")
#'
#' @references
#' Shibata, Yoshihiro, Hiroyuki Takahashi, and Yasuyuki Ishii. "A convenient in 
#' vitro screening method for predicting in vivo drug metabolic clearance using 
#' isolated hepatocytes suspended in serum." Drug metabolism and disposition 
#' 28.12 (2000): 1518-1523.
#'
#' @export format_clint
format_clint <- function(clint.data,
  FILENAME = "MYDATA",
  sample.col="Lab.Sample.Name",
  date.col="Date",
  compound.col="Compound.Name",
  dtxsid.col="DTXSID",
  lab.compound.col="Lab.Compound.Name",
  type.col="Sample.Type",
  dilution.col="Dilution.Factor",
  cal.col="Cal",
  istd.name.col="ISTD.Name",
  istd.conc.col="ISTD.Conc",
  istd.col="ISTD.Area",
  density.col="Hep.Density",
  conc.col="Conc", 
  time.col="Time", 
  area.col="Area"
  )
{
  clint.data <- as.data.frame(clint.data)

# We need all these columns in clint.data
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
  
  if (!(all(cols %in% colnames(clint.data))))
  {
    stop(paste("Missing columns named:",
      paste(cols[!(cols%in%colnames(clint.data))],collapse=", ")))
  }

  # Only include the data types used:
  clint.data <- subset(clint.data,clint.data[,type.col] %in% c(
    "Blank","Cvst"))
  
  # Organize the columns:
  clint.data <- clint.data[,cols]
    
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

  colnames(clint.data) <- c(
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
  
  # calculate the reponse:
  clint.data[,"Response"] <- clint.data[,area.col] /
     clint.data[,istd.col] * clint.data[,istd.conc.col]
  
# Write out a "level 1" file (data organized into a standard format):  
  write.table(clint.data, 
    file=paste(FILENAME,"-Clint-Level1.tsv",sep=""),
    sep="\t",
    row.names=F,
    quote=F)

  summarize_table(clint.data,
    req.types=c("Blank","Cvst"))

  return(clint.data)  
}


