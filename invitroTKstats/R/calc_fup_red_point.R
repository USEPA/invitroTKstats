#' Calculate Point Estimates of Fraction Unbound in Plasma (Fup) with Rapid Equilibrium Dialysis (RED) Data 
#'
#' This function calculates point estimates of fraction unbound in plasma (Fup) from 
#' mass spectrometry (MS) responses collected as part of in vitro measurement 
#' of chemical Fup using rapid equilibrium dialysis \insertCite{waters2008validation}{invitroTKstats}. 
#' Calculation formula is available in Details. 
#' 
#' The input to this function should be "Level-2" data. Level-2 data is Level-1,
#' data formatted with the \code{\link{format_fup_red}} function, and curated
#' with a verification column. "Y" in the verification column indicates the
#' data row is valid for analysis. 
#'
#' The data frame of observations should be annotated according to
#' of these types:
#' \tabular{rrrrr}{
#'   No Plasma Blank (no chemical, no plasma) \tab NoPlasma.Blank\cr
#'   Plasma Blank (no chemical, just plasma) \tab Plasma.Blank\cr
#'   Time zero chemical and plasma \tab T0\cr
#'   Equilibrium chemical in phosphate-buffered well (no plasma) \tab PBS\cr
#'   Equilibrium chemical in plasma well \tab Plasma\cr
#' }
#'
#' \eqn{f_{up}} is calculated from MS responses as:
#'
#'
#' \eqn{f_{up} = \frac{\max\left( 0, \frac{\sum_{i=1}^{n_P} (r_P * c_{DF})}{n_P} - \frac{\sum_{i=1}^{n_{NPB}} (r_{NPB}*c_{DF})}{n_{NPB}}\right)}
#' {\frac{\sum_{i=1}^{n_P} (r_P * c_{DF})}{n_P} - \frac{\sum_{i=1}^{n_B} (r_B * c_{DF})}{n_B}}}
#'
#' where \eqn{r_P} is PBS Response, \eqn{n_P} is the number of PBS Responses,
#' \eqn{c_{DF}} is Dilution Factor, \eqn{r_{NPB}} is No Plasma Blank Response,
#' \eqn{n_{NPB}} is the number of No Plasma Blank Responses, \eqn{r_{B}} is Plasma Blank Response,
#' and \eqn{n_B} is the number of Plasma Blank Responses.
#'
#' @param FILENAME (Character) A string used to identify the input Level-2 file.
#' "<FILENAME>-fup-RED-Level2.tsv".
#'
#' @param good.col (Character) Column name indicating which rows have been
#' verified, data rows valid for analysis are indicated with a "Y".
#' (Defaults to "Verified".)
#'
#' @return A data frame with one row per chemical, contains chemical identifiers 
#' such as preferred compound name, EPA's DSSTox Structure ID, calibration and a point estimate of 
#' fraction unbound in plasma (Fup) for all chemicals in the input data frame. 
#'
#' @author John Wambaugh
#'
#' @examples
#' red <- subset(wambaugh2019.red, Protein==100)
#' red$Date <- "2019"
#' red$Sample.Type <- "Blank"
#' red <- subset(red,!is.na(SampleName))
#' red[regexpr("PBS",red$SampleName)!=-1,"Sample.Type"] <- "PBS"
#' red[regexpr("Plasma",red$SampleName)!=-1,"Sample.Type"] <- "Plasma"
#' red$Dilution.Factor <- NA
#' red$Dilution.Factor <- as.numeric(red$Dilution.Factor)
#' red[red$Sample.Type=="PBS","Dilution.Factor"] <- 2
#' red[red$Sample.Type=="Plasma","Dilution.Factor"] <- 5
#' red[regexpr("T0",red$SampleName)!=-1,"Sample.Type"] <- "T0"
#' red$Analysis.Method <- "LC or GC"
#' red$Analysis.Instrument <- "No Idea"
#' red$Analysis.Parameters <- "None"
#'
#'
#' # Strip out protein conc information from compound names:
#' red$CompoundName <- gsub("-100P","",red$CompoundName)
#' red$CompoundName <- gsub("-30P","",red$CompoundName)
#' red$CompoundName <- gsub("-10P","",red$CompoundName)
#'
#' red$Test.Target.Conc <- 5
#' red$ISTD.Name <- "Bucetin and Diclofenac"
#' red$ISTD.Conc <- 1
#' red$Series <- 1
#'
#' level1 <- format_fup_red(red,
#'   FILENAME="Wambaugh2019",
#'   sample.col="SampleName",
#'   compound.col="Preferred.Name",
#'   lab.compound.col="CompoundName",
#'   cal.col="RawDataSet")
#'
#' level2 <- level1
#' level2$Verified <- "Y"
#'
#' write.table(level2,
#'   file="Wambaugh2019-fup-RED-Level2.tsv",
#'   sep="\t",
#'   row.names=F,
#'   quote=F)
#'
#' level3 <- calc_fup_red_point(FILENAME="Wambaugh2019")
#'
#' @references
#'  \insertRef{waters2008validation}{invitroTKstats}
#'
#' @export calc_fup_red_point
calc_fup_red_point <- function(FILENAME, good.col="Verified")
{
  MS.data <- read.csv(file=paste(FILENAME,"-fup-RED-Level2.tsv",sep=""),
    sep="\t",header=T)
  MS.data <- subset(MS.data,!is.na(Compound.Name))
  MS.data <- subset(MS.data,!is.na(Response))

  # Standardize the column names:
  sample.col <- "Lab.Sample.Name"
  date.col <- "Date"
  compound.col <- "Compound.Name"
  dtxsid.col <- "DTXSID"
  lab.compound.col <- "Lab.Compound.Name"
  type.col <- "Sample.Type"
  dilution.col <- "Dilution.Factor"
  replicate.col <- "Replicate"
  cal.col <- "Calibration"
  istd.name.col <- "ISTD.Name"
  istd.conc.col <- "ISTD.Conc"
  istd.col <- "ISTD.Area"
  std.conc.col <- "Std.Conc"
  test.nominal.conc.col <- "Test.Nominal.Conc"
  plasma.percent.col <- "Percent.Physiologic.Plasma"
  time.col <- "Time"
  area.col <- "Area"
  analysis.method.col <- "Analysis.Method"
  analysis.instrument.col <- "Analysis.Instrument"
  analysis.parameters.col <- "Analysis.Parameters"
  note.col <- "Note"
  level0.file.col <- "Level0.File"
  level0.sheet.col <- "Level0.Sheet"

# For a properly formatted level 2 file we should have all these columns:
  cols <-c(
    sample.col,
    date.col,
    compound.col,
    dtxsid.col,
    lab.compound.col,
    type.col,
    dilution.col,
    replicate.col,
    cal.col,
    istd.name.col,
    istd.conc.col,
    istd.col,
    std.conc.col,
    test.nominal.conc.col,
    plasma.percent.col,
    time.col,
    area.col,
    analysis.method.col,
    analysis.instrument.col,
    analysis.parameters.col,
    note.col,
    level0.file.col,
    level0.sheet.col,
    "Response",
    good.col)
# Throw error if not all columns present with expected names:
  if (!(all(cols %in% colnames(MS.data))))
  {
    warning("Run format_fup_red first (level 1) then curate to (level 2).")
    stop(paste("Missing columns named:",
      paste(cols[!(cols%in%colnames(MS.data))],collapse=", ")))
  }

  # Only include the data types used:
  MS.data <- subset(MS.data,MS.data[,type.col] %in% c(
    "Plasma","PBS","T0","Plasma.Blank","NoPlasma.Blank"))

  # Only used verfied data:
  MS.data <- subset(MS.data, MS.data[,good.col] == "Y")

  out.table <-NULL
  num.chem <- 0
  num.cal <- 0
  for (this.chem in unique(MS.data[,compound.col]))
  {
    this.subset <- subset(MS.data,MS.data[,compound.col]==this.chem)
    this.dtxsid <- this.subset$dtxsid[1]
    this.row <- cbind(this.subset[1,c(compound.col,dtxsid.col)],
      data.frame(Calibration="All Data",
        Fup=NaN))
    this.pbs <- subset(this.subset,Sample.Type=="PBS")
    if (dim(this.pbs)[1]==0) warning(paste0(
        "No PBS data for chemical ", this.chem))
    this.plasma <- subset(this.subset,Sample.Type=="Plasma")
    if (dim(this.plasma)[1]==0) warning(paste0(
        "No plasma data for chemical ", this.chem))
    this.plasma.blank <- subset(this.subset,Sample.Type=="Plasma.Blank")
    if (dim(this.plasma.blank)[1]==0) warning(paste0(
        "No plasma blank data for chemical ", this.chem))
    this.noplasma.blank <- subset(this.subset,Sample.Type=="NoPlasma.Blank")
    if (dim(this.noplasma.blank)[1]==0) warning(paste0(
        "No non-plasma blank data for chemical ", this.chem))
    if (length(unique(this.pbs$Dilution.Factor))>1) browser()
    df.pbs <- this.pbs$Dilution.Factor[1]
    if (length(unique(this.plasma$Dilution.Factor))>1) browser()
    df.plasma <- this.plasma$Dilution.Factor[1]
    if (length(unique(this.plasma.blank$Dilution.Factor))>1) browser()
    df.plasma.blank <- this.plasma.blank$Dilution.Factor[1]
    if (length(unique(this.noplasma.blank$Dilution.Factor))>1) browser()
    df.noplasma.blank <- this.noplasma.blank$Dilution.Factor[1]

  # Check to make sure there are data for PBS and plasma:
    if (dim(this.pbs)[1]> 0 &
        dim(this.plasma)[1] > 0 &
        dim(this.plasma.blank)[1] > 0)
    {
      num.chem <- num.chem + 1
      this.row$Fup <- signif(max(0,df.pbs*(mean(this.pbs$Response) -
        df.noplasma.blank*mean(this.noplasma.blank$Response))) /
        (df.plasma*(mean(this.plasma$Response) -
        df.plasma.blank*mean(this.plasma.blank$Response))),4)
      out.table <- rbind(out.table, this.row)
      print(paste(this.row$Compound.Name,"f_up =",signif(this.row$Fup,3)))
  # If fup is NA something is wrong, stop and figure it out:
      if(is.na(this.row$Fup)) browser()
  # If there are multiple measrument days, do separate calculations:
      if (length(unique(this.subset[,cal.col]))>1)
      {
        for (this.calibration in unique(this.subset[,cal.col]))
        {
          this.cal.subset <- subset(this.subset,
            this.subset[,cal.col]==this.calibration)
          this.row <- this.cal.subset[1,c(compound.col,dtxsid.col,cal.col)]
          this.pbs <- subset(this.cal.subset,Sample.Type=="PBS")
          this.plasma <- subset(this.cal.subset,Sample.Type=="Plasma")
          this.plasma.blank <- subset(this.cal.subset,Sample.Type=="Plasma.Blank")
          this.noplasma.blank <- subset(this.cal.subset,Sample.Type=="NoPlasma.Blank")
       # Check to make sure there are data for PBS and plasma:
          if (dim(this.pbs)[1]> 0 &
              dim(this.plasma)[1] > 0 &
              dim(this.plasma.blank)[1] > 0 &
              dim(this.noplasma.blank)[1] > 0)
          {
            this.row$Fup <- signif(max(0,df.pbs*(mean(this.pbs$Response) -
              df.noplasma.blank*mean(this.noplasma.blank$Response))) /
              (df.plasma*(mean(this.plasma$Response) -
              df.plasma.blank*mean(this.plasma.blank$Response))),4)
            out.table <- rbind(out.table, this.row)
            num.cal <- num.cal + 1
          }
        }
      } else num.cal <- num.cal + 1
    }
  }

  if (!is.null(out.table))
  {
    rownames(out.table) <- make.names(out.table$Compound.Name, unique=TRUE)
    out.table[,"Fup"] <- signif(as.numeric(out.table[,"Fup"]),3)
    out.table <- as.data.frame(out.table)
    out.table$Fup <- as.numeric(out.table$Fup)
  }

# Write out a "level 3" file (data organized into a standard format):
  write.table(out.table,
    file=paste(FILENAME,"-fup-RED-Level3.tsv",sep=""),
    sep="\t",
    row.names=F,
    quote=F)

  print(paste("Fraction unbound values calculated for",num.chem,"chemicals."))
  print(paste("Fraction unbound values calculated for",num.cal,"measurements."))

  return(out.table)
}


