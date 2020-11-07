#' Calculate a point estimate of Fraction Unbound in Plasma (UC)
#'
#' This function use describing mass spectrometry (MS) peak areas
#' from samples collected as part of in vitro measurement of chemical fraction
#' unbound in plasma using ultracentrifugation (Redgrave 1975?).
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
#'
#' F_up is calculated from MS responses as:
#'
#' f_up = mean(AF Response * Dilution.Factor) / mean(T5 Response * Dilution Factor)
#'
#' @param FILENAME A string used to identify the input file, whatever the
#' argument given, "-PPB-RED-Level2.tsv" is appended (defaults to "MYDATA")
#' 
#' @param good.col Name of a column indicating which rows have been verified for 
#' analysis, indicated by a "Y" (Defaults to "Verified")
#'
#' @return \item{data.frame}{A data.frame in standardized format} 
#'
#' @author John Wambaugh
#' 
#' @references
#' Redgrave, T. G., D. C. K. Roberts, and C. E. West. "Separation of plasma 
#' lipoproteins by density-gradient ultracentrifugation." Analytical 
#' Biochemistry 65.1-2 (1975): 42-49.#' 
#'
#' @export calc_fup_uc_point
calc_fup_red_point <- function(FILENAME, good.col="Verified")
{
  PPB.data <- read.csv(file=paste(FILENAME,"-PPB-UC-Level2.tsv",sep=""), 
    sep="\t",header=T)  
  PPB.data <- subset(PPB.data,!is.na(Compound.Name))
  PPB.data <- subset(PPB.data,!is.na(Response))

# Standardize the column names:
  sample.col <- "Lab.Sample.Name"
  date.col <- "Date"
  compound.col <- "Compound.Name"
  dtxsid.col <- "DTXSID"
  lab.compound.col <- "Lab.Compound.Name"
  type.col <- "Sample.Type"
  dilution.col <- "Dilution.Factor"
  compound.conc.col <- "Nominal.Conc"
  cal.col <- "Calibration"
  nominal.test.conc.col <- "Test.Target.Conc"
  istd.name.col <- "ISTD.Name"
  istd.conc.col <- "ISTD.Conc"
  istd.col <- "ISTD.Area"
  series.col <- "Series"
  area.col <- "Area"

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
    compound.conc.col,
    nominal.test.conc.col,
    istd.name.col,
    istd.conc.col,
    istd.col,
    series.col,
    area.col,
    "Response",
    good.col)
  if (!(all(cols %in% colnames(PPB.data))))
  {
    warning("Run format_fup_red first (level 1) then curate to (level 2).")
    stop(paste("Missing columns named:",
      paste(cols[!(cols%in%colnames(PPB.data))],collapse=", ")))
  }

  # Only include the data types used:
  PPB.data <- subset(PPB.data,PPB.data[,type.col] %in% c(
    "CC",
    "T1",
    "T5",
    "AF"))   
  
  # Only used verfied data:
  PPB.data <- subset(PPB.data, PPB.data[,good.col] == "Y")

  out.table <-NULL
  num.chem <- 0
  num.cal <- 0
  for (this.chem in unique(PPB.data[,compound.col]))
  {
    this.subset <- subset(PPB.data,PPB.data[,compound.col]==this.chem)
    this.dtxsid <- this.subset$dtxsid[1]
    this.row <- c(this.subset[1,c(compound.col,dtxsid.col)],
      data.frame(Calibration="All Data",
        Fup=NaN))
    this.af <- subset(this.subset,Sample.Type=="AF")
    this.t5 <- subset(this.subset,Sample.Type=="T5")
 # Check to make sure there are data for PBS and plasma: 
    if (dim(this.pbs)[1]> 0 & dim(this.plasma)[1] > 0 )
    {
      num.chem <- num.chem + 1
      this.row$Fup <- mean(this.af$Response*this.af$Dilution.Factor) /
        mean(this.t5$Response*this.r5$Dilution.Factor)
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
          this.af<- subset(this.cal.subset,Sample.Type=="AF")
          this.t5 <- subset(this.cal.subset,Sample.Type=="T5")
       # Check to make sure there are data for PBS and plasma: 
          if (dim(this.pbs)[1]> 0 & dim(this.plasma)[1] > 0 )
          {
            this.row$Fup <- mean(this.af$Response*this.af$Dilution.Factor) /
              mean(this.t5$Response*this.t5$Dilution.Factor)
            out.table <- rbind(out.table, this.row)
            num.cal <- num.cal + 1
          }
        }
      } else num.cal <- num.cal + 1
    }
  }

# Write out a "level 1" file (data organized into a standard format):  
  write.table(PPB.data, 
    file=paste(FILENAME,"-PPB-RED-Level3.tsv",sep=""),
    sep="\t",
    row.names=F,
    quote=F)
 
  print(paste("Fraction unbound values calculated for",num.chem,"chemicals."))
  print(paste("Fraction unbound values calculated for",num.cal,"measurements."))

  return(out.table)  
}


