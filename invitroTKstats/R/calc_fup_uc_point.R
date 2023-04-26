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
#' argument given, "-PPB-UC-Level2.tsv" is appended (defaults to "MYDATA")
#' 
#' @param good.col Name of a column indicating which rows have been verified for 
#' analysis, indicated by a "Y" (Defaults to "Verified")
#'
#' @return \item{data.frame}{A data.frame in standardized format} 
#'
#' @author John Wambaugh
#'
#' @examples
#' level0 <- kreutz2020
#' level0$Analysis.Method <- "GC" 
#' level0$Analysis.Instrument <- "No Idea"
#' level0$Analysis.Parameters <- "None"
#' level1 <- format_fup_uc(level0,
#'   FILENAME="Kreutz2020",
#'   compound.col="Name",
#'   compound.conc.col="Standard.Conc",
#'   area.col="Chem.Area"
#'   )
#' level2 <- level1
#' level2$Verified <- "Y"
#' 
#' write.table(level2,
#'   file="Kreutz2020-fup-UC-Level2.tsv",
#'   sep="\t",
#'   row.names=F,
#'   quote=F)
#'
#' level3 <- calc_fup_uc_point(FILENAME="Kreutz2020") 
#'
#' @references
#' Redgrave, T. G., D. C. K. Roberts, and C. E. West. "Separation of plasma 
#' lipoproteins by density-gradient ultracentrifugation." Analytical 
#' Biochemistry 65.1-2 (1975): 42-49.#' 
#'
#' @export calc_fup_uc_point
calc_fup_uc_point <- function(FILENAME, good.col="Verified")
{
  PPB.data <- read.csv(file=paste(FILENAME,"-fup-UC-Level2.tsv",sep=""), 
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
  if (!(all(cols %in% colnames(PPB.data))))
  {
    warning("Run format_fup_uc first (level 1) then curate to level 2.")
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
    this.dtxsid <- this.subset$DTXSID[1]
    this.row <- cbind(this.subset[1,c(compound.col,dtxsid.col,lab.compound.col)],
      data.frame(Calibration="All Data",
        Fup=NaN))
    this.af <- subset(this.subset,Sample.Type=="AF")
    this.t5 <- subset(this.subset,Sample.Type=="T5")
 # Check to make sure there are data for PBS and plasma: 
    if (dim(this.af)[1]> 0 & dim(this.t5)[1] > 0 )
    {
      num.chem <- num.chem + 1
      this.row$Fup <- mean(this.af$Response*this.af$Dilution.Factor) /
        mean(this.t5$Response*this.t5$Dilution.Factor)
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
          this.row <- this.cal.subset[1,c(compound.col,dtxsid.col,lab.compound.col,cal.col)]
          this.af<- subset(this.cal.subset,Sample.Type=="AF")
          this.t5 <- subset(this.cal.subset,Sample.Type=="T5")
       # Check to make sure there are data for PBS and plasma: 
          if (dim(this.af)[1]> 0 & dim(this.t5)[1] > 0 )
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

  rownames(out.table) <- make.names(out.table$Compound.Name, unique=TRUE)
  out.table[,"Fup"] <- signif(as.numeric(out.table[,"Fup"]),3) 
  out.table <- as.data.frame(out.table)
  out.table$Fup <- as.numeric(out.table$Fup)
 
# Write out a "level 3" file (data organized into a standard format):  
  write.table(out.table, 
    file=paste(FILENAME,"-fup-UC-Level3.tsv",sep=""),
    sep="\t",
    row.names=F,
    quote=F)
 
  print(paste("Fraction unbound values calculated for",num.chem,"chemicals."))
  print(paste("Fraction unbound values calculated for",num.cal,"measurements."))

  return(out.table)  
}


