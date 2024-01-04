#' Calculate Point Estimates of Fraction Unbound in Plasma with
#' Ultracentrifugation (UC) Data 
#'
#' This function calculates the point estimates for the fraction unbound in
#' plasma (Fup) using mass spectrometry (MS) peak areas from samples collected
#' as part of in vitro measurements of chemical Fup using ultracentrifugation
#' \insertCite{waters2008validation}{invitroTKstats}. See the Details section
#' for the equation(s) used in the point estimate.
#' 
#' The input to this function should be "Level-2" data. Level-2 data is Level-1,
#' data formatted with the \code{\link{format_fup_uc}} function, and curated
#' with a verification column. "Y" in the verification column indicates the
#' data row is valid for analysis. 
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
#' \eqn{f_{up}} is calculated from MS responses as:
#'
#' \eqn{f_{up} = \frac{\sum_{i = 1}^{n_A} (r_A * c_{DF}) / n_A}{\sum_{i = 1}^{n_{T5}} (r_{T5} * c_{DF}) / n_{T5}}}
#'
#' where \eqn{r_A} is Aqueous Fraction Response, \eqn{c_{DF}} is Dilution Factor,
#' \eqn{r_{T5}} is T5 Response, \eqn{n_A} is the number of Aqueous Fraction Responses,
#' and \eqn{n_{T5}} is the number of T5 Responses.
#'
#' @param FILENAME (Character) A string used to identify the input Level-2 file.
#' "<FILENAME>-fup-UC-Level2.tsv".
#' 
#' @param data.in (Data Frame) A Level-2 data frame containing
#' mass-spectrometry peak areas, indication of chemical identity,
#' and measurement type.
#'
#' @param good.col (Character) Column name indicating which rows have been
#' verified, data rows valid for analysis are indicated with a "Y".
#' (Defaults to "Verified".)
#' 
#' @param output.res (Logical) When set to \code{TRUE}, the result 
#' table (Level-3) will be exported the current directory as a .tsv file. 
#' (Defaults to \code{TRUE}.)
#' 
#' @param INPUT.DIR (Character) Path to the directory where the input level-2 file exists. 
#' If \code{NULL}, looking for the input level-2 file in the current working
#' directory. (Defaults to \code{NULL}.)
#' 
#' @param OUTPUT.DIR (Character) Path to the directory to save the output file. 
#' If \code{NULL}, the output file will be saved to the current working
#' directory or \code{INPUT.DIR} if specified. (Defaults to \code{NULL}.)
#' 
#' @return A data frame with one row per chemical, contains chemical identifiers 
#' such as preferred compound name, compound name used by the laboratory, 
#' EPA's DSSTox Structure ID, calibration, and point estimates for
#' the fraction unbound in plasma (Fup) for all chemicals in the input data frame. 
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
#' \insertRef{redgrave1975separation}{invitroTKstats}
#'
#' @import Rdpack
#'
#' @export calc_fup_uc_point
calc_fup_uc_point <- function(
    FILENAME, 
    data.in,
    good.col="Verified", 
    output.res=TRUE, 
    INPUT.DIR=NULL, 
    OUTPUT.DIR = NULL)
{
  
  if (!missing(data.in)) {
    PPB.data <- as.data.frame(data.in)
  } else if (missing(data.in)) {
    if (!is.null(INPUT.DIR)) {
    PPB.data <- read.csv(file=paste0(INPUT.DIR, "/", FILENAME,"-fup-UC-Level2.tsv"),
                         sep="\t",header=T)
    } else {
    PPB.data <- read.csv(file=paste0(FILENAME,"-fup-UC-Level2.tsv"),
                         sep="\t",header=T)
    }
  }
  
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
  
  if (output.res) {
    # Write out a "level 3" file (data organized into a standard format):
    # Determine the path for output
    
    if (!is.null(OUTPUT.DIR)) {
      file.path <- OUTPUT.DIR
    } else if (!is.null(INPUT.DIR)) {
      file.path <- INPUT.DIR
    } else {
      file.path <- getwd()
    }
    write.table(out.table,
                file=paste0(file.path, "/", FILENAME,"-fup-UC-Level3.tsv"),
                sep="\t",
                row.names=F,
                quote=F)
    
    # Print notification message stating where the file was output to
    cat(paste0("A Level-3 file named ",FILENAME,"-fup-UC-Level3.tsv", 
                " has been exported to the following directory: ", file.path), "\n")
  }
  

  print(paste("Fraction unbound values calculated for",num.chem,"chemicals."))
  print(paste("Fraction unbound values calculated for",num.cal,"measurements."))

  return(out.table)
}


