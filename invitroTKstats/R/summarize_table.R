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
#' @param input.table A data frame containing mass-spectrometry peak areas,
#' indication of chemical identiy, and measurment type. The data frame should
#' contain columns with names specified by the following arguments:
#' 
#' @param dtxsid.col Which column of input.table indicates EPA's DSSTox Structure 
#' ID (\url{http://comptox.epa.gov/dashboard}) (Defaults to "DTXSID")
#' 
#' @param compound.col Which column of input.table indicates the test compound
#' (Defaults to "Compound.Name")
#' 
#' @param cal.col Which column of input.table indicates the MS calibration -- for
#' instance different machines on the same day or different days with the same
#' MS analyzer (Defaults to "Calibration")
#'
#' @param type.col Which column of input.table indicates the sample type (see table
#' above)(Defaults to "Sample.Type")
#'
#' @param req.types If set (defaults to NULL) all of the measurment types 
#' included in this vecotor of character strings must
#' be available for each chemical-calibration pair to make an estimate.
#'
#' @return \item{list}{A list containing the summary counts from the table} 
#'
#' @author John Wambaugh
summarize_table <- function(input.table,
  dtxsid.col="DTXSID",
  compound.col="Compound.Name",
  cal.col="Calibration",
  type.col="Sample.Type",
  req.types=NULL)
{
  N.chems <- length(unique(input.table[,dtxsid.col]))
  N.obs <- dim(input.table)[1]
  N.meas <- length(unique(paste(input.table[,dtxsid.col],input.table[,cal.col])))
  
  cat(paste(N.obs,"observations of",N.chems,"chemicals based on",N.meas,"separate measrurements (calibrations).\n"))
  
  repeat.chems <- NULL
  N.complete <- 0
  incomplete.chems <- NULL
  for (this.chem in sort(unique(input.table[,dtxsid.col])))
  {
    this.subset <- subset(input.table,input.table[,dtxsid.col]==this.chem)
    if (length(unique(this.subset[,cal.col]))>1)
    {
      repeat.chems <- sort(unique(c(repeat.chems,this.subset[1,compound.col])))
      if (!is.null(req.types))
      {
        for (this.cal in sort(unique(this.subset[,cal.col])))
        {
          this.cal.subset <- subset(this.subset,this.subset[,cal.col]==this.cal)
          if (all(req.types %in% this.cal.subset[,type.col]))
          {
            N.complete <- N.complete + 1
          } else {
            incomplete.chems <- sort(unique(c(incomplete.chems,this.chem)))
          }
        }
      }
    } else {
      if (!is.null(req.types))
      {
        if (all(req.types %in% this.subset[,type.col]))
        {
          N.complete <- N.complete + 1
        } else {
          incomplete.chems <- sort(unique(c(incomplete.chems,this.chem)))
        }
      }
    }
  }
  
  if (!is.null(repeat.chems))
  {
    cat(paste("The following",
      length(repeat.chems),
      "chemicals have repeated observations:\n"))
    print(strwrap(paste(repeat.chems,collapse=", ")))
    cat("\n")
  }
  
  if (!is.null(incomplete.chems))
  {
    cat("The following",
      length(incomplete.chems),
      "chemicals have incomplete data sets:\n")
    print(strwrap(paste(incomplete.chems,collapse=", ")))
    cat("\n")
  }  
  
  return(list(
    N.chems=N.chems,
    N.obs=N.obs,
    N.meas=N.meas,
    N.complete=N.complete,
    repeat.chems=repeat.chems,
    incomplete.chems=incomplete.chems
    ))
}