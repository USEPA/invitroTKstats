#' Creates a Summary Table of MS Data Frame
#'
#' This function creates and returns a list containing summary counts of the input data table.
#' Summary counts include the counts of observations, unique chemicals and unique 
#' measurements in the input data table. If a vector of data types is specified in 
#' the argument \code{req.types}, the function also checks if each chemical has 
#' observations for every measurement type included the vector for each chemical-calibration pair. 
#' If it dose, the chemical is said to have a complete data set. Otherwise, it has an incomplete data set. 
#' The counts of chemicals with complete and incomplete data sets are returned in the output list. 
#' The input data frame can be Caco-2 data, ultracentrifugation (UC) data, rapid equilibrium dialysis (RED) data, 
#' or hepatocyte clearance (Clint) data. Tables of measurement types and annotations 
#' used in each assay are available in Details. 
#'
#' Sample types used in ultracentrifugation (UC) data should be 
#' annotated according to this table:
#' \tabular{rrrrr}{
#'   Calibration Curve \tab CC\cr
#'   Ultracentrifugation Aqueous Fraction \tab AF\cr
#'   Whole Plasma T1h Sample  \tab T1\cr
#'   Whole Plasma T5h Sample \tab T5\cr
#' }
#' 
#' Sample types in rapid equilibrium dialysis (RED) data should be annotated according to
#' this table:
#' \tabular{rrrrr}{
#'   Blank (ignored) \tab Blank\cr
#'   Plasma well concentration \tab Plasma\cr
#'   Phosphate-buffered well concentration\tab PBS\cr
#'   Time zero plasma concentration \tab T0\cr
#'   Plasma stability sample \tab Stability\cr
#'   Equilibrium Control Well 1 \tab EC1\cr
#'   Equilibrium Control Well 2 \tab EC2\cr
#'   Calibration Curve \tab CC\cr
#' }
#'
#' Sample types in hepatocyte clearance data should be annotated according to
#' this table:
#' \tabular{rrrrr}{
#'   Blank \tab Blank\cr
#'   Hepatocyte incubation concentration \tab Cvst\cr
#'   Inactivated Hepatocytes \tab Inactive\cr
#'   Calibration Curve \tab CC\cr
#' }
#' 
#' Sample types used in Caco-2 data should be annotated according to
#' this table:
#' \tabular{rr}{
#'   Blank with no chemical added \tab Blank \cr
#'   Dosing vehicle (C0) at target concentration \tab D0\cr
#'   Donor compartment at end of experiment \tab D2\cr
#'   Receiver compartment at end of experiment\tab R2\cr
#' }
#'
#' @param input.table (Data Frame) A data frame containing mass-spectrometry peak areas,
#' indication of chemical identity, and measurement type. The data frame should
#' contain columns with names specified by the following arguments:
#'
#' @param dtxsid.col (Character) Column name of input.table indicates EPA's DSSTox Structure
#' ID (\url{http://comptox.epa.gov/dashboard}). (Defaults to "DTXSID")
#'
#' @param compound.col (Character) Column name of input.table containing the test compound.
#' (Defaults to "Compound.Name".)
#'
#' @param cal.col (Character) Column name of input.table containing the MS calibration. 
#' Calibration typically uses indices or dates to represent if the analyses were done on 
#' different machines on the same day or on different days with the same MS analyzer. 
#' (Defaults to "Calibration".)
#'
#' @param type.col (Character) Column name of input.table containing the sample type (see tables
#' in Details). (Defaults to "Sample.Type".)
#'
#' @param req.types (Character Vector) A vector of character strings contains
#' measurement types. If used, the function checks if all of the measurement types included in this vector are
#' available for each chemical-calibration pair. (Defaults to \code{NULL}.)
#'
#' @return A list containing the summary counts from the input data table. The list includes: 
#' the number of observations, the number of unique chemicals, the number of unique measurements, 
#' the number of chemicals with complete data sets, the number of chemicals with incomplete data sets, 
#' and the number of chemicals with repeated observations.
#'
#' @author John Wambaugh
#' @export summarize_table
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

  cat(paste(N.obs,"observations of",N.chems,"chemicals based on",N.meas,"separate measurements (calibrations).\n"))

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
