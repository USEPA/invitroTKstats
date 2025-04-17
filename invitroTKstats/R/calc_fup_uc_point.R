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
#' @param data.in (Data Frame) A Level-2 data frame generated from the 
#' \code{format_fup_uc} function with a verification column added by 
#' \code{sample_verification}. Complement with manual verification if needed. 
#'
#' @param good.col (Character) Column name indicating which rows have been
#' verified, data rows valid for analysis are indicated with a "Y".
#' (Defaults to "Verified".)
#' 
#' @param output.res (Logical) When set to \code{TRUE}, the result 
#' table (Level-3) will be exported the current directory as a .tsv file. 
#' (Defaults to \code{TRUE}.)
#' 
#' @param sig.figs (Numeric) The number of significant figures to round the exported result table (Level-3). 
#' (Note: console print statements are also rounded to specified significant figures.) 
#' (Defaults to \code{3}.)
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
#' ## Load example level-2 data
#' level2 <- invitroTKstats::fup_uc_L2
#' 
#' ## scenario 1: 
#' ## input level-2 data from the R session and do not export the result table
#' level3 <- calc_fup_uc_point(data.in = level2, output.res = FALSE)
#'
#' ## scenario 2: 
#' ## import level-2 data from a 'tsv' file and export the result table
#' \dontrun{
#' ## Refer to sample_verification help file for how to export level-2 data to a directory.
#' ## Unless a different path is specified in OUTPUT.DIR,
#' ## the result table will be saved to the directory specified in INPUT.DIR.
#' ## Will need to replace FILENAME and INPUT.DIR with name prefix and location of level-2 'tsv'.
#' level3 <- calc_fup_uc_point(# e.g. replace with "Examples" from "Examples-fup-UC-Level2.tsv" 
#'                             FILENAME="<level-2 FILENAME prefix>", 
#'                             INPUT.DIR = "<level-2 FILE LOCATION>")
#' }
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
    sig.figs = 3, 
    INPUT.DIR=NULL, 
    OUTPUT.DIR = NULL)
{
  #assigning global variables
  Compound.Name <- Response <- Sample.Type <- NULL
  
  if (!missing(data.in)) {
    PPB.data <- as.data.frame(data.in)
  } else if (!is.null(INPUT.DIR)) {
    PPB.data <- read.csv(file=paste0(INPUT.DIR, "/", FILENAME,"-fup-UC-Level2.tsv"),
                         sep="\t",header=T)
    } else {
      PPB.data <- read.csv(file=paste0(FILENAME,"-fup-UC-Level2.tsv"),
                         sep="\t",header=T)
      }
  
  PPB.data <- subset(PPB.data,!is.na(Compound.Name))
  PPB.data <- subset(PPB.data,!is.na(Response))
  
  fup.uc.cols <- c(L1.common.cols,
                   test.conc.col = "Test.Compound.Conc",
                   test.nominal.conc.col = "Test.Nominal.Conc"
  )
  list2env(as.list(fup.uc.cols), envir = environment())
  cols <- c(unlist(mget(names(fup.uc.cols))), "Response", good.col)
  
  if (!any(c("Biological.Replicates", "Technical.Replicates") %in% colnames(PPB.data)))
    stop("Need at least one column representing replication, i.e. Biological.Replicates or Technical.Replicates. Run format_fup_uc first (level 1) then curate to (level 2).")
  
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
      if (!is.null(sig.figs)){
        print(paste(this.row$Compound.Name,"f_up =",signif(this.row$Fup,sig.figs)))
      } else {
        # If sig.figs = NULL, no rounding 
        print(paste(this.row$Compound.Name,"f_up =",this.row$Fup))
      }
  # If fup is NA something is wrong, stop and figure it out:
      if(is.na(this.row$Fup)){
        stop("calc_fup_uc_point - Fup value for `",this.chem,"` for the `All Data` Calibration is `NA`.")
        # browser()
      } 
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
  out.table[,"Fup"] <- as.numeric(out.table[,"Fup"])
  out.table <- as.data.frame(out.table)
  out.table$Fup <- as.numeric(out.table$Fup)
  
  if (output.res) {
    # Determine the path for output
    if (!is.null(OUTPUT.DIR)) {
      file.path <- OUTPUT.DIR
    } else if (!is.null(INPUT.DIR)) {
      file.path <- INPUT.DIR
    } else {
      file.path <- getwd()
    }
    
    rounded.out.table <- out.table 
    
    # Round results to desired number of sig figs
    if (!is.null(sig.figs)){
      rounded.out.table[,"Fup"] <- signif(rounded.out.table[,"Fup"],sig.figs)
      cat(paste0("\nData to export has been rounded to ", sig.figs, " significant figures.\n"))
    }
    
    # Write out a "level 3" file (data organized into a standard format):
    write.table(rounded.out.table,
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


