#' Calculate Point Estimates of Fraction Unbound in Plasma (Fup) with
#' Rapid Equilibrium Dialysis (RED) Data (Level-3)
#'
#' This function calculates the point estimates for the fraction unbound in plasma
#' (Fup) using mass spectrometry (MS) peak areas from samples collected as part of
#' \emph{in vitro} measurements of chemical Fup using rapid equilibrium dialysis
#' \insertCite{waters2008validation}{invitroTKstats}. See the Details section
#' for the equation(s) used in point estimation. 
#' 
#' The input to this function should be "level-2" data. Level-2 data is level-1,
#' data formatted with the \code{\link{format_fup_red}} function, and curated
#' with a verification column. "Y" in the verification column indicates the
#' data row is valid for analysis. 
#'
#' The data frame of observations should be annotated according to
#' these types:
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
#' {\frac{\sum_{i=1}^{n_{PL}} (r_{PL} * c_{DF})}{n_{PL}} - \frac{\sum_{i=1}^{n_B} (r_B * c_{DF})}{n_B}}}
#'
#' where \eqn{r_P} is PBS Response, \eqn{n_P} is the number of PBS Responses,
#' \eqn{c_{DF}} is the corresponding Dilution Factor, \eqn{r_{NPB}} is No Plasma Blank Response,
#' \eqn{n_{NPB}} is the number of No Plasma Blank Responses, \eqn{r_{PL}} is Plasma Response,
#' \eqn{n_{PL}} is the number of Plasma Responses, \eqn{r_{B}} is Plasma Blank Response,
#' and \eqn{n_B} is the number of Plasma Blank Responses.
#'
#' @param FILENAME (Character) A string used to identify the input level-2 file.
#' "<FILENAME>-fup-RED-Level2.tsv".
#' 
#' @param data.in (Data Frame) A level-2 data frame generated from the 
#' \code{format_fup_red} function with a verification column added by 
#' \code{sample_verification}. Complement with manual verification if needed.
#' 
#' @param good.col (Character) Column name indicating which rows have been
#' verified, data rows valid for analysis are indicated with a "Y".
#' (Defaults to "Verified".)
#' 
#' @param output.res (Logical) When set to \code{TRUE}, the result 
#' table (level-3) will be exported to the current directory as a .tsv file. 
#' (Defaults to \code{TRUE}.)
#' 
#' @param sig.figs (Numeric) The number of significant figures to round the exported result table (level-3). 
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
#' @return A level-3 data frame with one row per chemical, contains chemical identifiers 
#' such as preferred compound name, EPA's DSSTox Structure ID, calibration details,
#' and point estimates for the fraction unbound in plasma (Fup)
#' for all chemicals in the input data frame. 
#'
#' @author John Wambaugh
#'
#' @examples
#' ## Load example level-2 data
#' level2 <- invitroTKstats::fup_red_L2
#' 
#' ## scenario 1: 
#' ## input level-2 data from the R session and do not export the result table
#' level3 <- calc_fup_red_point(data.in = level2, output.res = FALSE)
#'
#' ## scenario 2: 
#' ## import level-2 data from a 'tsv' file and export the result table
#' \dontrun{
#' ## Refer to sample_verification help file for how to export level-2 data to a directory.
#' ## Unless a different path is specified in OUTPUT.DIR,
#' ## the result table will be saved to the directory specified in INPUT.DIR.
#' ## Will need to replace FILENAME and INPUT.DIR with name prefix and location of level-2 'tsv'.
#' level3 <- calc_fup_red_point(# e.g. replace with "Examples" from "Examples-fup-RED-Level2.tsv"
#'                              FILENAME="<level-2 FILENAME prefix>",
#'                              INPUT.DIR = "<level-2 FILE LOCATION>")
#' }
#'
#' @references
#'  \insertRef{waters2008validation}{invitroTKstats}
#'
#' @import Rdpack
#' @importFrom utils read.csv write.table
#'
#' @export calc_fup_red_point
calc_fup_red_point <- function(
    FILENAME, 
    data.in,
    good.col="Verified",
    output.res=TRUE, 
    sig.figs = 3, 
    INPUT.DIR=NULL, 
    OUTPUT.DIR = NULL)
{
  
  #assigning global variables
  Compound.Name <- Response <- Sample.Type <- Direction <- NULL
  
  if (!missing(data.in)) {
    MS.data <- as.data.frame(data.in)
  } else if (!is.null(INPUT.DIR)) {
    MS.data <- read.csv(file=paste0(INPUT.DIR, "/", FILENAME,"-fup-RED-Level2.tsv"),
                        sep="\t",header=T)
    } else {
      MS.data <- read.csv(file=paste0(FILENAME,"-fup-RED-Level2.tsv"),
                        sep="\t",header=T)
      }

  
  MS.data <- subset(MS.data,!is.na(Compound.Name))
  MS.data <- subset(MS.data,!is.na(Response))
  
  fup.red.cols <- c(L1.common.cols,
                    time.col = "Time",
                    test.conc.col = "Test.Compound.Conc",
                    test.nominal.conc.col = "Test.Nominal.Conc",
                    plasma.percent.col = "Percent.Physiologic.Plasma"
                    )
  list2env(as.list(fup.red.cols), envir = environment())
  cols <- c(unlist(mget(names(fup.red.cols))), "Response", good.col)
  
  # # Throw error if not all columns present with expected names:
  reps = c("Biological.Replicates", "Technical.Replicates")
  if (!(all(reps %in% colnames(MS.data))))
  {
    warning("Run format_fup_red first (level-1) then curate to (level-2).")
    stop(paste("Missing replication columns named:", 
               paste(reps[!(reps %in% colnames(MS.data))], collapse = ", ")))
  } else if (any(is.na(MS.data[,"Biological.Replicates"]))) 
    {
      warning("Run format_fup_red first (level-1) then curate to (level-2).")
      stop("NA values provided for Biological.Replicates")
    } 
      
  
  
  if (!(all(cols %in% colnames(MS.data))))
  {
    warning("Run format_fup_red first (level-1) then curate to (level-2).")
    stop(paste("Missing columns named:",
      paste(cols[!(cols%in%colnames(MS.data))],collapse=", ")))
  }

  # Only include the data types used:
  MS.data <- subset(MS.data,MS.data[,type.col] %in% c(
    "Plasma","PBS","T0","Plasma.Blank","NoPlasma.Blank"))

  # Only used verified data:
  MS.data <- subset(MS.data, MS.data[,good.col] == "Y")

  out.table <-NULL
  num.chem <- 0
  num.cal <- 0
  nonplasma.blanks.na.chem <- NULL
  plasma.blanks.na.chem <- NULL
  ignored.chem <- NULL
  for (this.chem in unique(MS.data[,compound.col]))
  {
    this.subset <- subset(MS.data,MS.data[,compound.col]==this.chem)
    this.dtxsid <- this.subset$dtxsid[1]
    this.row <- cbind(this.subset[1,c(compound.col,dtxsid.col)],
      data.frame(Calibration="All Data",
        Fup=NaN))
    
    this.pbs <- subset(this.subset,Sample.Type=="PBS")
    this.plasma <- subset(this.subset,Sample.Type=="Plasma")
    this.plasma.blank <- subset(this.subset,Sample.Type=="Plasma.Blank")
    this.noplasma.blank <- subset(this.subset,Sample.Type=="NoPlasma.Blank")
    
  # Check to make sure there are data for PBS and plasma:
    if (dim(this.pbs)[1]> 0 &
        dim(this.plasma)[1] > 0)
    {
      # Then check if there are any data for blanks
      if (dim(this.plasma.blank)[1]==0) {
        plasma.blank.mean <- 0
        df.plasma.blank <- 0
        # collect a list of name(s) of compounds missing non-plasma blanks
        plasma.blanks.na.chem <- c(plasma.blanks.na.chem, this.chem)
      } else {
        plasma.blank.mean <- mean(this.plasma.blank$Response)
        if (length(unique(this.plasma.blank$Dilution.Factor))>1){
          stop("calc_fup_red_point - Plasma.Blank samples for `",this.chem,"` have more than one `Dilution.Factor`.")
          # browser()
        } 
        df.plasma.blank <- this.plasma.blank$Dilution.Factor[1]
      }
      
      if (dim(this.noplasma.blank)[1]==0){
        noplasma.blank.mean <- 0
        df.noplasma.blank <- 0
        # collect a list of name(s) of compounds missing plasma blanks
        nonplasma.blanks.na.chem <- c(nonplasma.blanks.na.chem, this.chem)
      } else {
        noplasma.blank.mean <- mean(this.noplasma.blank$Response)
        if (length(unique(this.noplasma.blank$Dilution.Factor))>1){
          stop("calc_fup_red_point - No.Plasma.Blank samples for `",this.chem,"` have more than one `Dilution.Factor`.")
          # browser()
        } 
        df.noplasma.blank <- this.noplasma.blank$Dilution.Factor[1]
      }
      
      # Collect dilution factor for calculation
      if (length(unique(this.pbs$Dilution.Factor))>1){
        stop("calc_fup_red_point - PBS samples for `",this.chem,"` have more than one `Dilution.Factor`.")
        # browser()
      } 
      df.pbs <- this.pbs$Dilution.Factor[1]
      if (length(unique(this.plasma$Dilution.Factor))>1){
        stop("calc_fup_red_point - Plasma samples for `",this.chem,"` have more than one `Dilution.Factor`.")
        # browser()
      } 
      df.plasma <- this.plasma$Dilution.Factor[1]
      
      num.chem <- num.chem + 1
      fup.est <- max(0,df.pbs*mean(this.pbs$Response) - df.noplasma.blank*noplasma.blank.mean) /
        (df.plasma*mean(this.plasma$Response) - df.plasma.blank*plasma.blank.mean)
      this.row$Fup <- fup.est
      out.table <- rbind(out.table, this.row)
      if (!is.null(sig.figs)){
        print(paste(this.row$Compound.Name,"f_up =",signif(this.row$Fup,sig.figs)))
      } else {
        # If sig.figs = NULL, no rounding
        print(paste(this.row$Compound.Name,"f_up =",this.row$Fup))
      }
      # If fup is NA something is wrong, stop and figure it out:
      if(is.na(this.row$Fup)){
        stop("calc_fup_red_point - Fup value for `",this.chem,"` for the `All Data` Calibration is `NA`.")
        # browser()
      } 
      
      # If there are multiple measurement days, do separate calculations:
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
              dim(this.plasma)[1] > 0)
          {
            # Check to see if there are any blanks data
            if (dim(this.plasma.blank)[1]==0){
              plasma.blank.mean <- 0
              plasma.blanks.na.chem <- c(plasma.blanks.na.chem, paste(this.chem, "Calibration", this.calibration))
            } else {
              plasma.blank.mean <- mean(this.plasma.blank$Response)
            }
            
            if (dim(this.noplasma.blank)[1]==0){
              noplasma.blank.mean <- 0
              nonplasma.blanks.na.chem <- c(nonplasma.blanks.na.chem, paste(this.chem, "Calibration", this.calibration))
            } else {
              noplasma.blank.mean <- mean(this.noplasma.blank$Response)
            }
            fup.est <- max(0,df.pbs*mean(this.pbs$Response) - df.noplasma.blank*noplasma.blank.mean) /
              (df.plasma*mean(this.plasma$Response) - df.plasma.blank*plasma.blank.mean)
            this.row$Fup <- fup.est
            out.table <- rbind(out.table, this.row)
            if (!is.null(sig.figs)){
              print(paste(this.row$Compound.Name,"Calibration",this.calibration,"f_up =",signif(this.row$Fup,sig.figs)))
            } else {
              # If sig.figs = NULL, no rounding  
              print(paste(this.row$Compound.Name,"Calibration",this.calibration,"f_up =",this.row$Fup))
            }
            num.cal <- num.cal + 1
          } else ignored.chem <- c(ignored.chem, paste(this.chem, "Calibration", this.calibration))
        }
      } else num.cal <- num.cal + 1
    } else ignored.chem <- c(ignored.chem, this.chem)
  }
  
  ## issue notification messages
  if (!is.null(ignored.chem)) warning(paste0("The following chemical(s) were ignored due to missing PBS and/or Plasma data: ", paste(ignored.chem, collapse = ", "),"\n"))
  missingboth <- intersect(nonplasma.blanks.na.chem, plasma.blanks.na.chem)
  if (length(missingboth)!=0) warning(paste0("Plasma and non-plasma blank samples are missing for the following chemical(s): ", paste(missingboth, collapse = ", "), 
                                             ". Point estimations for these cases assume the blank adjustment is 0.\n"))
  plasma.blanks.na.chem <- setdiff(plasma.blanks.na.chem, missingboth)
  nonplasma.blanks.na.chem <- setdiff(nonplasma.blanks.na.chem, missingboth)
  if (length(nonplasma.blanks.na.chem)!=0) warning(paste0("Missing non-plasma blanks for chemical(s): ", paste(nonplasma.blanks.na.chem, collapse = ", "),
                                                         ". Point estimations for these cases assume the non-plasma blank adjustment is 0.\n"))
  if (length(plasma.blanks.na.chem)!=0) warning(paste0("Missing plasma blanks for chemical(s): ", paste(plasma.blanks.na.chem, collapse = ", "),
                                                      ". Point estimations for these cases assume the plasma blank adjustment is 0.\n"))

  if (!is.null(out.table))
  {
    rownames(out.table) <- make.names(out.table$Compound.Name, unique=TRUE)
    out.table[,"Fup"] <- as.numeric(out.table[,"Fup"])
    out.table <- as.data.frame(out.table)
    out.table$Fup <- as.numeric(out.table$Fup)
  }

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
    
    # Write out a "level-3" file (data organized into a standard format):
    write.table(rounded.out.table,
                file=paste0(file.path, "/", FILENAME,"-fup-RED-Level3.tsv"),
                sep="\t",
                row.names=F,
                quote=F)
    
  
    # Print notification message stating where the file was output to
    cat(paste0("A level-3 file named ",FILENAME,"-fup-RED-Level3.tsv", 
                " has been exported to the following directory: ", file.path), "\n")
  }

  print(paste("Fraction unbound values calculated for",num.chem,"chemicals."))
  print(paste("Fraction unbound values calculated for",num.cal,"measurements."))

  return(out.table)
}


