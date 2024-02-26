#' Extract level 1 ultracentrifugation data from wide level 0 file
#'
#' This function extracts data from a Microsoft Excel file containing many
#' columns corresponding to different types of data.
#'
#' The data frame of observations should be annotated according to
#' of these types:
#' \tabular{rrrrr}{
#'   Calibration Curve \tab CC\cr
#'   Ultracentrifugation Aqueous Fraction \tab UC\cr
#'   Whole Plasma T1h Sample  \tab T1\cr
#'   Whole Plasma T5h Sample \tab T5\cr
#' }
#'
#' @param data.set (Data Frame) A data frame containing a sheet of data for conversion.
#'
#' @param chem.name (Character) A string giving the lab name of the chemical analyzed. The
#' value provided is used for all rows in the output data frame. 
#'
#' @param area.col.num (Numeric) An integer indicating which column of data.set contains
#' the MS feature area for the chemical.
#'
#' @param ISTD.name (Character) A string indicating the internal standard used. The
#' value provided is used for all rows in the output data frame. 
#'
#' @param ISTD.offset (Numeric) An integer indicating how many columns difference there
#' is between the chemical of study MS area and the ISTD MS area. (Defaults to 2.)
#'
#' @param analysis.method (Character) A string describing the chemical analysis method. The
#' value provided is used for all rows in the output data frame. 
#' (Defaults to "GC", that is gas chromatography.)
#'
#' @param instrument (Character) A string describing the instrument used for chemical 
#' analysis. The value provided is used for all rows in the output data frame. 
#' (Defaults to "Something or Other 3000".)
#'
#' @param inst.param.offset (Numeric) An integer indicating the difference in the number
#' of columns between the MS peak area and the column giving the instrument
#' parameters. (Defaults to -3.)
#'
#' @param conc.offset (Numeric) An integer indicating the difference in the number
#' of columns between the MS peak area and the column giving the intended
#' concentration for calibration curves. (Defaults to -2.)
#'
#' @param area.base (Character) A character string used for forming the name of MS
#' feature area column names (used for both test chemical and ISTD). (Defaults to
#' "Area...".)
#'
#' @param inst.param.base (Character) A character string used for forming the name of the
#' chemical analysis instrument parameter column name. (Defaults to "RT...".)
#'
#' @param conc.base (Character) A character string used for forming the name of the
#' calibration curve intended concentration column name. (Defaults to
#' "Final Conc....".)
#'
#' @param id.cols (Character Vector) A vector of character strings used for identifying each
#' sample. (Defaults to c("Name", "Data File", "Type", "Acq. Date-Time").)
#' 
#' @param type.indicator.col (Character) A character string indicating which column of 
#' data.set contains the type of observation. (Defaults to "Name".)
#' 
#' @param AF.type.str (Character) String used to annotate observation of this type: Aqueous Fraction. 
#' (Defaults to "AF".)
#' 
#' @param T1.type.str (Character) String used to annotate observation of this type: Whole Plasma T1h Sample. 
#' (Defaults to "T1".)
#' 
#' @param T5.type.str (Character) String used to annotate observation of this type: Whole Plasma T5h Sample. 
#' (Defaults to "T5".)
#' 
#' @param CC.type.str (Character) String used to annotate observation of this type: Calibration Curve. 
#' (Defaults to "CC".)
#' 
#' @return \item{data.frame}{A data.frame in standardized "level1" format}
#'
#' @author John Wambaugh
#'
#' @examples
#' library(invitroTKstats)
#' level0 <-
#'
#' @references
#' \insertRef{redgrave1975separation}{invitroTKstats}
#'
# Lets make a function to automate annotation of Wetmore lab data:
extract_level1_fup_uc <- function(
  data.set,
  chem.name,
  area.col.num,
  ISTD.name,
  ISTD.offset=2,
  analysis.method = "GC",
  instrument="Something or Other 3000",
  inst.param.offset = -3,
  conc.offset = -2,
  area.base="Area...",
  inst.param.base="RT...",
  conc.base="Final Conc....",
  id.cols=c(
    "Name",
    "Data File",
    "Acq. Date-Time"
    ),
  type.indicator.col="Name",
  AF.type.str="AF",
  T1.type.str="T1",
  T5.type.str="T5",
  CC.type.str="CC"
  )
{
  #assigning global variables
  Type <- NULL
  
  area.col <- paste(area.base,area.col.num,sep="")
  instr.param.col <- paste(inst.param.base,area.col.num+inst.param.offset,sep="")
  conc.col <- paste(conc.base,area.col.num+conc.offset,sep="")
  istd.col <- paste(area.base,area.col.num+ISTD.offset,sep="")

  out <- data.set[,c(id.cols,conc.col,area.col,istd.col,instr.param.col)]
  out$Lab.Compound.Name <- chem.name
  out$ISTD.Name <- ISTD.name
  out$Analysis.Method <- analysis.method
  out$Instrument <- instrument

  colnames(out)[colnames(out)==area.col] <- "Area"
  colnames(out)[colnames(out)==istd.col] <- "ISTD.Area"
  colnames(out)[colnames(out)==conc.col] <- "Target.Conc"
  colnames(out)[colnames(out)==instr.param.col] <- "Analysis.Instrument.Param"

  out[,"Type"] <- ""
  out[regexpr(AF.type.str,unlist(out[,type.indicator.col]))!=-1,"Type"]<-"AF"
  out[regexpr(CC.type.str,unlist(out[,type.indicator.col]))!=-1,"Type"]<-"CC"
  out[regexpr(T1.type.str,unlist(out[,type.indicator.col]))!=-1,"Type"]<-"T1"
  out[regexpr(T5.type.str,unlist(out[,type.indicator.col]))!=-1,"Type"]<-"T5"

  out <- subset(out,Type!="")

  return(out)
}
