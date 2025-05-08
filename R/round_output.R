#' Round Numeric Data (Any Level and Assay)
#' 
#' This function rounds the numeric columns from any level of processing. Numeric
#' columns may include estimates of chemical-specific toxicokinetic (TK) parameters 
#' from the relevant \emph{in vitro} assays or numerical data measurements collected 
#' from the mass spectrometry experiments.  
#' 
#' For example, for level-3 or level-4 output results, estimates 
#' of intrinsic hepatic clearance (Cl~int~) from Hepatocyte Incubation data,
#' fraction unbound in plasma (F~up~) from Rapid Equilibrium Dialysis (RED) data, 
#' fraction unbound in plasma (F~up~) from Ultracentrifugation (UC) data, or
#' apparent membrane permeability from a Caco-2 assay can all be rounded to the 
#' desired number of significant figures. 
#' 
#' Note: Currently, for level-3 Caco-2 data, the "Frec_A2B.vec" and "Frec_B2A.vec" columns are 
#' not rounded. However, these columns can be rounded if the level-3 result table 
#' from \code{\link{calc_caco2_point}} is exported and the number of significant 
#' figures is specified. 
#' 
#' The input to this function can be any level of data (level-0 through level-4)
#' corresponding to any assay (Clint, Caco-2, Fup RED, Fup UC). The desired data object 
#' to be rounded can be a data.frame, specified with \code{data.in}, or a .tsv or 
#' .RData, specified with \code{FULL_FILENAME}. 
#' 
#' If the rounded output file is chosen to be exported and an output directory is
#' not specified, it will be exported to the user's R session temporary directory.
#' This temporary directory is a per-session directory whose path can be found 
#' with the following code: \code{tempdir()}. For more details, see 
#' \url{https://www.collinberke.com/til/posts/2023-10-24-temp-directories/}. 
#' 
#' As a best practice, \code{INPUT.DIR} (when importing a .tsv or .RData file) 
#' and/or \code{OUTPUT.DIR} should be specified to simplify the process of importing
#' and exporting files. This practice ensures that the exported files can easily 
#' be found and will not be exported to a temporary directory. 
#' 
#' @param FULL_FILENAME (Character) A string used to identify the full filename of
#' input .tsv or .RData file (i.e. "MYDATA-Clint-Level4.tsv" or "MYDATA-Clint-Level4Analysis-2025-04-23.RData").
#' The string is also used to name the exported data file (if chosen to be exported). 
#' (Note: \code{FULL_FILENAME} not required if \code{data.in} is provided.) 
#' (Defaults to \code{NULL}.)  
#' 
#' @param data.in (Data Frame) Any level data frame generated from \code{invitroTKstats}
#' package. 
#' (Note: \code{data.in} not required if \code{FULL_FILENAME} is provided.)
#' 
#' @param FILENAME (Character) A string used to name the start of the exported 
#' date file. Only required if input data is a data.frame and output file is being 
#' exported.
#' (Defaults to "MYDATA".)
#' 
#' @param assay (Character) A string used to name the assay used to generate the
#' input data. The string is appended to the name of the exported data file. Only
#' required if input data is a data.frame and output file is being exported. 
#' Must be one of the following assays: "Clint", "Caco-2", "fup-RED", or "fup-UC". 
#' (Defaults to \code{NULL}.)
#' 
#' @param level (Character) A string used the name the level of the input data. 
#' The string is appended to the name of the exported data file. Only required if 
#' input data is a data.frame and output file is being exported.  
#' Must be one of the following levels: "0", "1", "2", "3", "4". 
#' (Defaults to \code{NULL}.)
#' 
#' @param exclusion.cols (Character) Vector of column names to exclude from rounding. 
#' (Defaults to \code{NULL}.)
#' 
#' @param sig.figs (Numeric) The number of significant figures to round the desired 
#' numeric columns to. 
#' (Defaults to \code{3}.)
#' 
#' @param output.res (Logical) When set to \code{TRUE}, the rounded data file will 
#' be exported to the user's per-session temporary directory as a .tsv (if \code{data.in} is specified
#' or if \code{FULL_FILENAME} is a .tsv) or as an .RData (if \code{FULL_FILENAME}
#' is an .RData). 
#' (Defaults to \code{FALSE}.)
#' 
#' @param INPUT.DIR (Character) Path to the directory where the \code{FULL_FILENAME} exists. 
#' If \code{NULL}, looking for the input \code{FULL_FILENAME} in the current working
#' directory. 
#' (Defaults to \code{NULL}.)
#' 
#' @param OUTPUT.DIR (Character) Path to the directory to save the rounded data file. 
#' If \code{NULL}, the output file will be saved to the user's per-session temporary
#' directory or \code{INPUT.DIR} if specified. (Defaults to \code{NULL}.)
#' 
#' @return A rounded data frame 
#' 
#' @author Lindsay Knupp
#' 
#' @examples
#' ## Round Clint-L4 data, exclude p-value columns, and don't export results 
#' level4 <- invitroTKstats::clint_L4
#' round_output(data.in = level4,
#'              exclusion.cols = c("Clint.pValue", "Sat.pValue", "degrades.pValue"),
#'              output.res = FALSE)
#' 
#' ## Round Clint-L4 data and export results. 
#' ## Note: Will export as a .tsv file. 
#' \dontrun{
#' round_output(data.in = level4, assay = "Clint", level = "4")
#' }
#' 
#' ## Round Clint-L4 .tsv data and export to INPUT.DIR. 
#' ## Will need to replace FULL_FILENAME and INPUT.DIR with full filename and location of .tsv. 
#' \dontrun{
#' round_output(FULL_FILENAME = "Example-Clint-Level4.tsv", 
#'              INPUT.DIR = "<FULL_FILENAME FILE LOCATION>")
#' }
#' 
#' ## Round Clint-L4 .RData and export to OUTPUT.DIR 
#' ## Will need to replace FULL_FILENAME and INPUT.DIR with full filename and location
#' ## of .RData. Will also need to replace OUTPUT.DIR with desired location of rounded 
#' ## data file. 
#' \dontrun{
#' round_output(FULL_FILENAME = "Example-Clint-Level4Analysis-2025-04-17.RData",
#'              INPUT.DIR = "<FULL_FILENAME FILE LOCATION>",
#'              OUTPUT.DIR = "<DESIRED ROUNDED FILE LOCATION>")
#' }
#' 
#' @export round_output
round_output <- function(FULL_FILENAME = NULL,
                         data.in, 
                         FILENAME = "MYDATA", 
                         assay = NULL, 
                         level = NULL,
                         exclusion.cols = NULL,
                         sig.figs = 3,
                         output.res = FALSE,
                         INPUT.DIR = NULL, 
                         OUTPUT.DIR = NULL){
  
  # Extract file type from FULL_FILENAME. If FULL_FILENAME not provided (i.e. 
  # reading in a data.frame), assign .tsv in order to write to .tsv file 
  # (if output.res = TRUE)
  if (!is.null(FULL_FILENAME)) {
    split_FULL_FILENAME <- unlist(strsplit(FULL_FILENAME, split = "\\."))
    file_type <- split_FULL_FILENAME[[2]]
  } else{
    file_type <- "tsv"
  }

  # Read in input data. Supported file types include data.frame, .tsv, or .RData
  if (!missing(data.in)){
    input.table <- as.data.frame(data.in)
  } else if (!is.null(INPUT.DIR)) { # INPUT.DIR is specified 
      if (file_type == "tsv") input.table <- read.csv(file = paste(INPUT.DIR, FULL_FILENAME, sep = "/"), sep = "\t", header = TRUE)
      else input.table <- get(load(file = paste(INPUT.DIR, FULL_FILENAME, sep = "/")))
  } else { # INPUT.DIR is not specified 
      if (file_type == "tsv") input.table <- read.csv(file = FULL_FILENAME, sep = "\t", header = TRUE)
      else input.table <- get(load(file = FULL_FILENAME))
  }
  
  output.table <- input.table
  col.names <- colnames(output.table)
  
  # Numeric columns 
  numeric.cols <- col.names[sapply(output.table, is.numeric)]
  
  # Don't round columns provided in exclusion.cols 
  if (!is.null(exclusion.cols)){
    
    if (!all(exclusion.cols %in% numeric.cols)){
      stop(paste("\nExclusion columns: ", 
                 exclusion.cols[!exclusion.cols %in% numeric.cols], "not in numeric column names of input data." ))
    }
    # Remove exclusion.cols from all column names 
    indices <- which(numeric.cols %in% exclusion.cols)
    numeric.cols <- numeric.cols[-c(indices)]
  }
  # Round the numeric columns 
  rounded_cols <- signif(output.table[,numeric.cols], sig.figs)
  output.table[,numeric.cols] <- rounded_cols
  cat(paste0("\nData in ", paste(numeric.cols, collapse = ", "), " has been rounded to ", sig.figs, " significant figures."))
  
  # Export the data as the same file 
  if (output.res){
  
    # Create the exported filename 
    if (!is.null(FULL_FILENAME)) {
      full_filename <-  split_FULL_FILENAME[[1]]
      }
    else {
      approved_assays <- c("Clint", "Caco-2", "fup-RED", "fup-UC")
      if (is.null(assay) | !is.character(assay)) stop(paste0("If exporting data, must provide one of the approved assays: ", paste(approved_assays, collapse = ", "), " corresponding to input data.frame."))
      if (!is.null(assay)){
        if (!assay %in% approved_assays) stop(paste0("If exporting data, must provide one of the approved assays: ", paste(approved_assays, collapse = ", "), " corresponding to input data.frame."))
      }
      approved_levels <- c("0", "1", "2", "3", "4")
      if (is.null(level) | !is.character(level)) stop(paste0("If exporting data, must provide one of the approved levels: ", paste(approved_levels, collapse = ", "), " corresponding to input data.frame."))
      if (!is.null(level)){
        if (!level %in% approved_levels) stop(paste0("If exporting data, must provide one of the approved levels: ", paste(approved_levels, collapse = ", "), " corresponding to input data.frame."))
      }
      full_filename <- paste0(FILENAME, "-", assay, "-Level", level) 
    } 

    # Create the exported filepath 
    if (!is.null(OUTPUT.DIR)) file_path <- paste0(OUTPUT.DIR, "/", full_filename, "-rounded")  
    else if (!is.null(INPUT.DIR)) file_path <- paste0(INPUT.DIR, "/", full_filename, "-rounded")  
    else file_path <- paste0(tempdir(), "\\", full_filename, "-rounded")
    
    # Create the correct file extension
    if (file_type == "tsv") {
      write.table(output.table, file = paste0(file_path,".tsv"),
                                                      sep = "\t",
                                                      row.names = F,
                                                      quote = F)
      cat(paste0("\nData has been saved to ", paste0(file_path, ".tsv")))
    } else {
      save(output.table, file = paste0(file_path,".RData"))
      cat(paste0("\nData has been saved to ", paste0(file_path, ".RData")))
    }
  }
  return(output.table)
}
