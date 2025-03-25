#' Add Sample Verification Column (Level-2)
#'
#' This function takes in a Level-1 data frame and an exclusion list and 
#' returns a Level-2 data frame with a verification column. The
#' verification column contains either "Y", indicating the row is good for analysis,
#' or messages contained in the exclusion list for why the data rows are excluded. 
#' If an exclusion list is not provided, all rows are assumed to be good for use 
#' in further analyses and are verified with "Y".
#' 
#' The `exclusion.info` should be a data frame with the following columns:
#' \tabular{rr}{
#'   Variables \tab Level-1 variable(s) used to filter rows for exclusion\cr
#'   Values \tab Value(s) to exclude\cr
#'   Message \tab Simple explanation for the exclusion\cr
#' }
#' When filtering on multiple variable-value pairs, the character input for 
#' "Variables" and "Values" should be separated by a vertical bar "|" ,
#' and the variable-value pairs should match. See demonstration in Examples, Scenario 1. 
#'
#' @param FILENAME (Character) A string used to identify the output Level-1 file.
#' "<FILENAME>-<assay>-Level1.tsv". 
#' 
#' @param data.in (Data Frame) A Level-1 data frame from the format functions.
#'
#' @param exclusion.info (Data Frame) A data frame containing the variables and 
#' values of the corresponding variables to exclude rows. 
#' See details for full explanation.
#' 
#' @param assay (Character) A string indicating what assay data the input file is. Valid 
#' input is one of the following: "Clint", "fup-UC", "fup-RED", or "Caco-2". 
#' This argument only needs to be specified when importing input data set with \code{FILENAME} 
#' or exporting a data file.
#'
#' @param output.res (Logical) When set to \code{TRUE}, the result 
#' data frame (Level-2) will be exported as a .tsv file to the current directory. 
#' (Defaults to \code{TRUE}.)
#'
#' @param INPUT.DIR (Character) Path to the directory where the input level-1 file exists. 
#' If \code{NULL}, looking for the input level-1 file in the current working
#' directory. (Defaults to \code{NULL}.)
#' 
#' @param OUTPUT.DIR (Character) Path to the directory to save the output file. 
#' If \code{NULL}, the output file will be saved to the current working
#' directory or \code{INPUT.DIR} if specified. (Defaults to \code{NULL}.)
#' 
#' @return A Level-2 data frame with a verification column. 
#' 
#' @importFrom utils read.csv write.table
#' 
#' @export
#' 
#' @examples 
#' level2 <- invitroTKstats::kreutz2023.clint
#' 
#' # Data sets in the package are level-2 data. For demonstration purposes,
#' # remove the verification columns from the data and use it as level-1.
#' level1 <- subset(level2, select = -c(Verifed,Verified) )
#' 
#' # Scenario 1: Pass in data.in and exclusion.info data frame from R session 
#' 
#' # Create a exclusion criteria data frame
#' # If more than one variable is used to define a set of samples to be excluded,
#' # enter them as one string, separate the variables with a vertical bar, "|", and do the same for 
#' # values. 
#' exclusion_criteria <- data.frame(
#'   Variables = c(
#'     "DTXSID","Date","Compound.Name|Lab.Sample.Name"
#'   ),
#'   Values = c(
#'     "DTXSID00380798","21021",
#'     "Nonafluoropentanamide|Amide Hep120220 T0C2_020921"
#'   ),
#'   Message = c(
#'     "Exclude all samples of this compound.",
#'     "Exclude all samples from this date.",
#'     "These samples for this compound were contaminated."
#'     )
#' )
#' 
#' # Run the verification function.
#' my.level2 <- sample_verification(data.in=level1,
#'                                  exclusion.info = exclusion_criteria,
#'                                  output.res = FALSE)
#' 
#' # scenario 2: import 'tsv' as input data and do not pass in a exclusion.info data frame
#' 
#' \dontrun{
#' # Write the level-1 file to some folder
#' write.table(level1,
#' file="~/invitrotkstats/invitroTKstats/data-raw/kreutz-Clint-Level1.tsv",
#' sep="\t",
#' row.names=F,
#' quote=F)
#' 
#' # Run the verification function.
#' # Specify the path to import level-1 data with INPUT.DIR.
#' # No exclusion.info data frame used will label all samples as verified.
#' my.level2 <- sample_verification(FILENAME="kreutz", 
#' assay="Clint", INPUT.DIR = "~/invitrotkstats/invitroTKstats/data-raw")
#' }
#' 
#' @author Zhihui (Grace) Zhao
#' @import dplyr
#' 
sample_verification <- function(
    FILENAME, 
    data.in, 
    exclusion.info,
    assay,
    output.res = TRUE,
    INPUT.DIR = NULL,
    OUTPUT.DIR = NULL
    ){
  
  approved_assays <- c("Clint", "Caco-2", "fup-UC", "fup-RED")
  # if either importing or exporting data file, check if the assay given is valid.
  if (missing(data.in) | output.res)
    if (!(assay %in% approved_assays))
    stop("Invalid assay. ", "Use one of the approved assays: ", paste(approved_assays, collapse = ", "), ".")
  
  if (!missing(data.in)) {
    data.out <- as.data.frame(data.in)
    rm(data.in)
    } else if (!missing(assay) & !missing(FILENAME)) {
      if (!is.null(INPUT.DIR)) {
        data.out <- read.csv(file=paste0(INPUT.DIR, "/", FILENAME,"-", assay, "-Level1.tsv"),
                     sep="\t",header=T)
        } else {
          data.out <- read.csv(file=paste0(FILENAME,"-", assay, "-Level1.tsv"),
                          sep="\t",header=T)  
          }
      } else {
        stop(strwrap("A valid input data must be provided. If using a data frame, data.in is not specified. 
                           If importing a data file, missing either FILENAME and/or assay. 
                           Unable to import data from the 'tsv' without both FILENAME and assay."))
    }
  
  # Add the verification column with all "Y"
  data.out$Verified <- "Y"
  
  if (!missing(exclusion.info)) {
    exclusion.info <- as.data.frame(exclusion.info)
    for (i in 1:nrow(exclusion.info)) {
      ## split the list, consider all possible separators
      ## exclude period and underscore - period is used in column names and underscore can be used in sample names
      var.list <- strsplit(exclusion.info[i, "Variables"], "\\|", perl = TRUE)[[1]]
      var.list <- trimws(var.list)
      value.list <- strsplit(exclusion.info[i, "Values"], "\\|", perl = TRUE)[[1]]
      value.list <- trimws(value.list)
      ## check if every variable has a matched value 
      if (length(var.list) != length(value.list))
        stop("The lengths of variable list and value list do not match.")
      ## check if the variable names are valid
      if (!all(var.list %in% colnames(data.out))) 
        stop("Names of the variables use to determine the exclusion criteria do not match the column names of the level-1 data.")
      
      ## concatenate exclusion criteria into one string
      ec <- paste(paste(var.list, paste0("\"",value.list,"\""), sep = "=="), collapse = "&")
      ## filter the rows and attach exclusion message 
      ## suppressWarnings() is used here to suppress the warning from replace(). 
      ## The number of rows to replace does not match the replacement length, which should be the total number of rows in data.out
      data.out <- suppressWarnings(data.out %>% dplyr::mutate(Verified = replace(Verified, !! rlang::parse_expr(ec), 
                                                                          paste(exclusion.info[i, "Message"], Verified, sep = ", "))))
      

      }
  }
  
  ## clean up the format of messages
  data.out[,"Verified"] <- gsub(", Y", "", data.out[,"Verified"])
  
  if (output.res) {
    if (missing(assay) | missing(FILENAME)) stop("Missing either FILENAME and/or assay. Unable to export data to a 'tsv' without a FILENAME and assay.")
    
    if (!is.null(OUTPUT.DIR)) {
      file.path <- OUTPUT.DIR
    } else if (!is.null(INPUT.DIR)) {
      file.path <- INPUT.DIR
    } else {
      file.path <- getwd()
    }
    
    write.table(data.out,
                file=paste0(file.path, "/", FILENAME,"-", assay, "-Level2.tsv"),
                sep="\t",
                row.names=F,
                quote=F)
    cat(paste0("A Level-2 file named ",FILENAME,"-",assay,"-Level2.tsv", 
               " has been exported to the following directory: ", file.path), "\n")
      
  
  }
  
  return(data.out)
}

