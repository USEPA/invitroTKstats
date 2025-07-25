#' Add Sample Verification Column (Level-2)
#'
#' This function takes in a level-1 data frame and an exclusion list and 
#' returns a level-2 data frame with a verification column. The
#' verification column contains either "Y", indicating the row is good for analysis,
#' or messages contained in the exclusion list for why the data rows are excluded. 
#' If an exclusion list is not provided, all rows are assumed to be good for use 
#' in further analyses and are verified with "Y".
#' 
#' The `exclusion.info` should be a data frame with the following columns:
#' \tabular{rr}{
#'   Variables \tab level-1 variable(s) used to filter rows for exclusion\cr
#'   Values \tab Value(s) to exclude\cr
#'   Message \tab Simple explanation for the exclusion\cr
#' }
#' When filtering on multiple variable-value pairs, the character input for 
#' "Variables" and "Values" should be separated by a vertical bar "|" ,
#' and the variable-value pairs should match. See demonstration in Examples, Scenario 1. 
#' 
#' NOTE: Currently if NA's exist in a variable of interest for 'verification' assignments,
#' then that variable cannot be used for assigning verification. Thus, either alternative
#' variable-value pairs will need to be used in lieu of variable with missing values, or
#' (though less ideal) "manual coding" adjustments in the verification column may be necessary.
#' 
#' If the output level-2 data frame is chosen to be exported and an output directory 
#' is not specified, it will be exported to the user's R session temporary directory. 
#' This temporary directory is a per-session directory whose path can be found
#' with the following code: \code{tempdir()}. For more details, see 
#' \url{https://www.collinberke.com/til/posts/2023-10-24-temp-directories/}.
#' 
#' As a best practice, \code{INPUT.DIR} (when importing a .tsv file) and/or 
#' \code{OUTPUT.DIR} should be specified to simplify the process of importing 
#' and exporting files. This practice ensures that the exported files can easily 
#' be found and will not be exported to a temporary directory. 
#'
#' @param FILENAME (Character) A string used to identify the output level-1 file.
#' "<FILENAME>-<assay>-Level1.tsv". 
#' 
#' @param data.in (Data Frame) A level-1 data frame from the format functions.
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
#' @param output.res (Logical) When set to \code{TRUE}, the resulting
#' data frame (level-2) will be exported to the user's per-session temporary directory
#' or \code{OUTPUT.DIR} (if specified) as a .tsv file. 
#' (Defaults to \code{FALSE}.)
#'
#' @param INPUT.DIR (Character) Path to the directory where the input level-1 file exists. 
#' If \code{NULL}, looking for the input level-1 file in the current working
#' directory. (Defaults to \code{NULL}.)
#' 
#' @param OUTPUT.DIR (Character) Path to the directory to save the output file. 
#' If \code{NULL}, the output file will be saved to the user's per-session temporary
#' directory or \code{INPUT.DIR} if specified. (Defaults to \code{NULL}.)
#' 
#' @return A level-2 data frame with a verification column. 
#' 
#' @importFrom utils read.csv write.table
#' 
#' @export
#' 
#' @examples 
#' level1 <- invitroTKstats::clint_L1
#' 
#' # Scenario 1: Pass in data.in and exclusion.info data frame from R session 
#' 
#' # Create a exclusion criteria data frame
#' # Use the excluded samples found in \code{invitroTKstats::clint_L2_heldout}
#' # If more than one variable is used to define a set of samples to be excluded,
#' # enter them as one string, separate the Variables with a vertical bar, "|",
#' # and do the same for Values. 
#' 
#' excluded_level2 <- invitroTKstats::clint_L2_heldout
#' 
#' exclusion_criteria <- data.frame(
#'   Variables = paste("Compound.Name", "Lab.Sample.Name", sep = "|"), 
#'   Values = paste(excluded_level2[,"Compound.Name"], excluded_level2[,"Lab.Sample.Name"], sep = "|"),
#'   Message = excluded_level2[,"Verified"]
#'   )
#'   
#' # Run the verification function.
#' my.level2 <- sample_verification(data.in=level1,
#'                                  exclusion.info = exclusion_criteria,
#'                                  output.res = FALSE)
#' 
#' # Scenario 2: Import 'tsv' as input data and do not pass in an exclusion.info data frame
#' 
#' \dontrun{
#' # Write the level-1 file to some folder
#' # Will need to replace <desired level-1 FOLDER> with desired export folder location.
#' # The <desired level-1 FOLDER> needs to already exist.   
#' 
#' write.table(level1,
#' file=here::here("<desired level-1 FOLDER>/Smeltz-Clint-Level1.tsv"),
#' sep="\t",
#' row.names=F,
#' quote=F)
#' 
#' # Run the verification function.
#' # Specify the path to import level-1 data with INPUT.DIR.
#' # Will need to replace INPUT.DIR = <desired level-1 FOLDER> with chosen output
#' # folder location from above 
#' # If no exclusion.info data frame is used, will label all samples as verified.
#' # A level-2 file is also exported to INPUT.DIR when OUTPUT.DIR is not specified.
#' my.level2 <- sample_verification(FILENAME="Smeltz", 
#' assay="Clint", INPUT.DIR = here::here("<desired level-1 FOLDER>"))
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
    output.res = FALSE,
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
    
    if (!is.null(OUTPUT.DIR)) {  # export output file to OUTPUT.DIR (OUTPUT.DIR specified) 
      file.path <- OUTPUT.DIR
    } else if (!is.null(INPUT.DIR)) { # export output file to INPUT.DIR (OUTPUT.DIR not specified)
      file.path <- INPUT.DIR
    } else { # export output file to tempdir() (OUTPUT.DIR & INPUT.DIR not specified)
      file.path <- tempdir()
    }
    
    write.table(data.out,
                file=paste0(file.path, "/", FILENAME,"-", assay, "-Level2.tsv"),
                sep="\t",
                row.names=F,
                quote=F)
    cat(paste0("A level-2 file named ",FILENAME,"-",assay,"-Level2.tsv", 
               " has been exported to the following directory: ", file.path), "\n")
      
  
  }
  
  return(data.out)
}

