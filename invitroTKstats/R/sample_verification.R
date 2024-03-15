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
#'   Variable \tab Level-1 variable used to filter rows for exclusion\cr
#'   Value \tab Value to exclude\cr
#'   Message \tab Simple explanation for the exclusion\cr
#' }
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
#' @export
#' 
#' @examples 
#' level2 <- invitroTKstats::kreutz2023.clint
#' 
#' # Data sets in the package are level-2 data. For demonstration purposes,
#' # remove the verification columns from the data and use it as level-1.
#' level1 <- dplyr::select(level2, -c(Verifed, Verified))
#' 
#' # scenario 1: input data is in the R session
#' # create a exclusion criteria data frame
#' exclusion_criteria <- data.frame(
#' Variable = c("DTXSID", "Lab.Sample.Name", "Compound.Name"),
#' Value = c("DTXSID00380798","G4-Ametryn Hep052521 T0a", "Nonafluoropentanamide"), 
#' Message = c("bad id #1","bad sample", "bad compound")
#' )
#' 
#' my.level2 <- sample_verification(FILENAME="kreutz", 
#' data.in=level1, exclusion_criteria, assay="Clint")
#' 
#' # scenario 2: import 'tsv' as input data
#' 
#' # Write the level-1 file to some folder
#' write.table(data.out,
#' file="~/invitrotkstats/invitroTKstats/data-raw/kreutz-Clint-Level1.tsv",
#' sep="\t",
#' row.names=F,
#' quote=F)
#' 
#' # import level-1 with INPUT.DIR 
#' my.level2 <- sample_verification(FILENAME="kreutz", 
#' exclusion_criteria, assay="Clint", INPUT.DIR = "~/invitrotkstats/invitroTKstats/data-raw")
#' 
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
  if ((missing(data.in) |  output.res) & !(assay %in% approved_assays)) 
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
    # approach one
    
    # if (!all(unique(exclusion.info[, "primary_var"]) %in% colnames(data.out)) |
    #     !all(unique(exclusion.info[, "secondary_var"]) %in% colnames(data.out))) 
    #   stop("Name(s) of the list(s) do not match the column names of the level-1 data.")
    #    
    # for (i in 1:nrow(exclusion.info)) {
    #   this.prim.variable <- exclusion.info[i ,"primary_var"]
    #   this.prim.value <- exclusion.info[i, "primary_value"]
    #   this.sec.variable <- exclusion.info[i, "secondary_var"]
    #   this.sec.value <- exclusion.info[i, "secondary_value"]
    #   which.rows <- data.out[, this.prim.variable] == this.prim.value
    #   if (this.sec.variable != "")
    #     which.rows <- (which.rows & data.out[, this.sec.variable] == this.sec.value)
    #   
    #   which.rows[is.na(which.rows)] <- FALSE
    #   # append this message in case the sample is excluded for multiple reasons
    #   data.out[which.rows,"Verified"] <- paste(exclusion.info[i, "message"], data.out[which.rows,"Verified"], sep = ", ")
    # }
    
      ## approach two
      for (i in 1:nrow(exclusion.info)) {
        ## split the list, consider all possible separators
        ## exclude period - period is used in column names 
        var.list <- strsplit(exclusion.info[i, "primary_var"], "\\\\|[^.[:^punct:]]", perl = TRUE)[[1]]
        var.list <- trimws(var.list)
        value.list <- strsplit(exclusion.info[i, "primary_value"], "\\\\|[^.[:^punct:]]", perl = TRUE)[[1]]
        value.list <- trimws(value.list)
        ## check if every variable has a matched value 
        if (length(var.list) != length(value.list))
          stop("The lengths of variable list and value list do not match.")
        ## check if the variable names are valid
        if (!all(var.list %in% colnames(data.out))) 
          stop("Names of the variables use to determine the exclusion criteria do not match the column names of the level-1 data.")
        for (j in 1:length(var.list)) {
          this.variable <- var.list[j]
          this.value <- value.list[j]
          if (j == 1)
            which.rows <- data.out[, this.variable] == this.value
          else
            which.rows <- (which.rows & data.out[, this.variable] == this.value)
        }
        
        which.rows[is.na(which.rows)] <- FALSE
        # append this message in case the sample has multiple reasons to be excluded
        data.out[which.rows,"Verified"] <- paste(exclusion.info[i, "message"], data.out[which.rows,"Verified"], sep = ", ")
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

