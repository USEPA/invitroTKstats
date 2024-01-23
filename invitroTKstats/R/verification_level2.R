#' Creates a Data Frame with a Verification Column (Level-2)
#'
#' This function takes in a Level-1 data frame and an exclusion list and 
#' returns a Level-2 data frame with a verification column. The
#' verification column contains either "Y", indicating the row is good for analysis,
#' or messages contained in the exclusion list for why the data rows are excluded. 
#' If an exclusion list is not provided, all rows are assumed to be verified 
#' and are good to use in the analysis. 
#' 
#' The `exclusion.list` should be a data frame with the following columns:
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
#' @param exclusion.list (Data Frame) A data frame containing the variables and 
#' values of the corresponding variables to exclude rows. 
#' See details for full explanation.
#' 
#' @param assay (Character) A string indicating what assay data the input file is. Valid 
#' input is one of the following: "Clint", "fup-UC", "fup-RED", or "Caco2".
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
#' level1 <- dplyr::select(level2, -c(Verified))
#' 
#' exclusion_list <- data.frame(
#' Variable = c("DTXSID", "Lab.Sample.Name", "Compound.Name"),
#' Value = c("DTXSID00380798","G4-Ametryn Hep052521 T0a", "Nonafluoropentanamide"), 
#' Message = c("bad id #1","bad sample", "bad compound")
#' )
#' 
#' my.level2 <- verification_function(FILENAME="kreutz", 
#' data.in=level1, exclusion_list, assay="Clint")
#' 
verification_function <- function(
    FILENAME, 
    data.in, 
    exclusion.list,
    assay,
    output.res = TRUE,
    INPUT.DIR = NULL,
    OUTPUT.DIR = NULL
    ){
  
  if (!missing(data.in)) {
    data.out <- as.data.frame(data.in)
    rm(data.in)
    } else if (!is.null(INPUT.DIR)) {
    data.out <- read.csv(file=paste0(INPUT.DIR, "/", FILENAME,"-", assay, "-Level1.tsv"),
                     sep="\t",header=T)
    } else {
    data.out <- read.csv(file=paste0(FILENAME,"-", assay, "-Level1.tsv"),
                          sep="\t",header=T)  
    }
  
  # add a column with all "Y"
  data.out$Verified <- "Y"
  
  if (!missing(exclusion.list)) {
    if (!all(unique(exclusion.list[, "Variable"]) %in% colnames(data.out))) 
      stop("Name(s) of the list(s) do not match the column names of the level-1 data.")
       
    for (i in 1:nrow(exclusion.list)) {
      this.variable <- exclusion.list[i ,"Variable"]
      this.value <- exclusion.list[i, "Value"]
      which.rows <- data.out[, this.variable] == this.value
      which.rows[is.na(which.rows)] <- FALSE
      data.out[which.rows,"Verified"] <- exclusion.list[i, "Message"]
    }
  }
  
  if (output.res) {
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

