#' Function to create a catalog of level 0 files to be merged.
#' 
#' This function is meant for creating a catalog of all level 0 data
#' files listed that will be merged with the `merge_level0` function.
#' All arguments are required, with exception of `additional.info`.
#' 
#' @param file (\emph{character vector}) Vector of character strings with
#'             the file names of level 0 data.
#' @param sheet (\emph{character vector}) Vector of character strings containing
#'              the sheet name with MS data. 
#' @param skip.rows (\emph{numeric vector}) Numeric vector containing the
#'                  number of rows to skip in data file.
#' @param date (\emph{character vector}) Vector of character strings containing
#'             the date of data collection, format "MMDDYY".
#'             "MM" = 2 digit month, "DD" = 2 digit day, and "YY" = 2 digit year.
#' @param compound (\emph{character vector}) Vector of character strings with
#'                 the relevant chemical identifier.
#' @param istd (\emph{character vector}) Vector of character strings with the
#'             internal standard.
#' @param sample (\emph{character vector}) Vector of character strings with
#'               column names containing samples. 
#' @param type (\emph{character vector}) Vector of character strings with column
#'             names containing type information.
#' @param peak (\emph{character vector}) Vector of character strings with the
#'             column names containing mass spectrometry (MS) peak data.
#' @param istd.peak (\emph{character vector}) Vector of character strings with
#'                  column names containing internal standard (ITSD) peak data.
#' @param conc (\emph{character vector}) Vector of character strings with column
#'             names containing exposure concentration data.
#' @param analysis.param (\emph{character vector}) Vector of character strings
#'                       with column names containing analysis parameters.
#' @param num.rows (\emph{numeric vector}) Numeric vector containing the number
#'                 of rows with data to be pulled. (Default is \code{NULL}.)
#' @param additional.info (\emph{list} or \emph{data.frame}) Named list or
#'                        data.frame of additional columns to
#'                        include in the catalog. Additional columns should
#'                        follow the nomenclature of "<Fill-in>.ColName" if
#'                        indicating column names with information to pull,
#'                        otherwise a short name.  All spaces in additional
#'                        column names should be designated with a period, "." .
#'                        (Default is \code{NULL}, i.e. no additional columns.)
#' 
#' @seealso merge_level0
#' 
#' @example 
#' create_catalog(
#'   file = "testME.xlsx",sheet = 3,skip.rows = NA,
#'   date = "092723",compound = "80-05-7",
#'   istd = NA,sample = "ABC",type = "CC",
#'   peak = NA,istd.peak = NA,conc = 34,analysis.param = "clint"
#' )
#' 
#' @export
create_catalog <- function(
    file,sheet,skip.rows,date,compound,istd,sample,
    type,peak,istd.peak,conc,analysis.param,
    num.rows = NULL,
    additional.info = NULL){
  
  data.check <- c(file = missing(file),
                  sheet = missing(sheet),
                  skip.rows = missing(skip.rows),
                  date = missing(date),
                  compound = missing(compound),
                  istd = missing(istd),
                  sample = missing(sample),
                  type = missing(type),
                  peak = missing(peak),
                  istd.peak = missing(istd.peak),
                  conc = missing(conc),
                  analysis.param = missing(analysis.param))
  # check if any of the necessary arguments are not filled in
  if(any(data.check)){
    stop("The following arguments need to be specified:\n\t",
         paste(names(data.check)[which(data.check)],collapse = "\n\t"))
  }
  
  length.check <- c(file = length(file),
                    sheet = length(sheet),
                    skip.rows = length(skip.rows),
                    date = length(date),
                    compound = length(compound),
                    istd = length(istd),
                    sample = length(sample),
                    type = length(type),
                    peak = length(peak),
                    istd.peak = length(istd.peak),
                    conc = length(conc),
                    analysis.param = length(analysis.param))
  u.len.check <- unique(length.check)
  if(length(u.len.check) > 1 & !(1%in%u.len.check)|length(u.len.check) > 2){
    stop("The following columns have mis-matching lengths preventing data.frame creation:\n\t",
         paste(paste(names(length.check)[which(length.check!=1)],
                     length.check[which(length.check!=1)],sep = ": "),
               collapse = "\n\t"))
  }
  # build the base catalog
  catalog <- cbind.data.frame(
    file,sheet,skip.rows,
    date,compound,istd,sample,istd.peak,
    type,peak,conc,analysis.param
  )
  colnames(catalog) <- std.catcols
  
  # check if we need to add a column with the number of rows
  if(!is.null(num.rows)){
    if(length(num.rows)!=nrow(catalog) & length(num.rows)!=1){
      stop("Length of `num.rows` is greater than 1 and does not match the number of rows in the required catalog information.")
    }
    
    catalog <- cbind.data.frame(catalog,Number.Data.Rows = num.rows)
  }
  
  # check if we need to add columns with additional information
  if(!is.null(additional.info)){
    # check the class of the `additional.info` object
    stopifnot(is.data.frame(additional.info)|is.list(additional.info),
              "The 'additional.info' argument needs to be a data.frame or named list.")
    # if `additional.info` is a list make it a data.frame
    if(is.list(additional.info)){
      additional.info <- do.call("cbind.data.frame",additional.info)
    }
    # check the number of rows between `catalog`
    stopifnot(nrow(additional.info) == nrow(catalog),
              "Number of rows for the 'additional.info' does not match standard catalog column rows.")
    
    catalog <- cbind.data.frame(catalog,additional.info)
  }
  
  # Verify the catalog is in the appropriate format
  cat("##################################",
      "## Data Catalog Checks",
      "##################################",
      sep = "\n")
  check_catalog(catalog = catalog)
  cat("\n")
  cat("##################################")
  # output the catalog object
  return(catalog)
}