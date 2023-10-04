#' Function to create a catalog of level 0 files to be merged.
#' 
#' This function is meant for creating a catalog of all level 0 data
#' files listed that will be merged with the `merge_level0` function.
#' All arguments are required, with exception of `additional.info`.
#' 
#' @param file Vector of character strings with the file names of level 0 data.
#' @param sheet Vector of character strings containing the sheet with MS data. 
#' @param skip.rows Numeric vector containing the number of rows to skip in
#'                  data file.
#' @param date Vector of character strings containing the date of data
#'             collection, format "MMDDYY". "MM" = 2 digit month,
#'             "DD" = 2 digit month, and "YY" = 2 digit month.
#' @param compound Vector of character strings with the relevant chemical
#'                 identifier.
#' @param istd Vector of character strings with the internal standard.
#' @param sample Vector of character strings with column names containing
#'               samples. 
#' @param type Vector of character strings with column names containing type
#'             information.
#' @param peak Vector of character strings with the column names containing
#'             mass spectrometry (MS) peak data.
#' @param istd.peak Vector of character strings with column names containing
#'                  internal standard (ITSD) peak data.
#' @param conc Vector of character strings with column names containing
#'             exposure concentration data.
#' @param analysis.param Vector of character strings with column names
#'                       containing analysis parameters.
#' @param num.rows Numeric vector containing the number of rows with data to be
#'                 pulled. (Default is NULL.)
#' @param additional.info Named list or data.frame of additional columns to
#'                        include in the catalog. Additional columns should
#'                        follow the nomenclature of "<Fill-in>.ColName" if
#'                        indicating column names with information to pull,
#'                        otherwise a short name.  All spaces in additional
#'                        column names should be designated with a period, "." .
#'                        (Default is NULL, i.e. no additional columns.)
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
  # output the catalog object
  return(catalog)
}