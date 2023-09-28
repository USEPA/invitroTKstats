#' Function to create a catalog of level 0 files to be merged.
#' 
#' This function is meant for creating a catalog for all of the level 0 data
#' files that will be merged with the `merge_level0` function.
#' 
#' @param file Vector of character strings with the file names of level 0 data.
#' @param sheet Vector of character strings containing the sheet with MS data. 
#' @param skip.rows Numeric vector containing the number of rows to skip.
#' @param date Character vector containing the date of data collection,
#'             format "MMDDYY". "MM" = 2 digit month, "DD" = 2 digit month,
#'             and "YY" = 2 digit month.
#' @param compound Vector of character strings with the relevant chemical identifier.
#' @param istd
#' @param sample
#' @param type
#' @param peak
#' @param istd.peak
#' @param conc
#' @param analysis.param
#' @param additional.info Named list or data.frame of additional columns to
#'                        include in the catalog.  Additional columns should be
#'                        named with the following structure "<Fill-in>.ColName",
#'                        and all spaces should be designated by a period, "." .
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