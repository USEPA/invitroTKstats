#' Function to create a catalog of level 0 files to be merged.
#' 
#' DESCRIPTION NEEDED HERE
#' 
#' @param file Vector of character strings with the fileanames of level 0 data.
#' @param sheet Numeric vector containing the sheet where MS data is 
#' @param skip.rows Numeric vector containing the number of rows to skip
#' @param date 
#' @param compound 
#' @param istd
#' @param sample
#' @param type
#' @param peak
#' @param istd.peak
#' @param conc
#' @param analysis.param
#' @param additional.info Named list of additional columns to add to the
#' 
#' @example 
#' create_catalog(file = "testME.xlsx",sheet = 3,skip.rows = NA,date = "09-27-2023",compound = "80-05-7",istd = NA,sample = "ABC",type = "CC",peak = NA,istd.peak = NA,conc = 34,analysis.param = "clint")
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
    date,compound,istd,sample,
    type,peak,conc,analysis.param
  )
  # check if we need to add a column with the number of rows
  if(!is.null(additional.info)){
    stopifnot(is.data.frame(additional.info)|is.list(additional.info),
              "The 'additional.info' argument needs to be a data.frame or named list.")
    
    if(is.list(additional.info)){
      additional.info <- do.call("cbind.data.frame",additional.info)
    }
    catalog <- cbind.data.frame(catalog,additional.info)
  }
  # output the catalog object
  return(catalog)
}