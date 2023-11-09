#' Check the standard column names are in the data.
#' @param data Data frame to check.
#' @param std.colnames Vector of character strings with standard column names
#'                     to check for in the data.
.check_std_colnames_in_data <- function(data,std.colnames,data.name = NULL){
  # check if the standard catalog column names are in the catalog
  if(!all(std.colnames%in%colnames(data))){
    # if not then find which are not in then identify those to be flagged
    flag <- std.colnames[which(!(std.colnames%in%colnames(data)))]
    # STOP and print what is missing or mis-named
    stop(paste(flag,collapse=", ")," - missing or mis-named columns in the ",
         ifelse(is.null(data.name),"data.",paste0("data ",data.name,".")))
  }
}

#' Check the character columns are correctly of character class.
#' @param data Data frame to check.
#' @param char.cols Column names that should be of the character class.
.check_char_cols <- function(data,char.cols){
  # check 'character' class
  check.char <- sapply(data[,char.cols],is.character)
  if(!all(check.char)){
    stop("The following columns are not of class `character`: \n\t",
         paste(char.cols[-which(check.char)],collapse = "\n\t"))
  }
}

#' Check the numeric columns are correctly of numeric class.
#' @param data Data frame to check.
#' @param num.cols Column names that should be of the numeric class.
.check_num_cols <- function(data,num.cols){
  # check 'numeric' class
  check.num <- sapply(data[,num.cols],is.numeric)
  if(!all(check.num)){
    stop("The following columns are not of class `numeric`: \n\t",
         paste(num.cols[-which(check.num)],collapse = "\n\t"))
  }
}

#' Check if all the data is missing for specified columns.
#' @description This function checks for whether any of the specified columns
#'   are missing all of their data, either `NA` and/or `NULL`.
#' @param data Data frame to check.
#' @param req.cols Column names that should be checked for whether all data is missing.
.check_all_miss_cols <- function(data,req.cols){
  check <- sapply(data[,req.cols],function(x){all(is.na(x)|is.null(x))})
  if(any(check)){
    stop("The following column has no data entered, all values missing: \n\t",
         paste(req.cols[which(check)],collapse = "\n\t"))
  }
}

#' Check there is no missing data for specified columns.
#' @description This function checks for whether any of the required columns
#'   have a data entry of `NA` or `NULL`.
#' @param data Data frame to check.
#' @param req.cols Columns with required data.
#' @param return.missing Logical argument, if `TRUE` return rows missing data in
#'   column (list or vector by column name). (Default is `FALSE`.)
.check_no_miss_cols <- function(data,req.cols,return.missing = FALSE){
  # check 'required' columns for missing data
  check <- sapply(data[,req.cols],function(x){any(is.na(x)|is.null(x))})
  if(any(check)){
    # which columns have at least one cell (row) of missing data
    incomp.req.cols <- which(check)
    # which rows in the id'ed columns have missing data
    incomp.req.data <- sapply(data[,incomp.req.cols],
                              function(x){which(is.na(x)|is.null(x))})
    # if the user want know the cells missing data, then print.
    if(return.missing){print(incomp.req.data)}
    # throw an error and stop all actions
    stop("The following cells have missing data where values are required: \n\t",
         paste(missing_cells,collapse = ", "))
  }
}