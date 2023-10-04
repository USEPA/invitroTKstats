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

#' Check no missing data for specified required columns.
#' @description This function checks for whether any of the required columns have a data entry of `NA` or `NULL`.
#' @param data Data frame to check.
#' @param req.cols Columns with required data.
.check_req_cols <- function(data,req.cols){
  # check 'required' columns for missing data
  check <- sapply(data[,req.cols],function(x){any(is.na(x)|is.null(x))})
  if(any(check)){
    incomp.req.cols <- which(check)
    print(incomp.req.cols)
    incomp.req.data <- sapply(data[,incomp.req.cols],function(x)which(is.na(x)|is.null(x))) %>% 
      do.call("cbind",.)
    print(incomp.req.data)
    # for(i in 1:length()incomp.req.data)
    # data[incomp.req.data,incomp.req.cols]
  }
}