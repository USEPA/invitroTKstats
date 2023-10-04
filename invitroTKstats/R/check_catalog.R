#' Function to Check Level 0 Data Catalog
#' 
#' This function is meant to check whether the catalog file is in the anticipated
#' format with required information.
#' 
#' @param catalog The catalog to be checked, format `data.frame`.
#' 
#' @example
#' check_catalog(catalog = data.guide) # note the data.guide is not currently in `invitroTKstats`
#' 
#' @export
check_catalog <- function(catalog){
  ### Catalog Standard Column Names ###
  # set-up the standard catalog column name object
  std.cols <- std.catcols
  # check if the standard catalog column names are in the catalog
  .check_std_colnames_in_data(data = catalog,std.colnames = std.catcols,data.name = "catalog")
  # print passing message
  cat("All of the standard columns exist in the catalog. \n")
  
  ### Check that Required Columns have No Missing Data Entries ###
  # if(any(sapply(catalog,function(x){any(is.na(x))})) | any(sapply(catalog,function(x){any(is.null(x))}))){
  #   someNA <- names(which(sapply(catalog,function(x){any(is.na(x))})))
  #   someNULL <- names(which(sapply(catalog,function(x){any(is.null(x))})))
  #   
  #   if(length(someNA)>0){
  #     c("File",)
  #   }
  #   if(length(someNULL)>0){}
  # }
  ### Check Class of Standard Column Names ###
  # check if there are any columns with only missing (NA or NULL) data and remove from `catalog` just for further checks 
  if(any(sapply(catalog,function(x){all(is.na(x))})) | any(sapply(catalog,function(x){all(is.null(x))}))){
    allNA <- which(sapply(catalog,function(x){all(is.na(x))}))
    allNULL <- which(sapply(catalog,function(x){all(is.null(x))}))
    # remove NA and NULL columns from catalog for checks
    if(length(allNA)>0){
      catalog <- dplyr::select(catalog,-names(allNA)) # catalog[,names(!allNA)]
      std.cols <- std.cols[which(std.cols != names(allNA))]
    }
    if(length(allNULL)>0){
      catalog <- dplyr::select(catalog,-names(allNULL)) # catalog[,names(!allNULL)]
      std.cols <- std.cols[which(std.cols != names(allNULL))]
    }
  }
  
  # check if the standard catalog column names are the correct class
  std.cols.char <- std.cols[-which(std.cols == "Skip.Rows")]
  if("Number.Data.Rows"%in%colnames(catalog)){
    std.cols.num <- c(std.cols[which(std.cols == "Skip.Rows")],
                      num.rows = "Number.Data.Rows")
  }else{
    std.cols.num <- std.cols[which(std.cols == "Skip.Rows")]
  }
  # check 'character' class
  .check_char_cols(data = catalog,char.cols = std.cols.char)
  # check 'numeric' class
  .check_num_cols(data = catalog,num.cols = std.cols.num)
  
  cat("All of the standard columns in the catalog are of the correct class. \n")
  
  ### Final Check ###
  cat("Your data catalog is ready for merge_level0.")
}