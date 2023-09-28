#' Function to Check Level 0 Data Catalog
#' 
#' This function is meant to check whether the catalog file is in the anticipated
#' format with required information.
#' 
#' @param cat The catalog to be checked, format `data.frame`.
#' 
#' @example
#' check_catalog(cat = data.guide) # note the data.guide is not currently in `invitroTKstats`
#' 
#' @export
check_catalog <- function(cat){
  ### Catalog Standard Column Names ###
  # set-up the standard catalog column name object
  std.cols <- std.catcols
  # check if the standard catalog column names are in the catalog
  if(!all(std.cols%in%colnames(cat))){
    # if not then find which are not in then identify those to be flagged
    flag <- std.cols[which(!(std.cols%in%colnames(cat)))]
    # STOP and print what is missing or mis-named
    stop(paste(flag,collapse=", ")," - missing or mis-named columns in the",
         "data catalog.")
  }
  # cat("All of the standard columns in the catalog.") # only un-comment if more checks added
  
  ### Final Check ###
  cat("Your data catalog is ready for merge_level0.")
}