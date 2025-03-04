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
  # check if the standard catalog column names are in the catalog
  .check_std_colnames_in_data(data = catalog,std.colnames = std.catcols,data.name = "catalog")
  # print passing message
  cat("All of the standard columns exist in the catalog. \n") # <may need to check with John which standard columns can have some missing data but can't all be missing out of standard columns and or others>
  
  ### Check if there Columns with only Missing Data (and are Problematic) ###
  .check_all_miss_cols(data = catalog,req.cols = std.catcols) # <may need to check with John which standard columns require all data to be filled>
  
  ### Check that Required Columns have No Missing Data Entries ###
  .check_no_miss_cols(data = catalog,req.cols = std.catcols,return.missing = TRUE)
  cat("All standard columns are data complete.\n")
  
  ### Check Class of Standard Column Names ###
  # check if the standard catalog column names are the correct class
  std.cols.char <- std.catcols[-which(std.catcols %in% c("Skip.Rows","Col.Names.Loc"))]
  if("Number.Data.Rows"%in%colnames(catalog)){
    std.cols.num <- c(std.catcols[which(std.catcols %in% c("Skip.Rows","Col.Names.Loc"))],
                      num.rows = "Number.Data.Rows")
  }else{
    std.cols.num <- std.catcols[which(std.catcols == "Skip.Rows")]
  }
  # check 'character' class
  .check_char_cols(data = catalog,char.cols = std.cols.char)
  # check 'numeric' class
  .check_num_cols(data = catalog,num.cols = std.cols.num)
  
  cat("All of the standard columns in the catalog are of the correct class.\n")
  
  ### Final Check ###
  cat("Your data catalog is ready for merge_level0.")
}