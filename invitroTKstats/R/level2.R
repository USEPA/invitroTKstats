verification_fup_red <- function(
    FILENAME, 
    data.in, 
    exclusion.list,
    output.res = TRUE,
    INPUT.DIR = NULL,
    OUTPUT.DIR = NULL
    ){
  
  if (missing(data.in)) {
    if (!is.null(INPUT.DIR)) {
      data.in <- read.csv(file=paste0(INPUT.DIR, "/", FILENAME,"-fup-RED-Level1.tsv"),
                          sep="\t",header=T)
    } else {
      data.in <- read.csv(file=paste0(FILENAME,"-fup-RED-Level1.tsv"),
                        sep="\t",header=T)
    }
  } 
  
  # add a column with all "Y"
  data.in$Verified <- "Y"
  
  # option 1
  # Check if the naming of each list match the column name
  if (!(all(names(exclusion.list) %in% colnames(data.in)))) 
    cat("error.\n")
  
  for (i in length(exclusion.list)) {
    column <- names(exclusion.list)[i]
    criteria <- exclusion.list[[i]][["criteria"]]
    
    if (length(criteria) == 1) {
      which.rows <- data.in[, column] == exclusion.list[[i]][[column]][k]
      which.rows[is.na(which.rows)] <- FALSE
      data.in[which.rows,"Verified"] <- criteria[k]
    } else if (length(criteria) == exclusion.list[[i]][[column]]) {
      for (k in length(exclusion.list[[i]][[column]]))
      {
        which.rows <- data.in[, column] == exclusion.list[[i]][[column]][k]
        which.rows[is.na(which.rows)] <- FALSE
        data.in[which.rows,"Verified"] <- criteria[k]
      }
      
    }
    
    
    # which.rows <- data.in[, column] %in% exclusion.list[[i]][[column]]
    # which.rows[is.na(which.rows)] <- FALSE
    # data.in[which.rows,"Verified"] <- exclusion.list[[i]][["criteria"]]
    
    # which.rows <- grepl(pattern, data.in[, column])
    # data.in[which.rows,"Verified"] <- exclusion.list[[i]][["criteria"]]
  }
  
  # option 2
  # if (!missing(exclude_DTXSID)) {
  #   for (i in 1:nrow(exclude_DTXSID)) {
  #     this.id <- exclude_DTXSID[i, "DTXSID"]
  #     which.rows <- data.in[, "DTXSID"] == this.id 
  #     which.rows[is.na(which.rows)] <- FALSE
  #     data.in[which.rows,"Verified"] <- exclude_DTXSID[i, "criteria"]
  #   }
  # }
  
  # repeat for lab sample name and lab compound name 
  
  if (output.res) {
    if (!is.null(OUTPUT.DIR)) {
      file.path <- OUTPUT.DIR
    } else if (!is.null(INPUT.DIR)) {
      file.path <- INPUT.DIR
    } else {
      file.path <- getwd()
    }
    
    write.table(data.in,
                file=paste(file.pth, "/", FILENAME,"-fup_RED-Level2.tsv"),
                sep="\t",
                row.names=F,
                quote=F)
      
  
  }
  
  return(data.in)
}

