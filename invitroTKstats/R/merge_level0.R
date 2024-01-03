#' Merge multiple level0 files into a single table for processing
#'
#' This function reads multiple Excel files containing mass-spectrometry (MS) data
#' and extracts the chemical sample data from the specified
#' sheets. The argument `level0.catalog` is a table that
#' provides the necessary information to find the data for each chemical. The
#' primary data of interest are the analyte peak area, the internal standard
#' peak area, and the target concentration for calibration curve (CC) samples.
#' The argument `data.label` is used to annotate this particular mapping of level0
#' files into data ready to be organized into a level1 file.
#'
#' Unless specified to be a single value for all the files, for example sheet="Data",
#' the argument `level0.catalog` should be a data frame with the following columns:
#' \tabular{rr}{
#'   File \tab The Excel filename to be loaded\cr
#'   Sheet \tab The name of the Sheet to examine within in the Excel file\cr
#'   Skip.Rows \tab How many rows should be skipped on the sheet to get usable column names\cr
#'   Date \tab The date the measurements were made\cr
#'   Chemical.ID \tab The laboratory chemical identity\cr
#'   ISTD \tab The internal standard used\cr
#'   Sample.ColName \tab The column name on the sheet that contains sample identity\cr
#'   Type.ColName \tab The column name on the sheet that contains the type of sample\cr
#'   Peak.ColName \tab The column name on the sheet that contains the analyte MS peak area \cr
#'   ISTD.Peak.ColName \tab The column name on the sheet that contains the internal standard MS peak area\cr
#'   Conc.ColName \tab The column name on the sheet that contains the intended concentration for calibration curves \cr
#'   AnalysisParam.ColName \tab The column name on the sheet that contains the MS instrument parameters for the analyte\cr
#' }
#' Columns with names ending in ".ColName" indicate the columns to be extracted
#' from the specified Excel file and sheet containing level 0 data.
#'
#' @param FILENAME (Character) A string used to identify outputs of the function call.
#' (Default to "MYDATA")
#' 
#' @param level0.catalog A data frame describing which columns of which sheets
#' in which Excel files contain MS data for analysis. See details for full
#' explanation.
#' 
#' @param file.col (Character) Column name containing level 0 file names
#' to pull data from.
#' 
#' @param sheet (Character) Excel file sheet name/identifier containing
#' level 0 where data is to be pulled from. (Defaults to `NULL`.) (Note: Single
#' entry only, use only if all files have the same sheet identifier for
#' level 0 data.) 
#' 
#' @param sheet.col (Character) Catalog column name containing `sheet`
#' information. (Default to "Sheet")
#' 
#' @param skip.rows (Numeric) Number of rows to skip when extracting level 0
#' data from the specified Excel file(s). (Defaults to `NULL`.) (Note: Single
#' entry only, use only if all files need to skip the same number of rows
#' for extracting level 0 data.)
#' 
#' @param skip.rows.col (Character) Catalog column name containing `skip.rows`
#' information. (Default to "Skip.Rows")
#' 
#' @param num.rows (Numeric) Number of rows to pull when extracting level 0
#' data from the specified Excel file(s). (Defaults to `NULL`.) (Note: Single
#' entry only, use only if all files need to pull the same number of rows for
#' extracting level 0 data.)
#' 
#' @param num.rows.col (Character) Catalog column name containing `num.rows`
#' information. (Default to `NULL`)
#' 
#' @param date (Character) Date of laboratory measurements. Typical format
#' "MMDDYY" ("MM" = 2 digit month, "DD" = 2 digit day, and "YY" = 2 digit year).
#' (Defaults to `NULL`.) (Note: Single entry only, use only if all files have
#' the same laboratory measurement date.)
#' 
#' @param date.col (Character) Catalog column name containing `date`
#' information. (Defaults to "Date")
#' 
#' @param compound.col (Character) Column name of input.data (level0?) indicates the
#' test compound. (Defaults to "Chemical.ID")
#' 
#' @param istd.col (Character) Column name of input.data (level0?) indicating the
#' MS peak area for the internal standard. (Defaults to "ISTD")
#' 
#' @param sample.colname (Character) Column name of level 0 data containing
#' sample information. (Defaults to `NULL`.) (Note: Single entry only, use only
#' if all files use the same column name for sample names when extracting
#' level 0 data.)
#' 
#' @param sample.colname.col (Character) Catalog column name containing 
#' `sample.colname` information. (Defaults to "Sample.ColName") 
#' 
#' @param type.colname (Character) Column name of the level 0 data containing
#'  the type of sample. (Defaults to `NULL`.) (Note: Single entry only, use
#'  only if all files use the same column name for sample type information
#'  when extracting level 0 data.)
#' 
#' @param type.colname.col (Character) Catalog column name containing
#' `type.colname` information. (Defaults to "Type".)
#' 
#' @param peak.colname (Character) Column name of the level 0 data containing
#'  the analyte Mass Spectrometry peak area. (Defaults to `NULL`.)
#'  (Note: Single entry only, use only if all files use the same column name
#'  for analyte peak area information when extracting level 0 data.)
#' 
#' @param peak.colname.col (Character) Catalog column name containing
#' `peak.colname` information. (Defaults to "Peak.ColName")
#' 
#' @param istd.peak.colname (Character) Column name of the level 0 data
#'  containing the internal standard Mass Spectrometry peak area. (Note: Single
#'  entry only, use only if all files use the same column name for internal
#'  standard MS peak area information when extracting level 0 data.)
#' 
#' @param istd.peak.colname.col (Character) Catalog column name containing
#' `istd.peak.colname` information. (Defaults to "ISTD.Peak.ColName")
#' 
#' @param conc.colname (Character) Column name of the level 0 data containing
#'  intended concentrations for calibration curves. (Defaults to `NULL`.)
#'  (Note: Single entry only, use only if all files use the same column name
#'  for intended concentration information when extracting level 0 data.)
#'  
#' @param conc.colname.col (Character) Catalog column name containing 
#' `conc.colname` information. (Defaults to "Conc.ColName")
#' 
#' @param analysis.param.colname (Character) Column name of the level 0 data
#'  containing Mass Spectrometry instrument parameters for the analyte.
#'  (Defaults to `NULL`.) (Note: Single entry only, use only if all files use
#'  the same column name for analysis parameter information when extracting
#'  level 0 data.)
#' 
#' @param analysis.param.colname.col (Character) Catalog column name containing
#' `analysis.param.colname` information. (Defaults to "AnalysisParam.ColName")
#' 
#' @param additional.colnames Additional columns from the level 0 data files to
#'  pull information from when extracting level 0 data and include in the
#'  compiled level 0 returned from `merge_level0`. (Defaults to `NULL`.)
#' 
#' @param additional.colname.cols Catalog column name(s) containing 
#'  `additional.colnames` information, (Defaults to `NULL`.)
#' 
#' @param chem.ids (Data frame) A data frame containing basic chemical
#'  identification information for tested chemicals.
#' 
#' @param chem.lab.id.col (Character) Column in `chem.ids` containing
#'  the compound/chemical identifier used by the laboratory in level 0 measured
#'  data. (Defaults to "Chem.Lab.ID")
#' 
#' @param chem.name.col (Character) `chem.ids` column name containing the
#'  "standard" chemical name to use for annotation of the compiled level 0
#'  returned from `merge_level0`. (Defaults to "Compound")
#' 
#' @param chem.dtxsid.col (Character) `chem.ids` column name containing EPA's
#'  DSSTox Structure ID (\url{http://comptox.epa.gov/dashboard})
#'  (Defaults to "DTXSID") 
#'  
#' @param catalog.out (Logical) When set to \code{TRUE}, the data frame 
#' specified in \code{level0.catalog} will be exported the current directory 
#' or \code{OUTPUT.DIR} as a .tsv file. (Defaults to \code{TRUE}.)
#' 
#' @param output.res (Logical) When set to \code{TRUE}, the result 
#' table (Level-0) will be exported the current directory or \code{OUTPUT.DIR} 
#' as a .tsv file. (Defaults to \code{TRUE}.)
#' 
#' @param OUTPUT.DIR (Character) Path to the directory to save the output file. 
#' If \code{NULL}, the output file will be saved to the current working
#' directory. (Defaults to \code{NULL}.)
#' 
#' @return \item{data.frame}{A data.frame in standardized Level-0 format} 
#'
#' @author John Wambaugh
#' 
#' @import readxl
#' 
#' @export merge_level0
merge_level0 <- function(FILENAME="MYDATA",
  level0.catalog,
  file.col="File",
  sheet=NULL,
  sheet.col="Sheet",
  skip.rows=NULL,
  skip.rows.col="Skip.Rows",
  num.rows=NULL,
  num.rows.col=NULL,
  date=NULL,
  date.col="Date",
  compound.col="Chemical.ID",
  istd.col="ISTD",
  sample.colname=NULL,
  sample.colname.col="Sample.ColName",
  type.colname=NULL,
  type.colname.col="Type",
  peak.colname=NULL,
  peak.colname.col="Peak.ColName",
  istd.peak.colname=NULL,
  istd.peak.colname.col="ISTD.Peak.ColName",
  conc.colname=NULL,
  conc.colname.col="Conc.ColName",
  analysis.param.colname=NULL,
  analysis.param.colname.col="AnalysisParam.ColName",
  additional.colnames=NULL,
  additional.colname.cols=NULL,
  chem.ids,
  chem.lab.id.col="Chem.Lab.ID",
  chem.name.col="Compound",
  chem.dtxsid.col="DTXSID",
  catalog.out = TRUE,
  output.res = TRUE,
  OUTPUT.DIR = NULL
  )
{
  level0.catalog <- as.data.frame(level0.catalog)
  
# These arguments allow the user to specify a single value for every observation 
# in the table:  
  if (!is.null(sheet)) level0.catalog[,sheet.col] <- sheet
  if (!is.null(skip.rows)) level0.catalog[,skip.rows.col] <- skip.rows
  if (!is.null(num.rows)) level0.catalog[,num.rows.col] <- num.rows
  if (!is.null(date)) level0.catalog[,date.col] <- date
  if (!is.null(sample.colname)) level0.catalog[,sample.colname.col] <- 
    sample.colname
  if (!is.null(type.colname)) level0.catalog[,type.colname.col] <- 
    type.colname
  if (!is.null(peak.colname)) level0.catalog[,peak.colname.col] <- 
    std.conc
  if (!is.null(istd.peak.colname)) level0.catalog[,istd.peak.colname.col] <- 
    istd.peak.colname
  if (!is.null(conc.colname)) level0.catalog[,conc.colname.col] <- 
    conc.colname
  if (!is.null(analysis.param.colname)) level0.catalog[,analysis.param.colname.col] <- 
    analysis.param.colname
  
# We need all these columns in clint.data
  cols <-c(
    file.col,
    sheet.col,
    skip.rows.col,
    date.col,
    compound.col,
    istd.col,
    sample.colname.col,
    type.colname.col,
    peak.colname.col,
    istd.peak.colname.col,
    conc.colname.col,
    analysis.param.colname.col
    )
  
  if (!is.null(num.rows.col)) cols <- c(cols, num.rows.col)
  
  if (!is.null(additional.colname.cols)) cols <- c(cols,
      additional.colname.cols)

  if (!(all(cols %in% colnames(level0.catalog))))
  {
    stop(paste("Missing columns named:",
      paste(cols[!(cols%in%colnames(level0.catalog))],collapse=", ")))
  }

  # Organize the columns:
  level0.catalog <- level0.catalog[,cols]
    
# Standardize the column names:
  file.col <- "File"
  sheet.col <- "Sheet"
  skip.rows.col <- "Skip.Rows"
  date.col <- "Date"
  compound.col <- "Chemical.ID"
  istd.col <- "ISTD.Name"
  sample.colname.col <- "Sample.ColName"
  type.colname.col <- "Type.ColName"
  peak.colname.col <- "Peak.ColName"
  istd.peak.colname.col <- "ISTD.Peak.ColName"
  conc.colname.col <- "Conc.ColName"
  analysis.param.colname.col <- "AnalysisParam.ColName"
  
  std.colnames <- c(
    file.col,
    sheet.col,
    skip.rows.col,
    date.col,
    compound.col,
    istd.col,
    sample.colname.col,
    type.colname.col,
    peak.colname.col,
    istd.peak.colname.col,
    conc.colname.col,
    analysis.param.colname.col
    )

  if (!is.null(num.rows.col)) std.colnames <- c(std.colnames, num.rows.col)
  
  if (!is.null(additional.colname.cols)) std.colnames <- c(std.colnames,
    additional.colname.cols)
      
  colnames(level0.catalog) <- std.colnames

  out.data <- NULL
  for (this.row in 1:dim(level0.catalog)[1])
  {
    this.file <- as.character(level0.catalog[this.row,"File"])
    this.sheet <- as.character(level0.catalog[this.row, "Sheet"])
    this.skip <- as.numeric(level0.catalog[this.row, "Skip.Rows"]) 
    this.date <- as.character(level0.catalog[this.row, "Date"])
    this.chem <- as.character(level0.catalog[this.row, "Chemical.ID"])
    if (!(this.chem %in% chem.ids[,chem.lab.id.col]))
    {
      stop(paste0("Chem ID ",this.chem," not found in table chem.ids column ",
                  chem.lab.id.col))
    } else {
      id.index <- which(chem.ids[, chem.lab.id.col]==this.chem)
      this.name <- chem.ids[id.index, chem.name.col]
      this.dtxsid <- chem.ids[id.index, chem.dtxsid.col]
    }
    this.istd <- as.character(level0.catalog[this.row, "ISTD.Name"])
    this.sample.name.col <- as.character(level0.catalog[this.row, "Sample.ColName"])
    this.peak.col <- as.character(level0.catalog[this.row, "Peak.ColName"])
    this.istd.peak.col <- as.character(level0.catalog[this.row, "ISTD.Peak.ColName"])
    this.conc.col <- as.character(level0.catalog[this.row, "Conc.ColName"])
    this.type.col <- as.character(level0.catalog[this.row, "Type.ColName"])
    this.analysis.param.col <- as.character(level0.catalog[this.row, "AnalysisParam.ColName"])
# Read the data:
    this.data <- as.data.frame(read_excel(this.file, sheet=this.sheet, skip=this.skip))
# Trim the data if num.rows.col specified:
    if (!is.null(num.rows.col))
    {
      this.rows <- suppressWarnings(as.numeric(level0.catalog[this.row,num.rows.col]))
      if (!is.na(this.rows)) this.data <- this.data[1:this.rows,]
    }
# Annotate the data:
    this.data$Compound <- this.name
    this.data$DTXSID <- this.dtxsid
    this.data$Lab.Compound.ID <- this.chem
    this.data$Level0.File <- this.file
    this.data$level0.Sheet <- this.sheet
    this.data$Date <- this.date
    this.data$ISTD.Name <- this.istd
    needed.columns <- c("Compound",
                              "DTXSID",
                              "Lab.Compound.ID",
                              "Date",
                              this.sample.name.col,
                              this.type.col,
                              this.conc.col,
                              this.peak.col,
                              this.istd.peak.col,
                              "ISTD.Name",
                              this.analysis.param.col,
                              "Level0.File",
                              "level0.Sheet")

    if (!is.null(additional.colname.cols))
    {
      for (this.col in additional.colname.cols)
      {
        needed.columns <- c(needed.columns,
                            as.character(level0.catalog[this.row, this.col]))
      } 
    }
    cat(paste0(paste(this.file,this.sheet,this.chem,sep=", "),"\n"))
    reordered.data <- try(this.data[,needed.columns])
    if (is(reordered.data,"try-error")) 
    {
      print(paste("Columns needed:",paste(needed.columns,collapse=", ")))
      print(head(this.data))
      print(paste0("Missing columns: ",paste(needed.columns[!(needed.columns %in% colnames(this.data))],collapse=", ")))
      browser()
    }
    this.data <- reordered.data
    
    colnames(this.data)[1:13] <- c("Compound",
                             "DTXSID",
                             "Lab.Compound.ID",
                             "Date",
                             "Sample",
                             "Type",
                             "Compound.Conc",
                             "Peak.Area",
                             "ISTD.Peak.Area",
                             "ISTD.Name",
                             "Analysis.Params",
                             "Level0.File",
                             "Level0.Sheet"
                             )
    if (!is.null(additional.colname.cols))
    {
      colnames(this.data)[14:(13+length(additional.colname.cols))] <- 
        additional.colnames
    }
 
    out.data <- rbind(out.data, this.data[,1:(13+length(additional.colname.cols))])  
  }

  
  if (!is.null(OUTPUT.DIR)) {
    file.path <- OUTPUT.DIR
  } else {
    file.path <- getwd()
  }

  ## Keep outputting level-0 catalog for now but this functionality may be deprecated later
  if (catalog.out) {
    # Write out a "Catalog" file that explains how level 0 data were mapped to level 1
    write.table(level0.catalog, 
                file=paste(file.path, "/", FILENAME,"-level0-Catalog.tsv",sep=""),
                sep="\t",
                row.names=F,
                quote=F)
    cat(paste0("A Level-0 Catalog file named ",FILENAME,"-level0-Catalog.tsv", 
                " has been exported to the following directory: ", file.path),"\n")
  }

  if (output.res) {
    # Write out the merged Level-0 file
    write.table(out.data, 
                file=paste(file.path, "/", FILENAME,"-level0.tsv",sep=""),
                sep="\t",
                row.names=F,
                quote=F)
    cat(paste0("A Level-0 file named ",FILENAME,"-level0.tsv", 
                " has been exported to the following directory: ", file.path), "\n")
  }
 
  return(out.data)  
}


