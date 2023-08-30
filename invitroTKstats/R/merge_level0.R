#' Merge multiple level0 files into a single table for processing
#'
#' This function reads multiple Excel files containing mass-spectrometry (MS) data
#' and extracts the chemical sample data from the specified
#' sheets. The argument level0.catalog is a table that
#' provides the necessary information to find the data for each chemical. The
#' primary data of interest are the analyte peak area, the internal standard
#' peak area, and the target concentration for calibration curve (CC) samples.
#' The argument data.label is used to annotate this particular mapping of level0
#' files into data ready to be organized into a level1 file.
#'
#' Unless specified to be a single value for all the files, for example sheet="Data",
#' the argument level0.catalog should be a data frame with the following columns:
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
#' from the specified file and sheet.
#'
#' @param data.label A string used to identify outputs of the function call.
#' (defaults to "MYDATA")
#' 
#' @param level0.catalog A data frame describing which columns of which sheets
#' in which Excel files contain MS data for analysis. See details for full
#' explanation.
#' 
#' @param sample.col Which column of clint.data indicates the unique mass 
#' spectrometry (MS) sample name used by the laboratory. (Defaults to 
#' "Lab.Sample.Name")
#' 
#' @param lab.compound.col Which column of clint.data indicates The test compound 
#' name used by the laboratory (Defaults to "Lab.Compound.Name")
#' 
#' @param dtxsid.col Which column of clint.data indicates EPA's DSSTox Structure 
#' ID (\url{http://comptox.epa.gov/dashboard}) (Defaults to "DTXSID")
#' 
#' @param date.col Which column of clint.data indicates the laboratory measurement
#' date (Defaults to "Date")
#' 
#' @param compound.col Which column of clint.data indicates the test compound
#' (Defaults to "Compound.Name")
#' 
#' @param area.col Which column of clint.data indicates the target analyte (that 
#' is, the test compound) MS peak area (Defaults to "Area")
#' 
#' @param series.col Which column of clint.data indicates the "series", that is
#' a simultaneous replicate (Defaults to "Series")
#' 
#' @param type.col Which column of clint.data indicates the sample type (see table
#' above)(Defaults to "Sample.Type")
#' 
#' @param cal.col Which column of clint.data indicates the MS calibration -- for
#' instance different machines on the same day or different days with the same
#' MS analyzer (Defaults to "Cal")
#' 
#' @param dilution.col Which column of clint.data indicates how many times the
#' sample was diluted before MS analysis (Defaults to "Dilution.Factor")
#'
#' @param density.col Which column of clint.data indicates the density (units of
#' millions of hepatocytes per mL) hepatocytes in the in vitro incubation 
#' (Defaults to "Hep.Density" )
#' 
#' @param istd.col Which column of clint.data indicates the MS peak area for the
#' internal standard (Defaults to "ISTD.Area")
#' 
#' @param istd.name.col Which column of clint.data indicates identity of the 
#' internal standard (Defaults to "ISTD.Name")
#' 
#' @param istd.conc.col Which column of clint.data indicates the concentration 
#' (units if uM) of
#' the internal standard (Defaults to "ISTD.Conc")
#' 
#' @param conc.col Which column of clint.data indicates the intended
#' test chemical concentration 
#' (units if uM) of
#' at time zero (Defaults to "Conc") 
#'
#' @param time.col Which column of clint.data indicates the intended
#' time of the measurement (in minutes) since the test chemical was introduced
#' into the hepatocyte incubation (Defaults to "Time") 
#'
#' @param analysis.method.col Which column of PPB.data indicates the analytical
#' chemistry analysis method, typically "LCMS" or "GCMS" (Defaults to 
#' "Analysis.Method")
#'
#' @param analysis.instrument.col Which column of PPB.data indicates the 
#' instrument used for chemical analysis, for example 
#' "Agilent 6890 GC with model 5973 MS" (Defaults to 
#' "Analysis.Instrument")
#'
#' @param analysis.parameters.col Which column of PPB.data indicates the 
#' parameters used to identify the compound on the chemical analysis instrument,
#' for example 
#' "Negative Mode, 221.6/161.6, -DPb=26, FPc=-200, EPd=-10, CEe=-20, CXPf=-25.0"
#' (Defaulys to "Analysis.Parameters"). 
#'
#' @param FILENAME A string used to identify outputs of the function call.
#' (defaults to "MYDATA")
#' 
#' @param input.data A data frame containing mass-spectrometry peak areas,
#' indication of chemical identity, and measurement type. The data frame should
#' contain columns with names specified by the following arguments:
#' 
#' @param sample.col Which column of input.data indicates the unique mass 
#' spectrometry (MS) sample name used by the laboratory. (Defaults to 
#' "Lab.Sample.Name")
#' 
#' @param density.col Which column of input.data indicates the density of 
#' hepatocytes in suspension (10^6 hepatocytes / mL) (Defaults to "Hep.Density")
#' 
#' @param density.col A single value to be used for all samples indicating
#' the density of hepatocytes in suspension (10^6 hepatocytes / mL) 
#' (Defaults to NULL)
#' 
#' @param lab.compound.col Which column of input.data indicates The test compound 
#' name used by the laboratory (Defaults to "Lab.Compound.Name")
#' 
#' @param dtxsid.col Which column of input.data indicates EPA's DSSTox Structure 
#' ID (\url{http://comptox.epa.gov/dashboard}) (Defaults to "DTXSID")
#' 
#' @param date.col Which column of input.data indicates the laboratory measurement
#' date (Defaults to "Date")
#' 
#' @param series.col Which column of PPB.data indicates the "series", that is
#' a simultaneous replicate with the same analytical chemistry 
#' (Defaults to "Series")
#' 
#' @param series If this argument is used (defaults to NULL) every observation 
#' in the table is assigned the value of the argument and the corresponding
#' column in input.table (if present) is ignored.
#' 
#' @param compound.col Which column of input.data indicates the test compound
#' (Defaults to "Compound.Name")
#' 
#' @param area.col Which column of input.data indicates the target analyte (that 
#' is, the test compound) MS peak area (Defaults to "Area")
#' 
#' @param type.col Which column of input.data indicates the sample type (see table
#' above)(Defaults to "Type")
#' 
#' @param type.col Which column of input.data indicates the direction of the 
#' measurements (either "AtoB" for apical to basolateral or "BtoA" for vice 
#' versa) (Defaults to "Direction")
#' 
#' @param cal.col Which column of input.data indicates the MS calibration -- for
#' instance different machines on the same day or different days with the same
#' MS analyzer (Defaults to "Cal")
#' 
#' @param cal If this argument is used (defaults to NULL) every observation in
#' the table is assigned the value of the argument and the corresponding
#' column in input.table (if present) is ignored.
#' 
#' #param compound.conc.col Which column indicates the intended concentration 
#' of the test chemical for calibration curves (Defaults to "Standard.Conc")
#'
#' @param dilution If this argument is used (defaults to NULL) every 
#' observation in the table is assigned the value of the argument and the 
#' corresponding column in input.table (if present) is ignored.
#' 
#' 
#' @param istd.col Which column of input.data indicates the MS peak area for the
#' internal standard (Defaults to "ISTD.Area")
#' 
#' @param istd.name.col Which column of input.data indicates identity of the 
#' internal standard (Defaults to "ISTD.Name")
#' 
#' @param istd.name If this argument is used (defaults to NULL) every 
#' observation in the table is assigned the value of the argument and the 
#' corresponding column in input.table (if present) is ignored.
#' 
#' @param istd.conc.col Which column of input.data indicates the concentration of
#' the internal standard (Defaults to "ISTD.Conc")
#' 
#' @param istd.conc If this argument is used (defaults to NULL) every 
#' observation in the table is assigned the value of the argument and the 
#' corresponding column in input.table (if present) is ignored.
#' 
#' @param nominal.test.conc.col Which column of input.data indicates the intended
#' test chemical concentration at time zero in the dosing solution (added to the
#' donor side of the Caco-2 test well) (Defaults to "Test.Target.Conc") 
#' 
#' @param nominal.test.conc If this argument is used (defaults to NULL) every 
#' observation in the table is assigned the value of the argument and the 
#' corresponding column in input.table (if present) is ignored.
#'
#' @param analysis.method.col Which column of input.data indicates the analytical
#' chemistry analysis method, typically "LCMS" or "GCMS" (Defaults to 
#' "Analysis.Method")
#' 
#' @param analysis.method If this argument is used (defaults to NULL) every 
#' observation in the table is assigned the value of the argument and the 
#' corresponding column in input.table (if present) is ignored.
#'
#' @param analysis.instrument.col Which column of input.data indicates the 
#' instrument used for chemical analysis, for example 
#' "Agilent 6890 GC with model 5973 MS" (Defaults to 
#' "Analysis.Instrument")
#' 
#' @param analysis.instrument If this argument is used (defaults to NULL) every 
#' observation in the table is assigned the value of the argument and the 
#' corresponding column in input.table (if present) is ignored.
#'
#' @param analysis.parameters.col Which column of input.data indicates the 
#' parameters used to identify the compound on the chemical analysis instrument,
#' for example 
#' "Negative Mode, 221.6/161.6, -DPb=26, FPc=-200, EPd=-10, CEe=-20, CXPf=-25.0"
#' (Defaults to "Analysis.Parameters"). 
#' 
#' @param analysis.parameters If this argument is used (defaults to NULL) every 
#' observation in the table is assigned the value of the argument and the 
#' corresponding column in input.table (if present) is ignored.
#'
#' @param level0.file.col Which column of PPB.data indicates the file from
#' which the data were obtained (for example "MyWorkbook.xlsx").
#'  
#' @param level0.file If this argument is used (defaults to NULL) every
#' observation in the table is assigned the value of the argument and the
#' corresponding column in input.table (if present) is ignored.
#' 
#' @param level0.sheet.col Which column of PPB.data indicates the specific 
#' sheet containing the data if the file is an Excel workbook
#'  
#' @param level0.sheet If this argument is used (defaults to NULL) every
#' observation in the table is assigned the value of the argument and the
#' corresponding column in input.table (if present) is ignored.
#' 
#' @return \item{data.frame}{A data.frame in standardized "level1" format} 
#'
#' @author John Wambaugh
#' 
#' @export merge_level0
merge_level0 <- function(data.label="MYDATA",
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
  chem.dtxsid.col="DTXSID"
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

# Write out a "Catalog" file that explains how level 0 data were mapped to level 1
  write.table(level0.catalog, 
    file=paste(data.label,"-level0-Catalog.tsv",sep=""),
    sep="\t",
    row.names=F,
    quote=F)

  return(out.data)  
}


