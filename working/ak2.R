# There are multiple packages for loading Excel files, but I've been using this
# one lately:
library(readxl)
# This package allows us to convert rows and columns od data to different
# formats:
library(reshape2)

# Change to the directory that has your Excel file 
# (this is where it is on my computer):
setwd("c:/users/jwambaug/git/invitroTKstats/working/")

# load the data one sheet at a time from the Excel file, we skip the first row
# so that we get more of the column names loaded:
sheet1 <- read_excel(
  "021220_915_965_476_267_906_273_913_899_900_analysis_090920.xlsx",
  sheet=6,
  skip=1)

# For a level 1 file we need 1 row per observation, so we can use the "melt"
# function from the package reshape2 to simplify the more complicated structure 
# used in the Wetmore lab data sheets 
# (see https://www.r-bloggers.com/2012/04/melt/):
chem.data1.mt <- melt(
  sheet1, 
  id.vars = c(
    "Name",
    "Data File",
    "Type",
    "Acq. Date-Time"
    ) 
  )
  
  
# To figure out which columns we need for which variables, it helps to add
# the column numbers to the excel file (I used the =COLUMN() command:
  
# Now lets find the internal standard curves:
chem.data1.ISTDs <- list()

# BET is given by "Area...16" 
# (the sixteenth column of the sheet, which is named "Area")
chem.data1.ISTDs[["BET"]] <-subset(chem.data1.mt,variable=="Area...16")
chem.data1.ISTDs[["BET"]]$Internal.Standard <- "BET"

# HET is given by "Area...24" 
# (the twenty fourth column of the sheet, which is named "Area")
chem.data1.ISTDs[["HET"]] <-subset(chem.data1.mt,variable=="Area...24")
chem.data1.ISTDs[["HET"]]$Internal.Standard <- "HET"

# MFHET is given by "Area...32" 
# (the thirty second column of the sheet, which is named "Area")
chem.data1.ISTDs[["MFHET"]] <-subset(chem.data1.mt,variable=="Area...32")
chem.data1.ISTDs[["MFHET"]]$Internal.Standard <- "MFHET"

# MFHET is given by "Area...40" 
# (the thirty second column of the sheet, which is named "Area")
chem.data1.ISTDs[["DET"]] <-subset(chem.data1.mt,variable=="Area...40")
chem.data1.ISTDs[["DET"]]$Internal.Standard <- "DET"

# 4NT13C6 is given by "Area...80"
chem.data1.ISTDs[["4NT13C6"]] <-subset(chem.data1.mt,variable=="Area...80")
chem.data1.ISTDs[["4NT13C6"]]$Internal.Standard <- "4NT13C6"

# Lets make a function to automate anotation of Wetmore lab data:
convert_to_level1 <- function(data.set, ISTD.set,
  chem.name, area.col.num, ISTD, analysis.method,
  instrument, inst.param.offset = -3, conc.offset = -2, area.base="Area...",
  inst.param.base="RT...",conc.base="Final Conc....")
{
  area.col <- paste(area.base,area.col.num,sep="")
  instr.param.col <- paste(inst.param.base,area.col.num+inst.param.offset,sep="")
  conc.col <- paste(conc.base,area.col.num+conc.offset,sep="")
  
  out <- subset(data.set,variable==area.col)
  out$Lab.Compound.Name <- chem.name
# (the merge function matches columns with like data):
  out <- merge(out,ISTD.set[[ISTD]], by=c(
    "Name",
    "Data File",
    "Type",
    "Acq. Date-Time"
    ) 
  )
  # Change the value column names to reflect what we've got:
  colnames(out)[colnames(out)=="value.x"] <- "Area"
  colnames(out)[colnames(out)=="value.y"] <- "ISTD.Area"
  # We don't need the "variable columns anymore:'
  out <- out[,!(colnames(out) %in% c("variable.x","variable.y"))]
  # Add in instrument parameters:
  out <- merge(out,
  subset(data.set,variable==instr.param.col), 
  by=c(
    "Name",
    "Data File",
    "Type",
    "Acq. Date-Time"
    ) 
  )
  # Change the value column names to reflect what we've got:
  colnames(out)[colnames(out)=="value"] <- "Instrument.Param"
  # We don't need the "variable columns anymore:'
  out <- out[,!(colnames(out) %in% c("variable"))]
  # Add in intended conc for cal curves:
  out <- merge(out,
  subset(data.set,variable==conc.col), 
  by=c(
    "Name",
    "Data File",
    "Type",
    "Acq. Date-Time"
    ) 
  )
  # Change the value column names to reflect what we've got:
  colnames(out)[colnames(out)=="value"] <- "Conc"
  # We don't need the "variable columns anymore:'
  out <- out[,!(colnames(out) %in% c("variable"))]
  # Set instrument parameters:
  out$Analysis.Method <- analysis.method
  out$Instrument <- instrument

  return(out)
}

# Then do each test chemical:
# 915 is given by "Area...14"
chem.data1.915 <- convert_to_level1(
  chem.data1.mt,
  chem.data1.ISTDs,
  chem.name="915",
  area.col=14,
  ISTD="BET",
  analysis.method="GC",
  instrument="Some Machine 3000")
  
# 965 is given by "Area...22", we also need to change the ISTD, instrument
# param (retention time) column, and intended conc columns appropriately:
chem.data1.965 <- convert_to_level1(
  chem.data1.mt,
  chem.data1.ISTDs,
  chem.name="965",
  area.col="Area...22",
  ISTD="HET",
  instr.param="RT...19",
  intended.conc="Final Conc....20",
  analysis.method="GC",
  instrument="Some Machine 3000")

# 476 is given by "Area...30", we also need to change the ISTD, instrument
# param (retention time) column, and intended conc columns appropriately: 
chem.data1.476 <- convert_to_level1(
  chem.data1.mt,
  chem.data1.ISTDs,
  chem.name="476",
  area.col="Area...30",
  ISTD="MFHET",
  instr.param="RT...27",
  intended.conc="Final Conc....28",
  analysis.method="GC",
  instrument="Some Machine 3000")

# 267 is given by "Area...38", we also need to change the ISTD, instrument
# param (retention time) column, and intended conc columns appropriately: 
chem.data1.267 <- convert_to_level1(
  chem.data1.mt,
  chem.data1.ISTDs,
  chem.name="267",
  area.col="Area...38",
  ISTD="DET",
  instr.param="RT...27",
  intended.conc="Final Conc....28",
  analysis.method="GC",
  instrument="Some Machine 3000")
  
  



