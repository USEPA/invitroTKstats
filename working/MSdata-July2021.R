# Here's the new R package for analyzing these data:

library(invitroTKstats)

# There are multiple packages for loading Excel files, but I've been using this
# one lately:
library(readxl)

# Change to the directory that has your Excel file 
# (this is where it is on my computer):
setwd("c:/users/jwambaug/git/invitroTKstats/working/")

# load the data one sheet at a time from the Excel file, we skip the first row
# so that we get more of the column names loaded:
DTXSIDs <- c(
  "DTXSID4059833",
  "DTXSID3037707",
  "DTXSID1037303",
  "DTXSID3031860",
  "DXTSID80108992",
  "DTXSID50375114",
  "DTXSID90868151",
  "DTXSID3059921",
  "DTXSID3020209")

# Merge all data into single table:
UC.data <- NULL
for (this.sheet in 6:14)
{
  temp <- read_excel("PFAS LC-MS Plasma Protein Binding Report - 20201110 UC Assay_Reviewed.xlsx",
    sheet=this.sheet,
    skip=1)
# Annotate each data set with the chemical ID
  temp$DTXSID <- DTXSIDs[this.sheet-5]
# Separate the internal standard from the test chemical:
  sep.row <- which(is.na(unlist(temp[,11])))
# Find the name of the internal standard:
  temp$ISTD.Name<-as.character(temp[sep.row+1,11])
# Subset to just the test chemical:
  temp <- temp[1:(sep.row-1),]
  UC.data <- rbind(UC.data,temp)
}

# Assumed dilution factors:
CC.DILUTE <- 1
BLANK.DILUTE <- 1
AF.DILUTE <- 2*16
T5.DILUTE <- 5*16
T1.DILUTE <- 5*16

# Extract the sample type from column Sample Text:
UC.data[regexpr("AF",unlist(UC.data[,"Sample Text"]))!=-1,"Sample.Type"] <- "AF"
UC.data[regexpr("T1",unlist(UC.data[,"Sample Text"]))!=-1,"Sample.Type"] <- "T1"
UC.data[regexpr("T5",unlist(UC.data[,"Sample Text"]))!=-1,"Sample.Type"] <- "T5"
UC.data[regexpr("CC",unlist(UC.data[,"Sample Text"]))!=-1,"Sample.Type"] <- "CC"
UC.data[regexpr("Blank",unlist(UC.data[,"Sample Text"]))!=-1,"Sample.Type"] <- "Blank"

# Get rid of unused samples (QC samples):
UC.data <- subset(UC.data,!is.na(Sample.Type))

# Extract the standard concentration fro column Sample Text:
UC.data[,"Std.Conc"] <- unlist(lapply(
  strsplit(gsub(" pg/uL","",unlist(UC.data[,"Sample Text"]))," - "),
  function(x) ifelse(length(x)==2,as.numeric(x[2]),NA)))
UC.data[,"Std.Units"] <- "pg/uL" 

# No replicates:
UC.data$Series <- 1
# Extract the date from the time stamp
UC.data$Date <- unlist(lapply(strsplit(unlist(UC.data[,"Name"]),"_PFAS"),
  function(x) x[1]))
# Assume all samples analyzed on same date were the same calibration
UC.data$Cal <- UC.data$Date
# We'll need to go back and set this per sample type:
UC.data[UC.data[,"Sample.Type"]=="AF","Dilution.Factor"] <- AF.DILUTE
UC.data[UC.data[,"Sample.Type"]=="T1","Dilution.Factor"] <- T1.DILUTE
UC.data[UC.data[,"Sample.Type"]=="T5","Dilution.Factor"] <- T5.DILUTE
UC.data[UC.data[,"Sample.Type"]=="CC","Dilution.Factor"] <- CC.DILUTE
UC.data[UC.data[,"Sample.Type"]=="Blank","Dilution.Factor"] <- BLANK.DILUTE
 
# Treat the blanks as calibration data with concentration 0:
UC.data[UC.data[,"Sample.Type"]=="Blank","Std.Conc"] <- 0
UC.data[UC.data[,"Sample.Type"]=="Blank","Std. Conc (nM)"] <- "0"

UC.data[UC.data[,"Sample.Type"]=="Blank","Sample.Type"] <- "CC"

# Convert to uM:
UC.data[,"Std.Conc"] <- as.numeric(unlist(UC.data[,"Std. Conc (nM)"]))/1000



# Should update this with input from Marci:
UC.data$Analysis.Method <- "UPLC-MS/MS"
UC.data$Analysis.Instrument <- "Waters Xevo TQ-S micro (QEB0036)"
UC.data$Analysis.Parameters <- "None"

# ISTD cocn 1ppm:
UC.data$ISTD.Conc <- 1 #ppm

UC.data$Test.Target.Conc <- 0.01 # uM
# Make the numeric values numeric:
UC.data[,"Area"] <- as.numeric(unlist(UC.data[,"Area"]))
UC.data[,"IS Area"] <- as.numeric(unlist(UC.data[,"IS Area"]))

level1 <- format_fup_uc(UC.data,
  FILENAME="Smeltz2021",
  sample.col="Name",
  compound.col="DTXSID",
  compound.conc.col="Std.Conc", 
  lab.compound.col="DTXSID", 
  type.col="Sample.Type", 
  istd.col="IS Area"
  )
 
level2 <- level1
# Taking all data in spreadsheet as human verfied
level2$Verified <- "Y"

write.table(level2,
  file="Smeltz2021-PPB-UC-Level2.tsv",
  sep="\t",
  row.names=F,
  quote=F)

level3 <- calc_fup_uc_point(FILENAME="Smeltz2021") 

level4 <- calc_fup_uc(FILENAME="Smeltz2021") 
