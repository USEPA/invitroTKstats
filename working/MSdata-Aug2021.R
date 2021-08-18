# Here's the new R package for analyzing these data:

library(invitroTKstats)

# There are multiple packages for loading Excel files, but I've been using this
# one lately:
library(readxl)

# Change to the directory that has your Excel file 
# (this is where it is on my computer):
setwd("c:/users/jwambaug/git/invitroTKstats/working/")

# Assumed dilution factors:
CC.DILUTE <- 1
BLANK.DILUTE <- 1
AF.DILUTE <- 2*16
T5.DILUTE <- 5*16
T1.DILUTE <- 5*16

# There are five data files in this batch, store them in a list:
files <- list()
# This list stores the chemical ID's in each file
DTXSIDs <- list()           
# This list stores the internal standards for each chemical in each file:
ISTDs <- list()

files[[1]] <- read_excel(
           "SmeltzPFASAug2021/20210127_UCRepeat2_Data_MGS.xlsx",
           sheet=2,
           skip=5)

DTXSIDs[[1]] <- c(
  "DTXSID0060985",  # Compound 1:  Hexafluoroglutaryl Chloride
  "DTXSID20375106", # Compound 2:  Diox8-dioic
  "DTXSID60380390",  # Compound 3:  TFE-PFBS
  "DTXSID90315130", # Compound 4:  PFOS-Cl
  "DTXSID50469320",  # Compound 5:  PFHxSA
  "DTXSID3020209", # Compound 6:  n-Butylparaben
  "DTXSID30382104",  # Compound 7:  Cl-PFNA
  "DTXSID8051419",  # Compound 8:  PFOSAA
  "MPFBA", # Compound 9:  MPFBA
  "M3PFBS", # Compound 10:  M3PFBS
  "M9PFNA", # Compound 11:  M9PFNA
  "M8PFOS", # Compound 12:  M8PFOS
  "13C6 n-Butylparaben", # Compound 13:  13C6 n-Butylparaben
  "M8FOSA")  # Compound 14:  M8FOSA

ISTDs[[1]] <- c(
  "MPFBA",  # Compound 1:  Hexafluoroglutaryl Chloride
  "MPFBA", # Compound 2:  Diox8-dioic
  "M3PFBS",  # Compound 3:  TFE-PFBS
  "M8PFOS", # Compound 4:  PFOS-Cl
  "M8PFOS",  # Compound 5:  PFHxSA
  "13C6 n-Butylparaben", # Compound 6:  n-Butylparaben
  "M9PFNA",  # Compound 7:  Cl-PFNA
  "M8FOSA",  # Compound 8:  PFOSAA
  "Eponymous", # Compound 9:  MPFBA
  "Eponymous", # Compound 10:  M3PFBS
  "Eponymous", # Compound 11:  M9PFNA
  "Eponymous", # Compound 12:  M8PFOS
  "Eponymous", # Compound 13:  13C6 n-Butylparaben
  "Eponymous")  # Compound 14:  M8FOSA

files[[2]] <- read_excel(
           "SmeltzPFASAug2021/20210312_UC_PFAS36_Data_MGS.xlsx",
           sheet=2,
           skip=5)
           
DTXSIDs[[2]] <- c(
  "DTXSID80380256",  # Compound 1:  TFMFPA
  "DTXSID90558000",  # Compound 2:  6:2 monoPAP
  "DTXSID50379814",  # Compound 3:  PFEOES
  "DTXSID50226894",  # Compound 4:  9H-PFNA 
  "DTXSID40880025",  # Compound 5:  NaPFOA 
  "DTXSID50381073",  # Compound 6:  PFPE-5
  "DTXSID5061954",  # Compound 7:  11H-PFUnDA
  "DTXSID3020209",  # Compound 8:  n-Butylparaben
  "M5PFPeA",  # Compound 9:  M5PFPeA
  "M4PFHpA",  # Compound 10:  M4PFHpA
  "M8PFOA",  # Compound 11:  M8PFOA
  "M3PFHxS",  # Compound 12:  M3PFHxS
  "M9PFNA",  # Compound 13:  M9PFNA
  "M6PFDA",  # Compound 14:  M6PFDA
  "13C6 n-Butylparaben",  # Compound 15:  13C6 n-Butylparaben
  "M7PFUdA ")  # Compound 16:  M7PFUdA 

ISTDs[[2]] <- c(
  "M5PFPeA",  # Compound 1:  TFMFPA
  "M4PFHpA",  # Compound 2:  6:2 monoPAP
  "M3PFHxS",  # Compound 3:  PFEOES
  "M9PFNA",  # Compound 4:  9H-PFNA
  "M8PFOA",  # Compound 5:  NaPFOA  
  "M6PFDA",  # Compound 6:  PFPE-5
  "M7PFUdA",  # Compound 7:  11H-PFUnDA
  "13C6 n-Butylparaben",  # Compound 8:  n-Butylparaben
  "Eponymous",  # Compound 9:  M5PFPeA
  "Eponymous",  # Compound 10:  M4PFHpA
  "Eponymous",  # Compound 11:  M8PFOA
  "Eponymous",  # Compound 12:  M3PFHxS 
  "Eponymous",  # Compound 13:  M9PFNA      
  "Eponymous",  # Compound 14:  M6PFDA
  "Eponymous",  # Compound 15:  13C6 n-Butylparaben
  "Eponymous") # Compound 16:  M7PFUdA 
           
           
files[[3]] <- read_excel(
           "SmeltzPFASAug2021/20210503_UCRepeat1_Data_MGS.xlsx",
           sheet=2,
           skip=5)   

DTXSIDs[[3]] <- c(
  "DTXSID4059833",  # Compound 1:  octafluoroadipic
  "pyruvate")  # Compound 2:  pyruvate

ISTDs[[3]] <- c(
  "pyruvate",  # Compound 1:  octafluoroadipic
  "Eponymous") # Compound 2:  pyruvate

files[[4]] <- read_excel(
           "SmeltzPFASAug2021/20210504_UCRepeat2_Data_MGS.xlsx",
           sheet=2,
           skip=5)   

DTXSIDs[[4]] <- c(
  "DTXSID20375106",  # Compound 1:  diox8
  "DTXSID0060985",  # Compound 2:  hexafluoroCl
  "pyruvate")  # Compound 3:  pyruvate

ISTDs[[4]] <- c(
  "pyruvate",  # Compound 1:  diox8
  "pyruvate",  # Compound 2:  hexafluoroCl
  "Eponymous") # Compound 3:  pyruvate

files[[5]] <- read_excel(
           "SmeltzPFASAug2021/20210507_UC_PFAS36_EarlyEluter_Data_MGS.xlsx",
           sheet=2,
           skip=5)   

DTXSIDs[[5]] <- c(
  "DTXSID8059928",  # Compound 1:  tetrafluorosuccinate
  "DTXSID8059926",  # Compound 2:  hexafluoroglutarate
  "pyruvate")  # Compound 3:  pyruvate

ISTDs[[5]] <- c(
  "pyruvate",  # Compound 1:  tetrafluorosuccinate
  "pyruvate",  # Compound 2:  hexafluoroglutarate
  "Eponymous") # Compound 3:  pyruvate


# Now loop over the files and annotate the chemicals:
UC.data <- NULL
for (this.file in 1:length(files))
{
  total.rows <- dim(files[[this.file]])[1]
  this.compound <- 1
  
  files[[this.file]]$DTXSID <- DTXSIDs[[this.file]][this.compound] 
  files[[this.file]]$ISTD.Name <- ISTDs[[this.file]][this.compound] 
  for (this.row in 1:total.rows)
  {
    if (!is.na(as.character(files[[this.file]][this.row,1])))
      if (regexpr("Compound",as.character(files[[this.file]][this.row,1]))!=-1)
      {
        this.compound <- this.compound + 1
        files[[this.file]][this.row:total.rows,"DTXSID"] <- 
          DTXSIDs[[this.file]][this.compound]  
        files[[this.file]][this.row:total.rows,"ISTD.Name"] <- 
          ISTDs[[this.file]][this.compound] 
      }
  }
  UC.data <- rbind(UC.data,files[[this.file]])
}
dim(UC.data)

UC.data <- as.data.frame(UC.data)
UC.data <- subset(UC.data,!is.na(UC.data[,"Sample Text"]))
dim(UC.data)

# Extract the sample type from column Sample Text:
UC.data[regexpr("AF",unlist(UC.data[,"Sample Text"]))!=-1,"Sample.Type"] <- "AF"
UC.data[regexpr("T1",unlist(UC.data[,"Sample Text"]))!=-1,"Sample.Type"] <- "T1"
UC.data[regexpr("T5",unlist(UC.data[,"Sample Text"]))!=-1,"Sample.Type"] <- "T5"
UC.data[regexpr("CC",unlist(UC.data[,"Sample Text"]))!=-1,"Sample.Type"] <- "CC"
UC.data[regexpr("Blank",unlist(UC.data[,"Sample Text"]))!=-1,"Sample.Type"] <- "Blank"

# Get rid of unused samples (QC samples):
UC.data <- subset(UC.data,!is.na(Sample.Type))
dim(UC.data)

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
dim(UC.data)

# Convert to uM:
UC.data[,"Std.Conc"] <- as.numeric(unlist(UC.data[,"nM"]))/1000
# Get rid of CC samples that don't have a concentration:
UC.data <- subset(UC.data,Sample.Type=!"CC" | !is.na(Std.Conc))

dim(UC.data)



# Should update this with input from Marci:
UC.data$Analysis.Method <- "UPLC-MS/MS"
UC.data$Analysis.Instrument <- "Waters Xevo TQ-S micro (QEB0036)"
UC.data$Analysis.Parameters <- "None"

# ISTD cocn 1ppm:
UC.data$ISTD.Conc <- 1 #ppm

# Test chemical concentration:
UC.data$Test.Target.Conc <- 10 # uM

# Make the numeric values numeric:
UC.data[,"Area"] <- as.numeric(unlist(UC.data[,"Area"]))
UC.data[,"IS Area"] <- as.numeric(unlist(UC.data[,"IS Area"]))

level1 <- format_fup_uc(UC.data,
  FILENAME="SmeltzPFASAug2021/SmeltzAug2021",
  sample.col="Name",
  compound.col="DTXSID",
  compound.conc.col="Std.Conc", 
  lab.compound.col="DTXSID", 
  type.col="Sample.Type", 
  istd.col="IS Area"
  )
 
level2 <- level1
# Taking all data in spreadsheet as human verfied for starters
level2$Verified <- "Y"
# From Marci:
level2[level2$DTXSID=="DTXSID0060985" & level2$Date=="20201124","Verified"] <- 
  "MS-Drop"
level2[level2$DTXSID=="DTXSID20375106" & level2$Date=="20201124","Verified"] <- 
  "MS-Drop"
level2[level2$DTXSID=="DTXSID8059926" & level2$Date=="20210225","Verified"] <- 
  "MS-Drop"
level2[level2$DTXSID=="DTXSID8059928" & level2$Date=="20210225","Verified"] <- 
  "MS-Drop"

write.table(level2,
  file="SmeltzPFASAug2021/SmeltzAug2021-PPB-UC-Level2.tsv",
  sep="\t",
  row.names=F,
  quote=F)

level3 <- calc_fup_uc_point(FILENAME="SmeltzPFASAug2021/SmeltzAug2021") 

level4 <- calc_fup_uc(FILENAME="SmeltzPFASAug2021/SmeltzAug2021") 
