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

assayinfo <- read_excel(
           "20220201_PFAS-LC_FractionUnbound_MGS.xlsx",
           sheet=1)
assayinfo <- as.data.frame(assayinfo)

cheminfo <- read_excel(
           "20220201_PFAS-LC_FractionUnbound_MGS.xlsx",
           sheet=2)
cheminfo <- as.data.frame(cheminfo)



# Fill in the date column for all rows:
this.date <- assayinfo[1,"LCMS Analysis Date"]
this.row <- 1
while (this.row <= dim(assayinfo)[1])
{
  if (is.na(assayinfo[this.row,"LCMS Analysis Date"]))
  {
    assayinfo[this.row,"LCMS Analysis Date"] <- this.date 
  } else {
    this.date <- assayinfo[this.row,"LCMS Analysis Date"] 
  }
  this.row <- this.row + 1
}
# Get rid of the UTC:
assayinfo[,1] <- sapply(assayinfo[,1],function(x) gsub("  UTC","",x))

# Define the info for control chemicals:
controls <- data.frame(
  Compound="n-butylparaben",
  DTXSID="DTXSID3020209",
  ISTD="13C6-n-butylparaben")
                                     
# All the data (for now) are in different tabs (by date) of a single file 
# (still going with a list called files but each "file") is now a tab:
files <- list()
# This list stores the chemical ID's in each file
compounds <- list()           
DTXSIDs <- list()           
# This list stores the internal standards for each chemical in each file:
ISTDs <- list()

for (this.date in unique(assayinfo[,1]))
{
  this.sheet <- gsub("-","",this.date)
  this.index <- which(unique(assayinfo[,1])==this.date)
      
  files[[this.index]] <- as.data.frame(read_excel(
           "20220201_PFAS-LC_FractionUnbound_MGS.xlsx",
           sheet=this.sheet))
  files[[this.index]]$Date <- this.date        
  this.assayinfo <- subset(assayinfo,regexpr(this.date,assayinfo[,1])!=-1)
  
  compounds[[this.index]] <- this.assayinfo[,"Analyte"]
  DTXSIDs[[this.index]] <- this.assayinfo[,"DTXSID"]
  ISTDs[[this.index]] <- this.assayinfo[,"Int Std"]
}
         


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
        this.id <- NA
        # Try to get the id using two spaces after colon:
        this.id <- 
          strsplit(as.character(files[[this.file]][this.row,1]),":  ")[[1]][2]
        # If is.na, try with only one:
        if (is.na(this.id)) this.id <- 
          strsplit(as.character(files[[this.file]][this.row,1]),": ")[[1]][2]
        
        this.chem.index <- regexpr(paste("\\(",this.id,"\\)",sep=""), 
          cheminfo[,2])!=-1
        # check to see if this is one of the test chemicals:
        if (any(this.chem.index))
        {
          this.dtxsid <- cheminfo[this.chem.index,"DTXSID"]
          if (length(this.dtxsid)>1)
          {
            match.strings <- unlist(lapply(
              strsplit(cheminfo[this.chem.index,2],this.id),
              function(x) x[[1]]))
            match.strings <- gsub(" \\(","",match.strings)
            for (this.string in match.strings)
              if (any(regexpr(tolower(this.string),
                tolower(compounds[[this.file]]))!=-1))
              {
                this.assay.index <- which(regexpr(tolower(this.string),
                tolower(compounds[[this.file]]))!=-1)
                files[[this.file]][this.row:total.rows,"Compound.Name"] <- 
                  compounds[[this.file]][this.assay.index]
                files[[this.file]][this.row:total.rows,"DTXSID"] <- 
                  DTXSIDs[[this.file]][this.assay.index]
                this.dtxsid <- DTXSIDs[[this.file]][this.assay.index]
                files[[this.file]][this.row:total.rows,"ISTD.Name"] <- 
                  ISTDs[[this.file]][this.assay.index]
              }
          }
          this.assay.index <- which(DTXSIDs[[this.file]] == this.dtxsid)
          if (any(this.assay.index))
          {
            files[[this.file]][this.row:total.rows,"Compound.Name"] <- 
              compounds[[this.file]][this.assay.index]
            files[[this.file]][this.row:total.rows,"DTXSID"] <- 
              this.dtxsid
            files[[this.file]][this.row:total.rows,"ISTD.Name"] <- 
              ISTDs[[this.file]][this.assay.index]
          }
        # Check to see if we can find the chemical in the assay table:
        } else if (any(
          regexpr(tolower(this.id), tolower(assayinfo[,"Analyte"]))!=-1 &
          assayinfo[,"LCMS Analysis Date"] == files[[this.file]][this.row,"Date"]))
        {
          this.assay.index <- which(any(
            regexpr(tolower(this.id), tolower(assayinfo[,"Analyte"]))!=-1 &
            assayinfo[,"LCMS Analysis Date"] == files[[this.file]][this.row,"Date"]))
          files[[this.file]][this.row:total.rows,"Compound.Name"] <- 
            compounds[[this.file]][this.assay.index]
          files[[this.file]][this.row:total.rows,"DTXSID"] <- 
            DTXSIDs[[this.file]][this.assay.index]
          files[[this.file]][this.row:total.rows,"ISTD.Name"] <- 
            ISTDs[[this.file]][this.assay.index]          
        } else {
        # Maybe it's a control:
          this.index <- regexpr(tolower(this.id), tolower(controls[,"Compound"]))!=-1
          if (any(this.index))
          {
            this.compound <- this.compound + 1
            files[[this.file]][this.row:total.rows,"Compound.Name"] <- 
              controls[this.index,"Compound"]
            files[[this.file]][this.row:total.rows,"DTXSID"]  <- 
              controls[this.index,"DTXSID"]
            files[[this.file]][this.row:total.rows,"ISTD.Name"]  <- 
              controls[this.index,"ISTD"]
          } else {
          # Assume it's an internal standard:
            files[[this.file]][this.row:total.rows,"Compound.Name"] <- this.id
            files[[this.file]][this.row:total.rows,"DTXSID"] <- "ISTD"
            files[[this.file]][this.row:total.rows,"ISTD.Name"] <- "Eponymous"
          }
        } 
      }
  }
  UC.data <- rbind(UC.data,files[[this.file]])
}
dim(UC.data)
misses <- unique(subset(UC.data,DTXSID=="ISTD")$Compound.Name)
misses <- misses[!(misses %in% unique(assayinfo[,"Int Std"]))]

miss.table <- NULL
for (i in 1:length(misses))
{
  this.miss <- misses[i]
  dates <- unique(subset(UC.data,Compound.Name==this.miss)$Date)
  for (this.date in dates)
  {
    this.row <- data.frame(Compound = this.miss, Date = this.date)
    miss.table <- rbind(miss.table,this.row)
  }
}
miss.table[order(miss.table$Date),]

write.csv(miss.table[order(miss.table$Date),],file="chems-not-identified.csv",row.names=FALSE)







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
