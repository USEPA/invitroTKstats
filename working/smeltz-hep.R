library(invitroTKstats)
library(readxl)

setwd("c:/users/jwambaug/git/invitroTKstats/working/")

chem.ids <- read_excel("Hep12 Data for Uncertainty Feb2022.xlsx", sheet=1)
chem.ids <- as.data.frame(chem.ids)
 
smeltz.hep.ref <- read_excel("Hep12 Data for Uncertainty Feb2022.xlsx", sheet=2, skip=6)
smeltz.hep.ref <- as.data.frame(smeltz.hep.ref)

smeltz.hep.pfas <- read_excel("Hep12 Data for Uncertainty Feb2022.xlsx", sheet=3, skip=2)
smeltz.hep.pfas <- as.data.frame(smeltz.hep.pfas)
colnames(smeltz.hep.pfas) <- colnames(smeltz.hep.ref)

# Merge the ref chems and pfas:
smeltz.hep <- rbind(smeltz.hep.ref,smeltz.hep.pfas)

# get rid of blank rows:
smeltz.hep <- subset(smeltz.hep, !is.na(smeltz.hep[,1]))

# Annotate Compounds
this.compound <- "Phenacetin"
this.dtxsid <- chem.ids[
  tolower(chem.ids[,"Analyte Name"]) == tolower(this.compound),
  "DTXSID"]
this.istdname <- chem.ids[
  tolower(chem.ids[,"Analyte Name"]) == tolower(this.compound),
  "Internal Standard"]
  
this.row <- 1
while (this.row <= dim(smeltz.hep)[1])
{
  if (regexpr("Compound",smeltz.hep[this.row,1])!=-1)
  {
    this.compound <- strsplit(smeltz.hep[this.row,1],"  ")[[1]][2] 
    if (tolower(this.compound) %in% tolower(chem.ids[,"Internal Standard"]))
    {
      this.istd <- this.compound
      this.compound <- "ISTD"
      this.dtxsid <- NA
    } else {
      row.index <- which(regexpr(tolower(this.compound),
        tolower(chem.ids[,"Analyte Name"]))!=-1)
      this.dtxsid <- chem.ids[row.index, "DTXSID"]
      this.istdname <- chem.ids[row.index, "Internal Standard"]
    }
  }
  smeltz.hep[this.row, "Compound"] <- this.compound
  smeltz.hep[this.row, "DTXSID"] <- this.dtxsid
  smeltz.hep[this.row, "ISTD.Name"] <- this.istdname
  this.row <- this.row + 1
}

# Get rid of internal standard data:
smeltz.hep <- subset(smeltz.hep, Compound != "ISTD")

# Use invitroTKstats annotation of type:
smeltz.hep <- subset(smeltz.hep,!is.na(Type))
smeltz.hep[smeltz.hep$Type == "Analyte", "Type"] <- "Cvst"

# Add time of sample:
smeltz.hep <- subset(smeltz.hep,!is.na(smeltz.hep[,"Sample Text"]))
smeltz.hep[,"Time"] <- NA
smeltz.hep[regexpr("t240",tolower(smeltz.hep[,"Sample Text"]))!=-1,"Time"] <- 240/60
smeltz.hep[regexpr("t120",tolower(smeltz.hep[,"Sample Text"]))!=-1,"Time"] <- 120/60
smeltz.hep[regexpr("t60",tolower(smeltz.hep[,"Sample Text"]))!=-1,"Time"] <- 60/60
smeltz.hep[regexpr("t30",tolower(smeltz.hep[,"Sample Text"]))!=-1,"Time"] <- 30/60
smeltz.hep[regexpr("t15",tolower(smeltz.hep[,"Sample Text"]))!=-1,"Time"] <- 15/60
smeltz.hep[regexpr("t0",tolower(smeltz.hep[,"Sample Text"]))!=-1,"Time"] <- 0/60
smeltz.hep <- subset(smeltz.hep,!is.na(Time) | Type != "Cvst")

# Indicate whether hepatocytes have been inactivated:
smeltz.hep$Active.Hep <- NA

smeltz.hep[regexpr("living",tolower(smeltz.hep[,"Sample Text"]))!=-1,
  "Active.Hep"] <- 1
smeltz.hep[regexpr("inactive",tolower(smeltz.hep[,"Sample Text"]))!=-1,
  "Active.Hep"] <- 0
  
# For now we don't handle the inactives:
smeltz.hep <- subset(smeltz.hep, 
  regexpr("inactive",tolower(smeltz.hep[,"Sample Text"]))==-1)
  
# Make sure the areas are numeric:
smeltz.hep$Area <- as.numeric(smeltz.hep$Area)
smeltz.hep[,"IS Area"] <- as.numeric(smeltz.hep[,"IS Area"])
smeltz.hep[,"Std. Conc"] <- as.numeric(smeltz.hep[,"Std. Conc"])
smeltz.hep <- subset(smeltz.hep, !is.na(Area) & 
  !is.na(smeltz.hep[,"IS Area"]))
  
  
level1 <- format_clint(smeltz.hep,
  sample.col ="Name",
  date.col="Acq.Date",
  compound.col="Compound",
  lab.compound.col="Compound",
  type.col="Type",
  dilution=1,
  cal=1,
  istd.conc = 10/1000,
  istd.col= "IS Area",
  density = 0.5,
  conc = 1,
  time.col = "Time",
  analysis.method = "LCMS",
  analysis.instrument = "Unknown",
  analysis.parameters = "RT"
  )

level2 <- level1
level2$Verified <- "Y"
  
write.table(level2,
  file="Smeltz2022-Clint-Level2.tsv",
  sep="\t",
  row.names=F,
  quote=F)

level3 <- calc_clint_point(FILENAME="Smeltz2022")
   
# repeat these bits in case a markov chain crashes and we need to restart:
library(invitroTKstats)
setwd("c:/users/jwambaug/git/invitroTKstats/working/")

level4 <- calc_clint(FILENAME="Smeltz2022")   