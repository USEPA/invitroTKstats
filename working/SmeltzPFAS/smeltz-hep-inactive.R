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
    print(this.compound)
  }
  smeltz.hep[this.row, "Compound"] <- this.compound
  smeltz.hep[this.row, "DTXSID"] <- this.dtxsid
  smeltz.hep[this.row, "ISTD.Name"] <- this.istdname
  this.row <- this.row + 1
}
# Get rid of blank rows:
smeltz.hep <- subset(smeltz.hep,!is.na(smeltz.hep[,"Sample Text"]))

subset(smeltz.hep,DTXSID=="DTXSID5061954")


# Get rid of internal standard data:
smeltz.hep <- subset(smeltz.hep, Compound != "ISTD")
# Indicate whether hepatocytes have been inactivated:
smeltz.hep$Active.Hep <- NA




# Use invitroTKstats annotation of type:
smeltz.hep <- subset(smeltz.hep,!is.na(Type))
smeltz.hep[smeltz.hep$Type == "Analyte", "Type"] <- "Cvst"
smeltz.hep[smeltz.hep$Type == "Standard", "Type"] <- "CC"


# Were the hepatocytes alive?:
smeltz.hep[smeltz.hep[,"Type"]=="Cvst", "Active.Hep"] <- 1
smeltz.hep[regexpr("living",tolower(smeltz.hep[,"Sample Text"]))!=-1,
  "Active.Hep"] <- 1
smeltz.hep[regexpr("inactive",tolower(smeltz.hep[,"Sample Text"]))!=-1,
  "Active.Hep"] <- 0
smeltz.hep[
  sapply(smeltz.hep$Active.Hep,function(x) ifelse(is.na(x),FALSE,x==0)), 
  "Type"] <- "Inactive"
  
# Add time of sample:
smeltz.hep <- subset(smeltz.hep,!is.na(smeltz.hep[,"Sample Text"]))
smeltz.hep[,"Time"] <- NA
smeltz.hep[regexpr("t240",tolower(smeltz.hep[,"Sample Text"]))!=-1,"Time"] <- 240/60
smeltz.hep[regexpr("t120",tolower(smeltz.hep[,"Sample Text"]))!=-1,"Time"] <- 120/60
smeltz.hep[regexpr("t60",tolower(smeltz.hep[,"Sample Text"]))!=-1,"Time"] <- 60/60
smeltz.hep[regexpr("t30",tolower(smeltz.hep[,"Sample Text"]))!=-1,"Time"] <- 30/60
smeltz.hep[regexpr("t15",tolower(smeltz.hep[,"Sample Text"]))!=-1,"Time"] <- 15/60
smeltz.hep[regexpr("t0",tolower(smeltz.hep[,"Sample Text"]))!=-1,"Time"] <- 0/60

# Don't want media-only samples:
smeltz.hep[regexpr("wem",tolower(smeltz.hep[,"Sample Text"]))!=-1,"Time"] <- NA
smeltz.hep <- subset(smeltz.hep,!is.na(Time) | Type != "Cvst")

# Set the dilution factors:
# The dilution factor is indeed different across different sample types.
# -	The calibration curves (all, regardless of analyte) are diluted 240x
# -	The ref cmpds (may have label “HLB”) and PFAS labeled as “WAX1” are diluted 480x
# -	PFAS marked as WAX2 are diluted 720x
# -	It seems that when you asked for a pared down sheet this information was lost. We could add it to the Analytes page. I can see it in the more comprehensive sheet she#   provided to me.

smeltz.hep$Dilution.Factor <- 240
smeltz.hep[regexpr("hlb",tolower(smeltz.hep[,"Sample Text"]))!=-1,
  "Dilution.Factor"] <- 480
smeltz.hep[regexpr("wax1",tolower(smeltz.hep[,"Sample Text"]))!=-1,
  "Dilution.Factor"] <- 480
smeltz.hep[regexpr("wax2",tolower(smeltz.hep[,"Sample Text"]))!=-1,
  "Dilution.Factor"] <- 720
smeltz.hep[regexpr("cc",tolower(smeltz.hep[,"Sample Text"]))!=-1,
  "Dilution.Factor"] <- 240
    
## For now we don't handle the inactives:
#smeltz.hep <- subset(smeltz.hep,                           
#  regexpr("inactive",tolower(smeltz.hep[,"Sample Text"]))==-1)
  
# Make sure the areas are numeric:
smeltz.hep$Area <- as.numeric(smeltz.hep$Area)
smeltz.hep[,"IS Area"] <- as.numeric(smeltz.hep[,"IS Area"])
smeltz.hep[,"Std. Conc"] <- as.numeric(smeltz.hep[,"Std. Conc"])
smeltz.hep <- subset(smeltz.hep, !is.na(Area) & 
  !is.na(smeltz.hep[,"IS Area"]))
  
# Convert nM standard concs to uM:
smeltz.hep[,"nM"] <- as.numeric(smeltz.hep[,"nM"])/1000
# Only set a standard conc for calibration curve points:
smeltz.hep[smeltz.hep[,"Type"]!="CC","nM"] <- NA
# Concentrations calculated including dilution:
smeltz.hep[,"nM"] <- smeltz.hep[,"nM"]*smeltz.hep$Dilution.Factor
  
  
level1 <- format_clint(smeltz.hep,
  FILENAME="SmeltzPFAS/SmeltzPFAS",
  sample.col ="Name",
  date.col="Acq.Date",
  compound.col="Compound",
  lab.compound.col="Compound",
  type.col="Type",
  dilution.col="Dilution.Factor",
  cal=1,
  istd.conc = 10/1000,
  istd.col= "IS Area",
  density = 0.5,
  clint.assay.conc = 1,
  std.conc.col="nM",
  time.col = "Time",
  analysis.method = "LCMS",
  analysis.instrument = "Unknown",
  analysis.parameters.col = "RT",
  note="Sample Text"
  )

level2 <- level1
level2$Verified <- "Y"
  
for (this.id in chem.ids$DTXSID)
{
  this.mix <- chem.ids[chem.ids$DTXSID==this.id,"Mix"]
  if (regexpr("WAX1",this.mix)!=-1)
  {
    level2[level2$DTXSID==this.id & 
      (level2$Sample.Type %in% c("Cvst", "Inactive")) &
      regexpr("WAX2",level2$Note)!=-1,"Verified"] <- "Wrong Mix"
  } else if (regexpr("WAX2",this.mix)!=-1)
  {
    level2[level2$DTXSID==this.id & 
      (level2$Sample.Type %in% c("Cvst", "Inactive")) &
      regexpr("WAX1",level2$Note)!=-1,"Verified"] <- "Wrong Mix"
  }
}

level2[level2$Sample.Type=="CC" & is.na(level2$Std.Conc),"Verified"] <- 
  "Unknown Conc."
  

  
write.table(level2,
  file="SmeltzPFAS/SmeltzPFAS-Clint-Level2.tsv",
  sep="\t",
  row.names=F,
  quote=F)

for (this.id in unique(level2$DTXSID))
{
  this.subset <- subset(level2,DTXSID==this.id)
  plot(this.subset$Time,this.subset$Area,main=unique(this.subset$Compound.Name))
}
this.subset <- subset(level2,DTXSID==unique(level2$DTXSID)[3])
plot(this.subset$Time,this.subset$Area,main=unique(this.subset$Compound.Name))
# If units are right then we should have the same calibration for CC and Cvst:
subset(this.subset,Sample.Type=="CC")$Response/
  subset(this.subset,Sample.Type=="CC")$Std.Conc*
  subset(this.subset,Sample.Type=="CC")$Dilution.Factor
subset(this.subset,Time==0)$Response/1*subset(this.subset,Time==0)$Dilution.Factor



level3 <- calc_clint_point(FILENAME="SmeltzPFAS/SmeltzPFAS")
   
# repeat these bits in case a markov chain crashes and we need to restart:
library(invitroTKstats)
setwd("c:/users/jwambaug/git/invitroTKstats/working/")

level4 <- calc_clint(FILENAME="SmeltzPFAS/SmeltzPFAS")  
 