library(invitroTKstats)
library(readxl)

setwd("c:/users/jwambaug/git/invitroTKstats/working/SmeltzPFAS")

chem.ids <- read_excel("PFAS LC-MS RED Summary 20220709.xlsx", sheet=1, skip=1)[1:29,1:2]
chem.ids <- as.data.frame(chem.ids)
 
smeltz.red <- NULL
for (this.sheet in c(3,5,7,9))
{
  smeltz.red1 <- read_excel("PFAS LC-MS RED Summary 20220709.xlsx", sheet=this.sheet, skip=6)
  smeltz.red1 <- as.data.frame(smeltz.red1)
  this.compound <- "DTXSID00379268"
  this.row <- 1
  while(this.row <= dim(smeltz.red1)[1])
  {
    if (!is.na(smeltz.red1[this.row,1]))
    {
      if (regexpr("Compound",smeltz.red1[this.row,1])!=-1)
      {
        temp <- trimws(strsplit(smeltz.red1[this.row,1],": ")[[1]][2])
        this.compound <- unique(chem.ids[regexpr(paste0("\\(",temp,"\\)"),chem.ids[,2])!=-1,1])
        if (length(this.compound)==0) this.compound <- temp
      } else 
      {  
        smeltz.red1[this.row,"DTXSID"] <- this.compound
      }
   }
   this.row <- this.row + 1
  }
  smeltz.red <- rbind(smeltz.red, smeltz.red1)
}
smeltz.red <- subset(smeltz.red,!is.na(DTXSID))


smeltz.red2 <- read_excel("PFAS LC-MS RED Summary 20220709.xlsx", sheet=5, skip=6)
smeltz.red2 <- as.data.frame(smeltz.red2)
this.compound <- "DTXSID00379268"
this.row <- 1
while(this.row <= dim(smeltz.red1)[1])
{
  if (!is.na(smeltz.red1[this.row,1]))
  {
    if (regexpr("Compound",smeltz.red1[this.row,1])!=-1)
    {
      temp <- trimws(strsplit(smeltz.red2[this.row,1],":")[[1]][2])
      this.compound <- unique(chem.ids[regexpr(temp,chem.ids[,2])!=-1,1])
      if (length(this.compound)==0) this.compound <- temp
    } else 
    {  
      smeltz.red2[this.row,"DTXSID"] <- this.compound
    }
 }
 this.row <- this.row + 1
}

smeltz.red.pfas <- read_excel("Hep12 Data for Uncertainty Feb2022.xlsx", sheet=3, skip=2)
smeltz.red.pfas <- as.data.frame(smeltz.red.pfas)
colnames(smeltz.red.pfas) <- colnames(smeltz.red.ref)

# Merge the ref chems and pfas:
smeltz.red <- rbind(smeltz.red.ref,smeltz.red.pfas)

# get rid of blank rows:
smeltz.red <- subset(smeltz.red, !is.na(smeltz.red[,1]))

# Annotate Compounds
this.compound <- "Phenacetin"
this.dtxsid <- chem.ids[
  tolower(chem.ids[,"Analyte Name"]) == tolower(this.compound),
  "DTXSID"]
this.istdname <- chem.ids[
  tolower(chem.ids[,"Analyte Name"]) == tolower(this.compound),
  "Internal Standard"]
  
this.row <- 1
while (this.row <= dim(smeltz.red)[1])
{
  if (regexpr("Compound",smeltz.red[this.row,1])!=-1)
  {
    this.compound <- strsplit(smeltz.red[this.row,1],"  ")[[1]][2] 
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
  smeltz.red[this.row, "Compound"] <- this.compound
  smeltz.red[this.row, "DTXSID"] <- this.dtxsid
  smeltz.red[this.row, "ISTD.Name"] <- this.istdname
  this.row <- this.row + 1
}

subset(smeltz.red,DTXSID=="DTXSID5061954")

# Get rid of internal standard data:
smeltz.red <- subset(smeltz.red, Compound != "ISTD")

# Use invitroTKstats annotation of type:
smeltz.red <- subset(smeltz.red,!is.na(Type))
smeltz.red[smeltz.red$Type == "Analyte", "Type"] <- "Cvst"

# Add time of sample:
smeltz.red <- subset(smeltz.red,!is.na(smeltz.red[,"Sample Text"]))
smeltz.red[,"Time"] <- NA
smeltz.red[regexpr("t240",tolower(smeltz.red[,"Sample Text"]))!=-1,"Time"] <- 240/60
smeltz.red[regexpr("t120",tolower(smeltz.red[,"Sample Text"]))!=-1,"Time"] <- 120/60
smeltz.red[regexpr("t60",tolower(smeltz.red[,"Sample Text"]))!=-1,"Time"] <- 60/60
smeltz.red[regexpr("t30",tolower(smeltz.red[,"Sample Text"]))!=-1,"Time"] <- 30/60
smeltz.red[regexpr("t15",tolower(smeltz.red[,"Sample Text"]))!=-1,"Time"] <- 15/60
smeltz.red[regexpr("t0",tolower(smeltz.red[,"Sample Text"]))!=-1,"Time"] <- 0/60

# Don't want media-only samples:
 smeltz.red[regexpr("wem",tolower(smeltz.red[,"Sample Text"]))!=-1,"Time"] <- NA
smeltz.red <- subset(smeltz.red,!is.na(Time) | Type != "Cvst")

# Indicate whether hepatocytes have been inactivated:
smeltz.red$Active.Hep <- NA

smeltz.red[regexpr("living",tolower(smeltz.red[,"Sample Text"]))!=-1,
  "Active.Hep"] <- 1
smeltz.red[regexpr("inactive",tolower(smeltz.red[,"Sample Text"]))!=-1,
  "Active.Hep"] <- 0
  
# For now we don't handle the inactives:
smeltz.red <- subset(smeltz.red,                           
  regexpr("inactive",tolower(smeltz.red[,"Sample Text"]))==-1)
  
# Make sure the areas are numeric:
smeltz.red$Area <- as.numeric(smeltz.red$Area)
smeltz.red[,"IS Area"] <- as.numeric(smeltz.red[,"IS Area"])
smeltz.red[,"Std. Conc"] <- as.numeric(smeltz.red[,"Std. Conc"])
smeltz.red <- subset(smeltz.red, !is.na(Area) & 
  !is.na(smeltz.red[,"IS Area"]))
  
  
level1 <- format_clint(smeltz.red,
  FILENAME="SmeltzPFAS/SmeltzPFAS",
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
  clint.assay.conc = 1,
  time.col = "Time",
  analysis.method = "LCMS",
  analysis.instrument = "Unknown",
  analysis.parameters = "RT"
  )

level2 <- level1
level2$Verified <- "Y"
  
write.table(level2,
  file="SmeltzPFAS/SmeltzPFAS-nocc-noinactive-Clint-Level2.tsv",
  sep="\t",
  row.names=F,
  quote=F)

for (this.id in unique(level2$DTXSID))
{
  this.subset <- subset(level2,DTXSID==this.id)
  plot(this.subset$Time,this.subset$Area,main=unique(this.subset$Compound.Name))
}

level3 <- calc_clint_point(FILENAME="SmeltzPFAS/SmeltzPFAS-nocc-noinactive")
   
# repeat these bits in case a markov chain crashes and we need to restart:
library(invitroTKstats)
setwd("c:/users/jwambaug/git/invitroTKstats/working/")

level4 <- calc_clint(FILENAME="SmeltzPFAS/SmeltzPFAS-nocc-noinactive")  
 