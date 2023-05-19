library(invitroTKstats)
library(readxl)

# clear mermory:
rm(list=ls())

setwd("c:/users/jwambaug/git/invitroTKstats/working/CrizerPFAS")

chem.ids <- as.data.frame(read_excel("CrizerChemIDs.xlsx"))

data.guide <- as.data.frame(read_excel("dataguide-DC-hep.xlsx"))

dc.hep <- merge_level0(data.label="CrizerPFASHep",
# describe the MS data table:
             level0.catalog=data.guide,
             type.colname.col="Type.ColName",
             istd.col="ISTD.Name",
             analysis.param.colname.col="RT.ColName",
             additional.colname.cols="Note.ColName",
# describe the chemical ID table:
             chem.ids=chem.ids,
             chem.lab.id.col="Chemical",
             chem.name.col="Chemical")

# Set reasonable sig figs on retention time:
dc.hep$Analysis.Params <- signif(as.numeric(dc.hep$Analysis.Params), 6)

# Use invitroTKstats annotation of type:
dc.hep <- subset(dc.hep,!is.na(Type))
dc.hep[regexpr("std",tolower(dc.hep[,"Type"]))!=-1, "Sample.Type"] <- "CC"
dc.hep[(regexpr("blank", tolower(dc.hep[,"Sample"])) != -1 &
       regexpr("media", tolower(dc.hep[,"Sample"])) != -1), "Sample.Type"] <- 
       "Blank"
dc.hep[regexpr("unknown",tolower(dc.hep[,"Type"]))!=-1, 
               "Sample.Type"] <- "Cvst"
dc.hep[regexpr("no_cells",tolower(dc.hep[,"Sample"]))!=-1, 
               "Sample.Type"] <- "Inactive"
               
# Figure out from sample name whether the hepatocytes were alive:
# Indicate whether hepatocytes have been inactivated:
dc.hep$Active.Hep <- NA

 
# Add time of sample:
dc.hep <- subset(dc.hep,!is.na(dc.hep[,"Sample"]))
dc.hep[,"Time"] <- NA
dc.hep[regexpr("no_cells",tolower(dc.hep[,"Sample"]))!=-1,"Time"] <- 240/60
# Do this time first since zeros show up in other names:
dc.hep[regexpr("0min",tolower(dc.hep[,"Sample"]))!=-1,"Time"] <- 0/60
dc.hep[regexpr("240min",tolower(dc.hep[,"Sample"]))!=-1,"Time"] <- 240/60
dc.hep[regexpr("120min",tolower(dc.hep[,"Sample"]))!=-1,"Time"] <- 120/60
dc.hep[regexpr("90min",tolower(dc.hep[,"Sample"]))!=-1,"Time"] <- 90/60
dc.hep[regexpr("60min",tolower(dc.hep[,"Sample"]))!=-1,"Time"] <- 60/60
dc.hep[regexpr("30min",tolower(dc.hep[,"Sample"]))!=-1,"Time"] <- 30/60
dc.hep[regexpr("15min",tolower(dc.hep[,"Sample"]))!=-1,"Time"] <- 15/60

# Differentiate live hepatocytes from inactive:
dc.hep[dc.hep[,"Sample.Type"] %in% "Cvst", "Active.Hep"] <- 1
dc.hep[dc.hep[,"Sample.Type"] %in% "Inactive", "Active.Hep"] <- 0

# Set the dilution factors:
# The dilution factor is indeed different across different sample types.
# -	The calibration curves (all, regardless of analyte) are diluted 240x
# -	The ref cmpds (may have label “HLB”) and PFAS labeled as “WAX1” are diluted 480x
# -	PFAS marked as WAX2 are diluted 720x
# -	It seems that when you asked for a pared down sheet this information was lost. We could add it to the Analytes page. I can see it in the more comprehensive sheet she#   provided to me.

# Se the dilution facor for the various samples:
dc.hep$Dilution.Factor <- 1
dc.hep[dc.hep[,"Sample.Type"] %in% c("Cvst","Inactive"), "Dilution.Factor"] <- 4

# Make sure the areas are numeric:
dc.hep$Peak.Area <- as.numeric(dc.hep$Peak.Area)
dc.hep[,"Compound.Conc"] <- as.numeric(dc.hep[,"Compound.Conc"])

# No ISTD:
dc.hep[,"ISTD.Peak.Area"] <- 1
  
# Convert nM standard concs to uM:
dc.hep[,"Compound.Conc"] <- as.numeric(dc.hep[,"Compound.Conc"])/1000

level1 <- format_clint(dc.hep,
  FILENAME="CrizerPFASHep",
  sample.col ="Sample",
  date.col="Date",
  compound.col="Compound",
  lab.compound.col="Lab.Compound.ID",
  type.col="Sample.Type",
  dilution.col="Dilution.Factor",
  cal.col="Date",
  istd.conc = 10/1000,
  area.col="Peak.Area",
  istd.col= "ISTD.Peak.Area",
  density = 0.5,
  clint.assay.conc = 1,
  note.col="Peak Status",
  std.conc.col="Compound.Conc",
  time.col = "Time",
  analysis.method = "UPLC-MS/MS",
  analysis.instrument = "Waters Xevo TQ-S micro (QEB0036)",
  analysis.parameters.col = "Analysis.Params"
  )

level2 <- level1
level2$Verified <- "Y"
  
level2[sapply(level2[,"Note"], function(x) "Excluded" %in% x), 
              "Verified"] <- "Excluded"

write.table(level2,
  file="CrizerPFAS-Clint-Level2.tsv",
  sep="\t",
  row.names=F,
  quote=F)

for (this.id in unique(level2$DTXSID))
{
  this.subset <- subset(level2,DTXSID==this.id)
  plot(this.subset$Time,this.subset$Area,main=unique(this.subset$Compound.Name))
}



level3 <- calc_clint_point(FILENAME="CrizerPFAS")
   
# repeat these bits in case a markov chain crashes and we need to restart:
library(invitroTKstats)
setwd("c:/users/jwambaug/git/invitroTKstats/working/CrizerPFAS")

level4 <- calc_clint(FILENAME="CrizerPFAS",
                          NUM.CORES=8,
                          JAGS.PATH="C:/Users/jwambaug/AppData/Local/JAGS/JAGS-4.3.0/x64")  
 