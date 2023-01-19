library(invitroTKstats)
library(readxl)

setwd("c:/users/jwambaug/git/invitroTKstats/working/KreutzPFAS")

chem.ids <- read_excel("GCMS-PFAS_Data_Summary-toJFW.xlsx", sheet=1)
chem.ids <- as.data.frame(chem.ids)
# Deal with "or" statements in chemical ids:
for (this.row in 1:dim(chem.ids)[1])
{
  if (!is.na(chem.ids[this.row, "Abbrev_SPID"]))
    if (regexpr("or", chem.ids[this.row,"Abbrev_SPID"])!=-1)
    {
      ids <- strsplit(chem.ids[this.row,"Abbrev_SPID"]," or ")
      if (length(ids[[1]])>1)
      for (this.id in ids[[1]])
      {
        this.id <- gsub(";","",this.id)
        this.id <- gsub("\\?\\)","",this.id)
        new.row <- chem.ids[this.row,]
        new.row[,"Abbrev_SPID"] <- this.id
        chem.ids <- rbind(chem.ids,new.row) 
      }
      ids <- strsplit(chem.ids[this.row,"Abbrev_SPID"]," \\(or ")
      if (length(ids[[1]])>1)
      for (this.id in ids[[1]])
      {
        this.id <- gsub(";","",this.id)
        this.id <- gsub("\\?\\)","",this.id)
        new.row <- chem.ids[this.row,]
        new.row[,"Abbrev_SPID"] <- this.id
        chem.ids <- rbind(chem.ids,new.row) 
      }
    } 
}
new.row <- chem.ids[1,]
new.row[,] <- NA
new.row[,"DTXSID"] <- "DTXSID1023869"
new.row[,2]<- "Ametryn"
new.row[,4]<- "Ametryn"
chem.ids <- rbind(chem.ids,new.row)

data.guide <- as.data.frame(read_excel("dataguide-AK-hep.xlsx"))

ak.hep <- merge_level0(data.label="KreutzPFASHep",
# describe the MS data table:
             level0.catalog=data.guide,
             sample.colname="Name",
             type.colname="Type",
             istd.col="ISTD.Name",
             analysis.param.colname.col="Note.ColName",
# describe the chemical ID table:
             chem.ids=chem.ids,
             chem.lab.id.col="Abbrev_SPID",
             chem.name.col="Chemical Name (Common Abbreviation)")




# Use invitroTKstats annotation of type:
ak.hep <- subset(ak.hep,!is.na(Type))
ak.hep[ak.hep$Type == "Cal", "Type"] <- "CC"
ak.hep[ak.hep$Type == "MatrixBlank", "Type"] <- "Blank"

# Figure out from sample name whether the hepatocytes were alive:
# Indicate whether hepatocytes have been inactivated:
ak.hep$Active.Hep <- NA

ak.hep[ak.hep[,"Type"]=="Cvst", "Active.Hep"] <- 1
ak.hep[regexpr("HITC",tolower(ak.hep[,"Sample"]))!=-1,
  "Active.Hep"] <- 0
ak.hep[
  sapply(ak.hep$Active.Hep,function(x) ifelse(is.na(x),FALSE,x==0)), 
  "Type"] <- "Inactive"
  
# Add time of sample:
ak.hep <- subset(ak.hep,!is.na(ak.hep[,"Sample"]))
ak.hep[,"Time"] <- NA
ak.hep[regexpr("t240",tolower(ak.hep[,"Sample"]))!=-1,"Time"] <- 240/60
ak.hep[regexpr("t120",tolower(ak.hep[,"Sample"]))!=-1,"Time"] <- 120/60
ak.hep[regexpr("t60",tolower(ak.hep[,"Sample"]))!=-1,"Time"] <- 60/60
ak.hep[regexpr("t30",tolower(ak.hep[,"Sample"]))!=-1,"Time"] <- 30/60
ak.hep[regexpr("t15",tolower(ak.hep[,"Sample"]))!=-1,"Time"] <- 15/60
ak.hep[regexpr("t0",tolower(ak.hep[,"Sample"]))!=-1,"Time"] <- 0/60
ak.hep[!is.na(ak.hep$Time),"Type"]  <- "Cvst"

# Set the dilution factors:
# The dilution factor is indeed different across different sample types.
# -	The calibration curves (all, regardless of analyte) are diluted 240x
# -	The ref cmpds (may have label “HLB”) and PFAS labeled as “WAX1” are diluted 480x
# -	PFAS marked as WAX2 are diluted 720x
# -	It seems that when you asked for a pared down sheet this information was lost. We could add it to the Analytes page. I can see it in the more comprehensive sheet she#   provided to me.

ak.hep$Dilution.Factor <- 480
ak.hep[ak.hep$Type=="CC",
  "Dilution.Factor"] <- 240
  
# Make sure the areas are numeric:
ak.hep$Peak.Area <- as.numeric(ak.hep$Peak.Area)
ak.hep[,"ISTD.Peak.Area"] <- as.numeric(ak.hep[,"ISTD.Peak.Area"])
ak.hep[,"Compound.Conc"] <- as.numeric(ak.hep[,"Compound.Conc"])
ak.hep <- subset(ak.hep, !is.na(Peak.Area) & 
  !is.na(ak.hep[,"ISTD.Peak.Area"]))
  
# Convert nM standard concs to uM:
ak.hep[,"Compound.Conc"] <- as.numeric(ak.hep[,"Compound.Conc"])/1000

# Concentrations calculated including dilution:
#ak.hep[,"nM"] <- ak.hep[,"nM"]*ak.hep$Dilution.Factor

ak.hep$Note <-NA

# Get rid of samples with no ISTD measured:
ak.hep <- subset(ak.hep, ISTD.Peak.Area>0)
  
level1 <- format_clint(ak.hep,
  FILENAME="KreutzPFASHep",
  sample.col ="Sample",
  date.col="Date",
  compound.col="Compound",
  lab.compound.col="Lab.Compound.ID",
  type.col="Type",
  dilution.col="Dilution.Factor",
  cal.col="Date",
  istd.conc = 10/1000,
  area.col="Peak.Area",
  istd.col= "ISTD.Peak.Area",
  density = 0.5,
  clint.assay.conc = 1,
  std.conc.col="Compound.Conc",
  time.col = "Time",
  analysis.method = "UPLC-MS/MS",
  analysis.instrument = "Waters Xevo TQ-S micro (QEB0036)",
  analysis.parameters.col = "Analysis.Params"
  )

level2 <- level1
level2$Verified <- "Y"
  
level2[level2[,"Level0.File"]=="G6HC_474_3096_101621_final.xlsx" &
       level2[,"Lab.Compound.Name"]=="760", "Verified"] <- "Not present"
level2[level2[,"Level0.File"]=="G6HC_474_3096_101621_final.xlsx" &
       regexpr("474",level2[,"Lab.Sample.Name"])!=-1 &
       level2[,"Lab.Compound.Name"]!="474", "Verified"] <- "Not present"
level2[level2[,"Level0.File"]=="G6HC_474_3096_101621_final.xlsx" &
       regexpr("3096",level2[,"Lab.Sample.Name"])!=-1 &
       level2[,"Lab.Compound.Name"]!="3096", "Verified"] <- "Not present"       
  
  
write.table(level2,
  file="KreutzPFAS-Clint-Level2.tsv",
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



level3 <- calc_clint_point(FILENAME="KreutzPFAS")
   
# repeat these bits in case a markov chain crashes and we need to restart:
library(invitroTKstats)
setwd("c:/users/jwambaug/git/invitroTKstats/working/KreutzPFAS")

level4 <- calc_clint(FILENAME="KreutzPFAS",
                          NUM.CORES=8,
                          JAGS.PATH="C:/Users/jwambaug/AppData/Local/Programs/JAGS/JAGS-4.3.1/x64")  
 