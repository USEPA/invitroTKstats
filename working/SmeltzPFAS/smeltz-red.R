rm(list=ls())

library(invitroTKstats)
library(readxl)

setwd("c:/users/jwambaug/git/invitroTKstats/working/SmeltzPFAS")

chem.ids <- read_excel("PFAS LC-MS RED Summary 20220709.xlsx", sheet=1, skip=1)[1:29,1:2]
chem.ids <- as.data.frame(chem.ids)
chem.ids <- subset(chem.ids, !duplicated(chem.ids[,2]))
 
this.file <- "PFAS LC-MS RED Summary 20220709.xlsx"
smeltz.red <- NULL
for (this.sheet in c(3,5,7,9))
{
  smeltz.red1 <- read_excel(this.file, sheet=this.sheet, skip=6)

  this.sheet.name <- excel_sheets(this.file)[this.sheet]
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
  smeltz.red1$File <- this.file
  smeltz.red1$Sheet <- this.sheet.name
  smeltz.red <- rbind(smeltz.red, smeltz.red1)
}
smeltz.red <- subset(smeltz.red,!is.na(DTXSID))

# Convert date and time to strings:
smeltz.red$Acq.Date <- as.Date(as.numeric(smeltz.red$Acq.Date),origin = "1899-12-30")
h <- floor(as.numeric(smeltz.red$Acq.Time)*24)
m <- as.character(floor((as.numeric(smeltz.red$Acq.Time)*24-h)*60))
m[is.na(m)] <- "00"
m[sapply(m,nchar)==1] <- paste("0",m[sapply(m,nchar)==1],sep="")
smeltz.red$Acq.Time <- paste(h,m,sep=":")


  
# Set reasonable precision:
for (this.col in c("Area", "Height", "IS Area", "RT", "%Dev",
                   "Response", "Coeff. Of Determination", "Std. Conc", "nM"))
  smeltz.red[,this.col] <- signif(as.numeric(smeltz.red[,this.col]),6)
  
# Set Compound Identities:
chem.ids$Compound <- unlist(lapply(strsplit(chem.ids[,2]," \\("),function(x) x[[1]])) 
for (this.id in unique(smeltz.red$DTXSID))
{
  if (this.id %in% chem.ids$DTXSID)
  {
    smeltz.red[smeltz.red$DTXSID==this.id,"Compound"] <- 
      chem.ids[chem.ids$DTXSID==this.id,"Compound"] 
  } 
}  
smeltz.red[smeltz.red$DTXSID=="n-Butylparaben", "Compound"] <- "n-Butylparaben"
smeltz.red[smeltz.red$DTXSID=="n-Butylparaben", "DTXSID"] <- "DTXSID3020209"
smeltz.red[is.na(smeltz.red$Compound), "Compound"] <- "ISTD"

# Match chemicals to their internal standards:
test.chems <- unique(smeltz.red$DTXSID)
test.chems <- test.chems[regexpr("DTXSID",test.chems)!=-1]
for (this.chem in test.chems)
{
  this.data <- subset(smeltz.red,!is.na(smeltz.red[,"IS Area"]) &
                      DTXSID==this.chem)
  this.istd <- smeltz.red[sapply(smeltz.red$Area,
                                 function(x) this.data[1,"IS Area"]%in%x),
                                 "DTXSID"]    
  smeltz.red[smeltz.red$DTXSID==this.chem,"ISTD"] <- this.istd              
}

# Get rid of internal standard data:
smeltz.red <- subset(smeltz.red, Compound != "ISTD")

# Use invitroTKstats annotation of type:
smeltz.red <- subset(smeltz.red,!is.na(smeltz.red[,"Sample Text"]))
smeltz.red[regexpr("CC",smeltz.red[,"Sample Text"])!=-1,
                   "Sample.Type"] <- "CC"
smeltz.red[regexpr("-Pl-",smeltz.red[,"Sample Text"])!=-1,
                   "Sample.Type"] <- "Plasma"
smeltz.red[regexpr("-S-",smeltz.red[,"Sample Text"])!=-1,
                   "Sample.Type"] <- "PBS"
smeltz.red[regexpr("-A",smeltz.red[,"Sample Text"])!=-1,
                   "Replicate"] <- "A"
smeltz.red[regexpr("-B",smeltz.red[,"Sample Text"])!=-1,
                   "Replicate"] <- "B"
smeltz.red[regexpr("-C",smeltz.red[,"Sample Text"])!=-1,
                   "Replicate"] <- "C"
smeltz.red[regexpr("-EC1-",smeltz.red[,"Sample Text"])!=-1,
                   "Sample.Type"] <- "EC1"                   
smeltz.red[regexpr("-EC2-",smeltz.red[,"Sample Text"])!=-1,
                   "Sample.Type"] <- "EC2"
smeltz.red[regexpr("/T1",smeltz.red[,"Sample Text"])!=-1,
                   "Sample.Type"] <- "T0"
smeltz.red[regexpr("/T5",smeltz.red[,"Sample Text"])!=-1,
                   "Sample.Type"] <- "Stability"
smeltz.red[regexpr("/T1",smeltz.red[,"Sample Text"])!=-1,
                   "Time"] <- 1
smeltz.red[regexpr("/T5",smeltz.red[,"Sample Text"])!=-1,
                   "Time"] <- 5
smeltz.red[regexpr("Crash Blank",smeltz.red[,"Sample Text"])!=-1,
                  "Sample.Type"] <- "NoPlasma.Blank"
smeltz.red[regexpr("Matrix Blank",smeltz.red[,"Sample Text"])!=-1,
                  "Sample.Type"] <- "Plasma.Blank"
                                     
# Smeltz Crash Blanks weren't analyzed in MS, set to zero:
smeltz.red[smeltz.red[,"Sample.Type"]%in%"NoPlasma.Blank","Area"] <- 0
smeltz.red[smeltz.red[,"Sample.Type"]%in%"NoPlasma.Blank","IS Area"] <- 1
                                     
# Restrict to samples with peak areas for analyte ans internal standard
smeltz.red[,"Std. Conc"] <- as.numeric(smeltz.red[,"Std. Conc"])
smeltz.red <- subset(smeltz.red, !is.na(Area) & 
  !is.na(smeltz.red[,"IS Area"]))
  
# Annotate dilution factors:
smeltz.red[sapply(smeltz.red[,"Sample.Type"],function(x) "CC" %in% x),"Dilution.Factor"] <- 1
smeltz.red[sapply(smeltz.red[,"Sample.Type"],function(x) "PBS" %in% x),"Dilution.Factor"] <- 2
smeltz.red[sapply(smeltz.red[,"Sample.Type"],function(x) "EC1" %in% x),"Dilution.Factor"] <- 10
smeltz.red[sapply(smeltz.red[,"Sample.Type"],function(x) "EC2" %in% x),"Dilution.Factor"] <- 10
smeltz.red[sapply(smeltz.red[,"Sample.Type"],function(x) "Plasma" %in% x),"Dilution.Factor"] <- 20
smeltz.red[sapply(smeltz.red[,"Sample.Type"],function(x) "T0" %in% x),"Dilution.Factor"] <- 10
smeltz.red[sapply(smeltz.red[,"Sample.Type"],function(x) "Stability" %in% x),"Dilution.Factor"] <- 10
smeltz.red[sapply(smeltz.red[,"Sample.Type"],function(x) "Plasma.Blank" %in% x),"Dilution.Factor"] <- 1
smeltz.red[sapply(smeltz.red[,"Sample.Type"],function(x) "NoPlasma.Blank" %in% x),"Dilution.Factor"] <- 1

# Standard concentrations differ slightly from Protocol (CC17 = CC16, etc)
# Units should be uM:
smeltz.red[,"Std. Conc"] <- as.numeric(smeltz.red[,"Std. Conc"])/100


  
level1 <- format_fup_red(smeltz.red,
  FILENAME="SmeltzPFAS",
  sample.col ="Name",
  date.col="Acq.Date",
  compound.col="Compound",
  lab.compound.col="Compound",
  type.col="Sample.Type",
  dilution.col="Dilution.Factor",
  replicate.col="Replicate",
  cal=1,
  istd.conc = 10/1000,
  istd.col= "IS Area",
  istd.name.col = "ISTD", 
  std.conc.col = "Std. Conc", 
  level0.file.col = "File", 
  level0.sheet.col = "Sheet",
  nominal.test.conc = 10,
  plasma.percent = 100,
  time.col = "Time",
  analysis.method = "LCMS",
  analysis.instrument = "Waters ACQUITY I-Class UHPLC - Xevo TQ-S uTQMS",
  analysis.parameters = "RT",
  note.col=NULL
  )
# Produce errors if smeltz.red acccessed:
rm(smeltz.red)

level2 <- level1
level2$Verified <- "Y"
  

write.table(level2,
  file="SmeltzPFAS-PPB-RED-Level2.tsv",
  sep="\t",
  row.names=F,
  quote=F)


level3 <- calc_fup_red_point(FILENAME="SmeltzPFAS")
   
# repeat these bits in case a markov chain crashes and we need to restart:
library(invitroTKstats)
rm(list=ls())
setwd("c:/users/jwambaug/git/invitroTKstats/working/SmeltzPFAS")

level4 <- calc_fup_red(FILENAME="SmeltzPFAS",
                       NUM.CORES=8,
                       JAGS.PATH="C:/Users/jwambaug/AppData/Local/JAGS/JAGS-4.3.0/x64")  
 
 