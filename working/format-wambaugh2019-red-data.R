setwd("C:/Users/jwambaug/git/invitroTKstats/working")

library(readxl)
source("../invitroTKstats/R/format_fup_red.R")
source("../invitroTKstats/R/calc_fup_red_point.R")

# read from the Excel file using library(readxl)
wambaugh2019.red <- as.data.frame(read_excel("toxsci-19-0394-File011.xlsx"))
save(wambaugh2019.red,file="../invitroTKstats/data/wambaugh2019.RData")

red <- subset(wambaugh2019.red, Protein==100)
red$Date <- "2019"
red$Sample.Type <- "Blank"
red <- subset(red,!is.na(SampleName))
red[regexpr("PBS",red$SampleName)!=-1,"Sample.Type"] <- "PBS"
red[regexpr("Plasma",red$SampleName)!=-1,"Sample.Type"] <- "Plasma"
red$Dilution.Factor <- NA
red$Dilution.Factor <- as.numeric(red$Dilution.Factor)
red[red$Sample.Type=="PBS","Dilution.Factor"] <- 2
red[red$Sample.Type=="Plasma","Dilution.Factor"] <- 5
red[regexpr("T0",red$SampleName)!=-1,"Sample.Type"] <- "T0"

red$Test.Target.Conc <- 5
red$ISTD.Name <- "Bucetin and Diclofenac"
red$ISTD.Conc <- 1
red$Series <- 1

level1 <- format_fup_red(red,
  FILENAME="Wambaugh2019",
  sample.col="SampleName",
  compound.col="Preferred.Name",
  lab.compound.col="CompoundName",
  cal.col="RawDataSet")

level2 <- level1
level2$Verified <- "Y"

# Just use first 1000 observations for speed:
write.table(level2,
  file="Wambaugh2019-PPB-RED-Level2.tsv",
  sep="\t",
  row.names=F,
  quote=F)
  
level3 <- calc_fup_red_point(FILENAME="Wambaugh2019")



  

 