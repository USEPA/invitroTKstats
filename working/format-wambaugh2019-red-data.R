setwd("C:/Users/jwambaug/git/invitroTKstats/working")

library(readxl)
source("../invitroTKstats/R/format_fup_red.R")

# read from the Excel file using library(readxl)
red <- read_excel("toxsci-19-0394-File011.xlsx")
red <- as.data.frame(red)


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

wambaugh2019.red <- format_fup_red(red,
  FILENAME="Wambaugh2019",
  sample.col="SampleName",
  compound.col="Preferred.Name",
  lab.compound.col="CompoundName",
  cal.col="RawDataSet")
  
