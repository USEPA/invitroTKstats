setwd("C:/Users/jwambaug/git/invitroTKstats/working")

library(readxl)
library(invitroTKstats)

# read from the Excel file using library(readxl)
load("wambaugh2019.RData")
wambaugh2019.methods <- as.data.frame(read_excel("toxsci-19-0394-File010.xlsx"))
save(wambaugh2019.clint,wambaugh2019.red,wambaugh2019.methods,file="wambaugh2019.RData")

# Create the sort of method column invitroTKstats is expecting from the
# Wambaugh et al. (2019) supplemental table 10:
wambaugh2019.methods$Method <- ""
# Need to convert NA's to something to allow logic to work on whole column:
wambaugh2019.methods[is.na(wambaugh2019.methods$LC),"LC"] <- ""
# Describe the chemcials with LC as such:
wambaugh2019.methods[wambaugh2019.methods$LC=="Y","Method"] <- "LC"
# Need to convert NA's to something to allow logic to work on whole column:
wambaugh2019.methods[is.na(wambaugh2019.methods$GC),"GC"] <- ""
# Describe the chemcials with GC as such:
wambaugh2019.methods[wambaugh2019.methods$GC=="Y","Method"] <- "GC"
# Remove the non GC/LC-able chemicals:

# Set instruments and notes:
# Agilent QQQ:
# Need to convert NA's to something to allow logic to work on whole column:
wambaugh2019.methods[is.na(wambaugh2019.methods$Agilent.QQQ),"Agilent.QQQ"] <-
  "Failed"
wambaugh2019.methods[wambaugh2019.methods$Agilent.QQQ!="Failed","Instrument"] <-
  "Agilent QQQ"
wambaugh2019.methods[wambaugh2019.methods$Agilent.QQQ!="Failed",
  "Analysis.Notes"] <- 
  wambaugh2019.methods[wambaugh2019.methods$Agilent.QQQ!="Failed","Agilent.QQQ"]

# Water's Xevo:
# Need to convert NA's to something to allow logic to work on whole column:
wambaugh2019.methods[is.na(wambaugh2019.methods$Water.s.Xevo),"Water.s.Xevo"] <-
  "Failed"
wambaugh2019.methods[wambaugh2019.methods$Water.s.Xevo!="Failed",
  "Instrument"] <- "Waters Xevo"
wambaugh2019.methods[wambaugh2019.methods$Water.s.Xevo!="Failed",
  "Analysis.Notes"] <- 
  wambaugh2019.methods[wambaugh2019.methods$Water.s.Xevo!="Failed",
  "Water.s.Xevo"]
  
# AB.Sciex.Qtrap:
# Need to convert NA's to something to allow logic to work on whole column:
wambaugh2019.methods[is.na(wambaugh2019.methods$AB.Sciex.Qtrap),
  "AB.Sciex.Qtrap"] <-  "Failed"
wambaugh2019.methods[wambaugh2019.methods$AB.Sciex.Qtrap!="Failed",
  "Instrument"] <- "AB Sciex QTRAP"
wambaugh2019.methods[wambaugh2019.methods$AB.Sciex.Qtrap!="Failed",
  "Analysis.Notes"] <- 
  wambaugh2019.methods[wambaugh2019.methods$AB.Sciex.Qtrap!="Failed",
  "AB.Sciex.Qtrap"]  

# Agilent.GCMS:
# Need to convert NA's to something to allow logic to work on whole column:
wambaugh2019.methods[is.na(wambaugh2019.methods$Agilent.GCMS),
  "Agilent.GCMS"] <-  "Failed"
wambaugh2019.methods[wambaugh2019.methods$Agilent.GCMS!="Failed",
  "Instrument"] <- "Agilent GCMS"
wambaugh2019.methods[wambaugh2019.methods$Agilent.GCMS!="Failed",
  "Analysis.Notes"] <- 
  wambaugh2019.methods[wambaugh2019.methods$Agilent.GCMS!="Failed",
  "Agilent.GCMS"]  

# GCTOF:
# Need to convert NA's to something to allow logic to work on whole column:
wambaugh2019.methods[is.na(wambaugh2019.methods$GCTOF),
  "GCTOF"] <-  "Failed"
wambaugh2019.methods[wambaugh2019.methods$GCTOF!="Failed",
  "Instrument"] <- "GCTOF"
wambaugh2019.methods[wambaugh2019.methods$GCTOF!="Failed",
  "Analysis.Notes"] <- 
  wambaugh2019.methods[wambaugh2019.methods$GCTOF!="Failed",
  "GCTOF"]  
  
# Add other needed columns:
wambaugh2019.methods$ISTD.Name <- "Bucetin and Diclofenac"

create_method_table(wambaugh2019.methods,
  compound.col="PREFERRED_NAME",
  analysis.method.col="Method",
  analysis.instrument.col="Instrument", 
  analysis.parameters.col="Analysis.Notes")
  
