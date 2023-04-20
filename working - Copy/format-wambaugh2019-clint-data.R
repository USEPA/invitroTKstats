setwd("C:/Users/jwambaug/git/invitroTKstats/working")

library(readxl)
library(invitroTKstats)

# read from the Excel file using library(readxl)
load("wambaugh2019.RData")
wambaugh2019.clint <- as.data.frame(read_excel("toxsci-19-0394-File012.xlsx"))
save(wambaugh2019.clint,wambaugh2019.red,file="wambaugh2019.RData")

clint <- wambaugh2019.clint
clint$Date <- "2019"
clint$Sample.Type <- "Blank"
clint$Time..mins. <- as.numeric(clint$Time..mins.)
clint[!is.na(clint$Time..mins.),"Sample.Type"] <- "Cvst"
clint$ISTD.Name <- "Bucetin, Propranolol, and Diclofenac"
clint$ISTD.Conc <- 1
clint$Dilution.Factor <- 1
clint[is.na(clint$FileName),"FileName"]<-"Wambaugh2019"
clint$Hep.Density <- 0.5


level1 <- format_clint(clint,
  FILENAME="Wambaugh2019",
  sample.col="Sample.Name",
  compound.col="Preferred.Name",
  lab.compound.col="Name",
  time.col="Time..mins.",
  cal.col="FileName")

level2 <- level1
level2$Verified <- "Y"

# All data (allows test for saturation):
write.table(level2,
  file="Wambaugh2019-Clint-Level2.tsv",
  sep="\t",
  row.names=F,
  quote=F)

level3 <- calc_clint_point(FILENAME="Wambaugh2019")
 
# Just 1 uM data:
write.table(subset(level2,Conc==1),
  file="Wambaugh2019-1-Clint-Level2.tsv",
  sep="\t",
  row.names=F,
  quote=F)

level3.1 <- calc_clint_point(FILENAME="Wambaugh2019-1")

# Just 10 uM data:
write.table(subset(level2,Conc==10),
  file="Wambaugh2019-10-Clint-Level2.tsv",
  sep="\t",
  row.names=F,
  quote=F)
  
level3.10 <- calc_clint_point(FILENAME="Wambaugh2019-10")
  
 