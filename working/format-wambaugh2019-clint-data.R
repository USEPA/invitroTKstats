setwd("C:/Users/jwambaug/git/invitroTKstats/working")

library(readxl)
library(invitroTKstats)
           
# read from the Excel file using library(readxl)
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

# Just use first 1000 observations for speed:
write.table(level2,
  file="Wambaugh2019-PPB-RED-Level2.tsv",
  sep="\t",
  row.names=F,
  quote=F)

level3 <- calc_clint_point(FILENAME="Wambaugh2019")



  

 