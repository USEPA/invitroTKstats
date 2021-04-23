library(invitroTKstats)

# There are multiple packages for loading Excel files, but I've been using this
# one lately:
library(readxl)

# Change to the directory that has your Excel file 
# (this is where it is on my computer):
setwd("c:/users/jwambaug/git/invitroTKstats/working/")

# load the data one sheet at a time from the Excel file, we skip the first row
# so that we get more of the column names loaded:
TO1caco2 <- read_excel("HTTK2-TO1-caco2-data.xlsx")
TO1caco2$Type <- "Receiver"
TO1caco2[regexpr("Blank",TO1caco2$SampleName)!=-1,"Type"] <- "Blank"
TO1caco2[regexpr("A_B_dos",TO1caco2$SampleName)!=-1,"Type"] <- "Dosing"
TO1caco2[regexpr("A_B_don",TO1caco2$SampleName)!=-1,"Type"] <- "Donor"
TO1caco2[regexpr("B_A_rec",TO1caco2$SampleName)!=-1,"Type"] <- "Receiver"
TO1caco2[regexpr("B_A_dos",TO1caco2$SampleName)!=-1,"Type"] <- "Dosing"
TO1caco2[regexpr("B_A_don",TO1caco2$SampleName)!=-1,"Type"] <- "Donor"

TO1caco2$Direction <- "AtoB"
TO1caco2[regexpr("B_A",TO1caco2$SampleName)!=-1,"Direction"] <- "BtoA"

TO1caco2$ISTD.Area <- as.numeric(TO1caco2$ISTD.Area)
TO1caco2$Vol.Basal <- 0.25
TO1caco2$Vol.Apical <- 0.075
TO1caco2$Date <- "June2020-May2021"

TO1caco2$ISTD.Name <- "Diclofenac and Bucetin"
TO1caco2$ISTD.Conc <- 1

TO1caco2$Test.Target.Conc <- 10

TO1caco2$Test.Target.Conc <- 10



for (this.compound $ unique(TO1caco2$CompoundName))
TO1caco2$DQdt <- 



  
  

  