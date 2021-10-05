library(invitroTKstats)

# There are multiple packages for loading Excel files, but I've been using this
# one lately:
library(readxl)

# Change to the directory that has your Excel file 
# (this is where it is on my computer):
setwd("c:/users/jwambaug/git/invitroTKstats/working/")

# load the data one sheet at a time from the Excel file, we skip the first row
# so that we get more of the column names loaded:
TO1caco2 <- read.table("HTTK2TO1-Caco2-all.txt",header=TRUE)
TO1caco2$Type <- "R2"
TO1caco2[regexpr("Blank",TO1caco2$SampleName)!=-1,"Type"] <- "Blank"
TO1caco2[regexpr("A_B_dos",TO1caco2$SampleName)!=-1,"Type"] <- "D0"
TO1caco2[regexpr("A_B_don",TO1caco2$SampleName)!=-1,"Type"] <- "D2"
#TO1caco2[regexpr("B_A_rec",TO1caco2$SampleName)!=-1,"Type"] <- "R2"
TO1caco2[regexpr("B_A_dos",TO1caco2$SampleName)!=-1,"Type"] <- "D0"
TO1caco2[regexpr("B_A_don",TO1caco2$SampleName)!=-1,"Type"] <- "D2"

TO1caco2$Dilution.Factor <- 4
TO1caco2[regexpr("Blank",TO1caco2$SampleName)!=-1,"Dilution.Factor"] <- 1

TO1caco2$Direction <- "AtoB"
TO1caco2$Vol.Receiver <- 0.25
TO1caco2$Vol.Donor <- 0.075
TO1caco2[regexpr("B_A",TO1caco2$SampleName)!=-1,"Direction"] <- "BtoA"
TO1caco2[regexpr("B_A",TO1caco2$SampleName)!=-1,"Vol.Receiver"] <- 0.075
TO1caco2[regexpr("B_A",TO1caco2$SampleName)!=-1,"Vol.Donor"] <- 0.25

TO1caco2$ISTD.Area <- as.numeric(TO1caco2$ISTD.Area)
TO1caco2$Date <- "June2020-May2021"

TO1caco2$ISTD.Name <- "Diclofenac and Bucetin"
TO1caco2$ISTD.Conc <- 1

TO1caco2$Test.Target.Conc <- 10

save(TO1caco2,file="EPACyprotex2021.RData")




  
  

  