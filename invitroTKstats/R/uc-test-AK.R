setwd("C:/Users/jwambaug/git/invitroTKstats/invitroTKstats/R")

library(readxl)
library(dplyr)
library(parallel)
library(runjags)
source("calc_uc_fup.R")
source("plot_uc_results.R")

CC.DILUTE <- 1
BLANK.DILUTE <- 1
AF.DILUTE <- 2
T5.DILUTE <- 10
T1.DILUTE <- 10
T1.TARGET <- 10000

# read from the Excel file using library(readxl)
chems <- read_excel("120619_PFAS_PPB_Amides_UC_AK_072020.xlsx",sheet=2,skip=3)[1:6,1:3]

# info on what concentration each calibration sample was supposed to be:
sampleinfo <-  read_excel("120619_PFAS_PPB_Amides_UC_AK_072020.xlsx",sheet=8,skip=3)
sampleinfo <- sampleinfo[,c(3,8)]
sampleinfo <- subset(sampleinfo,!is.na(Name))
sampleinfo <- subset(sampleinfo,!duplicated(Name))


dat1 <- read_excel("120619_PFAS_PPB_Amides_UC_AK_072020.xlsx",sheet=12,skip=1)
dat2 <- read_excel("120619_PFAS_PPB_Amides_UC_AK_072020.xlsx",sheet=14,skip=1)
dat3 <- read_excel("120619_PFAS_PPB_Amides_UC_AK_072020.xlsx",sheet=15,skip=1)

# Set the date:
dat1$Date <- "102919"
dat2$Date <- "110119"
dat3$Date <- "121019"

#Label the samples:
dat1$Lab.Sample.Name <- dat1$Name
dat2$Lab.Sample.Name <- dat2$Name
dat3$Lab.Sample.Name <- dat3$Name

#Create a column to hold sample type variable:
dat1$Sample.Type <- ""
dat2$Sample.Type <- ""
dat3$Sample.Type <- ""


#Create a column to hold the concentration for standards:
dat1$Standard.Conc <- NaN
dat2$Standard.Conc <- NaN
dat3$Standard.Conc <- NaN
 
# All standards are sample type "Cal" for calibration curve:
dat1[dat1$Type=="Cal","Sample.Type"] <- "CC"
dat2[dat2$Type=="Cal","Sample.Type"] <- "CC"
dat3[dat3$Type=="Cal","Sample.Type"] <- "CC"

# Set the dilution factor:
dat1[dat1$Type=="Cal","Dilution.Factor"] <- CC.DILUTE
dat2[dat2$Type=="Cal","Dilution.Factor"] <- CC.DILUTE
dat3[dat3$Type=="Cal","Dilution.Factor"] <- CC.DILUTE

# Everything in each data set uses the same calibration:
dat1$Cal <- "102919"
dat2$Cal <- "110119"
dat3$Cal <- "121019"

# Add concentrations of standards
for (this.sample in unlist(dat1[dat1$Sample.Type=="CC",3]))
{
  dat1[dat1[,3]==this.sample,"Standard.Conc"] <- 
    as.numeric(sampleinfo[sampleinfo[,1]==this.sample,2])
}
for (this.sample in unlist(dat2[dat2$Sample.Type=="CC",3]))
{
  dat2[dat2[,3]==this.sample,"Standard.Conc"] <- 
    as.numeric(sampleinfo[sampleinfo[,1]==this.sample,2])
}
for (this.sample in unlist(dat3[dat3$Sample.Type=="CC",3]))
{
  dat3[dat3[,3]==this.sample,"Standard.Conc"] <- 
    as.numeric(sampleinfo[sampleinfo[,1]==this.sample,2])
}

# All blanks are treated as standards with concentration zero:
# Set type:
dat1[dat1$Type=="Blank","Sample.Type"] <- "CC"
dat2[dat2$Type=="Blank","Sample.Type"] <- "CC"
dat3[dat3$Type=="Blank","Sample.Type"] <- "CC"
# Set conc:
dat1[dat1$Type=="Blank","Standard.Conc"] <- 0
dat2[dat2$Type=="Blank","Standard.Conc"] <- 0
dat3[dat3$Type=="Blank","Standard.Conc"] <- 0
# Set dilution factor:
dat1[dat1$Type=="Blank","Dilution.Factor"] <- BLANK.DILUTE
dat2[dat2$Type=="Blank","Dilution.Factor"] <- BLANK.DILUTE
dat3[dat3$Type=="Blank","Dilution.Factor"] <- BLANK.DILUTE

# Identify the aqueous fraction samples:
# Set type:
dat1[regexpr("UCG1AF",unlist(dat1[,"Lab.Sample.Name"]))!=-1,"Sample.Type"] <- "AF"
dat2[regexpr("UCG1 AF",unlist(dat2[,"Lab.Sample.Name"]))!=-1,"Sample.Type"] <- "AF"
dat3[regexpr("AF",unlist(dat3[,"Lab.Sample.Name"]))!=-1,"Sample.Type"] <- "AF"

# Set dilution factor:
dat1[regexpr("UCG1AF",unlist(dat1[,"Lab.Sample.Name"]))!=-1,"Dilution.Factor"] <- AF.DILUTE
dat2[regexpr("UCG1 AF",unlist(dat2[,"Lab.Sample.Name"]))!=-1,"Dilution.Factor"] <- AF.DILUTE
dat3[regexpr("AF",unlist(dat3[,"Lab.Sample.Name"]))!=-1,"Dilution.Factor"] <- AF.DILUTE

# Identify the T1h and T5h whole samples:
# Set type:
dat1[regexpr("UCG1T1h",unlist(dat1[,"Lab.Sample.Name"]))!=-1,"Sample.Type"] <- "T1"
dat2[regexpr("UCG1 T1",unlist(dat2[,"Lab.Sample.Name"]))!=-1,"Sample.Type"] <- "T1"
dat3[regexpr("T1",unlist(dat3[,"Lab.Sample.Name"]))!=-1,"Sample.Type"] <- "T1"
# Set dilution factor:
dat1[regexpr("UCG1T1h",unlist(dat1[,"Lab.Sample.Name"]))!=-1,"Dilution.Factor"] <- T1.DILUTE
dat2[regexpr("UCG1 T1",unlist(dat2[,"Lab.Sample.Name"]))!=-1,"Dilution.Factor"] <- T1.DILUTE
dat3[regexpr("T1",unlist(dat3[,"Lab.Sample.Name"]))!=-1,"Dilution.Factor"] <- T1.DILUTE
# Set type:
dat1[regexpr("UCG1T5h",unlist(dat1[,"Lab.Sample.Name"]))!=-1,"Sample.Type"] <- "T5"
dat2[regexpr("UCG1 T5",unlist(dat2[,"Lab.Sample.Name"]))!=-1,"Sample.Type"] <- "T5"
dat3[regexpr("T5",unlist(dat3[,"Lab.Sample.Name"]))!=-1,"Sample.Type"] <- "T5"
# Set dilution factor:
dat1[regexpr("UCG1T5h",unlist(dat1[,"Lab.Sample.Name"]))!=-1,"Dilution.Factor"] <- T5.DILUTE
dat2[regexpr("UCG1 T5",unlist(dat2[,"Lab.Sample.Name"]))!=-1,"Dilution.Factor"] <- T5.DILUTE
dat3[regexpr("T5",unlist(dat3[,"Lab.Sample.Name"]))!=-1,"Dilution.Factor"] <- T5.DILUTE

# Set the target concentration:
dat1[dat1$Sample.Type=="T1","Test.Target.Conc"] <- T1.TARGET
dat2[dat2$Sample.Type=="T1","Test.Target.Conc"] <- T1.TARGET
dat3[dat3$Sample.Type=="T1","Test.Target.Conc"] <- T1.TARGET

# Identify the sets:
dat1$Set <- 0
dat1[regexpr("S1",unlist(dat1[,"Lab.Sample.Name"]))!=-1,"Set"] <- 1
dat1[regexpr("S2",unlist(dat1[,"Lab.Sample.Name"]))!=-1,"Set"] <- 2
dat1[regexpr("S3",unlist(dat1[,"Lab.Sample.Name"]))!=-1,"Set"] <- 3
dat2$Set <- 0
dat2[regexpr("S1",unlist(dat2[,"Lab.Sample.Name"]))!=-1,"Set"] <- 1
dat2[regexpr("S2",unlist(dat2[,"Lab.Sample.Name"]))!=-1,"Set"] <- 2
dat2[regexpr("S3",unlist(dat2[,"Lab.Sample.Name"]))!=-1,"Set"] <- 3



# Identify the internal standard used:
dat1$ISTD.Name <- "M8FOSA" 
dat1$ISTD.Conc <- 12000
# Be sure to check this!!
dat1$ISTD.Area <- unlist(dat1[,"Area...16"])
dat1 <- subset(dat1,!is.na(ISTD.Area))
dat1 <- subset(dat1,ISTD.Area>0)
dat2$ISTD.Name <- "M8FOSA" 
dat2$ISTD.Conc <- 12000
# BE SURE TO CHECK THIS:
dat2$ISTD.Area <- unlist(dat2[,"Area...16"])
dat2 <- subset(dat2,!is.na(ISTD.Area))
dat2 <- subset(dat2,ISTD.Area>0)
dat3$ISTD.Name <- "M8FOSA" 
dat3$ISTD.Conc <- 12000
# BE SURE TO CHECK THIS:
dat3$ISTD.Area <- unlist(dat3[,"Area...16"])
dat3 <- subset(dat3,!is.na(ISTD.Area))
dat3 <- subset(dat3,ISTD.Area>0)

# Now annotate data by chemical:
# 102919 has data on all five
# 110119 has data only on 916 and 923
# 120619 has data on 3117, 916, and 923

# Chemical 908:
WHICH.CHEM <- 1
dat1.chem1 <- dat1
dat1.chem1$DTXSID<- unlist(chems[WHICH.CHEM,"DTXSID"])
dat1.chem1$Name<- unlist(chems[WHICH.CHEM,"Name"])
dat1.chem1$Lab.Compound.Name <- unlist(chems[WHICH.CHEM,"Sample ID"])
# BE SURE TO CHECK THIS:
dat1.chem1$Chem.Area <- unlist(dat1[,"Area...14"])
#rownames(dat1.chem1) <- paste(chems[1,"DTXSID"],rownames(dat1.chem1),sep=".")


# Chemical 909:
WHICH.CHEM <- 2
dat1.chem2 <- dat1
dat1.chem2$DTXSID<- unlist(chems[WHICH.CHEM,"DTXSID"])
dat1.chem2$Name<- unlist(chems[WHICH.CHEM,"Name"])
dat1.chem2$Lab.Compound.Name <- unlist(chems[WHICH.CHEM,"Sample ID"])
# BE SURE TO CHECK THIS:
dat1.chem2$Chem.Area <- unlist(dat1[,"Area...22"])
#rownames(dat1.chem2) <- paste(chems[2,"DTXSID"],rownames(dat1.chem2),sep=".")


# Chemical 916
WHICH.CHEM <- 3
dat1.chem3 <- dat1
dat1.chem3$DTXSID<- unlist(chems[WHICH.CHEM,"DTXSID"])
dat1.chem3$Name<- unlist(chems[WHICH.CHEM,"Name"])
dat1.chem3$Lab.Compound.Name <- unlist(chems[WHICH.CHEM,"Sample ID"])
# BE SURE TO CHECK THIS:
dat1.chem3$Chem.Area <- unlist(dat1[,"Area...38"])
#rownames(dat1.chem3) <- paste(chems[3,"DTXSID"],rownames(dat1.chem3),sep=".")

dat2.chem3 <- dat2
dat2.chem3$DTXSID<- unlist(chems[WHICH.CHEM,"DTXSID"])
dat2.chem3$Name<- unlist(chems[WHICH.CHEM,"Name"])
dat2.chem3$Lab.Compound.Name <- unlist(chems[WHICH.CHEM,"Sample ID"])
# BE SURE TO CHECK THIS:
dat2.chem3$Chem.Area <- unlist(dat2[,"Area...46"])

dat3.chem3 <- dat3
dat3.chem3$DTXSID<- unlist(chems[WHICH.CHEM,"DTXSID"])
dat3.chem3$Name<- unlist(chems[WHICH.CHEM,"Name"])
dat3.chem3$Lab.Compound.Name <- unlist(chems[WHICH.CHEM,"Sample ID"])
# BE SURE TO CHECK THIS:
dat3.chem3$Chem.Area <- unlist(dat3[,"Area...30"])
dat3.chem3$Set <- 0
dat3.chem3[dat3.chem3$Sample.Type %in% c("T1","T5","AF"),"Set"] <- 2

# Chemical 923:
WHICH.CHEM <- 4
dat1.chem4 <- dat1
dat1.chem4$DTXSID<- unlist(chems[WHICH.CHEM,"DTXSID"])
dat1.chem4$Name<- unlist(chems[WHICH.CHEM,"Name"])
dat1.chem4$Lab.Compound.Name <- unlist(chems[WHICH.CHEM,"Sample ID"])
# BE SURE TO CHECK THIS:
dat1.chem4$Chem.Area <- unlist(dat1[,"Area...46"])
#rownames(dat1.chem4) <- paste(chems[4,"DTXSID"],rownames(dat1.chem4),sep=".")

dat2.chem4 <- dat2
dat2.chem4$DTXSID<- unlist(chems[WHICH.CHEM,"DTXSID"])
dat2.chem4$Name<- unlist(chems[WHICH.CHEM,"Name"])
dat2.chem4$Lab.Compound.Name <- unlist(chems[WHICH.CHEM,"Sample ID"])
# BE SURE TO CHECK THIS:
dat2.chem4$Chem.Area <- unlist(dat2[,"Area...54"])

dat3.chem4 <- dat3
dat3.chem4$DTXSID<- unlist(chems[WHICH.CHEM,"DTXSID"])
dat3.chem4$Name<- unlist(chems[WHICH.CHEM,"Name"])
dat3.chem4$Lab.Compound.Name <- unlist(chems[WHICH.CHEM,"Sample ID"])
# BE SURE TO CHECK THIS:
dat3.chem4$Chem.Area <- unlist(dat3[,"Area...38"])
dat3.chem4$Set <- 0
dat3.chem4[dat3.chem4$Sample.Type %in% c("T1","T5","AF"),"Set"] <- 2



# Chemical 3117:
WHICH.CHEM <- 5
dat1.chem5 <- dat1
dat1.chem5$DTXSID<- unlist(chems[WHICH.CHEM,"DTXSID"])
dat1.chem5$Name<- unlist(chems[WHICH.CHEM,"Name"])
dat1.chem5$Lab.Compound.Name <- unlist(chems[WHICH.CHEM,"Sample ID"])
# BE SURE TO CHECK THIS:
dat1.chem5$Chem.Area <- unlist(dat1[,"Area...30"])

dat3.chem5 <- dat3
dat3.chem5$DTXSID<- unlist(chems[WHICH.CHEM,"DTXSID"])
dat3.chem5$Name<- unlist(chems[WHICH.CHEM,"Name"])
dat3.chem5$Lab.Compound.Name <- unlist(chems[WHICH.CHEM,"Sample ID"])
# BE SURE TO CHECK THIS:
dat3.chem5$Chem.Area <- unlist(dat3[,"Area...14"])
dat3.chem5$Set <- 0
dat3.chem5[dat3.chem5$Sample.Type %in% c("T1","T5","AF"),"Set"] <- 3


all.data <- bind_rows(
  dat1.chem1,
  dat1.chem2,
  dat1.chem3,
  dat2.chem3,
  dat3.chem3,
  dat1.chem4,
  dat2.chem4,
  dat3.chem4,
  dat1.chem5,
  dat3.chem5)




# Only some chemicals are in each set:
all.data <- subset(all.data, Sample.Type=="CC" |
  (Set==1 & DTXSID %in% unlist(chems[1:2,"DTXSID"])) |
  (Set==3 & DTXSID %in% unlist(chems[5,"DTXSID"])) |
  (Set==2 & DTXSID %in% unlist(chems[3:4,"DTXSID"])))
  
# Some data did not pass QC:
all.data <- subset(all.data,
  sapply(all.data$Comment,
    function(x) ifelse(is.na(x),TRUE,regexpr("ropped",x)==-1)))
all.data <- subset(all.data,
  sapply(all.data$Lab.Sample.Name,
    function(x) ifelse(is.na(x),TRUE,regexpr("Cr10/11",x)==-1)))    
all.data <- subset(all.data,
  sapply(all.data$Lab.Sample.Name,
    function(x) ifelse(is.na(x),TRUE,regexpr("UCG1AFS2C_SE 10/11",x)==-1)))    
all.data <- subset(all.data,
  sapply(all.data$Lab.Sample.Name,
    function(x) ifelse(is.na(x),TRUE,regexpr("UCG1 AF S2C1009_SE 11/1",x)==-1)))    


# No replicate series in this data set:
all.data$Series <- 1
all.data[all.data$Sample.Type=="CC","Series"] <- NA



source("calc_uc_fup.R")
out <- calc_uc_fup(all.data,
  FILENAME = "AK_UC_Model_Results",
  compound.col="Name",
  area.col="Chem.Area",
  compound.conc.col="Standard.Conc")
  
  
p1 <- plot_uc_results(all.data,out$coda,"PFOA","041219",500)
p2 <- plot_uc_results(all.data,out$coda,"PFOA","100119",500)
p3 <- plot_uc_results(all.data,out$coda,"PFOS","072319",500)
p4 <- plot_uc_results(all.data,out$coda,"PFOS","010720",500,quad.cal=c(avar=0.105423,bvar=6.38662,cvar=0.002752060))

print(p4$plog)



