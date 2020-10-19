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

#Label the samples:
dat1$Sample.Name <- dat1$Name

#Create a column to hold sample type variable:
dat1$Sample.Type <- ""
#Create a column to hold the concentration for standards:
dat1$Standard.Conc <- NaN
 
# All standards are sample type "Cal" for calibration curve:
dat1[dat1$Type=="Cal","Sample.Type"] <- "CC"
# Set the dilution factor:
dat1[dat1$Type=="Cal","Dilution.Factor"] <- CC.DILUTE

# Everything in this data set using the same calibration:
dat1$Cal <- 1

# Add concentrations of standards
for (this.sample in unlist(dat1[dat1$Sample.Type=="CC",3]))
{
  dat1[dat1[,3]==this.sample,"Standard.Conc"] <- 
    as.numeric(sampleinfo[sampleinfo[,1]==this.sample,2])
}

# All blanks are treated as standards with concentration zero:
dat1[dat1$Type=="Blank","Sample.Type"] <- "CC"
dat1[dat1$Type=="Blank","Standard.Conc"] <- 0
dat1[dat1$Type=="Blank","Dilution.Factor"] <- BLANK.DILUTE

# Identify the aqueous fraction samples:
dat1[regexpr("UCG1AF",unlist(dat1[,3]))!=-1,"Sample.Type"] <- "AF"
dat1[regexpr("UCG1AF",unlist(dat1[,3]))!=-1,"Dilution.Factor"] <- AF.DILUTE
# Identify the T1h and T5h whole samples:
dat1[regexpr("UCG1T1h",unlist(dat1[,3]))!=-1,"Sample.Type"] <- "T1"
dat1[regexpr("UCG1T1h",unlist(dat1[,3]))!=-1,"Dilution.Factor"] <- T1.DILUTE
dat1[regexpr("UCG1T5h",unlist(dat1[,3]))!=-1,"Sample.Type"] <- "T5"
dat1[regexpr("UCG1T5h",unlist(dat1[,3]))!=-1,"Dilution.Factor"] <- T5.DILUTE

# Set the target concentration:
dat1[data1$Sample.Type=="T1","Test.Target.Conc"] <- T1.TARGET

# Identify the set:
dat1$Set <- 0
dat1[regexpr("S1",unlist(dat1[,3]))!=-1,"Set"] <- 1
dat1[regexpr("S2",unlist(dat1[,3]))!=-1,"Set"] <- 2
dat1[regexpr("S3",unlist(dat1[,3]))!=-1,"Set"] <- 3

# Identify the internal standard used:
dat1$ISTD.Name <- "M8FOSA" 
dat1$ISTD.Conc <- 12000
dat1$ISTD.Area <- unlist(dat1[,"Area...16"])
dat1 <- subset(dat1,!is.na(ISTD.Area))
dat1 <- subset(dat1,ISTD.Area>0)

# Now annotate data by chemical:
dat1.chem1 <- dat1
dat1.chem1$DTXSID<- unlist(chems[1,"DTXSID"])
dat1.chem1$Name<- unlist(chems[1,"Name"])
dat1.chem1$Chem.Area <- unlist(dat1[,"Area...14"])
rownames(dat1.chem1) <- paste(chems[1,"DTXSID"],rownames(dat1.chem1),sep=".")

dat1.chem2 <- dat1
dat1.chem2$DTXSID<- unlist(chems[2,"DTXSID"])
dat1.chem2$Name<- unlist(chems[2,"Name"])
dat1.chem2$Chem.Area <- unlist(dat1[,"Area...22"])
rownames(dat1.chem2) <- paste(chems[2,"DTXSID"],rownames(dat1.chem2),sep=".")

dat1.chem3 <- dat1
dat1.chem3$DTXSID<- unlist(chems[3,"DTXSID"])
dat1.chem3$Name<- unlist(chems[3,"Name"])
dat1.chem3$Chem.Area <- unlist(dat1[,"Area...38"])
rownames(dat1.chem3) <- paste(chems[3,"DTXSID"],rownames(dat1.chem3),sep=".")

dat1.chem4 <- dat1
dat1.chem4$DTXSID<- unlist(chems[4,"DTXSID"])
dat1.chem4$Name<- unlist(chems[4,"Name"])
dat1.chem4$Chem.Area <- unlist(dat1[,"Area...46"])
rownames(dat1.chem4) <- paste(chems[4,"DTXSID"],rownames(dat1.chem4),sep=".")

dat1.chem5 <- dat1
dat1.chem5$DTXSID<- unlist(chems[5,"DTXSID"])
dat1.chem5$Name<- unlist(chems[5,"Name"])
dat1.chem5$Chem.Area <- unlist(dat1[,"Area...30"])
rownames(dat1.chem5) <- paste(chems[5,"DTXSID"],rownames(dat1.chem5),sep=".")


dat1 <- bind_rows(dat1.chem1,dat1.chem2,dat1.chem3,dat1.chem4,dat1.chem5)

# Only some chemicals are in each set:
dat1 <- subset(dat1, Sample.Type=="CC" |
  (Set==1 & DTXSID %in% unlist(chems[1:2,"DTXSID"])) |
  (Set==3 & DTXSID %in% unlist(chems[5,"DTXSID"])) |
  (Set==2 & DTXSID %in% unlist(chems[3:4,"DTXSID"])))
  

# No replicate series in this data set:
dat1$Series <- 1
dat1[dat1$Sample.Type=="CC","Series"] <- NA

# Set the date:
dat1$Date <- "102919"
dat1$Lab.Sample.ID <- dat1$Name

   
out <- calc_uc_fup(dat1,
  compound.col="Name",
  area.col="Chem.Area",
  compound.conc.col="Standard.Conc")
  
  
p1 <- plot_uc_results(all.data,out$coda,"PFOA","041219",500)
p2 <- plot_uc_results(all.data,out$coda,"PFOA","100119",500)
p3 <- plot_uc_results(all.data,out$coda,"PFOS","072319",500)
p4 <- plot_uc_results(all.data,out$coda,"PFOS","010720",500,quad.cal=c(avar=0.105423,bvar=6.38662,cvar=0.002752060))

print(p4$plog)



