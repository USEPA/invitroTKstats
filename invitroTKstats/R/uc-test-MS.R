setwd("C:/Users/jwambaug/git/invitroTKstats/invitroTKstats/R")

library(gdata)
library(parallel)
library(runjags)
source("calc_uc_fup.R")
source("plot_uc_results.R")

CC.DILUTE <- 1
BLANK.DILUTE <- 1
AF.DILUTE <- 2*16
T5.DILUTE <- 5*16
T1.DILUTE <- 5*16

# read from the Excel file using library(gdata)
dat <- read.xls("20200402_PFAS_UC_PFOA_PFOS.xlsx",stringsAsFactors=F,sheet=3,skip=12)
#Create a column to hold sample type variable:
dat$Sample.Type <- ""
# All standards are sample type "CC" for calibration curve:
dat[dat$Type=="Standard","Sample.Type"] <- "CC"
# All mixed matrix blanks are treated as standards with concentration zero:
dat[regexpr("Mixed Matrix Blank",dat$Sample.Text)!=-1,"Sample.Type"] <- "CC"
dat[regexpr("Mixed Matrix Blank",dat$Sample.Text)!=-1,"Std..Conc"] <- 0
# Identify the aqueous fraction samples:
dat[regexpr("UC-CR",dat$Sample.Text)!=-1,"Sample.Type"] <- "AF"
# Identify the T1h and T5h whole samples:
dat[regexpr("UC-T1",dat$Sample.Text)!=-1,"Sample.Type"] <- "T1"
dat[regexpr("UC-T5",dat$Sample.Text)!=-1,"Sample.Type"] <- "T5"
# Identify the series (note, these are different chemicals)
dat[regexpr("-S2",dat$Sample.Text)!=-1,"Series"] <- 2
dat <- subset(dat,Sample.Type=="CC"|Series==2)
# CC is already diluted, everything else is diluted at least 16 times:
dat$Dilution.Factor <- CC.DILUTE
# Blanks are treated as part of the CC but need to be diluted:
dat[regexpr("Mixed Matrix Blank",dat$Sample.Text)!=-1,"Dilution.Factor"] <- BLANK.DILUTE
# Additional dilutions for AF and T1/T4 samples:
dat[dat[,"Sample.Type"]=="AF","Dilution.Factor"] <- AF.DILUTE
dat[dat[,"Sample.Type"]=="T5","Dilution.Factor"] <- T5.DILUTE
dat[dat[,"Sample.Type"]=="T1","Dilution.Factor"] <- T1.DILUTE
# Get rid of data that does not have a sample type:
dat <- subset(dat,Sample.Type!="")
# Get rid of non PFOA data
dat <- subset(dat,Series %in% c(NA,"2"))
dat$Compound.Name <- "PFOA"




# Give the calibration a name:
dat$Cal <- "041219"
# Set a date for the samples:
dat$Date <- "041219"
# Set info for chemical:
dat$Lab.Compound.Name <- "PFOA"
dat$DTXSID <- "DTXSID8031865"


colnames(dat)[colnames(dat)=="Name"] <- "Lab.Sample.Name"
colnames(dat)[colnames(dat)=="Std..Conc"] <- "Standard.Conc"
# concentration units (pg/uL -> ug/L -> uM):
dat$Standard.Conc <- dat$Standard.Conc / 414.07

# Recalculate the response column:
dat$ISTD.Name <-""
dat$ISTD.Conc <- 3 / 414.07

# What concentration were we trying for:
dat[dat$Sample.Type=="T1","Test.Target.Conc"] <- 10 / 414.07


# Get rid of data that for whatever reason doesn't have a response value:
dat <- subset(dat,!is.na(Response))



all.data <- dat



# read from the Excel file using library(gdata)
dat <- read.xls("20200402_PFAS_UC_PFOA_PFOS.xlsx",stringsAsFactors=F,sheet=4,skip=11)
#Create a column to hold sample type variable:
dat$Sample.Type <- ""
# All standards are sample type "CC" for calibration curve:
dat[dat$Type=="Standard","Sample.Type"] <- "CC"
# All mixed matrix blanks are treated as standards with concentration zero:
dat[regexpr("Mixed Matrix Blank",dat$Sample.Text)!=-1,"Sample.Type"] <- "CC"
dat[regexpr("Mixed Matrix Blank",dat$Sample.Text)!=-1,"Std..Conc"] <- 0
# Identify the aqueous fraction samples:
dat[regexpr("UC_UF_Mix1",dat$Sample.Text)!=-1,"Sample.Type"] <- "AF"
# Identify the T1h and T5h whole samples:
dat[regexpr("UC_T1hr_Mix1",dat$Sample.Text)!=-1,"Sample.Type"] <- "T1"
dat[regexpr("UC_T5hr_Mix1",dat$Sample.Text)!=-1,"Sample.Type"] <- "T5"
# Identify the series (note, these are different chemicals)
dat[regexpr("Mix1",dat$Sample.Text)!=-1,"Series"] <- 1
# CC is already diluted, everything else is diluted at least 16 times:
dat$Dilution.Factor <- CC.DILUTE
# Blanks are treated as part of the CC but need to be diluted:
dat[regexpr("Mixed Matrix Blank",dat$Sample.Text)!=-1,"Dilution.Factor"] <- BLANK.DILUTE
# Need to adjust the study conc column to reflect target concentration, not diluted:
dat[dat[,"Sample.Type"]=="CC","Std..Conc"] <- dat[dat[,"Sample.Type"]=="CC","Std..Conc"]
# Additional dilutions for AF and T1/T4 samples:
dat[dat[,"Sample.Type"]=="AF","Dilution.Factor"] <- AF.DILUTE
dat[dat[,"Sample.Type"]=="T5","Dilution.Factor"] <- T5.DILUTE
dat[dat[,"Sample.Type"]=="T1","Dilution.Factor"] <- T1.DILUTE
# Get rid of data that does not have a sample type:
dat <- subset(dat,Sample.Type!="")
# Get rid of non-PFOA data:
dat <- subset(dat,Series %in% c(NA,"1"))
dat$Compound.Name <- "PFOA"

# Give the calibration a name:
dat$Cal <- "100119"
# Set a date for the samples:
dat$Date <- "100119"

# Set info for chemical:
dat$Lab.Compound.Name <- "PFOA"
dat$DTXSID <- "DTXSID8031865"

colnames(dat)[colnames(dat)=="Name"] <- "Lab.Sample.Name"
colnames(dat)[colnames(dat)=="Std..Conc"] <- "Standard.Conc"
# adjust to actual molecular weight:
dat$Standard.Conc <- dat$Standard.Conc * 305 * 414.07

# Information on the internal standard:
dat$ISTD.Name <-""
dat$ISTD.Conc <- .01

# What concentration were we trying for:
dat[dat$Sample.Type=="T1","Test.Target.Conc"] <- 10 / 414.07


# Get rid of data that for whatever reason doesn't have a response value:
dat <- subset(dat,!is.na(Response))

# add these data to data object:
all.data <- all.data[,colnames(all.data)[colnames(all.data)%in%colnames(dat)]]
all.data <- rbind(all.data,dat[,colnames(all.data)])






  
# read from the Excel file using library(gdata)
dat <- read.xls("20200402_PFAS_UC_PFOA_PFOS.xlsx",stringsAsFactors=F,sheet=5,skip=11)
#Create a column to hold sample type variable:
dat$Sample.Type <- ""
# All standards are sample type "CC" for calibration curve:
dat[dat$Type=="Standard","Sample.Type"] <- "CC"
# All mixed matrix blanks are treated as standards with concentration zero:
dat[regexpr("Mixed Matrix Blank",dat$Sample.Text)!=-1,"Sample.Type"] <- "CC"
dat[regexpr("Mixed Matrix Blank",dat$Sample.Text)!=-1,"Std..Conc"] <- 0
# Identify the aqueous fraction samples:
dat[regexpr("UC_UF_Mix2",dat$Sample.Text)!=-1,"Sample.Type"] <- "AF"
# Identify the T1h and T5h whole samples:
dat[regexpr("UC_T1hr_Mix2",dat$Sample.Text)!=-1,"Sample.Type"] <- "T1"
dat[regexpr("UC_T5hr_Mix2",dat$Sample.Text)!=-1,"Sample.Type"] <- "T5"
# Identify the series (note, these are different chemicals)
dat[regexpr("Mix2",dat$Sample.Text)!=-1,"Series"] <- 1
# CC is already diluted, everything else is diluted at least 16 times:
dat$Dilution.Factor <- CC.DILUTE
# Blanks are treated as part of the CC but need to be diluted:
dat[regexpr("Mixed Matrix Blank",dat$Sample.Text)!=-1,"Dilution.Factor"] <- BLANK.DILUTE
# Need to adjust the study conc column to reflect target concentration, not diluted:
dat[dat[,"Sample.Type"]=="CC","Std..Conc"] <- dat[dat[,"Sample.Type"]=="CC","Std..Conc"]
# Additional dilutions for AF and T1/T4 samples:
dat[dat[,"Sample.Type"]=="AF","Dilution.Factor"] <- AF.DILUTE
dat[dat[,"Sample.Type"]=="T5","Dilution.Factor"] <- T5.DILUTE
dat[dat[,"Sample.Type"]=="T1","Dilution.Factor"] <- T1.DILUTE
# Get rid of data that does not have a sample type:
dat <- subset(dat,Sample.Type!="")
# Get rid of non-PFOA data:
dat <- subset(dat,Series %in% c(NA,"1"))
dat$Compound.Name <- "PFOS"

# Give the calibration a name:
dat$Cal <- "072319"
# Set the date for the samples:
dat$Date <- "072319"
# Get rid of data that for whatever reason doesn't have a response value:
dat <- subset(dat,!is.na(Response))
# adjust to actual molecular weight:
dat$Nominal.Conc <- dat$Nominal.Conc / 305 * 500.13 

# Set info for chemical:
dat$Lab.Compound.Name <- "PFOS"
dat$DTXSID <- "DTXSID3031864"

colnames(dat)[colnames(dat)=="Name"] <- "Lab.Sample.Name"
colnames(dat)[colnames(dat)=="Std..Conc"] <- "Standard.Conc"
# adjust to actual molecular weight:
dat$Standard.Conc <- dat$Standard.Conc * 305 * 500.13 

# Information on the internal standard:
dat$ISTD.Name <-""
dat$ISTD.Conc <- .01

# What concentration were we trying for:
dat[dat$Sample.Type=="T1","Test.Target.Conc"] <- 10 / 500.13 

# add these data to data object:
all.data <- rbind(all.data,dat[,colnames(all.data)])

  
  
  
# read from the Excel file using library(gdata)
dat <- read.xls("20200402_PFAS_UC_PFOA_PFOS.xlsx",stringsAsFactors=F,sheet=6,skip=11)
#Create a column to hold sample type variable:
dat$Sample.Type <- ""
# All standards are sample type "CC" for calibration curve:
dat[dat$Type=="Standard","Sample.Type"] <- "CC"
# All mixed matrix blanks are treated as standards with concentration zero:
dat[regexpr("Mixed Matrix Blank",dat$Sample.Text)!=-1,"Sample.Type"] <- "CC"
dat[regexpr("Mixed Matrix Blank",dat$Sample.Text)!=-1,"Std..Conc"] <- 0
  dat[regexpr("Mixed Matrix Blank",dat$Sample.Text)!=-1,"Std..Conc"] <- dat[regexpr("Mixed Matrix Blank",dat$Sample.Text)!=-1,"Std..Conc"] 
# Identify the aqueous fraction samples:
dat[regexpr("UC_UF_Mix3",dat$Sample.Text)!=-1,"Sample.Type"] <- "AF"
# Identify the T1h and T5h whole samples:
dat[regexpr("UC_T1hr_Mix3",dat$Sample.Text)!=-1,"Sample.Type"] <- "T1"
dat[regexpr("UC_T5hr_Mix3",dat$Sample.Text)!=-1,"Sample.Type"] <- "T5"
# Identify the series (note, these are different chemicals)
dat[regexpr("Mix3",dat$Sample.Text)!=-1,"Series"] <- 1
# CC is already diluted, everything else is diluted at least 16 times:
dat$Dilution.Factor <- CC.DILUTE
# Blanks are treated as part of the CC but need to be diluted:
dat[regexpr("Mixed Matrix Blank",dat$Sample.Text)!=-1,"Dilution.Factor"] <- BLANK.DILUTE
# Need to adjust the study conc column to reflect target concentration, not diluted:
dat[dat[,"Sample.Type"]=="CC","Std..Conc"] <- dat[dat[,"Sample.Type"]=="CC","Std..Conc"]
# Additional dilutions for AF and T1/T4 samples:
dat[dat[,"Sample.Type"]=="AF","Dilution.Factor"] <- AF.DILUTE
dat[dat[,"Sample.Type"]=="T5","Dilution.Factor"] <- T5.DILUTE
dat[dat[,"Sample.Type"]=="T1","Dilution.Factor"] <- T1.DILUTE
# Get rid of data that does not have a sample type:
dat <- subset(dat,Sample.Type!="")
# Get rid of non-PFOA data:
dat <- subset(dat,Series %in% c(NA,"1"))
dat$Compound.Name <- "PFOS"
# Recalculate the response:
IS.conc <- .01
dat$Response2 <- dat$Area/dat$IS.Area*IS.conc
# Compare the two response columns:
dat[,c("Response","Response2")]




# Subset to just the columns the Bayesian analysis needs: 
dat <- dat[,c("Compound.Name","Name","Sample.Type","Series","Std..Conc","Dilution.Factor","Response2")]
# Rename the columns:
colnames(dat) <- c("Compound.Name","Sample.Name","Sample.Type","Series","Nominal.Conc","Dilution.Factor","Response")
# Give the calibration a name:
dat$Cal <- "010720"
# Get rid of data that for whatever reason doesn't have a response value:
dat <- subset(dat,!is.na(Response))
# adjust to actual molecular weight and uM:
dat$Nominal.Conc <- dat$Nominal.Conc / 305 * 500.13
# add these data to data object:
all.data <- rbind(all.data,dat)

   
out <- calc_uc_fup(all.data)
  
p1 <- plot_uc_results(all.data,out$coda,"PFOA","041219",500)
p2 <- plot_uc_results(all.data,out$coda,"PFOA","100119",500)
p3 <- plot_uc_results(all.data,out$coda,"PFOS","072319",500)
p4 <- plot_uc_results(all.data,out$coda,"PFOS","010720",500,quad.cal=c(avar=0.105423,bvar=6.38662,cvar=0.002752060))

print(p4$plog)



