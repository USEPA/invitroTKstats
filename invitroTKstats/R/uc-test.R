setwd("C:/Users/jwambaug/git/invitroTKstats/invitroTKstats/R")

library(gdata)
library(parallel)
library(runjags)
source("calc_uc_fup.R")
source("plot_uc_results.R")

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
dat$Dilution.Factor <- 1
# Blanks are treated as part of the CC but need to be diluted:
dat[regexpr("Mixed Matrix Blank",dat$Sample.Text)!=-1,"Dilution.Factor"] <- 4*4
# Additional dilutions for AF and T1/T4 samples:
dat[dat[,"Sample.Type"]=="AF","Dilution.Factor"] <- 4*4*2
dat[dat[,"Sample.Type"]=="T5","Dilution.Factor"] <- 4*4*5
dat[dat[,"Sample.Type"]=="T1","Dilution.Factor"] <- 4*4*5
# Get rid of data that does not have a sample type:
dat <- subset(dat,Sample.Type!="")
# Get rid of non PFOA data
dat <- subset(dat,Series %in% c(NA,"2"))
dat$Compound.Name <- "PFOA"
# Recalculate the response column:
IS.conc <- 3
dat$Response2 <- dat$Area/dat$IS.Area*IS.conc
# Compare the two response columns:
dat[,c("Response","Response2")]
# Subset to just the columns the Bayesian analysis needs: 
dat <- dat[,c("Compound.Name","Name","Sample.Type","Series","Std..Conc","Dilution.Factor","Response2")]
# Rename the columns:
colnames(dat) <- c("Compound.Name","Sample.Name","Sample.Type","Series","Nominal.Conc","Dilution.Factor","Response")
# Give the calibration a name:
dat$Cal <- "041219"
# Get rid of data that for whatever reason doesn't have a response value:
dat <- subset(dat,!is.na(Response))
# concentration units (pg/uL -> ug/L -> uM):
dat$Nominal.Conc <- dat$Nominal.Conc / 414.07



# Save the data:
write.csv(dat,file="PFOA041219.csv",row.names=F)

#out <- calc_uc_fup(dat)

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
dat$Dilution.Factor <- 1
# Blanks are treated as part of the CC but need to be diluted:
dat[regexpr("Mixed Matrix Blank",dat$Sample.Text)!=-1,"Dilution.Factor"] <- 4*4
# Need to adjust the study conc column to reflect target concentration, not diluted:
dat[dat[,"Sample.Type"]=="CC","Std..Conc"] <- dat[dat[,"Sample.Type"]=="CC","Std..Conc"]
# Additional dilutions for AF and T1/T4 samples:
dat[dat[,"Sample.Type"]=="AF","Dilution.Factor"] <- 4*4*2
dat[dat[,"Sample.Type"]=="T5","Dilution.Factor"] <- 4*4*5
dat[dat[,"Sample.Type"]=="T1","Dilution.Factor"] <- 4*4*5
# Get rid of data that does not have a sample type:
dat <- subset(dat,Sample.Type!="")
# Get rid of non-PFOA data:
dat <- subset(dat,Series %in% c(NA,"1"))
dat$Compound.Name <- "PFOA"
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
dat$Cal <- "100119"
# Get rid of data that for whatever reason doesn't have a response value:
dat <- subset(dat,!is.na(Response))
# adjust to actual molecular weight:
dat$Nominal.Conc <- dat$Nominal.Conc / 305 * 414.07
# Save the data:
write.csv(dat,file="PFOA100119.csv",row.names=F)
# add these data to data object:
all.data <- rbind(all.data,dat)
  
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
dat$Dilution.Factor <- 1
# Blanks are treated as part of the CC but need to be diluted:
dat[regexpr("Mixed Matrix Blank",dat$Sample.Text)!=-1,"Dilution.Factor"] <- 4*4
# Need to adjust the study conc column to reflect target concentration, not diluted:
dat[dat[,"Sample.Type"]=="CC","Std..Conc"] <- dat[dat[,"Sample.Type"]=="CC","Std..Conc"]
# Additional dilutions for AF and T1/T4 samples:
dat[dat[,"Sample.Type"]=="AF","Dilution.Factor"] <- 4*4*2
dat[dat[,"Sample.Type"]=="T5","Dilution.Factor"] <- 4*4*5
dat[dat[,"Sample.Type"]=="T1","Dilution.Factor"] <- 4*4*5
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
dat$Cal <- "072319"
# Get rid of data that for whatever reason doesn't have a response value:
dat <- subset(dat,!is.na(Response))
# adjust to actual molecular weight:
dat$Nominal.Conc <- dat$Nominal.Conc / 305 * 500.13 
# Save the data:
write.csv(dat,file="PFOS072319.csv",row.names=F)
# add these data to data object:
all.data <- rbind(all.data,dat)
  
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
dat$Dilution.Factor <- 1
# Blanks are treated as part of the CC but need to be diluted:
dat[regexpr("Mixed Matrix Blank",dat$Sample.Text)!=-1,"Dilution.Factor"] <- 4*4
# Need to adjust the study conc column to reflect target concentration, not diluted:
dat[dat[,"Sample.Type"]=="CC","Std..Conc"] <- dat[dat[,"Sample.Type"]=="CC","Std..Conc"]
# Additional dilutions for AF and T1/T4 samples:
dat[dat[,"Sample.Type"]=="AF","Dilution.Factor"] <- 4*4*2
dat[dat[,"Sample.Type"]=="T5","Dilution.Factor"] <- 4*4*5
dat[dat[,"Sample.Type"]=="T1","Dilution.Factor"] <- 4*4*5
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
# Save the data:
write.csv(dat,file="PFOS010720.csv",row.names=F)
# add these data to data object:
all.data <- rbind(all.data,dat)

   
out <- calc_uc_fup(all.data)
  

plots <- plot_uc_results(all.data,out$coda,"PFOS","010720",500,quad.cal=c(avar=0.105423,bvar=6.38662,cvar=0.002752060))
print(plots$plinear)
print(plots$plog)


plots <- plot_uc_results(all.data,out$coda,"PFOS","072319",500)
plots <- plot_uc_results(all.data,out$coda,"PFOA","100119",500)
plots <- plot_uc_results(all.data,out$coda,"PFOA","041219",500)
plots <- plot_uc_results(all.data,out$coda,"PFOS","010720",500,quad.cal=c(avar=0.105423,bvar=6.38662,cvar=0.002752060))




