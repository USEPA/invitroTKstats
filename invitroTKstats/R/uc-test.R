setwd("C:/Users/jwambaug/git/invitroTKstats/invitroTKstats/R")

library(gdata)
library(parallel)
library(runjags)
source("calc_uc_fup.R")

 
dat <- read.xls("20200402_PFAS_UC_PFOA_PFOS.xlsx",stringsAsFactors=F,sheet=3,skip=12)

dat$Sample.Type <- ""
dat[dat$Type=="Standard","Sample.Type"] <- "CC"
dat[regexpr("Mixed Matrix Blank",dat$Sample.Text)!=-1,"Sample.Type"] <- "CC"
dat[regexpr("Mixed Matrix Blank",dat$Sample.Text)!=-1,"Std..Conc"] <- 0
dat[regexpr("UC-CR",dat$Sample.Text)!=-1,"Sample.Type"] <- "AF"
dat[regexpr("UC-T1",dat$Sample.Text)!=-1,"Sample.Type"] <- "T1"
  dat[regexpr("UC-T5",dat$Sample.Text)!=-1,"Sample.Type"] <- "T5"
dat[regexpr("-S1",dat$Sample.Text)!=-1,"Series"] <- 1
dat[regexpr("-S2",dat$Sample.Text)!=-1,"Series"] <- 2
dat[regexpr("-S3",dat$Sample.Text)!=-1,"Series"] <- 3
dat$Dilution.Factor <- 4*4
dat[dat[,"Sample.Type"]=="AF","Dilution.Factor"] <- 4*4*2
dat[dat[,"Sample.Type"]=="T5","Dilution.Factor"] <- 4*4*5
dat[dat[,"Sample.Type"]=="T1","Dilution.Factor"] <- 4*4*5

                                                    
dat <- subset(dat,Sample.Type!="")
dat <- subset(dat,Series %in% c(NA,"2"))

# Recalculate the response:
IS.conc <- 3
dat$Response2 <- dat$Area/dat$IS.Area*IS.conc
dat[,c("Response","Response2")] 
dat <- dat[,c("Name","Sample.Type","Series","Std..Conc","Dilution.Factor","Response2")]
colnames(dat) <- c("Sample.Name","Sample.Type","Series","Nominal.Conc","Dilution.Factor","Response")
dat$Compound.Name <- "PFOA"
dat$Cal <- 1
write.csv(dat,file="PFOA041219.csv",row.names=F)

all.data <- dat

dat <- read.xls("20200402_PFAS_UC_PFOA_PFOS.xlsx",stringsAsFactors=F,sheet=4,skip=11)

dat$Sample.Type <- ""
dat[dat$Type=="Standard","Sample.Type"] <- "CC"
dat[regexpr("Mixed Matrix Blank",dat$Sample.Text)!=-1,"Sample.Type"] <- "CC"
dat[regexpr("Mixed Matrix Blank",dat$Sample.Text)!=-1,"Std..Conc"] <- 0
dat[regexpr("UC_UF_Mix1",dat$Sample.Text)!=-1,"Sample.Type"] <- "AF"
dat[regexpr("UC_T1hr_Mix1",dat$Sample.Text)!=-1,"Sample.Type"] <- "T1"
dat[regexpr("UC_T5hr_Mix1",dat$Sample.Text)!=-1,"Sample.Type"] <- "T5"
dat[regexpr("Mix1",dat$Sample.Text)!=-1,"Series"] <- 1
dat$Dilution.Factor <- 4*4
dat[dat[,"Sample.Type"]=="CC","Std..Conc"] <- dat[dat[,"Sample.Type"]=="CC","Std..Conc"]*16
dat[dat[,"Sample.Type"]=="AF","Dilution.Factor"] <- 4*4*2
dat[dat[,"Sample.Type"]=="T5","Dilution.Factor"] <- 4*4*5
dat[dat[,"Sample.Type"]=="T1","Dilution.Factor"] <- 4*4*5
dat <- subset(dat,Sample.Type!="")
dat <- subset(dat,Series %in% c(NA,"1"))
# Recalculate the response:
IS.conc <- .01
dat$Response2 <- dat$Area/dat$IS.Area*IS.conc
dat[,c("Response","Response2")]
dat <- dat[,c("Name","Sample.Type","Series","Std..Conc","Dilution.Factor","Response2")]
colnames(dat) <- c("Sample.Name","Sample.Type","Series","Nominal.Conc","Dilution.Factor","Response")
dat$Compound.Name <- "PFOA"
dat$Cal <- 2
dat <- subset(dat,!is.na(Response))
write.csv(dat,file="PFOA100119.csv",row.names=F)


all.data <- rbind(all.data,dat)
  
  
  
   
out <- calc_uc_fup(all.data)
  


