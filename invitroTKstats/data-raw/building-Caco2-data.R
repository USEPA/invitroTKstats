## This R script should be ran from the command line using
## R CMD BATCH data-raw/building-Caco2-data.R

## Script used to create data examples for Caco2 to use for function documentations and vignette
## The Excel file containing the raw data can be tracked under the Jira ticket IVTKS-75

## load necessary package
library(readxl)
library(invitroTKstats)

## Build the data guide 
level0.catalog <- data.frame(File = rep("Edited_EPA_Task 10_13_Caco-2 Compiled_LCMSGC_10032017_Data Summary_GZ.xlsm", 3),
                             Sheet = rep("Raw Data", 3),
                             Skip.Rows = c(0, 26, 52),
                             Date = rep("02 Mar 2017", 3),
                             Chemical.ID = c("BF175258","BF175270","EV0000613"),
                             ISTD.Name = rep("ISTD Name", 3),
                             Sample.ColName = rep("SampleName", 3),
                             Type.ColName = rep("Type", 3),
                             Peak.ColName = rep("Area",3),
                             ISTD.Peak.ColName = rep("ISTD Area",3),
                             Conc.ColName = rep("Test Concentration", 3),
                             AnalysisParam.ColName = rep("Transition", 3)
)

## Build chem.ids table 
chem.ids <- data.frame(Compound = c("compound 1","compound 2","compound 3"),
                       DTXSID = c("ID 1", "ID 2", "ID 3"),
                       Chem.Lab.ID = c("BF175258","BF175270","EV0000613"))

## Specify the path to the Excel file 
path <- "~/invitrotkstats/invitroTKstats/data-raw/Caco2"
## Compile level-0 data 
level0 <- merge_level0(level0.catalog = level0.catalog,
                       ## each compound has 16 samples 
                       num.rows=16,
                       num.rows.col="Num.rows",
                       istd.col="ISTD.Name",
                       type.colname.col="Type.ColName",
                       chem.ids = chem.ids,
                       chem.lab.id.col="Chem.Lab.ID",
                       chem.name.col="Compound",
                       chem.dtxsid.col="DTXSID",
                       INPUT.DIR = path,
                       catalog.out = FALSE,
                       output.res = FALSE)


# Add direction column
level0[, "Direction"] <- NA
level0[grep("_A_B", level0$Sample), "Direction"] <- "AtoB"
level0[grep("_B_A", level0$Sample), "Direction"] <- "BtoA"

## Add Vol.Receiver and Vol.Donor columns 
level0[, "Vol.Receiver"] <- NA
level0[level0$Direction == "AtoB","Vol.Receiver"] <- 0.25
level0[level0$Direction == "BtoA","Vol.Receiver"] <- 0.075

level0[, "Vol.Donor"] <- NA
level0[level0$Direction == "AtoB","Vol.Donor"] <- 0.075
level0[level0$Direction == "BtoA","Vol.Donor"] <- 0.25

## Add Dilution column
level0[, "Dilution.Factor"] <- NA
level0[level0$Type == "Blank", "Dilution.Factor"] <- 1
level0[is.na(level0$Dilution.Factor), "Dilution.Factor"] <- 4

## Run through the format function 
level1 <- format_caco2(data.in = level0,
                       sample.col="Sample",
                       lab.compound.col = "Lab.Compound.ID",
                       compound.col = "Compound",
                       series=1,
                       area.col="Peak.Area",
                       istd.col="ISTD.Peak.Area",
                       type.col="Type",
                       direction.col="Direction",
                       membrane.area=0.11,
                       receiver.vol.col="Vol.Receiver",
                       donor.vol.col="Vol.Donor",
                       compound.conc.col="Compound.Conc",
                       cal=1,
                       dilution.col="Dilution.Factor", 
                       meas.time = 2,
                       istd.name="Some ISTD Name",
                       istd.conc=1,
                       nominal.test.conc=10,
                       analysis.method = "Mass Spec",
                       analysis.instrument.col="Analysis.Params",
                       analysis.parameters="TBD",
                       output.res = FALSE
)

## Confirmed that all rows with missing responses are Blank samples
## Replace NA with 0 for blank adjustment
level1[is.na(level1$Response), "Response"] <- 0

## Save level-0 and level-1 data to use for function demo/example documentation 
save(level0, level1, file = "~/invitrotkstats/invitroTKstats/data/Caco2-example.RData")

## Include session info
utils::sessionInfo()