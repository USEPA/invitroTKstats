## This R script should be ran from the command line using
## R CMD BATCH data-raw/building-Caco2-data.R

## Script used to create data examples for Caco2 to use for function documentations and vignette
## The original Excel file containing the raw data can be found in the following directory: 
## "<...>\CCTE_ExpoCast\ExpoCast2019\HTTKNewData\EPA_Task 10_13_Caco-2 Compiled_LCMSGC_10032017_Data Summary.xlsm"

## The actual file used for this script is an edited version of the original Excel file. 
## Edits were made to add information required for merge_level0.
## Edits include: added headers (column names) to the beginning of each compound chunk 
## and added a 'Test Concentration' column and a 'Type' column.

## load necessary package
library(readxl)
library(invitroTKstats)

## Build the data guide 
level0.catalog <- data.frame(File = rep("Edited_EPA_Task 10_13_Caco-2 Compiled_LCMSGC_10032017_Data Summary_GZ.xlsm", 3),
                             Sheet = rep("Raw Data", 3),
                             Skip.Rows = c(0, 26, 52),
                             Date = rep("03 Oct 2017", 3),
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
## Use placeholders chemical identifiers before finding the true identifiers for the chemical lab IDs
chem.ids <- data.frame(Compound = c("Compound 1","Compound 2","Compound 3"),
                       DTXSID = c("ID 1", "ID 2", "ID 3"),
                       Chem.Lab.ID = c("BF175258","BF175270","EV0000613"))

## Specify the path to the Excel file 
## The Excel file is not tracked with the package. When re-creating the data,
## retrieve the file from the directory mentioned above and save it to the path below.
## Make necessary adjustments if needed. 
path <- "~/invitrotkstats/invitroTKstats/data-raw/Caco2"
## Compile level-0 data 
caco2_L0 <- merge_level0(level0.catalog = level0.catalog,
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

## There are some additional columns needed for caco2_L0 to go to level-1.
## But currently these columns cannot be handled/added by merge_level0 or 
## additional utility functions. Need to manually add them in. 

## Add direction column
caco2_L0[, "Direction"] <- NA
caco2_L0[grep("_A_B", caco2_L0$Sample), "Direction"] <- "AtoB"
caco2_L0[grep("_B_A", caco2_L0$Sample), "Direction"] <- "BtoA"

## Add Vol.Receiver and Vol.Donor columns 
## The amounts can be found in the Volume column (column K)
## in the original Excel file, in sheet 'Raw Data'.
caco2_L0[, "Vol.Receiver"] <- NA
caco2_L0[caco2_L0$Direction == "AtoB","Vol.Receiver"] <- 0.25
caco2_L0[caco2_L0$Direction == "BtoA","Vol.Receiver"] <- 0.075

caco2_L0[, "Vol.Donor"] <- NA
caco2_L0[caco2_L0$Direction == "AtoB","Vol.Donor"] <- 0.075
caco2_L0[caco2_L0$Direction == "BtoA","Vol.Donor"] <- 0.25

## Add Dilution column
## Dilution information can be found in the original Excel file, in sheet 'Raw Data'
## in a row says 'Dilution' at the end of each compound section.
caco2_L0[, "Dilution.Factor"] <- NA
caco2_L0[caco2_L0$Type == "Blank", "Dilution.Factor"] <- 1
caco2_L0[is.na(caco2_L0$Dilution.Factor), "Dilution.Factor"] <- 4

## Run through the format function 
caco2_L1 <- format_caco2(data.in = caco2_L0,
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
all(caco2_L1[is.na(caco2_L1$Response), "Sample.Type"] == "Blank")

## Rows with NA responses will be dropped in calc_caco2_point (level-3 function), 
## in this example data, it is equivalent to dropping all Blank samples.
## Blank samples are needed for the calculation, so 
## need to replace NA with 0 for blank adjustment.
caco2_L1[is.na(caco2_L1$Response), "Response"] <- 0

## Save level-0 and level-1 data to use for function demo/example documentation 
save(caco2_L0, caco2_L1, file = "~/invitrotkstats/invitroTKstats/data/Caco2-example.RData")

## Include session info
utils::sessionInfo()