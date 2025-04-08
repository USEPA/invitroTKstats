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
# library(invitroTKstats) ## use when installed package is up-to-date
devtools::load_all("~/Git/invitrotkstats/invitroTKstats") ## use when installed package is not up-to-date, but branch is up-to-date with 'dev' branch
library(here)
library(dplyr)

## Build chem.ids table 
## The Excel file ("SupTable1-AnalyticMethods.xlsx" and "Caco2.xlsx") containing the chemical ID mapping is not tracked with the package. 
## When re-creating the data, retrieve the file from the directory "<...>\CCTE_ExpoCast\ExpoCast2019\HTTKNewData\Summary"
## and save it to the path: "data-raw/Honda-Caco2".
## Create the folder if need to.
suptab1 <- readxl::read_xlsx(
  path = here::here("data-raw/Honda-Caco2/SupTable1-AnalyticalMethods.xlsx"),
  sheet = "SupTable1-AnalyticalMethods"
)
suptab1 <- as.data.frame(suptab1)

caco2_compdata <- readxl::read_xlsx(
  path = here::here("data-raw/Honda-Caco2/Caco2.xlsx"),
  sheet = "Caco2"
)
caco2_compdata <- as.data.frame(caco2_compdata)
caco2_compdata <- dplyr::select(caco2_compdata,c("test_article","dtxsid","casrn")) # keep only the columns we need for chemical ID mapping
caco2_compdata <- dplyr::mutate(caco2_compdata,CASRN = `casrn`) # create a new column to enable merging `suptab1` and `caco2_compdata`

# combine the two sources of chemical ID mapping
chem.ids <- dplyr::left_join(suptab1,caco2_compdata,by = "CASRN")
chem.ids <- dplyr::filter(chem.ids,duplicated(`DTXSID`))

## Save the caco2 chemical ID mapping information for the package - remove columns not needed
caco2_cheminfo <- dplyr::select(chem.ids,-c("LC","Agilent.QQQ","Water.s.Xevo","AB.Sciex.Qtrap","GC","Agilent.GCMS","GCTOF","Comment","dtxsid","casrn"))

# check that the number of rows in the chem information matches the number of unique DTXSID's
length(unique(caco2_cheminfo$DTXSID))==nrow(caco2_cheminfo)

# create chem ID mapping table for level-0 compilation - we can overwrite previous `chem.ids`
chem.ids <- create_chem_table(input.table = caco2_cheminfo,
                              dtxsid.col = "DTXSID",
                              compound.col = "PREFERRED_NAME",
                              lab.compound.col = "test_article")

## Build the data guide 
level0.catalog <- create_catalog(file = "Edited_EPA_Task 10_13_Caco-2 Compiled_LCMSGC_10032017_Data Summary_GZ.xlsm",
                                 sheet = "Raw Data",
                                 skip.rows = c(1, 27, 53),
                                 col.names.loc = c(1, 27, 53),
                                 date = "100317",
                                 compound = c("BF175258","BF175270","EV0000613"),
                                 istd = "ISTD Name",
                                 num.rows = 16, 
                                 
                                 # column names 
                                 sample = "SampleName",
                                 type = "Type", 
                                 peak = "Area",
                                 istd.peak = "ISTD Area",
                                 conc = "Test Concentration",
                                 analysis.param = "Transition"
                                 )

## Specify the path to the Excel file 
## The Excel file is not tracked with the package. When re-creating the data,
## retrieve the file from the directory mentioned above and save it to the path below.
## Make necessary adjustments if needed. 
path <- "~/Git/invitrotkstats/invitroTKstats/data-raw/Caco2"
## Compile level-0 data 
caco2_L0 <- merge_level0(level0.catalog = level0.catalog,
                       num.rows.col="Number.Data.Rows",
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
                       biological.replicates = 1,
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
                       time = 2, 
                       istd.name.col = "ISTD.Name",
                       istd.conc=1,
                       nominal.test.conc=10,
                       analysis.method.col = "Analysis.Params",
                       # These data was collected in close time proximity to the 
                       # Wambaugh2019 work which used Agilent.GCMS for GC, 
                       # thus we assume the same instrument was used for these data.
                       analysis.instrument="Agilent.GCMS",
                       analysis.parameters="TBD",
                       output.res = FALSE,
                       note.col = NULL
)

## Confirmed that all rows with missing responses are Blank samples
all(caco2_L1[is.na(caco2_L1$Response), "Sample.Type"] == "Blank")

## Rows with NA responses will be dropped in calc_caco2_point (level-3 function), 
## in this example data, it is equivalent to dropping all Blank samples.
## Blank samples are needed for the calculation, so 
## need to replace NA with 0 for blank adjustment.
caco2_L1[is.na(caco2_L1$Response), "Response"] <- 0

## Set every sample in the data as verified and export level-2 as a .tsv to 
## the vignette directory for use of function example and vignette.
caco2_L2 <- sample_verification(FILENAME = "Examples", 
                                data.in = caco2_L1, 
                                assay = "Caco-2", 
                                OUTPUT.DIR="~/Git/invitrotkstats/invitroTKstats/vignettes")

## Verify that exported level-2 is unrounded 
level2_TSV <- read.delim("~/Git/invitrotkstats/invitroTKstats/vignettes/Examples-Caco-2-Level2.tsv",
           sep = "\t")
## Compare Responses to caco2_L2 Responses 
isTRUE(all.equal(level2_TSV$Response,caco2_L2$Response))

## Save level-0 and level-1 data to use for function demo/example documentation 
save(caco2_cheminfo,caco2_L0, caco2_L1, caco2_L2, file = "~/Git/invitrotkstats/invitroTKstats/data/Caco2-example.RData")

## Include session info
utils::sessionInfo()
