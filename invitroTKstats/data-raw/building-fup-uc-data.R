## This R script should be ran from the command line using
## R CMD BATCH data-raw/building-fup-uc-data.R

## Script used to create data examples for fup-UC to use for function documentation
## and vignettes. The data examples are small subsets of three compounds.

## load necessary packages
library(readxl)
# library(invitroTKstats) ## use when installed package is up-to-date
devtools::load_all(here::here()) ## use when installed package is not up-to-date, but branch is up-to-date with 'dev' branch
library(here)
library(dplyr)

## Chose three compounds that have all samples verified
uc.list <- c("DTXSID00192353", "DTXSID0059829", "DTXSID3037707")
unique(smeltz2023.uc[smeltz2023.uc$DTXSID %in% uc.list, "Verified"])

## Prepare Level-0
## The Excel file containing the level-0 samples is not tracked with the package. 
## When re-creating the data, retrieve the file from the 'invitrotkstats' repository 
## under directory: "working/SmeltzPFAS" and save it to the path: "data-raw/Smeltz-UC".
## Create the folder if need to.

## Read in the assay summary table from the Excel file. 
## This table records which date each compound was tested on and which mix was used.
## This table will be used later in the process to remove samples with the wrong mix.
assayinfo <- read_excel(here::here(
  "data-raw/Smeltz-UC/20220201_PFAS-LC_FractionUnbound_MGS.xlsx"),
  sheet=1)
assayinfo <- as.data.frame(assayinfo)
## Fill in the date column for all rows:
this.date <- assayinfo[1,"LCMS Analysis Date"]
this.row <- 1
while (this.row <= dim(assayinfo)[1])
{
  if (is.na(assayinfo[this.row,"LCMS Analysis Date"]))
  {
    assayinfo[this.row,"LCMS Analysis Date"] <- this.date 
  } else {
    this.date <- assayinfo[this.row,"LCMS Analysis Date"] 
  }
  this.row <- this.row + 1
}
## Remove "UTC" from the dates
assayinfo[,1] <- sapply(assayinfo[,1],function(x) gsub("  UTC","",x))

## Read in chem.ids
chem.ids <- readxl::read_xlsx(
  path = here::here("data-raw/Smeltz-UC/20220201_PFAS-LC_FractionUnbound_MGS.xlsx"),
  sheet = "Summarized Wetmore Fu Values"
)
chem.ids <- as.data.frame(chem.ids)
chem.ids <- subset(chem.ids, !duplicated(chem.ids[,"DTXSID"]))
## In this table, the chemical names and their lab IDs are in the same column 
## Extract them into two separate columns
chem.ids$Compound <- unlist(lapply(strsplit(chem.ids[,2]," \\("),function(x) x[[1]])) 
chem.ids$Chem.Lab.ID <- gsub(")", "", unlist(lapply(strsplit(chem.ids[,2]," \\("),function(x) if (length(x) != 1) x[[2]] else NA)))

## Save the fup uc chemical ID mapping information for the package - remove columns not needed
fup_uc_cheminfo <- dplyr::select(chem.ids,-c("Mean fu","SD fu","CV fu","Category","...7"))

# check that the number of rows in the chem information matches the number of unique DTXSID's
length(unique(fup_uc_cheminfo$DTXSID))==nrow(fup_uc_cheminfo)

# create chem ID mapping table for level-0 compilation - we can overwrite previous `chem.ids`
chem.ids <- create_chem_table(input.table = fup_uc_cheminfo,
                              dtxsid.col = "DTXSID",
                              compound.col = "Compound",
                              lab.compound.col = "Chem.Lab.ID")

## Prepare a data guide for merge_level0 
this.file <- "20220201_PFAS-LC_FractionUnbound_MGS.xlsx"

data.guide <- create_catalog(
  file = this.file,
  sheet = c("20200103","20210308","20201123"),
  skip.rows = c(572,7,138),
  col.names.loc = c(572, 7, 138),
  date = c("010320","030821","112320"),
  compound = c("8:2 FTS", "PFOA-F", "K-PFBS"),
  istd = c("M2-8:2FTS", "M8PFOA", "M3PFBS"),
  num.rows = c(109, 106, 127),
  
  # column names 
  sample = "Name",
  type = "Type",
  peak = "Area",
  istd.peak = "IS Area",
  conc = c("uM", "nM", "nM"),
  analysis.param = "RT",
  ## Need this sample text column to create new columns later
  additional.info = list(SampleText.ColName = rep("Sample Text",3))
)

## Pull in level-0 data
## In the merge_level0 function, specify the path to the level-0 Excel file 
## with the argument INPUT.DIR. 
fup_uc_L0 <- merge_level0(level0.catalog  = data.guide,
                           num.rows.col="Number.Data.Rows",
                           istd.col="ISTD.Name",
                           type.colname.col="Type.ColName",
                           additional.colnames = "Sample Text",
                           additional.colname.cols = "SampleText.ColName",
                           chem.ids = chem.ids,
                           chem.lab.id.col = "Lab.Compound.Name",
                           chem.name.col = "Compound.Name",
                           output.res = FALSE,
                           catalog.out = FALSE,
                           INPUT.DIR = here::here("data-raw/Smeltz-UC"))

## There are some additional columns needed for fup_uc_L0 to go to level-1.
## But these columns do not exist in the original data file and  
## currently cannot be handled/added by additional utility functions. 
## Need to manually add them in. These follow the steps found in MSdata-Mar2022.R 
## (starting at line 227), which can be found in the invitrotkstats repo under 'working/SmeltzPFAS'.

## Remove rows with empty Sample.Text
fup_uc_L0 <- subset(fup_uc_L0,!is.na(fup_uc_L0[,"Sample Text"]))

# Extract sample types from the column Sample.Text:
fup_uc_L0[regexpr("AF",unlist(fup_uc_L0[,"Sample Text"]))!=-1,"Sample.Type"] <- "AF"
fup_uc_L0[regexpr("UF",unlist(fup_uc_L0[,"Sample Text"]))!=-1,"Sample.Type"] <- "AF"
fup_uc_L0[regexpr("T1",unlist(fup_uc_L0[,"Sample Text"]))!=-1,"Sample.Type"] <- "T1"
fup_uc_L0[regexpr("T5",unlist(fup_uc_L0[,"Sample Text"]))!=-1,"Sample.Type"] <- "T5"
fup_uc_L0[fup_uc_L0$Type == "Standard","Sample.Type"] <- "CC"
fup_uc_L0[fup_uc_L0$Type == "Blank","Sample.Type"] <- "Blank"

## Remove unused samples (QC samples):
fup_uc_L0 <- subset(fup_uc_L0,!is.na(Sample.Type))

## Make sure numeric columns are in the correct class
fup_uc_L0[,"Compound.Conc"] <- as.numeric(fup_uc_L0[,"Compound.Conc"])
fup_uc_L0[,"Peak.Area"] <- as.numeric(fup_uc_L0[,"Peak.Area"])
fup_uc_L0[,"ISTD.Peak.Area"] <- as.numeric(fup_uc_L0[,"ISTD.Peak.Area"])
## The unit of compound concentrations varies from compound to compound.
## Concentrations of compounds other than 8:2 FTS are in nM. 
## Convert them to be in the unit of uM, which is what the package uses for the calculations.
fup_uc_L0[fup_uc_L0$Lab.Compound.ID != "8:2 FTS", "Compound.Conc"] <- fup_uc_L0[fup_uc_L0$Lab.Compound.ID != "8:2 FTS", "Compound.Conc"]/1000

fup_uc_L0[fup_uc_L0[,"Sample.Type"]!="CC","Compound.Conc"] <- NA 
## Create the Dilution.Factor column
## The information comes from lines 17 to 21 in MSdata-Mar2022.R:
fup_uc_L0[fup_uc_L0[,"Sample.Type"]=="AF","Dilution.Factor"] <- 2*16
fup_uc_L0[fup_uc_L0[,"Sample.Type"]=="T1","Dilution.Factor"] <- 5*16
fup_uc_L0[fup_uc_L0[,"Sample.Type"]=="T5","Dilution.Factor"] <- 5*16
fup_uc_L0[fup_uc_L0[,"Sample.Type"]=="CC","Dilution.Factor"] <- 1
fup_uc_L0[fup_uc_L0[,"Sample.Type"]=="Blank","Dilution.Factor"] <- 1

## Treat the blanks as calibration data with concentration 0:
fup_uc_L0[fup_uc_L0[,"Sample.Type"]=="Blank","Compound.Conc"] <- 0
fup_uc_L0[fup_uc_L0[,"Sample.Type"]=="Blank","Sample.Type"] <- "CC"

## Remove CC samples that don't have a concentration:
fup_uc_L0 <- subset(fup_uc_L0,!(Sample.Type=="CC" & is.na(Compound.Conc)))

## Create a replicate column
fup_uc_L0[,"Replicate"] <- ""
fup_uc_L0[regexpr("_A",fup_uc_L0[,"Sample Text"])!=-1,"Replicate"] <- "A"
fup_uc_L0[regexpr("_B",fup_uc_L0[,"Sample Text"])!=-1,"Replicate"] <- "B"
fup_uc_L0[regexpr("_C",fup_uc_L0[,"Sample Text"])!=-1,"Replicate"] <- "C"

## Identify and remove mixes that do not match what should be used for the analytes 
## according to the assay information table
bad.mix <- rep(FALSE,dim(fup_uc_L0)[1])
for (this.id in uc.list){
  this.date <- fup_uc_L0[fup_uc_L0$DTXSID == this.id, "Date"][1]
  this.assay.index <- which(
    assayinfo[,"DTXSID"] == this.id &
      assayinfo[,"LCMS Analysis Date"] == as.Date(this.date, "%m%d%y"))
  if (any(this.assay.index))
  {
    this.mix <- assayinfo[this.assay.index,"Mix"]
    for (other.mix in 1:3)
      if (other.mix != this.mix)
      {
        bad.mix <- bad.mix |
          (fup_uc_L0[,"DTXSID"] == this.id &
             fup_uc_L0[,"Date"] == this.date &
             regexpr(paste("Mix",other.mix,sep=""),fup_uc_L0[,"Sample Text"])!=-1)
      }
  }
  
}

fup_uc_L0 <- subset(fup_uc_L0,!bad.mix)

## Prepare level-1 data
fup_uc_L1 <- format_fup_uc(data.in = fup_uc_L0,
                           sample.col="Sample",
                           compound.col="Compound",
                           test.conc.col ="Compound.Conc", 
                           lab.compound.col="Lab.Compound.ID", 
                           type.col="Sample.Type", 
                           istd.col="ISTD.Peak.Area",
                           cal.col = "Date",
                           area.col = "Peak.Area",
                           istd.conc = 1,
                           note.col = NULL,
                           test.nominal.conc = 10,
                           analysis.method = "UPLC-MS/MS",
                           analysis.instrument = "Waters Xevo TQ-S micro (QEB0036)",
                           analysis.parameters.col = "Analysis.Params",
                           biological.replicates = 1,
                           technical.replicates.col = "Replicate",
                           output.res = FALSE
                          )

## Verify all samples
fup_uc_L2 <- sample_verification(data.in = fup_uc_L1,
                                 assay = "fup-UC",
                                 # don't export the output TSV file
                                 output.res = FALSE)

## Compare with smeltz2023.uc to make sure the subsets match the original datasets.
uc.sub <- smeltz2023.uc[smeltz2023.uc$DTXSID %in% uc.list, ]

## Check the dimensions
dim(uc.sub)
dim(fup_uc_L2)

colnames(uc.sub)
colnames(fup_uc_L2)
common.cols <- intersect(colnames(uc.sub),colnames(fup_uc_L2))
all(colnames(uc.sub[,common.cols]) == colnames(fup_uc_L2[,common.cols]))

## The order of the samples seems to be different
## Sort the data by lab sample name and DTXSID (sample name alone is not a unique identifier)
og_level2 <- uc.sub[with(uc.sub, order(Lab.Sample.Name, DTXSID)), ]
ex_level2 <- fup_uc_L2[with(fup_uc_L2, order(Lab.Sample.Name, DTXSID)), ]

## Compare values in the columns with the same names
all(og_level2[,common.cols] == ex_level2[,common.cols])
for(i in common.cols){
  test <- all(og_level2[,i] == ex_level2[,i])
  if(test == FALSE | is.na(test)){
    print(i)
  }
}

## Address the discrepancies one by one

## Discrepancies found in columns Date and Calibration:
## These two columns have the same information because we are using the dates as calibration indices.
## The two datasets should have matched dates but in different formats. 
## The original dataset has the dates in format 'yyyy-mm-dd' and the example dataset has the dates
## in format 'mmddyy'.

unique(ex_level2$Date)
unique(og_level2$Date)

## Discrepancies found in the Lab.Compound.Name column:
## The original dataset maps the 'Compound' column from Level-0 to this column (line 355 in MSdata-Mar2022.R),
## which is incorrect. I used the 'Lab.Compound.ID' column from Level-0. I will keep my input.

## Discrepancies found in the Area: 
## Some Areas were NA which causes all() to also return NA. Remove the NA entries and 
## compare again. 
all.equal(og_level2[!is.na(og_level2$Area), "Area"], ex_level2[!is.na(ex_level2$Area), "Area"])

## Discrepancies found in the ISTD.Area column:
## Some minor rounding discrepancies in the ISTD.Area column 
mismatches <- which(og_level2[,"ISTD.Area"] != ex_level2[,"ISTD.Area"])
diffs <- og_level2[mismatches,"ISTD.Area"]  - ex_level2[mismatches,"ISTD.Area"]
summary(diffs)
## The max difference in ISTD.Area is 0.0477439. 

## Discrepancies found in the Analysis.Parameters column:
## Analysis.Parameters is a numeric column in smeltz2023.uc while it is a character
## column in the example data.
## Convert Analysis.Parameters to numeric column and check again:
ex_level2$Analysis.Parameters <- as.numeric(ex_level2$Analysis.Parameters)
all.equal(ex_level2$Analysis.Parameters, og_level2$Analysis.Parameters)

## Discrepancies found in the Note column:
## The original data maps the Replicate column to the note column (line 358 in MSdata-Mar2022.R).
## I will keep my input (i.e. filled with "").

## Discrepancies found in the Response column:
## Differences in this column are most likely caused by rounding errors since 
## the same issue happened with area and istd.area, and responses are calculated from them.
diffs <- og_level2[,"Response"]- ex_level2[,"Response"]
summary(diffs[!is.na(diffs)])

## The maximum difference is 4.305e-03 and the minimum difference is -4.943e-03.
## Difference in the third decimal place is not reasonable.
## The sample with the largest difference has a peak area of 9925.596 and an 
## ISTD.Area of 581.952. 
ex_level2[which.max(diffs), "Area"]
ex_level2[which.max(diffs),"ISTD.Area"]

## To preserve the precision, use the area and ISTD area from Level-0 to calculate the responses
## and compare with responses from the original data set. 
ex_level0 <- fup_uc_L0[with(fup_uc_L0, order(Sample, DTXSID)), ]
## ISTD concentration is set to 1 (line 330 in MSdata-Mar2022.R)
## Set to 4 significant digits because that's that format_fup_uc does (lines 410-411)
ex_level0[,"Response"] <- signif(as.numeric(ex_level0[,"Peak.Area"]) /
                                  as.numeric(ex_level0[,"ISTD.Peak.Area"]) * 1,4)

diffs <- og_level2[,"Response"] - ex_level0[,"Response"]
## Replace missing values with 0
diffs[is.na(diffs)] <- 0
## The differences, if any, should be smaller than a reasonable tolerance. Here uses four decimal places.
all(abs(diffs) <= 1e-4)
## Additional check with calc_fup_uc_point is included below to ensure these rounding errors
## will not cause significant differences in later calculations. 

## Compare the columns with different column names  
colnames(og_level2)[which(!(colnames(og_level2) %in% common.cols))]
colnames(ex_level2)[which(!(colnames(ex_level2) %in% common.cols))]

## As mentioned above, the original data maps replicate to the Note column,
## and creates a Series column filled with 1 (line 261 in MSdata-Mar2022.R).
## Should compare Technical.Replicates from the example data to Note in the original data,  
## compare Biological.Replicates to Series in the original data, and compare Test.Nominal.Conc
## to UC.Assay.T1.Conc in the origindal data. 
all.equal(ex_level2[, "Biological.Replicates"], og_level2[, "Series"])
all.equal(ex_level2[, "Technical.Replicates"], og_level2[, "Note"])
all.equal(ex_level2[, "Test.Nominal.Conc"], og_level2[, "UC.Assay.T1.Conc"])


## Check the concentration columns:
all.equal(ex_level2[, "Test.Compound.Conc"], og_level2[, "Standard.Conc"])
## Standard.Conc was rounded to have 4 significant figures when the original dataset
## was processed, we don't want to do that for the example dataset.
## Compare the two columns in the same precision.
all.equal(signif(ex_level2[, "Test.Compound.Conc"],4), og_level2[, "Standard.Conc"])

## All columns are checked and any differences are documented

##---------------------------------------------##
## Run level-3 calculations with the example dataset.
fup_uc_L3 <- calc_fup_uc_point(data.in = fup_uc_L2, output.res = FALSE)

## The original dataset needs to update some column names.
colnames(og_level2)[which(names(og_level2) == "Series")] <- "Biological.Replicates"
colnames(og_level2)[which(names(og_level2) == "Standard.Conc")] <- "Test.Compound.Conc"
## Run level-3 calculations with the original dataset.
## Rename UC.Assay.T1.Conc as Test.Nominal.Conc
og_level2 <- rename(og_level2, "Test.Nominal.Conc" = "UC.Assay.T1.Conc")
og_level3 <- calc_fup_uc_point(data.in = og_level2, output.res = FALSE)

## Compare the results.
## The last two compounds are switched in the og_level2 data frame so need to 
## sort the Fup values. 
## They are matched with slight rounding discrepancies. 
## Note that calc_fup_uc_point rounds Fup to the fourth significant figures.
all.equal(sort(fup_uc_L3$Fup), sort(og_level3$Fup))
## A mean relative difference of 1.8e-4. 
##---------------------------------------------##

## Run level-4 calculations with the example dataset
path.to.jags <- runjags::findjags()
## runjags::findjags() not working as argument to calc_fup_uc
## manually remove trailing path
path.to.jags <- gsub("/bin/jags-terminal.exe", "", path.to.jags)

tictoc::tic()
fup_uc_L4 <- calc_fup_uc(FILENAME = "Example",
                         data.in = fup_uc_L2, 
                         JAGS.PATH = path.to.jags,
                         TEMP.DIR = here::here("data-raw/Smeltz-UC"),
                         OUTPUT.DIR = here::here("data-raw/Smeltz-UC")
                         )
tictoc::toc()
## Initial processing took 617.86 seconds 

## Load Results dataframe
## To recreate, will need to change FILENAME as date will be different 
load(here::here("data-raw/Smeltz-UC/Example-fup-UC-Level4Analysis-2025-04-15.RData"))
fup_uc_L4 <- Results 

## Load L2 heldout dataframe 
fup_uc_L2_heldout <- read.delim(here::here("data-raw/Smeltz-UC/Example-fup-UC-Level2-heldout.tsv"),
                                 sep = "\t")

## Load PREJAGS dataframe 
load(here::here("data-raw/Smeltz-UC/Example-fup-UC-PREJAGS.RData"))
fup_uc_PREJAGS <- mydata

## Save all levels to use for function demo/example documentation 
save(fup_uc_cheminfo,fup_uc_L0, fup_uc_L1, fup_uc_L2, fup_uc_L3,
     fup_uc_L4, fup_uc_L2_heldout, fup_uc_PREJAGS, file = here::here("data/Fup-UC-example.RData"))

## Include session info
utils::sessionInfo()
