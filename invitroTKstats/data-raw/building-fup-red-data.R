## This R script should be ran from the command line using
## R CMD BATCH data-raw/building-fup-red-data.R

## Script used to create data examples for fup-RED to use for function documentation
## and vignettes. The data examples should be small subsets of three compounds.

## load necessary package
library(readxl)
# library(invitroTKstats) ## use when installed package is up-to-date
devtools::load_all(here::here()) ## use when installed package is not up-to-date, but branch is up-to-date with 'dev' branch
library(here)
library(dplyr)

## Picked three compounds from the data that we know all samples are verified (with "Y").
red.list <- c("DTXSID6062599","DTXSID5030030","DTXSID8031865")

## Check all compounds have all samples verified with "Y".
## They are also from the same sheet in the original file. 
## Make them easier to work with later.
unique(smeltz2023.red[smeltz2023.red$DTXSID %in% red.list, "Verified"])
unique(smeltz2023.red[smeltz2023.red$DTXSID %in% red.list, "Level0.File"])
unique(smeltz2023.red[smeltz2023.red$DTXSID %in% red.list, "Level0.Sheet"])

## Prepare Level-0
## Read in chem.ids, which is the summary table from the Excel file containing the level-0 samples.
## The Excel file is not tracked with the package. When re-creating the data,
## retrieve the file from the 'invitrotkstats' repository under directory: "working/SmeltzPFAS"
## and save it to the path: "data-raw/Smeltz-RED". Create the folder if need to. 
chem.ids <- readxl::read_xlsx(
  path = here::here("data-raw/Smeltz-RED/PFAS LC-MS RED Summary 20220709.xlsx"),
  sheet = "Summary",skip = 1,n_max = 29,col_names = TRUE
)
chem.ids <- as.data.frame(chem.ids)
# Duplicated DTXSIDs have the same Lab Compound ID, therefore can remove duplicates
chem.ids <- subset(chem.ids, !duplicated(chem.ids[,"DTXSID"]))

## In this table, the chemical names and their lab IDs are in the same column 
## Extract them into two separate columns
chem.ids$Compound <- unlist(lapply(strsplit(chem.ids[,2]," \\("),function(x) x[[1]])) 
chem.ids$Chem.Lab.ID <- gsub(")", "", unlist(lapply(strsplit(chem.ids[,2]," \\("),function(x) x[[2]])))

## Save the fup red chemical ID mapping information for the package - remove columns not needed
fup_red_cheminfo <- dplyr::select(chem.ids,-c(3:64))

# check that the number of rows in the chem information matches the number of unique DTXSID's
length(unique(fup_red_cheminfo$DTXSID))==nrow(fup_red_cheminfo)


# create chem ID mapping table for level-0 compilation - we can overwrite previous `chem.ids`
chem.ids <- create_chem_table(input.table = fup_red_cheminfo,
                              dtxsid.col = "DTXSID",
                              compound.col = "Compound",
                              lab.compound.col = "Chem.Lab.ID")

## Read in level-0 file
## Prepare a data guide for merge_level0 
this.file <- "PFAS LC-MS RED Summary 20220709.xlsx"

data.guide <- create_catalog(
  file = this.file,
  sheet = "RED1 Raw",
  skip.rows = c(279,551,823),
  date = "022322",
  compound = c("PFPeA", "PFBS", "PFOA"),
  istd = c("M5PFPeA", "M3PFBS", "M8PFOA"),
  col.names.loc = c(279, 551, 823),
  num.rows = 268,
  
  # column names 
  sample = "Name",
  type = "Type",
  peak = "Area",
  istd.peak = "IS Area",
  conc = "Std. Conc",
  analysis.param = "RT",
  ## Need this sample text column to create new columns later
  additional.info = list(SampleText.ColName = rep("Sample Text",3))
)

## Pull in level-0 data
## In the merge_level0 function, specify the path to the level-0 Excel file 
## with the argument INPUT.DIR.
fup_red_L0 <- merge_level0(level0.catalog  = data.guide,
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
             INPUT.DIR = here::here("data-raw/Smeltz-RED/"))

## There are some additional columns needed for fup_red_L0 to go to level-1.
## But these columns do not exist in the original data file and  
## currently cannot be handled/added by additional utility functions. 
## Need to manually add them in. These follow the steps found in MS-REDdata-Apr2023.Rmd 
## (starting at line 125), which can be found in the invitrotkstats repo under 'working/SmeltzPFAS'.

## Convert the area columns and concentration column to numeric 
## The RMD I'm following sets the precision for these columns to have 6 significant figures,
## however that step is skipped in this script as we are trying to avoid rounding 
## before the end analysis.
for (this.col in c("Peak.Area", "Compound.Conc", "ISTD.Peak.Area"))
  fup_red_L0[,this.col] <- as.numeric(fup_red_L0[,this.col])

## Create the Sample Type column, use the package annotations
fup_red_L0 <- subset(fup_red_L0,!is.na(fup_red_L0[,"Sample Text"]))
fup_red_L0[regexpr("CC",fup_red_L0[,"Sample Text"])!=-1,
           "Sample.Type"] <- "CC"
fup_red_L0[regexpr("-Pl-",fup_red_L0[,"Sample Text"])!=-1,
           "Sample.Type"] <- "Plasma"
fup_red_L0[regexpr("-S-",fup_red_L0[,"Sample Text"])!=-1,
           "Sample.Type"] <- "PBS"
fup_red_L0[regexpr("-EC1-",fup_red_L0[,"Sample Text"])!=-1,
           "Sample.Type"] <- "EC_acceptor"                   
fup_red_L0[regexpr("-EC2-",fup_red_L0[,"Sample Text"])!=-1,
           "Sample.Type"] <- "EC_donor"
fup_red_L0[regexpr("/T1",fup_red_L0[,"Sample Text"])!=-1,
           "Sample.Type"] <- "T0"
fup_red_L0[regexpr("/T5",fup_red_L0[,"Sample Text"])!=-1,
           "Sample.Type"] <- "Stability"
fup_red_L0[regexpr("Crash Blank",fup_red_L0[,"Sample Text"])!=-1,
           "Sample.Type"] <- "NoPlasma.Blank"
fup_red_L0[regexpr("Matrix Blank",fup_red_L0[,"Sample Text"])!=-1,
           "Sample.Type"] <- "Plasma.Blank"

## Create the Replicate column
fup_red_L0[regexpr("-A",fup_red_L0[,"Sample Text"])!=-1,
           "Replicate"] <- "A"
fup_red_L0[regexpr("-B",fup_red_L0[,"Sample Text"])!=-1,
           "Replicate"] <- "B"
fup_red_L0[regexpr("-C",fup_red_L0[,"Sample Text"])!=-1,
           "Replicate"] <- "C"

## Create the Time column
fup_red_L0[regexpr("/T1",fup_red_L0[,"Sample Text"])!=-1,
           "Time"] <- 1
fup_red_L0[regexpr("/T5",fup_red_L0[,"Sample Text"])!=-1,
           "Time"] <- 5

## Set Area of blank samples to 0
fup_red_L0[fup_red_L0[,"Sample.Type"]%in%"NoPlasma.Blank","Peak.Area"] <- 0
fup_red_L0[fup_red_L0[,"Sample.Type"]%in%"NoPlasma.Blank","ISTD.Peak.Area"] <- 1

## Remove samples with missing peak areas for analytes and internal standards
fup_red_L0 <- subset(fup_red_L0, !is.na(Peak.Area) & 
                       !is.na(fup_red_L0[,"ISTD.Peak.Area"]))

## Create the Dilution.Factor column
## Information found in lines 177 to 188 in MS-REDdata-Apr2023.Rmd
fup_red_L0[sapply(fup_red_L0[,"Sample.Type"],function(x) "CC" %in% x),"Dilution.Factor"] <- 1
fup_red_L0[sapply(fup_red_L0[,"Sample.Type"],function(x) "PBS" %in% x),"Dilution.Factor"] <- 2
fup_red_L0[sapply(fup_red_L0[,"Sample.Type"],function(x) "EC_acceptor" %in% x),"Dilution.Factor"] <- 10
fup_red_L0[sapply(fup_red_L0[,"Sample.Type"],function(x) "EC_donor" %in% x),"Dilution.Factor"] <- 10
fup_red_L0[sapply(fup_red_L0[,"Sample.Type"],function(x) "Plasma" %in% x),"Dilution.Factor"] <- 20
fup_red_L0[sapply(fup_red_L0[,"Sample.Type"],function(x) "T0" %in% x),"Dilution.Factor"] <- 10
fup_red_L0[sapply(fup_red_L0[,"Sample.Type"],function(x) "Stability" %in% x),"Dilution.Factor"] <- 10
fup_red_L0[sapply(fup_red_L0[,"Sample.Type"],function(x) "Plasma.Blank" %in% x),"Dilution.Factor"] <- 1
fup_red_L0[sapply(fup_red_L0[,"Sample.Type"],function(x) "NoPlasma.Blank" %in% x),"Dilution.Factor"] <- 1

## Convert the unit to what the package uses: uM
fup_red_L0[,"Compound.Conc"] <- as.numeric(fup_red_L0[,"Compound.Conc"])/100

## Prepare level-1 data
fup_red_L1 <- format_fup_red(data.in = fup_red_L0,
                             sample.col ="Sample",
                             date.col="Date",
                             compound.col="Compound",
                             lab.compound.col="Lab.Compound.ID",
                             type.col="Sample.Type",
                             dilution.col="Dilution.Factor",
                             technical.replicates.col ="Replicate",
                             biological.replicates = 1, 
                             cal=1,
                             area.col = "Peak.Area",
                             istd.conc = 10/1000,
                             istd.col= "ISTD.Peak.Area",
                             test.conc.col = "Compound.Conc", 
                             test.nominal.conc = 10,
                             plasma.percent = 100,
                             time.col = "Time",
                             analysis.method = "LCMS",
                             analysis.instrument = "Waters ACQUITY I-Class UHPLC - Xevo TQ-S uTQMS",
                             analysis.parameters = "RT",
                             note.col=NULL,
                             output.res = FALSE
)

## Prepare Level-2 data
## All samples are verified 
fup_red_L2 <- sample_verification(data.in = fup_red_L1, 
                                  assay = "fup-RED",
                                  # don't export the output TSV file 
                                  output.res = FALSE
                                  )

# fup_red_L2 <- fup_red_L1
# fup_red_L2$Verified <- "Y"

## Compare with smeltz2023.red to make sure the subsets match the original datasets.
red.sub <- smeltz2023.red[smeltz2023.red$DTXSID %in% red.list, ]

## Compare the the dimensions
dim(fup_red_L2)
dim(red.sub)

## red.sub erroneously removed EC1 and EC2 sample types because of a error in 
## format_fup_red 

## Remove EC_donor and EC_acceptor types from fup_red_L2 and compare again 
subfup_red_L2 <- fup_red_L2 %>% 
  dplyr::filter(!Sample.Type %in% c("EC_donor", "EC_acceptor"))

dim(subfup_red_L2)
dim(red.sub)

## Check all the columns with the same names have matching values
common_cols <- intersect(colnames(subfup_red_L2),colnames(red.sub))
## Subset to only common columns
subred.sub <- red.sub[,common_cols]
subfup_red_L2 <- subfup_red_L2[,common_cols]
## check to see if the common columns have the same information
for(i in common_cols){
  test <- all(subred.sub[,i] == subfup_red_L2[,i])
  if(test == FALSE | is.na(test)){
    print(i)
  }
}

## Discrepancies found in columns "Date", "Note" and "Time":
## The Date column in the original (smeltz2023.red) dataset uses the 
## sample acquired dates from the Excel file, while I use the dates in the lab sample names. 
## The format of the dates are also different. In smeltz2023.red the dates are in 
## "yyyy-mm-dd" format while in my subsets the dates are in the "mmddyy" format.

## For the Note column, smeltz2023.red uses 'NA' (defaults to logical class type) 
## to fill the column while I used "" (defaults to character class type) 

## Time column contains missing values so the check with all() will return NA 
## Use other method to check for this column
table(subred.sub[,"Time"] == subfup_red_L2[,"Time"])

## Discrepancies found in column Lab.Compound.Name:
## The original dataset uses the compound names ('Compound' column from Level-0) to 
## fill this column while I use the Lab.Compound.ID column from Level-0.

## Discrepancies found in columns "ISTD.Area", "Area" and "Response":
## As noted above (lines 76-79), area columns were rounded when the original 
## dataset was processed, but I did not do that for the example data files.
## Discrepancies found in the comparisons are likely rounding errors. 

## ISTD.Area
diffs <- signif(subfup_red_L2[,"ISTD.Area"],5) - signif(red.sub[,"ISTD.Area"],5)
summary(diffs[!is.na(diffs)])
## The largest difference in ISTD.Area is 1.
## Examine the differences more closely - summary of the ISTD.Area column shows that  
## at least 75% of the values are larger than 15,000, and samples 
## with the largest differences have ISTD areas that are at least ten thousands.
## These samples having differences of 1 is not a concern.
summary(subfup_red_L2[,"ISTD.Area"])
## Samples with differences larger or equal to 1
subfup_red_L2[,"ISTD.Area"][which(diffs >= 1)]

## Area
diffs <- signif(subfup_red_L2[,"Area"],5) - signif(red.sub[,"Area"],5)
summary(diffs[!is.na(diffs)])
## The largest difference in Area is 10.
## Examine the differences more closely - samples with the largest differences have
## areas more than a hundred thousand. These samples having differences of 10 is not a concern.
summary(subfup_red_L2[,"Area"])
subfup_red_L2[,"ISTD.Area"][which(diffs >= 10)]

## Response
diffs <- signif(subfup_red_L2[,"Response"],4) - signif(red.sub[,"Response"],4)
summary(diffs[!is.na(diffs)])
## The largest difference in Response is 1.000e-04.
## Examine the differences more closely 
## Check the range of the Response values - zeros are from the blank samples with area = 0
summary(subfup_red_L2[,"Response"])
## Check the smallest non-zero response 
min(subfup_red_L2[,"Response"][subfup_red_L2[,"Response"] > 0])
## Given this range, it's reasonable to assume any difference smaller than
## 1e-07 is not a concern. 
## Samples with differences larger than 1e-07
diffs[which(diffs > 1e-07)]
subfup_red_L2[,"Response"][which(diffs > 1e-07)]
## There only eleven samples that have differences larger than 1e-07, and the 
## discrepancies, when compared to the actual values, happen at the last digit of the decimals.

## Check with calc_fup_red_point to ensure these rounding errors do not cause 
## significant differences in the later calculations included below following the comparisons.

## Compare the columns with different column names 
## Remove EC_donor and EC_acceptor types from fup_red_L2 and compare again 
subfup_red_L2 <- fup_red_L2 %>% 
  dplyr::filter(!Sample.Type %in% c("EC_donor", "EC_acceptor"))

colnames(subfup_red_L2)[which(!(colnames(subfup_red_L2) %in% common_cols))]
colnames(red.sub)[which(!(colnames(red.sub) %in% common_cols))]
## Replicate columns contain missing values too, check with all() returns NA
all(subfup_red_L2[,"Technical.Replicates"] == red.sub[,"Replicate"])
table(subfup_red_L2[,"Technical.Replicates"] == red.sub[,"Replicate"])

## This check returns FALSE
all(subfup_red_L2[,"Test.Compound.Conc"] == red.sub[,"Std.Conc"])
## Investigate the differences
diffs <- subfup_red_L2[,"Test.Compound.Conc"] - red.sub[,"Std.Conc"]
summary(diffs[!is.na(diffs)])

## The maximum difference is 3e-07 and the minimum difference is -5e-7
## These differences could be due to different systems or are rounding errors

## Set to the same precision and check again
diffs <- signif(subfup_red_L2[,"Test.Compound.Conc"],6) - signif(red.sub[,"Std.Conc"],6)
## Replace missing values with 0
diffs[is.na(diffs)] <- 0
## The differences, if any, should be smaller than a reasonable tolerance. Here uses four decimal places.
all(abs(diffs) <= 1e-4)

## All columns are checked and any differences are documented

##---------------------------------------------##
## Run level-3 calculations with the example dataset.
fup_red_L3 <- calc_fup_red_point(data.in = fup_red_L2, output.res = FALSE)

## The original dataset needs to update two columns to use new names.
colnames(red.sub)[which(names(red.sub) == "Replicate")] <- "Technical.Replicates"
colnames(red.sub)[which(names(red.sub) == "Std.Conc")] <- "Test.Compound.Conc"
## Run level-3 calculations with the original dataset.
## Need to add populate biological replicates column with 1
red.sub[,"Biological.Replicates"] <-  1
og_level3 <- calc_fup_red_point(data.in = red.sub, output.res = FALSE)

## Compare the results. Minor differences. The max difference is 3.269e-6.
summary(fup_red_L3$Fup - og_level3$Fup)
all.equal(fup_red_L3$Fup,og_level3$Fup)
##---------------------------------------------##

## Run level-4 calculations with the example dataset.
tictoc::tic()
fup_red_L4 <- calc_fup_red(FILENAME = "Example",
                           data.in = fup_red_L2, 
                           JAGS.PATH = runjags::findjags(),
                           TEMP.DIR = here::here("data-raw/Smeltz-RED"),
                           OUTPUT.DIR = here::here("data-raw/Smeltz-RED"))
tictoc::toc()
## Initial processing took 480.82 seconds.

## Load Results dataframe
## To recreate, will need to change FILENAME as date will be different 
load(here::here("data-raw/Smeltz-RED/Example-fup-RED-Level4Analysis-2025-04-17.RData"))
fup_red_L4 <- Results 

## Load L2 heldout dataframe 
fup_red_L2_heldout <- read.delim(here::here("data-raw/Smeltz-RED/Example-fup-RED-Level2-heldout.tsv"),
                               sep = "\t")

## Load PREJAGS dataframe 
load(here::here("data-raw/Smeltz-RED/Example-fup-RED-PREJAGS.RData"))
fup_red_PREJAGS <- mydata

## Save level-0 and level-1 data to use for function demo/example documentations 
save(fup_red_cheminfo,fup_red_L0, fup_red_L1, fup_red_L2, fup_red_L3, fup_red_L4,
     fup_red_L2_heldout, fup_red_PREJAGS, file = here::here("data/Fup-RED-example.RData"))

## Include session info
utils::sessionInfo()
