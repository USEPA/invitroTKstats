## This R script should be ran from the command line using
## R CMD BATCH data-raw/building-Clint-data.R

## Script used to create data examples for Clint data to use for function documentation
## and vignettes. The data examples are small subsets of three compounds.

## Load necessary package
library(invitroTKstats)
library(readxl)

## smeltz2023.clint only has data for seven compounds.
## Unfortunately there's only one compound which has all samples verified with a "Y",
## other compounds all have some samples excluded from the analysis. 
## Need to go through couple verification steps from level-1 to level-2.

## Choose three compounds for the subset
clint.list <- c("DTXSID1021116", "DTXSID6023525", "DTXSID80380256")

## Prepare Level-0
## The Excel file containing the level-0 samples is not tracked with the package. 
## When re-creating the data, retrieve the file from the 'invitrotkstats' repository 
## under directory: "working/SmeltzPFAS" and save it to the path: "data-raw/Smeltz-Clint".
## Create the folder if need to. 

## Read in chem.ids
chem.ids <- read_excel("~/invitrotkstats/invitroTKstats/data-raw/Smeltz-Clint/Hep12 Data for Uncertainty Feb2022.xlsx", 
                       sheet=1)
chem.ids <- as.data.frame(chem.ids)
## In this table, the chemical names and their lab IDs are in the same column 
## Extract them into two separate columns
chem.ids$Compound <- unlist(lapply(strsplit(chem.ids[,2]," \\("),function(x) x[[1]])) 
chem.ids$Chem.Lab.ID <- gsub(")", "", unlist(lapply(strsplit(chem.ids[,2]," \\("),function(x) if (length(x)!= 1) x[[2]] else tolower(x[[1]]))))

## Read in level-0 file
## Prepare a data guide for merge_level0 
this.file <- "Hep12 Data for Uncertainty Feb2022.xlsx"

data.guide <- create_catalog(
  file = rep(this.file, 3),
  sheet = c("Ref Chem Data","Ref Chem Data", "PFAS Data"),
  skip.rows = c(6, 137, 6),
  date = rep("012822", 3),
  compound = c("phenacetin", "propranolol", "TFMFPA"),
  istd = c("propranolol-d7", "propranolol-d7", "M5PFPeA"),
  sample = rep("Name",3),
  type = rep("Type", 3),
  peak = rep("Area", 3),
  istd.peak = rep("IS Area", 3),
  conc = rep("nM", 3),
  analysis.param = rep("RT",3),
  num.rows = c(127, 127, 165),
  ## Need this sample text column to create new columns later
  additional.info = list(SampleText.ColName = rep("Sample Text",3))
)

## Pull in level-0 data
## In the merge_level0 function, specify the path to the level-0 Excel file 
## with the argument INPUT.DIR. Make necessary adjustments if needed.
clint_L0 <- merge_level0(level0.catalog  = data.guide,
                           num.rows.col="Number.Data.Rows",
                           istd.col="ISTD.Name",
                           type.colname.col="Type.ColName",
                           additional.colnames = "Sample Text", 
                           additional.colname.cols = "SampleText.ColName",
                           chem.ids = chem.ids,
                           output.res = FALSE,
                           catalog.out = FALSE,
                           INPUT.DIR = "~/invitrotkstats/invitroTKstats/data-raw/Smeltz-Clint")

## There are some additional columns needed for clint_L0 to go to level-1.
## But these columns do not exist in the original data file and  
## currently cannot be handled/added by additional utility functions. 
## Need to manually add them in. Following the steps in smeltz-hep-inactive.R,
## this script can also be found under "working/SmeltzPFAS".

## Remove rows with blank sample text
clint_L0 <- subset(clint_L0,!is.na(clint_L0[,"Sample Text"]))

## Create sample type column
## Use the package annotation of type:
clint_L0 <- subset(clint_L0,!is.na(Type))
clint_L0[clint_L0$Type == "Analyte", "Type"] <- "Cvst"
clint_L0[clint_L0$Type == "Standard", "Type"] <- "CC"
clint_L0[regexpr("inactive",tolower(clint_L0[,"Sample Text"]))!=-1,
         "Type"] <- "Inactive"

## Create time column
clint_L0[,"Time"] <- NA
clint_L0[regexpr("t240",tolower(clint_L0[,"Sample Text"]))!=-1,"Time"] <- 240/60
clint_L0[regexpr("t120",tolower(clint_L0[,"Sample Text"]))!=-1,"Time"] <- 120/60
clint_L0[regexpr("t60",tolower(clint_L0[,"Sample Text"]))!=-1,"Time"] <- 60/60
clint_L0[regexpr("t30",tolower(clint_L0[,"Sample Text"]))!=-1,"Time"] <- 30/60
clint_L0[regexpr("t15",tolower(clint_L0[,"Sample Text"]))!=-1,"Time"] <- 15/60
clint_L0[regexpr("t0",tolower(clint_L0[,"Sample Text"]))!=-1,"Time"] <- 0/60

## Remove media-only samples
clint_L0[regexpr("wem",tolower(clint_L0[,"Sample Text"]))!=-1,"Time"] <- NA
clint_L0 <- subset(clint_L0,!is.na(Time) | Type != "Cvst")

## Create a column for the dilution factors:
## Information found in lines 100 to 115 in smeltz-hep-inactive.R
## The dilution factor is indeed different across different sample types.
## -	The calibration curves (all, regardless of analyte) are diluted 240x
## -	The ref cmpds (may have label “HLB”) and PFAS labeled as “WAX1” are diluted 480x
## -	PFAS marked as WAX2 are diluted 720x

clint_L0$Dilution.Factor <- 240
clint_L0[regexpr("hlb",tolower(clint_L0[,"Sample Text"]))!=-1,
           "Dilution.Factor"] <- 480
clint_L0[regexpr("wax1",tolower(clint_L0[,"Sample Text"]))!=-1,
           "Dilution.Factor"] <- 480
clint_L0[regexpr("wax2",tolower(clint_L0[,"Sample Text"]))!=-1,
           "Dilution.Factor"] <- 720
clint_L0[regexpr("cc",tolower(clint_L0[,"Sample Text"]))!=-1,
           "Dilution.Factor"] <- 240

## Make sure concentration and area columns are numeric
clint_L0$Peak.Area <- as.numeric(clint_L0$Peak.Area)
clint_L0[,"ISTD.Peak.Area"] <- as.numeric(clint_L0[,"ISTD.Peak.Area"])
clint_L0[,"Compound.Conc"] <- as.numeric(clint_L0[,"Compound.Conc"])
clint_L0 <- subset(clint_L0, !is.na(Peak.Area) & 
                       !is.na(clint_L0[,"ISTD.Peak.Area"]))

## The 'Compound.Conc' column maps to the 'nM' column in the original Excel file
## The unit is in nM, need to convert to uM which is what the package uses
clint_L0[,"Compound.Conc"] <- as.numeric(clint_L0[,"Compound.Conc"])/1000
## Only set a nominal/expected conc for calibration curve points
clint_L0[clint_L0[,"Type"]!="CC","Compound.Conc"] <- NA
## Concentrations calculated including dilution
clint_L0[,"Compound.Conc"] <- clint_L0[,"Compound.Conc"]*clint_L0$Dilution.Factor


## Prepare Level-1 
clint_L1 <- format_clint(data.in = clint_L0,
                       sample.col ="Sample",
                       date.col="Date",
                       compound.col="Compound",
                       lab.compound.col="Lab.Compound.ID",
                       type.col="Type",
                       dilution.col="Dilution.Factor",
                       cal=1,
                       istd.conc = 10/1000,
                       istd.col= "ISTD.Peak.Area",
                       area.col = "Peak.Area",
                       density = 0.5,
                       clint.assay.conc = 1,
                       biological.replicates = 1,
                       test.conc.col="Compound.Conc",
                       time.col = "Time",
                       analysis.method = "LCMS",
                       analysis.instrument = "Unknown",
                       analysis.parameters.col = "Analysis.Params",
                       note="Sample Text",
                       output.res = FALSE
)


## Verification Steps (Level-2)

## The first step: check for correct mix
## Column "Note" records which mix was actually used for each sample of the test compounds. 
## If the mix used for a sample is not matched with the expected mix recorded on the
## summary table (chem.ids), exclude the sample due to wrong mix.

## Check the correct mix for each compound:
## Two of the selected compounds are reference chemicals, so only to check 
## DTXSID80380256. Exclude any samples with a note of WAX2 being used.
chem.ids[chem.ids$DTXSID %in% clint.list, c("DTXSID", "Mix")]
## Look for notes contain the word "WAX2"
unique(clint_L1[clint_L1$DTXSID == "DTXSID80380256", "Note"])
## Put the information into an exclusion criteria
EC <- data.frame(
  Variables = c("DTXSID|Note"),
  Values = paste("DTXSID80380256", c("WAX2 T15-B", "WAX2 T240-A Living", "WAX2 T240-B Living","WAX2 T240-C Living",
                                     "WAX2 T240-B Inactive","WAX2 T240-C Inactive"), sep = "|"),
  Message = c("Wrong Mix")
)

## Run it through the verification function
clint_L2 <- sample_verification(data.in = clint_L1, assay = "Clint",
                                exclusion.info = EC,
                                output.res = FALSE)

## The second step: check for unknown concentration
## Exclude calibration curve (CC) samples if the concentration is unknown.
## Filtering on a numeric column or detecting if a numeric value is missing are not 
## supported by the current sample_verification function. 
## We will do this step manually. 
clint_L2[clint_L2$Sample.Type=="CC" & is.na(clint_L2$Test.Compound.Conc),"Verified"] <- 
  "Unknown Conc."

## Compare with the level-2 data that's already in the package 
## see if the subset matches what's in the full data.
clint.sub <- smeltz2023.clint[smeltz2023.clint$DTXSID %in% clint.list, ]

## Check the dimension
## The example data set is expected to have one more column because the original 
## data set does not have a replicate or a series column. With new update 
## to the level-1 format function (IVTKS-4) the data are required to have either
## a biological replicate column or a technical replicate column.
dim(clint.sub)
dim(clint_L2)

colnames(clint_L2)
colnames(clint.sub)
common.cols <- intersect(colnames(clint.sub),colnames(clint_L2))
## Check if the common columns have the same names
colnames(clint.sub[,common.cols]) == colnames(clint_L2[,common.cols])
## Check all the columns with the same name have matching values
all(clint.sub[,common.cols] == clint_L2[,common.cols])
for(i in common.cols){
  test <- all(clint.sub[,i] == clint_L2[,i])
  if(test == FALSE | is.na(test)){
    print(i)
  }
}

## Address the discrepancies one by one

## Discrepancies found in the Date column:
## The original (smeltz2023.clint) data set uses the sample acquired date
## from the Excel file and the dates were converted into numbers when read in.
## I use dates extracted from the lab sample names and enter them in format mmddyy.

## Discrepancies found in the Compound.Name column:
## The original data set uses lab compound IDs/abbreviations for the Compound.Name 
## column instead of the full name. When I created the subsets the full name is 
## used for Compound.Name and lab compound IDs/abbreviations are used for the 
## Lab.Compound.ID column. I will keep my input. 

## Discrepancies found in the Lab.Compound.Name column:
## The original data set was complied without the use of merge_level0 and the data guide.
## It was complied by reading a entire sheet from the Excel file in and then filling in the 
## compound names. The first compound name, "propranolol", had its first letter capitalized
## when being filled in (line 23 in smeltz-hep-inactive.R). 
## However, compound names in the raw Excel file are all in lower case.
## Convert to all lower case and check again:
all(tolower(clint.sub[,"Lab.Compound.Name"]) == tolower(clint_L2[,"Lab.Compound.Name"]))

## Discrepancies found in the Area and ISTD.Area columns:
## Area and ISTD.Area were rounded to 5 significant digits in format_clint (lines 442-443)
## when the example data set was complied. Rounding was not used when the 
## original data set was complied.
## In the example data set, areas larger than 10,000 were rounded to the nearest whole numbers.
## Compare the columns in 5 significant digits: 
all.equal(signif(clint.sub[,"Area"], 5), clint_L2[,"Area"])
all.equal(signif(clint.sub[,"ISTD.Area"], 5), clint_L2[,"ISTD.Area"])

## Discrepancies found in the Time column:
## Time column contains missing values so the check with all() will return NA 
## Use other method to check for this column
table(clint.sub[,"Time"] == clint_L2[,"Time"])

## Discrepancies found in the Analysis.Parameters column:
## Analysis.Parameters is a numeric column in smeltz2023.clint while it is a character
## column in the example data.
## Compare them as if Analysis.Parameters in clint_L2 is a numeric column:
all.equal(clint_L2$Analysis.Parameters, clint.sub$Analysis.Parameters)
all.equal(as.numeric(clint_L2$Analysis.Parameters), clint.sub$Analysis.Parameters)

## Discrepancies found in the Level0.Sheet column:
## The script smeltz-hep-inactive.R created this column with a single value for all rows 
## (in line 137) while it in fact pulls data from two different sheets from the original 
## Excel workbook (line 9). I will keep my input since they are accurate.

## Discrepancies found in the Response column:
## Differences in this column are most likely caused by rounding errors since 
## the same issue happened with area and istd.area, and responses are calculated from them.
diffs <- clint.sub[,"Response"]- clint_L2[,"Response"]
summary(diffs[!is.na(diffs)])
## The maximum difference is 1.000e-05 and the minimum difference is -1.000e-06

## Set to the same precision and check again
clint.og.resp <- signif(as.numeric(clint.sub[,"Response"]),6)
clint.ex.resp <- signif(as.numeric(clint_L2[,"Response"]),6)
diffs <- clint.ex.resp - clint.og.resp
## Replace missing values with 0
diffs[is.na(diffs)] <- 0
## The differences, if any, should be smaller than a reasonable tolerance. Here uses four decimal places.
all(abs(diffs) <= 1e-4)
## Additional check with calc_clint_point is included below to ensure these rounding errors
## in Area, ISTD.Area and Response will not cause significant differences in later calculations. 

## Compare the columns with different column names  
colnames(clint_L2)[which(!(colnames(clint_L2) %in% common.cols))]
colnames(clint.sub)[which(!(colnames(clint.sub) %in% common.cols))]

## Biological.Replicates is a new required column filled with 1.
## Compare the concentration columns:
all(clint_L2[,"Test.Compound.Conc"] == clint.sub[,"Std.Conc"])

clint.ex.conc <- signif(as.numeric(clint_L2[,"Test.Compound.Conc"]),6)
clint.og.conc <- signif(as.numeric(clint.sub[,"Std.Conc"]),6)
diffs <- clint.ex.conc - clint.og.conc
## Replace missing values with 0
diffs[is.na(diffs)] <- 0
## The differences, if any, should be smaller than a reasonable tolerance. Here uses four decimal places.
all(abs(diffs) <= 1e-4)

## All columns are checked and any differences are documented

##---------------------------------------------##
## Run level-3 calculations with the example dataset.
ex_level3 <- calc_clint_point(data.in = clint_L2, output.res = FALSE)

## The original dataset needs to update some column names.
og_level2 <- clint.sub
og_level2$Biological.Replicates <- 1
colnames(og_level2)[which(names(og_level2) == "Std.Conc")] <- "Test.Compound.Conc"
## Run level-3 calculations with the original dataset.
og_level3 <- calc_clint_point(data.in = og_level2, output.res = FALSE)

## Compare the results. They are matched.
## Note that calc_fup_uc_point rounds Fup to the fourth significant figures.
all.equal(ex_level3$Fup,og_level3$Fup)
##---------------------------------------------##

## Save level-0 to level-2 data to use for function demo/example documentation 
save(clint_L0, clint_L1, clint_L2, file = "~/invitrotkstats/invitroTKstats/data/Clint-example.RData")

## Include session info
utils::sessionInfo()

