## This R script should be ran from the command line using
## R CMD BATCH data-raw/building-fup-uc-data.R

## Script used to create data examples for fup-UC to use for function documentations 
## and vignette. The data examples are small subsets of three compounds.

## load necessary package
library(readxl)
library(invitroTKstats)

## Chose three compounds that have all samples verified
uc.list <- c("DTXSID00192353", "DTXSID0059829", "DTXSID3037707")
unique(smeltz2023.uc[smeltz2023.uc$DTXSID %in% uc.list, "Verified"])

## Prepare Level-0
## The Excel file containing the level-0 samples is not tracked with the package. 
## When re-creating the data, retrieve the file from the 'invitrotkstats' repository 
## under directory: "working/SmeltzPFAS" and save it to the path: "data-raw/Smeltz-UC".
## Create the folder if need to.

## Read in the assay summary table from the Excel file. 
## This table records what date each compound was tested on and what mix was used.
## This table will be used later in the process to remove samples with the wrong mix.
assayinfo <- read_excel(
  "~/invitrotkstats/invitroTKstats/data-raw/Smeltz-UC/20220201_PFAS-LC_FractionUnbound_MGS.xlsx",
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
## Remove "UTC" from the date
assayinfo[,1] <- sapply(assayinfo[,1],function(x) gsub("  UTC","",x))

## Read in chem.ids
cheminfo <- read_excel(
  "~/invitrotkstats/invitroTKstats/data-raw/Smeltz-UC/20220201_PFAS-LC_FractionUnbound_MGS.xlsx",
  sheet=2)[, 1:2]
cheminfo <- as.data.frame(cheminfo)
## In this table, the chemical names and their lab IDs are in the same column 
## Extract them into two separate columns
cheminfo$Compound <- unlist(lapply(strsplit(cheminfo[,2]," \\("),function(x) x[[1]])) 
cheminfo$Chem.Lab.ID <- gsub(")", "", unlist(lapply(strsplit(cheminfo[,2]," \\("),function(x) if (length(x) != 1) x[[2]] else NA)))

## Prepare a data guide for merge_level0 
this.file <- "20220201_PFAS-LC_FractionUnbound_MGS.xlsx"
data.guide <- create_catalog(
  file = rep(this.file, 3),
  sheet = c("20200103","20210308","20201123"),
  skip.rows = c(571,6,137),
  date = c("2020-01-03","2021-03-08","2020-11-23"),
  compound = c("8:2 FTS", "PFOA-F", "K-PFBS"),
  istd = c("M2-8:2FTS", "M8PFOA", "M3PFBS"),
  sample = rep("Name",3),
  type = rep("Type", 3),
  peak = rep("Area", 3),
  istd.peak = rep("IS Area", 3),
  conc = c("uM", "nM", "nM"),
  analysis.param = rep("RT",3),
  num.rows = c(109,106,127),
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
                           additional.colnames = "Sample Text", "nM",
                           additional.colname.cols = "SampleText.ColName",
                           chem.ids = cheminfo,
                           output.res = FALSE,
                           catalog.out = FALSE,
                           INPUT.DIR = "~/invitrotkstats/invitroTKstats/data-raw/Smeltz-UC")

## There are some additional columns needed for fup_uc_L0 to go to level-1.
## But these columns do not exist in the original data file and  
## currently cannot be handled/added by additional utility functions. 
## Need to manually add them in. Following the steps in MSdata-Mar2022.R,
## this script can also be found under "working/SmeltzPFAS".

## Remove rows with empty sample text
fup_uc_L0 <- subset(fup_uc_L0,!is.na(fup_uc_L0[,"Sample Text"]))

# Extract the sample type from column Sample Text:
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
fup_uc_L0[,"Peak.Area"] <- as.numeric(unlist(fup_uc_L0[,"Peak.Area"]))
fup_uc_L0[,"ISTD.Peak.Area"] <- as.numeric(unlist(fup_uc_L0[,"ISTD.Peak.Area"]))
## The unit of compound concentration varies from compound to compound
## concentration of compound other than 8:2 FTS are in nM, need to convert the unit to uM which is what the package uses
fup_uc_L0[fup_uc_L0$Lab.Compound.ID != "8:2 FTS", "Compound.Conc"] <- fup_uc_L0[fup_uc_L0$Lab.Compound.ID != "8:2 FTS", "Compound.Conc"]/1000

fup_uc_L0[fup_uc_L0[,"Sample.Type"]!="CC","Compound.Conc"] <- NA 
## Create the dilution factor column
## The information come from lines 17 to 21 in MSdata-Mar2022.R:
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

## Identify and remove mixes that do not match what's should be used for the analyte 
## according to the assay information table
bad.mix <- rep(FALSE,dim(fup_uc_L0)[1])
for (this.id in uc.list){
  this.date <- fup_uc_L0[fup_uc_L0$DTXSID == this.id, "Date"][1]
  this.assay.index <- which(
    assayinfo[,"DTXSID"] == this.id &
      assayinfo[,"LCMS Analysis Date"] == this.date)
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
                           uc.assay.conc = 10,
                           analysis.method = "UPLC-MS/MS",
                           analysis.instrument = "Waters Xevo TQ-S micro (QEB0036)",
                           analysis.parameters.col = "Analysis.Params",
                           biological.replicates.col = "Replicate",
                           output.res = FALSE
                          )

## Verify all samples
fup_uc_L2 <- sample_verification(data.in = fup_uc_L1, output.res = FALSE)

## Compare with smeltz2023.uc to make sure the subset of three compounds 
## matches what's in the full data
uc.sub <- smeltz2023.uc[smeltz2023.uc$DTXSID %in% uc.list, ]

## Check the dimension
dim(uc.sub)
dim(fup_uc_L2)

colnames(uc.sub)
colnames(fup_uc_L2)
common.cols <- intersect(colnames(uc.sub),colnames(fup_uc_L2))
all(colnames(uc.sub[,common.cols]) == colnames(fup_uc_L2[,common.cols]))

## The order of the samples seems to be different
## Sort the data by lab sample name and DTXSID (sample name along is not a unique identifier)
og_level2 <- uc.sub[with(uc.sub, order(Lab.Sample.Name, DTXSID)), ]
ex_level2 <- fup_uc_L2[with(fup_uc_L2, order(Lab.Sample.Name, DTXSID)), ]

all(og_level2[,common.cols] == ex_level2[,common.cols])
for(i in common.cols){
  test <- all(og_level2[,i] == ex_level2[,i])
  if(test == FALSE | is.na(test)){
    print(i)
  }
}

## Address the discrepancies one by one

## Discrepancies found in the Lab.Compound.Name column:
## The original data uses compound name for lab compound id (line 355 in MSdata-Mar2022.R).
## I will keep my input.

## Discrepancies found in the Area and ISTD.Area columns:
## Area and ISTD.Area were rounded to 5 significant digits in format_fup_uc (lines 406-407).
## Compare the columns in 5 significant digits: 
og_level2$Area <- signif(og_level2[,"Area"], 5)
## There are some rows with missing peak area, replace NA with 0
og_level2$Area[is.na(og_level2$Area)] <- 0
ex_level2$Area[is.na(ex_level2$Area)] <- 0
all(og_level2[,"Area"] == ex_level2[,"Area"])

og_level2$ISTD.Area <- signif(og_level2[,"ISTD.Area"], 5)
all(og_level2[,"ISTD.Area"] == ex_level2[,"ISTD.Area"])

## Discrepancies found in the Analysis.Parameters column:
## Analysis.Parameters is a numeric column in smeltz2023.uc while it is a character
## column in the example data.
## Convert Analysis.Parameters to numeric column and check again:
ex_level2$Analysis.Parameters <- as.numeric(ex_level2$Analysis.Parameters)
## There are some rows with missing values, replace NA with 0
og_level2$Analysis.Parameters[is.na(og_level2$Analysis.Parameters)] <- 0
ex_level2$Analysis.Parameters[is.na(ex_level2$Analysis.Parameters)] <- 0
all(ex_level2$Analysis.Parameters == og_level2$Analysis.Parameters)

## Discrepancies found in the Note column:
## The original data maps the replicate column to the note column (line 358 in MSdata-Mar2022.R).
## I will keep my input (i.e. filled with "").

## Discrepancies found in the Response column:
## Differences in this column are most likely caused by rounding errors since 
## the same issue happened with area and istd.area, and responses are calculated from them.
diffs <- og_level2[,"Response"]- ex_level2[,"Response"]
summary(diffs[!is.na(diffs)])

## The maximum difference is 1.000e-03 and the minimum difference is -1.000e-02.
## The difference in two decimal places is not reasonable.
## After further investigation, the sample that have the largest difference has 
## a peak area of 174260 and ISTD area of 17279. With 5 significant digits the 
## numbers are rounded to whole number and even the tens. 

## To preserve the precision, use the area and ISTD area from Level-0 to calculate response
## and compare with response from the original data set 
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

## Compare the columns with different column names  
colnames(og_level2)[which(!(colnames(og_level2) %in% common.cols))]
colnames(ex_level2)[which(!(colnames(ex_level2) %in% common.cols))]

## As mentioned above, the original data maps replicate to the Note column,
## and creates a series column filled with 1 (line 261 in MSdata-Mar2022.R).
## Should compare Biological.Replicates from the example data to Note in the original data 
all(ex_level2[, "Biological.Replicates"] == og_level2[, "Note"])

## Check the concentration columns:
## Replace missing values with 0
ex_level2[is.na(ex_level2$Test.Compound.Conc), "Test.Compound.Conc"] <- 0
og_level2[is.na(og_level2$Standard.Conc), "Standard.Conc"] <- 0
## The differences, if any, should be smaller than a reasonable tolerance. Here uses four decimal places.
diffs <- ex_level2[, "Test.Compound.Conc"] - og_level2[, "Standard.Conc"]
all(abs(diffs) <= 1e-4)

## All columns are checked and any differences are documented

## Save level-0 to level-2 data to use for function demo/example documentation 
save(fup_uc_L0, fup_uc_L1, fup_uc_L2, file = "~/invitrotkstats/invitroTKstats/data/Fup-UC-example.RData")

## Include session info
utils::sessionInfo()
