## This R script should be ran from the command line using
## R CMD BATCH data-raw/building-fup-red-data.R

## Script used to create data examples for fup-RED to use for function documentations
## and vignette. The data examples should be small subsets of three compounds.

## load necessary package
library(readxl)
library(invitroTKstats)

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
chem.ids <- read_excel("~/invitrotkstats/invitroTKstats/data-raw/Smeltz-RED/PFAS LC-MS RED Summary 20220709.xlsx", sheet=1, skip=1)[1:29,1:2]
chem.ids <- as.data.frame(chem.ids)
chem.ids <- subset(chem.ids, !duplicated(chem.ids[,2]))
## In this table, the chemical names and their lab IDs are in the same column 
## Extract them into two separate columns
chem.ids$Compound <- unlist(lapply(strsplit(chem.ids[,2]," \\("),function(x) x[[1]])) 
chem.ids$Chem.Lab.ID <- gsub(")", "", unlist(lapply(strsplit(chem.ids[,2]," \\("),function(x) x[[2]])))

## Read in level-0 file
## Prepare a data guide for merge_level0 
this.file <- "PFAS LC-MS RED Summary 20220709.xlsx"

data.guide <- create_catalog(
  file = rep(this.file, 3),
  sheet = rep("RED1 Raw", 3),
  skip.rows = c(278,550,822),
  date = rep("022322",3),
  compound = c("PFPeA", "PFBS", "PFOA"),
  istd = c("M5PFPeA", "M3PFBS", "M8PFOA"),
  sample = rep("Name",3),
  type = rep("Type", 3),
  peak = rep("Area", 3),
  istd.peak = rep("IS Area", 3),
  conc = rep("Std. Conc", 3),
  analysis.param = rep("RT",3),
  num.rows = rep(268, 3),
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
             additional.colnames = "Sample Text", "nM",
             additional.colname.cols = "SampleText.ColName",
             chem.ids = chem.ids,
             output.res = FALSE,
             catalog.out = FALSE,
             INPUT.DIR = "~/invitrotkstats/invitroTKstats/data-raw/Smeltz-RED/")

## There are some additional columns needed for fup_red_L0 to go to level-1.
## But these columns do not exist in the original data file and  
## currently cannot be handled/added by additional utility functions. 
## Need to manually add them in. Following the steps in MS-REDdata-Apr2023.Rmd,
## this markdown can also be found under "working/SmeltzPFAS".

## Set reasonable precision for numeric columns
for (this.col in c("Peak.Area", "Compound.Conc", "ISTD.Peak.Area"))
  fup_red_L0[,this.col] <- signif(as.numeric(fup_red_L0[,this.col]),6)

## Create the Sample Type column, use the package annotations
fup_red_L0 <- subset(fup_red_L0,!is.na(fup_red_L0[,"Sample Text"]))
fup_red_L0[regexpr("CC",fup_red_L0[,"Sample Text"])!=-1,
           "Sample.Type"] <- "CC"
fup_red_L0[regexpr("-Pl-",fup_red_L0[,"Sample Text"])!=-1,
           "Sample.Type"] <- "Plasma"
fup_red_L0[regexpr("-S-",fup_red_L0[,"Sample Text"])!=-1,
           "Sample.Type"] <- "PBS"
fup_red_L0[regexpr("-EC1-",fup_red_L0[,"Sample Text"])!=-1,
           "Sample.Type"] <- "EC1"                   
fup_red_L0[regexpr("-EC2-",fup_red_L0[,"Sample Text"])!=-1,
           "Sample.Type"] <- "EC2"
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

## Remove samples with missing peak areas for analyte and internal standard
fup_red_L0 <- subset(fup_red_L0, !is.na(Peak.Area) & 
                       !is.na(fup_red_L0[,"ISTD.Peak.Area"]))

## Create the dilution factors column
## Information found in lines 177 to 188 in MS-REDdata-Apr2023.Rmd
fup_red_L0[sapply(fup_red_L0[,"Sample.Type"],function(x) "CC" %in% x),"Dilution.Factor"] <- 1
fup_red_L0[sapply(fup_red_L0[,"Sample.Type"],function(x) "PBS" %in% x),"Dilution.Factor"] <- 2
fup_red_L0[sapply(fup_red_L0[,"Sample.Type"],function(x) "EC1" %in% x),"Dilution.Factor"] <- 10
fup_red_L0[sapply(fup_red_L0[,"Sample.Type"],function(x) "EC2" %in% x),"Dilution.Factor"] <- 10
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
                             lab.compound.col="Compound",
                             type.col="Sample.Type",
                             dilution.col="Dilution.Factor",
                             biological.replicates.col ="Replicate",
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
fup_red_L2 <- fup_red_L1
fup_red_L2$Verified <- "Y"

## Compare with smeltz2023.red to make sure the subset of three compounds 
## matches what's in the full data
red.sub <- smeltz2023.red[smeltz2023.red$DTXSID %in% red.list, ]

## Compare the the dimensions
dim(fup_red_L2)
dim(red.sub)

## check all the columns with the same name have matching values
common_cols <- intersect(colnames(fup_red_L2),colnames(red.sub))
## subset to only common columns
subred.sub <- red.sub[,common_cols]
subfup_red_L2 <- fup_red_L2[,common_cols]
## check to see if the subset data has all the same information
for(i in common_cols){
  test <- all(subred.sub[,i] == subfup_red_L2[,i])
  if(test == FALSE | is.na(test)){
    print(i)
  }
}

## Discrepancies in columns "Date", "Note" and "Time":
## The Date column in the original (smeltz2023.red) data set might use the 
## acquired date from the Excel file, while I use the date in the lab sample names. 
## The format of the date is also different. In smeltz2023.red the date is in 
## "yyyy-mm-dd" format while I use "mmddyy" format.

## For the Note column, smeltz2023.red uses 'NA' to fill the column while I used "".

## Time column contains missing values so the check with all() will return NA 
## Use other method to check for this column
table(subred.sub[,"Time"] == subfup_red_L2[,"Time"])

## Compare the columns with different column names  
colnames(fup_red_L2)[which(!(colnames(fup_red_L2) %in% common_cols))]
colnames(red.sub)[which(!(colnames(red.sub) %in% common_cols))]
## Replicate columns contain missing values too, check with all() returns NA
all(fup_red_L2[,"Biological.Replicates"] == red.sub[,"Replicate"])
table(fup_red_L2[,"Biological.Replicates"] == red.sub[,"Replicate"])

## This check returns FALSE
all(fup_red_L2[,"Test.Compound.Conc"] == red.sub[,"Std.Conc"])
## Investigate the differences
diffs <- fup_red_L2[,"Test.Compound.Conc"] - red.sub[,"Std.Conc"]
summary(diffs[!is.na(diffs)])

## The maximum difference is 5.551e-17 and the minimum difference is -5.551e-17
## These differences could be due to different systems or are rounding errors

## Set to the same precision and check again
fup_red_L2[,"Test.Compound.Conc"] <- signif(as.numeric(fup_red_L2[,"Test.Compound.Conc"]),6)
red.sub[,"Std.Conc"] <- signif(as.numeric(red.sub[,"Std.Conc"]),6)
diffs <- fup_red_L2[,"Test.Compound.Conc"] - red.sub[,"Std.Conc"]
## Replace missing values with 0
diffs[is.na(diffs)] <- 0
## The differences, if any, should be smaller than a reasonable tolerance. Here uses four decimal places.
all(abs(diffs) <= 1e-4)

## All columns are checked and any differences are documented

## Save level-0 and level-1 data to use for function demo/example documentation 
save(fup_red_L0, fup_red_L1, fup_red_L2, file = "~/invitrotkstats/invitroTKstats/data/Fup-RED-example.RData")

## Include session info
utils::sessionInfo()
