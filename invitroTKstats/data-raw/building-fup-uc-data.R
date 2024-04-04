## This R script should be ran from the command line using
## R CMD BATCH data-raw/building-fup-uc-data.R

## Script used to create data examples from the full Smeltz 2023
## fup UC data to use for function documentations and vignette.
## The data examples are smaller subsets of three compounds.

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
## Make necessary adjustments if needed. 

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
## with the argument INPUT.DIR. Make necessary adjustments if needed.
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

## Make sure numeric columns are in the correct class and set a reasonable precison
fup_uc_L0[,"Compound.Conc"] <- signif(as.numeric(fup_uc_L0[,"Compound.Conc"]),4)
fup_uc_L0[,"Peak.Area"] <- as.numeric(unlist(fup_uc_L0[,"Peak.Area"]))
fup_uc_L0[,"ISTD.Peak.Area"] <- as.numeric(unlist(fup_uc_L0[,"ISTD.Peak.Area"]))
## The unit of compound concentration varies
## concentration of compound other than 8:2 FTS are in nM, need to convert the unit to uM which is what the package uses
fup_uc_L0[fup_uc_L0$Lab.Compound.ID != "8:2 FTS", "Compound.Conc"] <- fup_uc_L0[fup_uc_L0$Lab.Compound.ID != "8:2 FTS", "Compound.Conc"]/1000

fup_uc_L0[fup_uc_L0[,"Sample.Type"]!="CC","Compound.Conc"] <- NA 
## Create the dilution factor column
fup_uc_L0[fup_uc_L0[,"Sample.Type"]=="AF","Dilution.Factor"] <- 2*16
fup_uc_L0[fup_uc_L0[,"Sample.Type"]=="T1","Dilution.Factor"] <- 5*16
fup_uc_L0[fup_uc_L0[,"Sample.Type"]=="T5","Dilution.Factor"] <- 5*16
fup_uc_L0[fup_uc_L0[,"Sample.Type"]=="CC","Dilution.Factor"] <- 1
fup_uc_L0[fup_uc_L0[,"Sample.Type"]=="Blank","Dilution.Factor"] <- 1

# Treat the blanks as calibration data with concentration 0:
fup_uc_L0[fup_uc_L0[,"Sample.Type"]=="Blank","Compound.Conc"] <- 0
fup_uc_L0[fup_uc_L0[,"Sample.Type"]=="Blank","Sample.Type"] <- "CC"

# Remove CC samples that don't have a concentration:
fup_uc_L0 <- subset(fup_uc_L0,!(Sample.Type=="CC" & is.na(Compound.Conc)))

# Create a replicate column
fup_uc_L0[,"Replicate"] <- ""
fup_uc_L0[regexpr("_A",fup_uc_L0[,"Sample Text"])!=-1,"Replicate"] <- "A"
fup_uc_L0[regexpr("_B",fup_uc_L0[,"Sample Text"])!=-1,"Replicate"] <- "B"
fup_uc_L0[regexpr("_C",fup_uc_L0[,"Sample Text"])!=-1,"Replicate"] <- "C"

# get rid of mixes that don't contain analyte:
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

## Compare some key parameters 
nrow(uc.sub) == nrow(fup_uc_L2)
all(unique(uc.sub$Compound.Name) %in% unique(fup_uc_L2$Compound.Name))

summary(uc.sub$Response)
summary(fup_uc_L2$Response)

summary(uc.sub$Standard.Conc)
summary(fup_uc_L2$Test.Compound.Conc)

table(uc.sub$Sample.Type)
table(fup_uc_L2$Sample.Type)

## Save level-0 to level-2 data to use for function demo/example documentation 
save(fup_uc_L0, fup_uc_L1, fup_uc_L2, file = "~/invitrotkstats/invitroTKstats/data/Fup-UC-example.RData")

## Include session info
utils::sessionInfo()
