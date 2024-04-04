## This R script should be ran from the command line using
## R CMD BATCH data-raw/building-Clint-data.R

## Script used to create data examples from the full Smeltz 2023
## Clint data to use for function documentations and vignette.
## The data examples are smaller subsets of three compounds.

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
## Make necessary adjustments if needed. 

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
  date = rep("2022-01-28", 3),
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
## Need to manually add them in. Following the steps in smeltz-hep-cc.R,
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
## The dilution factor is indeed different across different sample types.
## -	The calibration curves (all, regardless of analyte) are diluted 240x
## -	The ref cmpds (may have label “HLB”) and PFAS labeled as “WAX1” are diluted 480x
## -	PFAS marked as WAX2 are diluted 720x
## -	It seems that when you asked for a pared down sheet this information was lost. We could add it to the Analytes page. I can see it in the more comprehensive sheet she#   provided to me.

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
## Only set a standard/expected conc for calibration curve points
clint_L0[clint_L0[,"Type"]!="CC","Compound.Conc"] <- NA
## Concentrations calculated including dilution
clint_L0[,"Compound.Conc"] <- clint_L0[,"Compound.Conc"]*clint_L0$Dilution.Factor


## Prepare Level-1 
clint_L1 <- format_clint(data.in = clint_L0,
                       sample.col ="Sample",
                       date.col="Date",
                       compound.col="Compound",
                       lab.compound.col="Compound",
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

## First assign "Y" to all samples
clint_L2 <- sample_verification(data.in = clint_L1, assay = "Clint",output.res = FALSE)

## If the mix recorded in column "Note" does not match what's recorded on the 
## summary table (chem.ids), exclude the sample due to wrong mix.
for (this.id in clint.list)
{
  this.mix <- chem.ids[chem.ids$DTXSID==this.id,"Mix"]
  if (regexpr("WAX1",this.mix)!=-1)
  {
    clint_L2[clint_L2$DTXSID==this.id & clint_L2$Sample.Type =="Cvst" &
             regexpr("WAX2",clint_L2$Note)!=-1,"Verified"] <- "Wrong Mix"
  } else if (regexpr("WAX2",this.mix)!=-1)
  {
    clint_L2[clint_L2$DTXSID==this.id & clint_L2$Sample.Type =="Cvst" &
             regexpr("WAX1",clint_L2$Note)!=-1,"Verified"] <- "Wrong Mix"
  }
}
## Exclude calibration curve samples if the concentration is unknown.
clint_L2[clint_L2$Sample.Type=="CC" & is.na(clint_L2$Std.Conc),"Verified"] <- 
  "Unknown Conc."

## Compare with the level-2 data that's already in the package 
## see if the subset matches what's in the full data.
clint.sub <- smeltz2023.clint[smeltz2023.clint$DTXSID %in% clint.list, ]

## check if the columns are correct
names(clint.sub)
names(clint_L2)

## Compare some key parameters 
nrow(clint_L2) == nrow(clint.sub)
## This is FALSE because smeltz.clint 
## somehow uses lab compound ID/abbreviations for compound name
## instead of the full name, which is what I did when I create the subset
all(tolower(unique(clint.sub$Compound.Name)) %in% tolower(unique(clint_L2$Compound.Name)))
unique(clint.sub$Compound.Name)
unique(clint_L2$Compound.Name)

summary(clint.sub$Response)
summary(clint_L2$Response)

table(clint.sub$Sample.Type)
table(clint_L2$Sample.Type)

summary(clint.sub$Std.Conc)
summary(clint_L2$Test.Compound.Conc)

## Save level-0 to level-2 data to use for function demo/example documentation 
save(clint_L0, clint_L1, clint_L2, file = "~/invitrotkstats/invitroTKstats/data/Clint-example.RData")

## Include session info
utils::sessionInfo()

