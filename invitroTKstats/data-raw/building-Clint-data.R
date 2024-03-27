## This R script should be ran from the command line using
## R CMD BATCH data-raw/building-clint-data.R

## Load necessary package
library(invitroTKstats)
library(readxl)

## smeltz2023.clint only has data for seven compounds.
## Unfortunately there's only one compound which has all samples verified with a "Y",
## other compounds all have some samples excluded from the analysis. 
## Need to go through couple verification steps from level-1 to level-2.

## Choose three compounds
clint.list <- c("DTXSID1021116", "DTXSID6023525", "DTXSID80380256")

## Prepare Level-0
## read in chem.ids
chem.ids <- read_excel("~/invitrotkstats/working/SmeltzPFAS/Hep12 Data for Uncertainty Feb2022.xlsx", sheet=1)
chem.ids <- as.data.frame(chem.ids)

## merge_level0 is not needed here because all the data are in the same Excel file.
## Reference chemicals data is in the second sheet and test chemicals data is in the third sheet. 
smeltz.hep.ref <- read_excel("~/invitrotkstats/working/SmeltzPFAS/Hep12 Data for Uncertainty Feb2022.xlsx", sheet=2, skip=6)
smeltz.hep.ref <- as.data.frame(smeltz.hep.ref)

smeltz.hep.pfas <- read_excel("~/invitrotkstats/working/SmeltzPFAS/Hep12 Data for Uncertainty Feb2022.xlsx", sheet=3, skip=2)
smeltz.hep.pfas <- as.data.frame(smeltz.hep.pfas)
## Match the column names
colnames(smeltz.hep.pfas) <- colnames(smeltz.hep.ref)

## Merge the reference chemicals and PFAS
clint_L0 <- rbind(smeltz.hep.ref,smeltz.hep.pfas)

## Remove blank rows
clint_L0 <- subset(clint_L0, !is.na(clint_L0[,1]))

## There are some additional columns needed for clint_L0 to go to level-1.
## But these columns cannot do not exist in the original data file and  
## currently cannot be handled/added by additional utility functions. 
## Need to manually add them in. 

## Record the level-0 file name and sheet name
clint_L0$Level0.File <- "Hep12 Data for Uncertainty Feb2022.xlsx"
clint_L0$Level0.Sheet <- "PFAS Data"

## Annotate Compounds
## The first compound name is in the rows that were skipped. 
## Manually annotate this one.
this.compound <- "Phenacetin"
this.dtxsid <- chem.ids[
  tolower(chem.ids[,"Analyte Name"]) == tolower(this.compound),
  "DTXSID"]
this.istdname <- chem.ids[
  tolower(chem.ids[,"Analyte Name"]) == tolower(this.compound),
  "Internal Standard"]

## Loop through the first column to identify compound names and find the matching 
## DTXSID from chem.ids
this.row <- 1
while (this.row <= dim(clint_L0)[1])
{
  if (regexpr("Compound",clint_L0[this.row,1])!=-1)
  {
    this.compound <- strsplit(clint_L0[this.row,1],"  ")[[1]][2] 
    if (tolower(this.compound) %in% tolower(chem.ids[,"Internal Standard"]))
    {
      this.istd <- this.compound
      this.compound <- "ISTD"
      this.dtxsid <- NA
    } else {
      row.index <- which(regexpr(tolower(this.compound),
                                 tolower(chem.ids[,"Analyte Name"]))!=-1)
      this.dtxsid <- chem.ids[row.index, "DTXSID"]
      this.istdname <- chem.ids[row.index, "Internal Standard"]
    }
  }
  clint_L0[this.row, "Compound"] <- this.compound
  clint_L0[this.row, "DTXSID"] <- this.dtxsid
  clint_L0[this.row, "ISTD.Name"] <- this.istdname
  this.row <- this.row + 1
}

## Remove rows with blank sample text
clint_L0 <- subset(clint_L0,!is.na(clint_L0[,"Sample Text"]))

## Remove internal standard data
clint_L0 <- subset(clint_L0, Compound != "ISTD")

## Only keep data for the three chosen compounds
clint_L0 <- clint_L0[clint_L0$DTXSID %in% clint.list,]

## Indicate whether hepatocytes have been inactivated:
clint_L0$Active.Hep <- NA
clint_L0[regexpr("living",tolower(clint_L0[,"Sample Text"]))!=-1,
           "Active.Hep"] <- 1
clint_L0[regexpr("inactive",tolower(clint_L0[,"Sample Text"]))!=-1,
           "Active.Hep"] <- 0

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
clint_L0$Area <- as.numeric(clint_L0$Area)
clint_L0[,"IS Area"] <- as.numeric(clint_L0[,"IS Area"])
clint_L0[,"Std. Conc"] <- as.numeric(clint_L0[,"Std. Conc"])
clint_L0 <- subset(clint_L0, !is.na(Area) & 
                       !is.na(clint_L0[,"IS Area"]))

## Convert nM standard concs to uM
clint_L0[,"nM"] <- as.numeric(clint_L0[,"nM"])/1000
## Only set a standard conc for calibration curve points
clint_L0[clint_L0[,"Type"]!="CC","nM"] <- NA
## Concentrations calculated including dilution
clint_L0[,"nM"] <- clint_L0[,"nM"]*clint_L0$Dilution.Factor

## Prepare Level-1 
clint_L1 <- format_clint(data.in = clint_L0,
                       sample.col ="Name",
                       date.col="Acq.Date",
                       compound.col="Compound",
                       lab.compound.col="Compound",
                       type.col="Type",
                       dilution.col="Dilution.Factor",
                       cal=1,
                       istd.conc = 10/1000,
                       istd.col= "IS Area",
                       density = 0.5,
                       clint.assay.conc = 1,
                       biological.replicates = 1,
                       test.conc.col="nM",
                       time.col = "Time",
                       analysis.method = "LCMS",
                       analysis.instrument = "Unknown",
                       analysis.parameters.col = "RT",
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
all(unique(clint.sub$Compound.Name) %in% unique(clint_L2$Compound.Name))

summary(clint.sub$Response)
summary(clint_L2$Response)

table(clint.sub$Sample.Type)
table(clint_L2$Sample.Type)

## Save level-0 and level-1 data to use for function demo/example documentation 
save(clint_L0, clint_L1, clint_L2, file = "~/invitrotkstats/invitroTKstats/data/Clint-example.RData")

## Include session info
utils::sessionInfo()

