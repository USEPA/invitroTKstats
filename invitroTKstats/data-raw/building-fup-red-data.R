## This R script should be ran from the command line using
## R CMD BATCH data-raw/building-fup-red-data.R

## Script used to create data examples from the full Smeltz 2023 
## fup red data (level-2) to use for function documentations
## and vignette. The data examples should be smaller subsets of three compounds.

## load necessary package
library(readxl)
library(invitroTKstats)

## Picked three compounds from the data that we know all samples are verified (with "Y").
red.list <- c("DTXSID6062599","DTXSID5030030","DTXSID8031865")

## Check all compounds have all samples verified with "Y".
## They are also from the same sheet in the original file. 
## Make them easier to work with later.
unique(smeltz2023.red[smeltz2023.red$DTXSID %in% red.list, "Verified"])
unique(smeltz2023.red[smeltz2023.red$DTXSID %in% red.list, "Level0.Sheet"])

## Prepare Level-0
## read in chem.ids
chem.ids <- read_excel("~/invitrotkstats/invitroTKstats/data-raw/Smeltz-RED/PFAS LC-MS RED Summary 20220709.xlsx", sheet=1, skip=1)[1:29,1:2]
chem.ids <- as.data.frame(chem.ids)
chem.ids <- subset(chem.ids, !duplicated(chem.ids[,2]))

## read in level-0 file
## merge_level0 is not needed here because all the relevant the data exist in one file and one sheet
this.file <- "~/invitrotkstats/invitroTKstats/data-raw/Smeltz-RED/PFAS LC-MS RED Summary 20220709.xlsx"
fup_red_L0 <- read_excel(this.file, sheet=3, skip=6)
this.sheet.name <- excel_sheets(this.file)[3]
fup_red_L0 <- as.data.frame(fup_red_L0)

## Identify the compound names and the corresponding DTXSID's
## The first compound was skipped. Added back in. 
this.compound <- "DTXSID00379268"
this.row <- 1
while(this.row <= dim(fup_red_L0)[1])
  {
  if (!is.na(fup_red_L0[this.row,1]))
    {
    if (regexpr("Compound",fup_red_L0[this.row,1])!=-1)
      {
      temp <- trimws(strsplit(fup_red_L0[this.row,1],": ")[[1]][2])
      this.compound <- unique(chem.ids[regexpr(paste0("\\(",temp,"\\)"),chem.ids[,2])!=-1,1])
      if (length(this.compound)==0) this.compound <- temp
      } else 
        {  
          fup_red_L0[this.row,"DTXSID"] <- this.compound
        }
    }
  this.row <- this.row + 1
}

## There are some additional columns needed for fup_red_L0 to go to level-1.
## But these columns cannot do not exist in the original data file and  
## currently cannot be handled/added by additional utility functions. 
## Need to manually add them in. 

## Remove the first two columns that are just row numbers, created from reading in from Excel
fup_red_L0 <- fup_red_L0[, -c(1,2)]

## Record the level-0 file name and sheet name
fup_red_L0$File <- this.file
fup_red_L0$Sheet <- this.sheet.name

## Remove samples that are not able to identify which compounds they belong 
fup_red_L0 <- subset(fup_red_L0,!is.na(DTXSID))

## Convert date and time to strings
fup_red_L0$Acq.Date <- as.Date(as.numeric(fup_red_L0$Acq.Date),origin = "1899-12-30")
h <- floor(as.numeric(fup_red_L0$Acq.Time)*24)
m <- as.character(floor((as.numeric(fup_red_L0$Acq.Time)*24-h)*60))
m[is.na(m)] <- "00"
m[sapply(m,nchar)==1] <- paste("0",m[sapply(m,nchar)==1],sep="")
fup_red_L0$Acq.Time <- paste(h,m,sep=":")

## Set reasonable precision for numeric columns
for (this.col in c("Area", "Height", "IS Area", "RT", "%Dev",
                   "Response", "Coeff. Of Determination", "Std. Conc", "nM"))
  fup_red_L0[,this.col] <- signif(as.numeric(fup_red_L0[,this.col]),6)

## Find matching compound names from chem.ids table
chem.ids$Compound <- unlist(lapply(strsplit(chem.ids[,2]," \\("),function(x) x[[1]])) 
for (this.id in unique(fup_red_L0$DTXSID))
{
  if (this.id %in% chem.ids$DTXSID)
  {
    fup_red_L0[fup_red_L0$DTXSID==this.id,"Compound"] <- 
      chem.ids[chem.ids$DTXSID==this.id,"Compound"] 
  } 
}  

## Match chemicals to their internal standards
for (this.chem in red.list)
{
  this.data <- subset(fup_red_L0,!is.na(fup_red_L0[,"IS Area"]) &
                        DTXSID==this.chem)
  this.istd <- fup_red_L0[sapply(fup_red_L0$Area,
                                 function(x) this.data[1,"IS Area"]%in%x), "DTXSID"]    
  fup_red_L0[fup_red_L0$DTXSID==this.chem,"ISTD"] <- this.istd              
}
## Remove internal standard data
fup_red_L0 <- subset(fup_red_L0, Compound != "ISTD")
## Only keep data for the three chosen compounds
fup_red_L0 <- fup_red_L0[fup_red_L0$DTXSID %in% red.list,]

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
fup_red_L0[fup_red_L0[,"Sample.Type"]%in%"NoPlasma.Blank","Area"] <- 0
fup_red_L0[fup_red_L0[,"Sample.Type"]%in%"NoPlasma.Blank","IS Area"] <- 1

# Remove samples with missing peak areas for analyte and internal standard
fup_red_L0[,"Std. Conc"] <- as.numeric(fup_red_L0[,"Std. Conc"])
fup_red_L0 <- subset(fup_red_L0, !is.na(Area) & 
                       !is.na(fup_red_L0[,"IS Area"]))

# Create the dilution factors column
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
fup_red_L0[,"Std. Conc"] <- as.numeric(fup_red_L0[,"Std. Conc"])/100

## Prepare level-1 data
fup_red_L1 <- format_fup_red(data.in = fup_red_L0,
                             sample.col ="Name",
                             date.col="Acq.Date",
                             compound.col="Compound",
                             lab.compound.col="Compound",
                             type.col="Sample.Type",
                             dilution.col="Dilution.Factor",
                             technical.replicates.col ="Replicate",
                             cal=1,
                             istd.conc = 10/1000,
                             istd.col= "IS Area",
                             istd.name.col = "ISTD", 
                             test.conc.col = "Std. Conc", 
                             level0.file.col = "File", 
                             level0.sheet.col = "Sheet",
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

## Compare some key parameters 
nrow(red.sub) == nrow(fup_red_L1)
all(unique(red.sub$Compound.Name) %in% unique(fup_red_L1$Compound.Name))

summary(red.sub$Response)
summary(fup_red_L1$Response)

table(red.sub$Sample.Type)
table(fup_red_L1$Sample.Type)

## Save level-0 and level-1 data to use for function demo/example documentation 
save(fup_red_L0, fup_red_L1, fup_red_L2, file = "~/invitrotkstats/invitroTKstats/data/Fup-RED-example.RData")

## Include session info
utils::sessionInfo()
