## This R script should be ran from the command line using
## R CMD BATCH data-raw/building-fup-uc-data.R

## Script used to create data examples for fup-uc to use for function documentations and vignette

library(invitroTKstats)
## Load a previously complied Level-0 file of Smeltz 2023 work
fup_uc_L0 <- read.table("~/invitrotkstats/invitroTKstats/data-raw/Smeltz-UC/SmeltzPFAS-PPB-UC-Level0.tsv",
                 sep = '\t', header = TRUE)

## Chose three compounds that have all samples verified
uc.list <- c("DTXSID5061954", "DTXSID50892417", "DTXSID8037708")
## Keep only the data of the chosen compounds
fup_uc_L0 <- fup_uc_L0[fup_uc_L0$DTXSID %in% uc.list, ]

## Run through level-1 processing function 
fup_uc_L1 <- format_fup_uc(data.in = fup_uc_L0,
                           sample.col="Name",
                           compound.col="Compound.Name",
                           test.conc.col ="Std.Conc", 
                           lab.compound.col="Compound.Name", 
                           type.col="Sample.Type", 
                           istd.col="IS.Area",
                           note.col="Replicate",
                           uc.assay.conc.col="Test.Target.Conc",
                           technical.replicates.col = "Replicate",
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

table(uc.sub$Sample.Type)
table(fup_uc_L2$Sample.Type)

## Save level-0 to level-2 data to use for function demo/example documentation 
save(fup_uc_L0, fup_uc_L1, fup_uc_L2, file = "~/invitrotkstats/invitroTKstats/data/Fup-UC-example.RData")

## Include session info
utils::sessionInfo()
