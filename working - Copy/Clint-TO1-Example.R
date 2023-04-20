library(invitroTKstats)

level0 <- TO1caco2
level1 <- format_caco2(level0,
   FILENAME="EPACyprotex2021",
   sample.col="SampleName",
   dtxsid.col="CompoundName",
   lab.compound.col="CompoundName",
   cal=1,
   istd.conc.col="ISTD.Conc",
   compound.col="CompoundName",
   compound.conc.col="Test.Target.Conc",
   membrane.area=0.11,
   series=1,
   analysis.parameters="Feature",
   analysis.instrument="GC or LC",
   analysis.method="Mass Spec"
   )
#'    
#' level2 <- level1
#' level2$Verified <- "Y"
#'  
#' write.table(level2,
#'   file="EPACyprotex2021-Caco-2-Level2.tsv",
#'   sep="\t",
#'   row.names=F,
#'   quote=F)
#'    
#' level3 <- calc_caco2_point(FILENAME="EPACyprotex2021")
#' 
#' write.table(level3,
#'   file="EPACyprotex2021-Caco-2-Level3.tsv",
#'   sep="\t",
#'   row.names=F,
#'   quote=F)