library(invitroTKstats)

level0 <- kreutz2020
level1 <- format_fup_uc(level0,
  FILENAME="Kreutz2020",
  compound.col="Name",
  compound.conc.col="Standard.Conc",
  area.col="Chem.Area"
  )
level2 <- level1
level2$Verified <- "Y"

write.table(level2,
  file="Kreutz2020-PPB-UC-Level2.tsv",
  sep="\t",
  row.names=F,
  quote=F)

level3 <- calc_fup_uc_point(FILENAME="Kreutz2020")
