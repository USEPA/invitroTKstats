## code to prepare `std.catcols` dataset goes here
# standard catalog column names
# use: creating a catalog and/or checking a catalog structure.
std.catcols <- c(
  file.col = "File",
  sheet.col = "Sheet",
  skip.rows.col = "Skip.Rows",
  date.col = "Date",
  compound.col = "Chemical.ID",
  istd.col = "ISTD.Name",
  col.names.loc = "Col.Names.Loc",
  sample.colname.col = "Sample.ColName",
  type.colname.col = "Type.ColName",
  peak.colname.col = "Peak.ColName",
  istd.peak.colname.col = "ISTD.Peak.ColName",
  conc.colname.col = "Conc.ColName",
  analysis.param.colname.col = "AnalysisParam.ColName"
)
# create the .rda for the package
usethis::use_data(std.catcols, overwrite = TRUE)
# session information
sessionInfo()
