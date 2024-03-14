## code to prepare `L1.common.cols` dataset
## Standard names for columns that are expected for all assays.
## Used in format_<assay> functions to check & rename required columns.
L1.common.cols <- c(
  sample.col = "Lab.Sample.Name",
  date.col = "Date",
  compound.col = "Compound.Name",
  dtxsid.col = "DTXSID",
  lab.compound.col = "Lab.Compound.Name",
  type.col = "Sample.Type",
  dilution.col = "Dilution.Factor",
  cal.col = "Calibration",
  istd.name.col = "ISTD.Name",
  istd.conc.col = "ISTD.Conc",
  istd.col = "ISTD.Area",
  area.col = "Area",
  analysis.method.col = "Analysis.Method",
  analysis.instrument.col = "Analysis.Instrument",
  analysis.parameters.col = "Analysis.Parameters",
  note.col = "Note",
  level0.file.col = "Level0.File",
  level0.sheet.col = "Level0.Sheet"
)
# create the .rda for the package
usethis::use_data(L1.common.cols, overwrite = TRUE)
# session information
sessionInfo()