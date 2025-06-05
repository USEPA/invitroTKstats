##===========================##
## Global Requirement Script ##
##===========================##

## Deal with "no visible binding for global variable"
utils::globalVariables(
  c(
    "missing_cells","Test.Compound.Conc","L1.common.cols",
    "type.col","compound.col","dtxsid.col","time.col",
    "membrane.area.col","lab.compound.col","density.col",
    "cal.col","area.col","std.catcols","Verified"
  )
)