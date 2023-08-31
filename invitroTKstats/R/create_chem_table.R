#' Creates a standardized data table of chemical identities
#'
#' This function extracts the chemical analysis methods from a set of MS data.
#'
#' @param input.table A data frame containing mass-spectrometry peak areas,
#' indication of chemical identity, and analytical chemistry methods.
#' It should contain columns with names specified by the following arguments:
#' 
#' @param dtxsid.col Which column of input.table indicates EPA's DSSTox Structure 
#' ID (\url{http://comptox.epa.gov/dashboard}) (Defaults to "DTXSID")
#' 
#' @param compound.col Which column of input.table indicates the test compound
#' (Defaults to "Compound.Name")
#' 
#' @param lab.compound.col Which column of PPB.data indicates The test compound 
#' name used by the laboratory (Defaults to "Lab.Compound.Name")
#' 
#' @return A data.frame with one row per method 
#'
#' @author John Wambaugh
#'
#' @examples
#'
#' library(invitroTKstats)
#'
#' # Smeltz et al. (2020) data:
#' create_chem_table(smeltz2020)
#'
#' # Kreutz et al. (2020) data:
#' create_chem_table(kreutz2020,compound.col="Name")
#'
#' # Wambaugh et al. (2019) data:
#' # Strip out protein conc information from compound names:
#' wambaugh2019.red$CompoundName <- gsub("-100P","",wambaugh2019.red$CompoundName)
#' wambaugh2019.red$CompoundName <- gsub("-30P","",wambaugh2019.red$CompoundName)
#' wambaugh2019.red$CompoundName <- gsub("-10P","",wambaugh2019.red$CompoundName)
#'
#' wambaugh2019.redchems <- create_chem_table(wambaugh2019.red,
#'    compound.col="Preferred.Name",
#'    lab.compound.col="CompoundName")
#'
#' @export create_chem_table
create_chem_table <- function(input.table,
  dtxsid.col="DTXSID",
  compound.col="Compound.Name",
  lab.compound.col="Lab.Compound.Name"
  )
{
# We need all these columns in input.table
  cols <-c(
    compound.col,
    dtxsid.col,
    lab.compound.col
    )
  
  if (!(all(cols %in% colnames(input.table))))
  {
    stop(paste("Missing columns named:",
      paste(cols[!(cols%in%colnames(input.table))],collapse=", ")))
  }

  N.chems <- length(unique(input.table[,dtxsid.col]))
  
  out.table <- NULL
  for (this.chem in sort(unique(input.table[,dtxsid.col])))
  {
    this.subset <- subset(input.table,input.table[,dtxsid.col]==this.chem)
    this.row <- data.frame(
      Compound.Name=paste(unique(this.subset[,compound.col]),
        collapse=", "),
      DTXSID=this.chem,
      Lab.Compound.Name=paste(unique(this.subset[,lab.compound.col]),
        collapse=", "),
      stringsAsFactors=F
      )
    out.table <- rbind(out.table,this.row)
  }

  cat(paste(N.chems,"chemicals.\n"))
  
  return(out.table)
}