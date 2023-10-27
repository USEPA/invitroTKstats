#' Creates a Standardized Data Table of Chemical Identities
#'
#' This function creates a data frame summarizing chemical identifiers used for 
#' each tested chemical in MS data. Each row in the resulting data frame provides 
#' EPA's DSSTox Structure ID (dtxsid), preferred compound name, and the name used 
#' by the laboratory.
#'
#' @param input.table (Data Frame) A data frame containing mass-spectrometry peak areas,
#' indication of chemical identity, and analytical chemistry methods.
#' It should contain columns with names specified by the following arguments:
#' 
#' @param dtxsid.col (Character) Column name of input.table containing EPA's DSSTox Structure 
#' ID (\url{http://comptox.epa.gov/dashboard}). (Defaults to "DTXSID".)
#' 
#' @param compound.col (Character) Column name of input.table containing the test compound.
#' (Defaults to "Compound.Name".)
#' 
#' @param lab.compound.col (Character) Column name of input.table containing the test compound 
#' name used by the laboratory. (Defaults to "Lab.Compound.Name".)
#' 
#' @return A data frame with one row per unique chemical, contains chemical 
#' identifiers including DTXSID, standard/preferred compound name, 
#' and the chemical name used by the laboratory for all chemicals in the input data frame.  
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