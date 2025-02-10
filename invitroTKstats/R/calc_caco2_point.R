#' Calculate a point estimate of apparent membrane permeability from Caco-2 data
#'
#' This function uses mass spectrometry (MS) peak areas
#' from samples collected as part of in vitro measurement of membrane
#' permeability using Caco-2 cells \insertCite{hubatsch2007determination}{invitroTKstats}.
#' Data are read from a "Level2" text file that should have been formatted and
#' created
#' by \code{\link{format_caco2}} (this is the "Level1" file). The Level1 file
#' should have been curated by adding a column with the value "Y" indicating
#' that each row is verified as usable for analysis (that is, the Level2 file).
#'
#' The data frame of observations should be annotated according to direction
#' (either apical to basal -- "AtoB" -- or basal to apical -- "BtoA") and type
#' of concentration measured:
#' \tabular{rr}{
#'   Blank with no chemical added \tab Blank \cr
#'   Dosing vehicle (C0) at target concentration \tab D0\cr
#'   Donor compartment at end of experiment \tab D2\cr
#'   Receiver compartment at end of experiment\tab R2\cr
#' }
#'
#' Apparent membrane permeability (\eqn{P_{app}}) is calculated from MS responses as:
#'
#'
#' \eqn{P_{app} = \frac{dQ/dt}{c_0*A}}
#'
#' The rate of permeation, \eqn{\frac{dQ}{dt}}\eqn{\left(\frac{\text{peak area}}{\text{time (s)}} \right)} is calculated as:
#'
#' \eqn{\frac{dQ}{dt} = \max\left(0, \frac{\sum_{i=1}^{n_P} (r_P * c_{DF})}{n_P} - \frac{\sum_{i=1}^{n_B} (r_B * c_{DF})}{n_B}\right)}
#'
#' where \eqn{r_P} is PBS Response, \eqn{c_{DF}} is Dilution Factor, \eqn{r_B} is Blank Response,
#' \eqn{n_P} is the number of PBS Responses, and \eqn{n_B} is the number of Blank Responses.
#'
#' @param FILENAME (Character) A string used to identify the input Level-2 file.
#' "<FILENAME>-Caco-2-Level2.tsv".
#' 
#' @param data.in (Data Frame) A Level-2 data frame generated from the 
#' \code{format_caco2} function with a verification column added by 
#' \code{sample_verification}. Complement with manual verification if needed.
#' 
#' @param good.col (Character) Column name indicating which rows have been
#' verified, data rows valid for analysis are indicated with a "Y".
#' (Defaults to "Verified".)
#' 
#' @param output.res (Logical) When set to \code{TRUE}, the result 
#' table (Level-3) will be exported the current directory as a .tsv file. 
#' (Defaults to \code{TRUE}.)
#' 
#' @param INPUT.DIR (Character) Path to the directory where the input level-2 file exists. 
#' If \code{NULL}, looking for the input level-2 file in the current working
#' directory. (Defaults to \code{NULL}.)
#' 
#' @param OUTPUT.DIR (Character) Path to the directory to save the output file. 
#' If \code{NULL}, the output file will be saved to the current working
#' directory or \code{INPUT.DIR} if specified. (Defaults to \code{NULL}.)
#' 
#' @return \item{data.frame}{A data.frame in standardized format}
#' \tabular{rrr}{
#'   C0_A2B \tab Time zero donor concentration \tab Mass Spec Response Ratio (RR) \cr
#'   dQdt_A2B \tab Estimated rate of mass movement through membrane \tab RR*cm^3/s \cr
#'   Papp_A2B \tab Apparent membrane permeability \tab 10^-6 cm/s\cr
#'   C0_B2A \tab Time zero donor concentration \tab Mass Spec Response Ratio (RR) \cr
#'   dQdt_B2A \tab Estimated rate of mass movement through membrane \tab RR*cm^3/s \cr
#'   Papp_B2A \tab Apparent membrane permeability \tab 10^-6 cm/s\cr
#'   Refflux \tab Efflux ratio \tab unitless\cr
#'   Frec_A2B.vec \tab Fraction recovered for the apical-basal direction, calculated as the fraction of the initial donor amount recovered in the receiver compartment \tab unitless \cr
#'   Frec_B2B.vec \tab Fraction recovered for the basal-apical direction, calculated in the same way as Frec_A2B.vec but in the opposite transport direction \tab unitless \cr 
#'   Recovr_Class_A2B \tab Recovery classification for apical-to-basal permeability("Low Recovery" or "High Recovery") \tb qualitative category \cr
#'   Recovr_Class_B2A \tab Recovery classification for basal-to-apical permeability("Low Recovery" or "High Recovery") \tb qualitative category \cr
#' }
#'
#' @author John Wambaugh
#'
#' @examples
#' ## Load example level-2 data
#' level2 <- invitroTKstats::caco2_L2
#' 
#' ## scenario 1: 
#' ## input level-2 data from the R session and do not export the result table
#' level3 <- calc_caco2_point(data.in = level2, output.res = FALSE)
#'
#' ## scenario 2: 
#' ## import level-2 data from a 'tsv' file and export the result table
#' \dontrun{
#' ## Refer to sample_verification help file for how to export level-2 data to a directory.
#' ## Unless a different path is specified in OUTPUT.DIR,
#' ## the result table will be saved to the directory specified in INPUT.DIR.
#' level3 <- calc_caco2_point(FILENAME="Examples", 
#'                            INPUT.DIR = "invitroTKstats/vignettes")
#' }
#'
#' @references
#' \insertRef{hubatsch2007determination}{invitroTKstats}
#'
#' @import Rdpack
#'
#' @export calc_caco2_point
calc_caco2_point <- function(
    FILENAME, 
    data.in,
    good.col="Verified", 
    output.res=TRUE, 
    INPUT.DIR=NULL,
    OUTPUT.DIR = NULL)
{
  # These are the required data types as indicated by type.col.
  # In order to calculate the parameter a chemical must have peak areas for each
  # of these measurements:
  req.types=c("Blank","D0","D2","R2")
  
  if (!missing(data.in)) {
    input.table <- as.data.frame(data.in)
  } else if (!is.null(INPUT.DIR)) {
    input.table <- read.csv(file=paste0(INPUT.DIR, "/", FILENAME,"-Caco-2-Level2.tsv"),
                            sep="\t",header=T)
  } else {
    input.table <- read.csv(file=paste0(FILENAME,"-Caco-2-Level2.tsv"),
                            sep="\t",header=T)
  }
  
  input.table <- subset(input.table,!is.na(Compound.Name))
  input.table <- subset(input.table,!is.na(Response))
  
  caco2.cols <- c(L1.common.cols, 
                  time.col = "Time",
                  direction.col="Direction",
                  compound.conc.col="Nominal.Conc",
                  nominal.test.conc.col="Test.Target.Conc",
                  membrane.area.col="Membrane.Area",
                  receiver.vol.col="Vol.Receiver",
                  donor.vol.col="Vol.Donor"
  )
  
  list2env(as.list(caco2.cols), envir = environment())
  cols <- c(unlist(mget(names(caco2.cols))), "Response", good.col)
  
  if (!any(c("Biological.Replicates", "Technical.Replicates") %in% colnames(input.table)))
    stop("Need at least one column representing replication, i.e. Biological.Replicates or Technical.Replicates. Run format_caco2 first (level 1) then curate to (level 2).")
  
  if (!(all(cols %in% colnames(input.table))))
  {
    warning("Run format_fup_red first (level 1) then curate to (level 2).")
    stop(paste("Missing columns named:",
               paste(cols[!(cols%in%colnames(input.table))],collapse=", ")))
  }
  
  # Only include the data types used:
  input.table <- subset(input.table,input.table[,type.col] %in% req.types)
  
  # Only used verfied data:
  input.table <- subset(input.table, input.table[,good.col] == "Y")
  
  out.table <-NULL
  num.a2b <- 0
  num.b2a <- 0
  num.efflux <- 0
  for (this.chem in unique(input.table[,compound.col]))
  {
    this.subset <- subset(input.table, input.table[,compound.col]==this.chem)
    this.dtxsid <- this.subset$dtxsid[1]
    this.row <- cbind(this.subset[1,
                                  c(compound.col, dtxsid.col, time.col, membrane.area.col)],
                      data.frame(Calibration="All Data",
                                 C0_A2B = NaN, dQdt_A2B=NaN, Papp_A2B=NaN, Frec_A2B=NaN,
                                 C0_B2A = NaN, dQdt_B2A=NaN, Papp_B2A=NaN, Frec_B2A=NaN, Refflux=NaN))
    for (this.direction in c("AtoB","BtoA"))
    {
      this.blank <- subset(this.subset, Sample.Type=="Blank" &
                             Direction==this.direction)
      this.dosing <- subset(this.subset, Sample.Type=="D0" &
                              Direction==this.direction)
      this.donor <- subset(this.subset,Sample.Type=="D2" &
                             Direction==this.direction)
      this.receiver <- subset(this.subset,Sample.Type=="R2" &
                                Direction==this.direction)
      
      # Check to make sure there are data for PBS and plasma:
      if (dim(this.blank)[1]> 0 &
          dim(this.dosing)[1] > 0 &
          dim(this.donor)[1] > 0 &
          dim(this.receiver)[1] > 0)
      {
        if (this.direction == "AtoB")
        {
          dir.string <- "A2B"
          num.a2b <- num.a2b + 1
        } else {
          dir.string <- "B2A"
          num.b2a <- num.b2a+1
        }
        
        # Calculate C0
        # only can handle one dilution factor right now:
        if (length(unique(this.dosing$Dilution.Factor))>1) browser()
        this.row[paste("C0",dir.string,sep="_")] <- max(0,
                                                        unique(this.dosing$Dilution.Factor)*(mean(this.dosing$Response) -
                                                                                               mean(this.blank$Response))) # [C0] = Peak area (RR) 
        
        # Calculate dQ/dt
        # only can handle one dilution factor and one receiver volume right now:
        if (length(unique(this.receiver$Dilution.Factor))>1 |
            length(unique(this.receiver$Vol.Receiver))>1 |
            length(unique(this.dosing$Time))>1) browser()
        this.row[paste("dQdt",dir.string,sep="_")] <- max(0,
                                                          (unique(this.receiver$Dilution.Factor)*
                                                             mean(this.receiver$Response) -
                                                             mean(this.blank$Response)) * # Peak area (RR)
                                                            unique(this.receiver$Vol.Receiver) / # cm^3
                                                            unique(this.receiver$Time) / 3600 #  1/h -> 1/s
        ) # [dQdt] = Peak area (RR) * cm^3 / s 
        
        # Calculate Papp
        this.row[paste("Papp",dir.string,sep="_")] <- max(0,
                                                          as.numeric(this.row[paste("dQdt",dir.string,sep="_")]) /  # Peak area (RR) * cm^3 / s 
                                                            as.numeric(this.row[paste("C0",dir.string,sep="_")]) / # Peak area (RR)
                                                            as.numeric(this.row["Membrane.Area"]) * # cm^ 2
                                                            1e6 # cm -> 10-6 cm 
        ) # [Papp] = cm^2/s
        
        #Calculate Recovery
        if (length(unique(this.donor$Dilution.Factor))>1 |
            length(unique(this.dosing$Dilution.Factor))>1 |
            length(unique(this.receiver$Dilution.Factor))>1 |
            length(unique(this.receiver$Vol.Receiver))>1) browser()
        this.row[paste("Frec",dir.string,sep="_")] <- max(0,
                                                          (this.donor$Vol.Donor*(this.donor$Dilution.Factor)*(this.donor$Response-rep(mean(this.blank$Response),
                                                           length(this.donor$Response)))+this.receiver$Vol.Receiver*(this.receiver$Dilution.Factor)*
                                                          (this.receiver$Response-rep(mean(this.blank$Response),
                                                           length(this.receiver$Response))))/(this.dosing$Vol.Donor*(this.dosing$Dilution.Factor)*
                                                          (this.dosing$Response-rep(mean(this.blank$Response),
                                                           length(this.dosing$Response)))))
      }
    }
    
    if (!is.nan(unlist(this.row["Papp_A2B"])) &
        !is.nan(unlist(this.row["Papp_B2A"])))
    {
      num.efflux <- num.efflux + 1
      this.row["Refflux"] <- as.numeric(this.row["Papp_B2A"]) /
        as.numeric(this.row["Papp_A2B"])
    }
    out.table <- rbind(out.table, this.row)
    print(paste(this.row$Compound.Name,"Refflux =",
                signif(this.row$Refflux,3)))
  }
  
  rownames(out.table) <- make.names(out.table$Compound.Name, unique=TRUE)
  out.table[,"C0_A2B"] <- signif(as.numeric(out.table[,"C0_A2B"]),3)
  out.table[,"C0_B2A"] <- signif(as.numeric(out.table[,"C0_B2A"]),3)
  out.table[,"dQdt_A2B"] <- signif(as.numeric(out.table[,"dQdt_A2B"]),3)
  out.table[,"dQdt_B2A"] <- signif(as.numeric(out.table[,"dQdt_B2A"]),3)
  out.table[,"Papp_A2B"] <- signif(as.numeric(out.table[,"Papp_A2B"]),3)
  out.table[,"Papp_B2A"] <- signif(as.numeric(out.table[,"Papp_B2A"]),3)
  out.table[,"Refflux"] <- signif(as.numeric(out.table[,"Refflux"]),3)
  out.table <- as.data.frame(out.table)
  
  # Create new columns to store recovery classification separately
  out_table$Recovr_Class_A2B=NA
  out_table$Recovr_Class_B2A=NA

  # Assign recovery classifications without changing Papp values
  out_table$Recovr_Class_A2B[out.table$Frec_A2B.vec < 0.4] <- "Low Recovery"
  out_table$Recovr_Class_A2B[out.table$Frec_A2B.vec > 2.0] <- "High Recovery"
  out_table$Recovr_Class_B2A[out.table$Frec_B2A.vec < 0.4] <- "Low Recovery"
  out_table$Recovr_Class_B2A[out.table$Frec_B2A.vec > 2.0] <- "High Recovery"
  
  # Calculate efflux ratio:
  out.table[,"Refflux"] <- signif(as.numeric(out.table[,"Refflux"]),3)
  
  if (output.res) {
    # Write out a "level 3" file (data organized into a standard format):
    # Determine the path for output
    
    if (!is.null(OUTPUT.DIR)) {
      file.path <- OUTPUT.DIR
    } else if (!is.null(INPUT.DIR)) {
      file.path <- INPUT.DIR
    } else {
      file.path <- getwd()
    }
    write.table(out.table,
                file=paste0(file.path, "/", FILENAME,"-Caco-2-Level3.tsv"),
                sep="\t",
                row.names=F,
                quote=F)
    
    # Print notification message stating where the file was output to
    cat(paste0("A Level-3 file named ",FILENAME,"-Caco-2-Level3.tsv", 
               " has been exported to the following directory: ", file.path), "\n")
  }
  
  print(paste("Apical to basal permeability calculated for",num.a2b,"chemicals."))
  print(paste("Basal to apical permeability calculated for",num.b2a,"chemicals."))
  print(paste("Efflux ratio calculated for",num.efflux,"chemicals."))
  
  return(out.table)
}

