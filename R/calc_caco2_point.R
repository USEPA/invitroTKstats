#' Calculate a Point Estimate of Apparent Membrane Permeability (Papp) from Caco-2 data (Level-3)
#'
#' This function calculates a point estimate of apparent membrane permeability (Papp)
#' using mass spectrometry (MS) peak areas from samples collected as part of in 
#' vitro measurements of membrane permeability using Caco-2 cells \insertCite{hubatsch2007determination}{invitroTKstats}.
#' 
#' The input to this function should be "level-2" data. Level-2 data is level-1,
#' data formatted with the \code{\link{format_caco2}} function, and curated
#' with a verification column. "Y" in the verification column indicates the
#' data row is valid for analysis. 
#'
#' The data frame of observations should be annotated according to direction
#' (either apical to basolateral -- "AtoB" -- or basolateral to apical -- "BtoA") and type
#' of concentration measured:
#' \tabular{rr}{
#'   Blank with no chemical added \tab Blank \cr
#'   Target concentration added to donor compartment at time 0 (C0) \tab D0\cr
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
#' \eqn{\frac{dQ}{dt} = \max\left(0, \frac{\sum_{i=1}^{n_{R2}} (r_{R2} * c_{DF})}{n_{R2}} - \frac{\sum_{i=1}^{n_{BL}} (r_{BL} * c_{DF})}{n_{BL}}\right)}
#'
#' where \eqn{r_{R2}} is Receiver Response, \eqn{c_{DF}} is the corresponding Dilution Factor, \eqn{r_{BL}} is Blank Response,
#' \eqn{n_{R2}} is the number of Receiver Responses, and \eqn{n_{BL}} is the number of Blank Responses.
#' 
#' If the output level-3 result table is chosen to be exported and an output 
#' directory is not specified, it will be exported to the user's R session
#' temporary directory. This temporary directory is a per-session directory 
#' whose path can be found with the following code: \code{tempdir()}. For more 
#' details, see \url{https://www.collinberke.com/til/posts/2023-10-24-temp-directories/}.
#' 
#' As a best practice, \code{INPUT.DIR} (when importing a .tsv file) and/or \code{OUTPUT.DIR} should be 
#' specified to simplify the process of importing and exporting files. This 
#' practice ensures that the exported files can easily be found and will not be 
#' exported to a temporary directory. 
#'
#' @param FILENAME (Character) A string used to identify the input level-2 file,
#' "<FILENAME>-Caco-2-Level2.tsv" (if importing from a .tsv file), and/or used 
#' to identify the output level-3 file, "<FILENAME>-Caco-2-Level3.tsv" (if exporting).
#' 
#' @param data.in (Data Frame) A level-2 data frame generated from the 
#' \code{format_caco2} function with a verification column added by 
#' \code{sample_verification}. Complement with manual verification if needed.
#' 
#' @param good.col (Character) Column name indicating which rows have been
#' verified, data rows valid for analysis are indicated with a "Y".
#' (Defaults to "Verified".)
#' 
#' @param output.res (Logical) When set to \code{TRUE}, the result 
#' table (level-3) will be exported to the user's per-session temporary 
#' directory or \code{OUTPUT.DIR} (if specified) as a .tsv file. 
#' (Defaults to \code{FALSE}.)
#' 
#' @param sig.figs (Numeric) The number of significant figures to round the exported result table (level-3). 
#' (Note: console print statements are also rounded to specified significant figures.) 
#' (Defaults to \code{3}.)
#' 
#' @param INPUT.DIR (Character) Path to the directory where the input level-2 file exists. 
#' If \code{NULL}, looking for the input level-2 file in the current working
#' directory. (Defaults to \code{NULL}.)
#' 
#' @param OUTPUT.DIR (Character) Path to the directory to save the output file. 
#' If \code{NULL}, the output file will be saved to the user's per-session temporary
#' directory or \code{INPUT.DIR} if specified. (Defaults to \code{NULL}.)
#' 
#' @return \item{data.frame}{A level-3 data.frame in standardized format}
#' \tabular{rrr}{
#'   C0_A2B \tab Time zero donor concentration \tab Mass Spec Response Ratio (RR) \cr
#'   dQdt_A2B \tab Estimated rate of mass movement through membrane \tab RR*cm^3/s \cr
#'   Papp_A2B \tab Apparent membrane permeability \tab 10^-6 cm/s\cr
#'   C0_B2A \tab Time zero donor concentration \tab Mass Spec Response Ratio (RR) \cr
#'   dQdt_B2A \tab Estimated rate of mass movement through membrane \tab RR*cm^3/s \cr
#'   Papp_B2A \tab Apparent membrane permeability \tab 10^-6 cm/s\cr
#'   Refflux \tab Efflux ratio \tab unitless\cr
#'   Frec_A2B.vec \tab Fraction recovered for the apical-basolateral direction, calculated as the fraction of the initial donor amount recovered in the receiver compartment (collapsed numeric vector, values for replicates separated by a "|") \tab unitless \cr
#'   Frec_A2B.mean \tab Mean of the fraction recovered for the apical-basolateral direction \tab unitless \cr
#'   Frec_B2A.vec \tab Fraction recovered for the basolateral-apical direction, calculated in the same way as Frec_A2B.vec but in the opposite transport direction (collapsed numeric vector, values for replicates separated by a "|") \tab unitless \cr 
#'   Frec_B2A.mean \tab Mean of the fraction recovered for the basolateral-apical direction \tab unitless \cr
#'   Recovery_Class_A2B.vec \tab Recovery classification for apical-to-basolateral permeability("Low Recovery" if Frec_A2B.vec < 0.4 or "High Recovery" if Frec_A2B.vec > 2.0) (collapsed character vector, values for replicates separated by a "|") \tab qualitative category \cr
#'   Recovery_Class_A2B.mean \tab Recovery classification for the mean apical-to-basolateral permeability("Low Recovery" if Frec_A2B.mean < 0.4 or "High Recovery" if Frec_A2B.mean > 2.0) \tab qualitative category \cr
#'   Recovery_Class_B2A.vec \tab Recovery classification for basolateral-to-apical permeability("Low Recovery" if Frec_B2A.vec < 0.4 or "High Recovery" if Frec_B2A.vec > 2.0) (collapsed character vector, values for replicates separated by a "|") \tab qualitative category \cr
#'   Recovery_Class_B2A.mean \tab Recovery classification for the mean basolateral-to-apical permeability("Low Recovery" if Frec_B2A.mean < 0.4 or "High Recovery" if Frec_B2A.mean > 2.0) \tab qualitative category \cr
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
#' ## import level-2 data from a 'tsv' file and export the result table to 
#' ## same location as INPUT.DIR 
#' \dontrun{
#' ## Refer to sample_verification help file for how to export level-2 data to a directory.
#' ## Unless a different path is specified in OUTPUT.DIR,
#' ## the result table will be saved to the directory specified in INPUT.DIR.
#' ## Will need to replace FILENAME and INPUT.DIR with name prefix and location of level-2 'tsv'.
#' level3 <- calc_caco2_point(# e.g. replace with "Examples" from "Examples-Caco-2-Level2.tsv" 
#'                            FILENAME="<level-2 FILENAME prefix>", 
#'                            INPUT.DIR = "<level-2 FILE LOCATION>",
#'                            output.res = TRUE)
#' }
#' 
#' ## scenario 3: 
#' ## input level-2 data from the R session and export the result table to the 
#' ## user's temporary directory
#' ## Will need to replace FILENAME with desired level-2 filename prefix. 
#' \dontrun{
#' level3 <- calc_caco2_point(# e.g. replace with "MYDATA"
#'                            FILENAME = "<desired level-2 FILENAME prefix>",
#'                            data.in = level2,
#'                            output.res = TRUE)
#' # To delete, use the following code. For more details, see the link in the 
#' # "Details" section. 
#' file.remove(list.files(tempdir(), full.names = TRUE, 
#' pattern = "<desired level-2 FILENAME prefix>-Caco-2-Level3.tsv"))
#' }
#' 
#' @references
#' \insertRef{hubatsch2007determination}{invitroTKstats}
#'
#' @import Rdpack
#' @importFrom utils read.csv write.table
#'
#' @export calc_caco2_point
calc_caco2_point <- function(
    FILENAME, 
    data.in,
    good.col="Verified", 
    output.res=FALSE, 
    sig.figs = 3,
    INPUT.DIR=NULL,
    OUTPUT.DIR = NULL)
{
  # These are the required data types as indicated by type.col.
  # In order to calculate the parameter a chemical must have peak areas for each
  # of these measurements:
  req.types=c("Blank","D0","D2","R2")
  
  #assigning global variables
  Compound.Name <- Response <- Sample.Type <- Direction <- NULL

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
                  test.conc.col="Test.Compound.Conc",
                  test.nominal.conc.col="Test.Nominal.Conc",
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
                                 C0_A2B = NaN, dQdt_A2B=NaN, Papp_A2B=NaN, Frec_A2B.vec=NaN, Frec_A2B.mean = NaN,
                                 Recovery_Class_A2B.vec = NA, Recovery_Class_A2B.mean = NA,
                                 C0_B2A = NaN, dQdt_B2A=NaN, Papp_B2A=NaN, Frec_B2A.vec=NaN, Frec_B2A.mean = NaN,
                                 Recovery_Class_B2A.vec = NA, Recovery_Class_B2A.mean = NA,
                                 Refflux=NaN))
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
        if (length(unique(this.dosing$Dilution.Factor))>1){
          stop("calc_caco2_point - There is more than one `Dilution.Factor` for `D0` samples of `",this.chem,"` in direction ",this.direction,".")
          # browser()
        } 
        this.row[paste("C0",dir.string,sep="_")] <- max(0,
                                                        mean(this.dosing$Response * this.dosing$Dilution.Factor) -
                                                        mean(this.blank$Response * this.blank$Dilution.Factor)) # [C0] = Peak area (RR) 
        
        # Calculate dQ/dt
        # only can handle one dilution factor and one receiver volume right now:
        if (length(unique(this.receiver$Dilution.Factor))>1 | 
            length(unique(this.receiver$Vol.Receiver))>1 |
            length(unique(this.dosing$Time))>1){
          stop("calc_caco2_point - `Dilution.Factor`, `Vol.Receiver`, and/or `Time` has more than one unique value for `",this.chem,"` in direction ",this.direction,".")
          # browser()
        } 
        this.row[paste("dQdt",dir.string,sep="_")] <- max(0,
                                                          (mean(this.receiver$Response * this.receiver$Dilution.Factor) -
                                                             mean(this.blank$Response * this.blank$Dilution.Factor)) * # Peak area (RR)
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
            length(unique(this.receiver$Vol.Receiver))>1){
          
          stop("calc_caco2_point - There is more than one `Dilution.Factor` for 'D0', 'D2', or 'R2' samples, or more than one `Vol.Receiver` for 'R2' samples for `",this.chem,"` in direction ",this.direction,".")
          # browser()
        }
        this.Frec.vec <- base::pmax(0,
                             (this.donor$Vol.Donor*(this.donor$Dilution.Factor)*(this.donor$Response-rep(mean(this.blank$Response),
                                                                                                         length(this.donor$Response)))+this.receiver$Vol.Receiver*(this.receiver$Dilution.Factor)*
                                (this.receiver$Response-rep(mean(this.blank$Response),
                                                            length(this.receiver$Response))))/(this.dosing$Vol.Donor*(this.dosing$Dilution.Factor)*
                                                                                                 (this.dosing$Response-rep(mean(this.blank$Response),
                                                                                                                           length(this.dosing$Response)))))
        this.mean.Frec.vec <- mean(this.Frec.vec)
        this.row[paste0("Frec_",dir.string,".vec")] <- paste(this.Frec.vec, collapse = "|")
        this.row[paste0("Frec_",dir.string,".mean")] <- as.numeric(this.mean.Frec.vec)
        
        #Calculate Recovery Class
        this.Recovery.vec = rep(NA, length(this.Frec.vec))
        for (i in 1:length(this.Frec.vec)){
          if (this.Frec.vec[i] < 0.4) this.Recovery.vec[i] <- "Low Recovery"
          else if (this.Frec.vec[i] > 2) this.Recovery.vec[i] <- "High Recovery"
        }
        if (any(!is.na(this.Recovery.vec))) this.row[paste0("Recovery_Class_",dir.string,".vec")] <- paste(this.Recovery.vec, collapse = "|")
        
        if (this.mean.Frec.vec < 0.4) this.row[paste0("Recovery_Class_",dir.string,".mean")] <- "Low Recovery"
        else if (this.mean.Frec.vec > 2) this.row[paste0("Recovery_Class_",dir.string,".mean")] <- "High Recovery"
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
      if (!is.null(sig.figs)) {
        print(paste(this.row$Compound.Name,"Refflux =",
                    signif(this.row$Refflux,sig.figs)))
      } else {
        print(paste(this.row$Compound.Name,"Refflux =",
                    this.row$Refflux))
      }
  }
  
  rownames(out.table) <- make.names(out.table$Compound.Name, unique=TRUE)
  out.table[,"C0_A2B"] <- as.numeric(out.table[,"C0_A2B"])
  out.table[,"C0_B2A"] <- as.numeric(out.table[,"C0_B2A"])
  out.table[,"dQdt_A2B"] <- as.numeric(out.table[,"dQdt_A2B"])
  out.table[,"dQdt_B2A"] <- as.numeric(out.table[,"dQdt_B2A"])
  out.table[,"Papp_A2B"] <- as.numeric(out.table[,"Papp_A2B"])
  out.table[,"Papp_B2A"] <- as.numeric(out.table[,"Papp_B2A"])
  out.table[,"Refflux"] <- as.numeric(out.table[,"Refflux"])
  out.table[,"Frec_A2B.mean"] <- as.numeric(out.table[,"Frec_A2B.mean"])
  out.table[,"Frec_B2A.mean"] <- as.numeric(out.table[,"Frec_B2A.mean"])
  
  out.table[,"Frec_A2B.vec"] <- as.character(out.table[,"Frec_A2B.vec"])
  out.table[,"Frec_B2A.vec"] <- as.character(out.table[,"Frec_B2A.vec"])
  out.table[,"Recovery_Class_A2B.vec"] <- as.character(out.table[,"Recovery_Class_A2B.vec"])
  out.table[,"Recovery_Class_B2A.vec"] <- as.character(out.table[,"Recovery_Class_B2A.vec"])
  out.table[,"Recovery_Class_A2B.mean"] <- as.character(out.table[,"Recovery_Class_A2B.mean"])
  out.table[,"Recovery_Class_B2A.mean"] <- as.character(out.table[,"Recovery_Class_B2A.mean"])
  out.table <- as.data.frame(out.table)
  
  if (output.res) {
    # Determine the path for output
    if (!is.null(OUTPUT.DIR)) {
      file.path <- OUTPUT.DIR
    } else if (!is.null(INPUT.DIR)) {
      file.path <- INPUT.DIR
    } else {
      file.path <- tempdir() 
    }
    
    rounded.out.table <- out.table 
    
    # Round results to desired number of sig figs
    if (!is.null(sig.figs)){
      rounded.out.table[,"C0_A2B"] <- signif(rounded.out.table[,"C0_A2B"], sig.figs)
      rounded.out.table[,"C0_B2A"] <- signif(rounded.out.table[,"C0_B2A"], sig.figs)
      rounded.out.table[,"dQdt_A2B"] <- signif(rounded.out.table[,"dQdt_A2B"], sig.figs)
      rounded.out.table[,"dQdt_B2A"] <- signif(rounded.out.table[,"dQdt_B2A"], sig.figs)
      rounded.out.table[,"Papp_A2B"] <- signif(rounded.out.table[,"Papp_A2B"], sig.figs)
      rounded.out.table[,"Papp_B2A"] <- signif(rounded.out.table[,"Papp_B2A"], sig.figs)
      rounded.out.table[,"Refflux"] <- signif(rounded.out.table[,"Refflux"], sig.figs)
      rounded.out.table[,"Frec_A2B.mean"] <- signif(rounded.out.table[,"Frec_A2B.mean"], sig.figs)
      rounded.out.table[,"Frec_B2A.mean"] <- signif(rounded.out.table[,"Frec_B2A.mean"], sig.figs)
      # split collapsed numeric vector and round according to sig.figs
      split_round = function(val){
        val = val %>% 
          base::strsplit(split = "\\|") %>% 
          unlist() %>% 
          as.numeric() %>% 
          signif(sig.figs) %>% 
          paste(collapse = "|")
      }
      rounded.out.table[,"Frec_A2B.vec"] <- base::sapply(rounded.out.table[,"Frec_A2B.vec"], split_round, USE.NAMES = F)
      rounded.out.table[,"Frec_B2A.vec"] <- base::sapply(rounded.out.table[,"Frec_B2A.vec"], split_round, USE.NAMES = F)
      cat(paste0("\nData to export has been rounded to ", sig.figs, " significant figures.\n"))
    }
    
    # Write out a "level 3" file (data organized into a standard format):
    write.table(rounded.out.table,
      file=paste0(file.path, "/", FILENAME,"-Caco-2-Level3.tsv"),
      sep="\t",
      row.names=F,
      quote=F)
   
    # Print notification message stating where the file was output to
    cat(paste0("A Level-3 file named ",FILENAME,"-Caco-2-Level3.tsv", 
               " has been exported to the following directory: ", file.path), "\n")
    
    
  }
  
  print(paste("Apical to basolateral permeability calculated for",num.a2b,"chemicals."))
  print(paste("Basolateral to apical permeability calculated for",num.b2a,"chemicals."))
  print(paste("Efflux ratio calculated for",num.efflux,"chemicals."))
  
  return(out.table)
}

