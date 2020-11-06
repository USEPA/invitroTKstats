UC_PPB_model <- "
model {

# Measurement Model:
  # Mass-spec calibration:
  for (i in 1:Num.cal)
  {
    # Priors:
    log.const.analytic.sd[i] ~ dunif(-10,-0.5)
    log.hetero.analytic.slope[i] ~ dunif(-5,0)
    C.thresh[i] ~ dunif(0,Test.Nominal.Conc[i]/10)
    log.calibration[i] ~ dunif(-2, 2)
    # Scale conversions:
    const.analytic.sd[i] <- 10^log.const.analytic.sd[i]
    hetero.analytic.slope[i] <- 10^log.hetero.analytic.slope[i]
    calibration[i] <- 10^log.calibration[i]
    # Concentrations below this value are not detectable:
    background[i] <- C.thresh[i]/calibration[i]
  }
  
  # Mass-spec observations:  
  for (i in 1:Num.obs) 
  {
# The parameters for calibration curve
    slope[i] <- calibration[obs.cal[i]]
    intercept[i] <- background[obs.cal[i]]
# Mass spec response as a function of diluted concentration:        
    Response.pred[i] <- 
      slope[i] * 
#      (Conc[obs.conc[i]]/Dilution.Factor[i] +
      (Conc[obs.conc[i]]/Dilution.Factor[i] - C.thresh[obs.cal[i]]) *
      step(Conc[obs.conc[i]]/Dilution.Factor[i] - C.thresh[obs.cal[i]]) +
      intercept[i] 
# Heteroskedastic precision:
    Response.prec[i] <- (const.analytic.sd[obs.cal[i]] +
      hetero.analytic.slope[obs.cal[i]] * Response.pred[i])^(-2)
# Model for the observation:
    Response.obs[i] ~ dnorm(Response.pred[i],Response.prec[i])
  }

# Binding Model:
  # Prior on Fup: 
  log.Fup ~ dunif(-15, 0)
  # Scale conversion:
  Fup <- 10^log.Fup


  for (i in (Num.cc.obs +1):(Num.cc.obs + Num.series)) 
  {
  # Priors for whole samples for ultra centrigugation UC):
    Conc[i] ~ dnorm(Test.Nominal.Conc[obs.cal[i]],
      Test.Nominal.Conc[obs.cal[i]]^-2)
  # Aqueous fraction concentrations for UC samples:
    Conc[i+Num.series] <- Fup * Conc[i]
  }   
}
"

#' Calculate fraction unbound in plasma from ultracentrifugation data
#'
#' The data frame of observations should be annotated according to
#' of these types:
#' \tabular{rrrrr}{
#'   Calibration Curve \tab CC\cr
#'   Ultra-centrifugation Aqueous Fraction \tab UC\cr
#'   Whole Plasma T1h Sample  \tab T1\cr
#'   Whole Plasma T5h Sample \tab T5\cr
#' }
#' We don't currently use the T1 data, but you must have CC, AF, and T5 data.
#'
#' @param PPB.data A data frame containing mass-spectrometry peak areas,
#' indication of chemical identiy, and measurment type.
#' @param this.conc The plasma protein concentration relative to physiologic
#' levels (default 100%)
#' @param FILENAME A string used to identify outputs of the function call.
#' (defaults to "BASE_Model_Results")
#' @param TEMP.DIR An optional directory where file writing may be faster.
#' @param JAGS.PATH The file path to JAGS.
#' @param NUM.CHAINS The number of Markov Chains to use. This allows evaluation
#' of convergence according to Gelman and Rubin diagnostic.
#' @param NUM.CORES The number of processors to use (default 2)
#' @param RANDOM.SEED The seed used by the random number generator 
#' (default 1111)
#'
#' @return \item{data.frame}{A data.frame containing quunantiles of the 
#' Bayesian posteriors} 
#'
#' @author John Wambaugh and Chantel Nicolas
#' 
#' @import parallel 
#' @import runjags
#' 
#' @export format_fup_uc
format_fup_uc <- function(PPB.data,
  FILENAME = "UC_Model_Results",
  TEMP.DIR = NULL,
  NUM.CHAINS=5, 
  NUM.CORES=2,
  RANDOM.SEED=1111,
  sample.col="Lab.Sample.Name",
  lab.compound.col="Lab.Compound.Name",
  dtxsid.col="DTXSID",
  date.col="Date",
  compound.col="Compound.Name",
  area.col="Area",
  series.col="Series",
  type.col="Sample.Type",
  compound.conc.col="Nominal.Conc",
  cal.col="Cal",
  dilution.col="Dilution.Factor",
  istd.col="ISTD.Area",
  istd.name.col="ISTD.Name",
  istd.conc.col="ISTD.Conc",
  nominal.test.conc.col="Test.Target.Conc" 
  )
{
  if (!is.null(TEMP.DIR)) 
  {
    current.dir <- getwd()
    setwd(TEMP.DIR)
  }
  
  PPB.data <- as.data.frame(PPB.data)
  all.blanks <- subset(PPB.data,!is.na(eval(area.col)))

  # Only include the data types used:
  PPB.data <- subset(PPB.data,PPB.data[,type.col] %in% c("CC","T1","T5","AF"))
  
  # Organize the columns:
  PPB.data <- PPB.data[,c(
    sample.col,
    date.col,
    compound.col,
    dtxsid.col,
    lab.compound.col,
    type.col,
    dilution.col,
    cal.col,
    compound.conc.col,
    nominal.test.conc.col,
    istd.name.col,
    istd.conc.col,
    istd.col,
    series.col,
    area.col)]
    
  # Rename the columns:
    sample.col <- "Lab.Sample.Name"
    date.col <- "Date"
    compound.col <- "Compound.Name"
    dtxsid.col <- "DTXSID"
    lab.compound.col <- "Lab.Compound.Name"
    type.col <- "Sample.Type"
    dilution.col <- "Dilution.Factor"
    cal.col <- "Calibration"
    compound.conc.col <- "Standard.Conc"
    nominal.test.conc.col <- "Test.Target.Conc"
    istd.name.col <- "ISTD.Name"
    istd.conc.col <- "ISTD.Conc"
    istd.col <- "ISTD.Area"
    series.col <- "Series"
    area.col <- "Area"
  colnames(PPB.data) <- c(
    sample.col,
    date.col,
    compound.col,
    dtxsid.col,
    lab.compound.col,
    type.col,
    dilution.col,
    cal.col,
    compound.conc.col,
    nominal.test.conc.col,
    istd.name.col,
    istd.conc.col,
    istd.col,
    series.col,
    area.col)
  
  # calculate the reponse:
  PPB.data[,"Response"] <- PPB.data[,area.col] /
     PPB.data[,istd.col] *  PPB.data[,istd.conc.col]
  
# Write out a "level 1" file (data organized into a standard format):  
  write.table(PPB.data, 
    file=paste(FILENAME,"-Level1.tsv",sep=""),
    sep="\t",
    row.names=F,
    quote=F)

  return(PPB.data)  
}


