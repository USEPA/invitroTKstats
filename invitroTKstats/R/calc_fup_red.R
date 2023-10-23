fup_RED_model <- "
model {
# Measurement Model:
  # Mass-spec calibration:
  for (i in 1:Num.cal)
  {
    # Priors:
    # (Note that a uniform prior on the log variable is weighted toward lower
    # values)
    log.const.analytic.sd[i] ~ dunif(-6, 1)
    log.hetero.analytic.slope[i] ~ dunif(-6, 1)
    C.thresh[i] ~ dunif(0,Test.Nominal.Conc[i]/10)
    log.calibration[i] ~ dnorm(0,0.01)
    background[i] ~ dexp(100)
    # Scale conversions:
    const.analytic.sd[i] <- 10^log.const.analytic.sd[i]
    hetero.analytic.slope[i] <- 10^log.hetero.analytic.slope[i]
    calibration[i] <- 10^log.calibration[i]
  }

# Likelihood for blanks without plasma observations:
  for (i in 1:Num.NoPlasma.Blank.obs)
  {
    NoPlasma.Blank.pred[i] <-
      background[NoPlasma.Blank.cal[i]]/NoPlasma.Blank.df
    NoPlasma.Blank.prec[i] <- (const.analytic.sd[NoPlasma.Blank.cal[i]] +
                             hetero.analytic.slope[NoPlasma.Blank.cal[i]]*(NoPlasma.Blank.pred[i]))^(-2)
# Model for the observation:
    NoPlasma.Blank.obs[i] ~ dnorm(NoPlasma.Blank.pred[i], NoPlasma.Blank.prec[i])
  }

# Extent of interference from plasma:
  log.Plasma.Interference ~ dunif(-6, log(10/Test.Nominal.Conc)/log(10))
  Plasma.Interference <- 10^log.Plasma.Interference
# Likelihood for the plasma blank observations:
  for (i in 1:Num.Plasma.Blank.obs)
  {
    Plasma.Blank.pred[i] <-
      (calibration[Plasma.Blank.cal[i]] *
      (Plasma.Interference*Assay.Protein.Percent[Plasma.Blank.rep[i]]/100 -
      C.thresh[Plasma.Blank.cal[i]]) *
      step(Plasma.Interference*Assay.Protein.Percent[Plasma.Blank.rep[i]]/100 -
      C.thresh[Plasma.Blank.cal[i]]) +
      background[Plasma.Blank.cal[i]]) /
      Plasma.Blank.df
    Plasma.Blank.prec[i] <- (const.analytic.sd[Plasma.Blank.cal[i]] +
                             hetero.analytic.slope[Plasma.Blank.cal[i]]*(Plasma.Blank.pred[i]))^(-2)
# Model for the observation:
    Plasma.Blank.obs[i] ~ dnorm(Plasma.Blank.pred[i], Plasma.Blank.prec[i])
  }

# Likelihood for the T0 observations:
  for (i in 1:Num.T0.obs)
  {
# Mass spec response as a function of diluted concentration:
    T0.pred[i] <-
      (calibration[T0.cal[i]] *
      (Test.Nominal.Conc - Plasma.Interference - C.thresh[T0.cal[i]]) *
      step(Test.Nominal.Conc - Plasma.Interference - C.thresh[T0.cal[i]]) +
      calibration[T0.cal[i]] * Plasma.Interference +
      background[T0.cal[i]]) /
      T0.df
# Heteroskedastic precision:
    T0.prec[i] <- (const.analytic.sd[T0.cal[i]] +
      hetero.analytic.slope[T0.cal[i]] * T0.pred[i])^(-2)
# Model for the observation:
    T0.obs[i] ~ dnorm(T0.pred[i], T0.prec[i])
  }

# Likelihood for the calibration curve observations:
  for (i in 1:Num.CC.obs)
  {
# Mass spec response as a function of diluted concentration:
    CC.pred[i] <-
      (calibration[CC.cal[i]] *
      (CC.conc[i] - Plasma.Interference - C.thresh[CC.cal[i]]) *
      step(CC.conc[i] - Plasma.Interference - C.thresh[CC.cal[i]]) +
      calibration[CC.cal[i]] * Plasma.Interference +
      background[CC.cal[i]]) /
      CC.df
# Heteroskedastic precision:
    CC.prec[i] <- (const.analytic.sd[CC.cal[i]] +
      hetero.analytic.slope[CC.cal[i]] * CC.pred[i])^(-2)
# Model for the observation:
    CC.obs[i] ~ dnorm(CC.pred[i], CC.prec[i])
  }

# Likelihood for the RED plasma protein binding assay::
  log.Kd ~ dunif(-10,5)
  Kd <- 10^log.Kd
  Fup <- Kd /
         (Kd +
         Physiological.Protein.Conc)

# Concentrtion in each replicate:
  for (i in 1:Num.rep)
  {
# Calculate protein concentration for observation:
    C.protein[i] <- Physiological.Protein.Conc * Assay.Protein.Percent[i] / 100
# Missing (bound to walls/membrane) chemical:
    C.missing[i] ~ dunif(0, Test.Nominal.Conc)
# Unbound concentration in both wells:
    C.u[i] <- (Test.Nominal.Conc-C.missing[i])*(Kd/(2*Kd+C.protein[i]))
# Bound concentration in plasma well:
    C.b[i] <- (Test.Nominal.Conc-C.missing[i])*C.protein[i]/(2*Kd+C.protein[i])
# Toal concentration in plasma well:
    C.total[i] <- C.b[i] + C.u[i]
  }

# Likelihood for the PBS observations:
  for (i in 1:Num.PBS.obs)
  {
# Mass spec response as a function of diluted concentration:
    PBS.conc[i] <- C.u[PBS.rep[i]]
    PBS.pred[i] <-
      (calibration[PBS.cal[i]] *
      (PBS.conc[i] - C.thresh[PBS.cal[i]]) *
      step(PBS.conc[i]  - C.thresh[PBS.cal[i]]) +
      background[PBS.cal[i]]) /
      PBS.df
# Heteroskedastic precision:
    PBS.prec[i] <- (const.analytic.sd[PBS.cal[i]] +
      hetero.analytic.slope[PBS.cal[i]] * PBS.pred[i])^(-2)
# Model for the observation:
    PBS.obs[i] ~ dnorm(PBS.pred[i], PBS.prec[i])
  }

# Likelihood for the Plasma observations:
  for (i in 1:Num.Plasma.obs)
  {
# Mass spec response as a function of diluted concentration:
    Plasma.conc[i] <- C.total[Plasma.rep[i]]
    Plasma.pred[i] <-
      (calibration[Plasma.cal[i]] *
      (Plasma.conc[i] -
      Plasma.Interference*Assay.Protein.Percent[Plasma.rep[i]]/100 -
      C.thresh[Plasma.cal[i]]) *
      step(Plasma.conc[i] -
      Plasma.Interference*Assay.Protein.Percent[Plasma.rep[i]]/100 -
      C.thresh[Plasma.cal[i]]) +
      calibration[Plasma.cal[i]] *
      Plasma.Interference * Assay.Protein.Percent[Plasma.rep[i]]/100 +
      background[Plasma.cal[i]]) /
      Plasma.df
# Heteroskedastic precision:
    Plasma.prec[i] <- (const.analytic.sd[Plasma.cal[i]] +
      hetero.analytic.slope[Plasma.cal[i]] * Plasma.pred[i])^(-2)
# Model for the observation:
    Plasma.obs[i] ~ dnorm(Plasma.pred[i], Plasma.prec[i])
  }
}
"

#' Calculate Fraction Unbound in Plasma From Rapid Equilibrium Dialysis Data
#'
#' This function uses MCMC simulation to calculate fraction unbound in plasma (Fup) 
#' from rapid equilibrium dialysis (RED) data. This function then calculates quantiles 
#' of posterior samples and returns a summary table along with all MCMC results.
#' 
#' The input data should be a "Level-2" data that have been formatted and created
#' by \code{\link{format_fup_red}}, then have been curated and had a column added
#' with the value "Y" indicating that each row is verified as usable for analysis. 
#' 
#' Be aware that this function will write couple files to user's current working directory
#' unless other path is specified in the argument. Files being saved include the 
#' summary table (.RData), JAG argument (.RData), and any "unverified" data that 
#' were held out from the analysis (.tsv).  
#'
#' The data frame of observations should be annotated according to
#' of these types:
#' \tabular{rrrrr}{
#'   No Plasma Blank (no chemical, no plasma) \tab NoPlasma.Blank\cr
#'   Plasma Blank (no chemical, just plasma) \tab Plasma.Blank\cr
#'   Time zero chemical and plasma \tab T0\cr
#'   Equilibrium chemical in phosphate-buffered well (no plasma) \tab PBS\cr
#'   Equilibrium chemical in plasma well \tab Plasma\cr
#' }
#'
#'
#' @param FILENAME A string used to identify the input Level-2 file, whatever the
#' argument given, "-fup-RED-Level2.tsv" is appended.
#'
#' @param TEMP.DIR Alternative directory to output file
#' (Defaults to NULL, all files will be written to user's current working directory).
#'
#' @param JAGS.PATH Specific path to JAGS on user's computer (Defaults to NA).
#'
#' @param NUM.CHAINS The number of Markov Chains to use (Defaults to 5).
#'
#' @param NUM.CORES The number of processors to use (Defaults to 2).
#'
#' @param RANDOM.SEED The seed used by the random number generator
#' (Defaults to 1111). 
#'
#' @param good.col Name of a column indicating which rows have been verified for
#' analysis, indicated by a "Y" (Defaults to "Verified"). 
#'
#' @param Physiological.Protein.Conc Assumed physiological protein concentration 
#' for plasma protein binding calculations. Defaults to 70/(66.5*1000)*1000000 
#' per Berg and Lane (2011): 60-80 mg/mL, albumin is 66.5 kDa, pretend all protein is albumin to get uM. 
#'
#' @return A list of two objects: 
#' \enumerate{
#'    \item{A data frame containing quantiles of the Bayesian posteriors of 
#'    fraction unbound in plasma (Fup) of all compounds in input file. Column includes:
#'    Compound.Name - compound name, Lab.Compound.Name - compound name used by 
#'    the laboratory, DTXSID - EPA's DSSTox Structure ID, Fup.point - point estimate of Fup,
#'    Fup.Med - median of posteriors, Fup.Low - 2.5th quantile, and Fup.High - 97.5th quantile}
#'    \item{A runjags-class object containing results from JAGS model.}
#' }
#'
#' @references
#' \insertRef{waters2008validation}{invitroTKstats}
#'
#' \insertRef{wambaugh2019assessing}{invitroTKstats}
#'
#' @author John Wambaugh and Chantel Nicolas
#'
#' @examples
#' # Level-2 file
#' write.table(smeltz2023.red,
#'   file="SmeltzPFAS-fup-RED-Level2.tsv",
#'   sep="\t",
#'   row.names=F,
#'   quote=F)
#' 
#' # JAGS.PATH should be changed to the specific path to JAGS on user's computer
#' level4 <- calc_fup_red(FILENAME="SmeltzPFAS",
#'                        NUM.CORES=8,
#'                        JAGS.PATH="C:/Users/jwambaug/AppData/Local/JAGS/JAGS-4.3.0/x64")
#'
#' @export calc_fup_red
calc_fup_red <- function(
  FILENAME,
  TEMP.DIR = NULL,
  NUM.CHAINS=5,
  NUM.CORES=2,
  RANDOM.SEED=1111,
  good.col="Verified",
  JAGS.PATH = NA,
  Physiological.Protein.Conc = 70/(66.5*1000)*1000000 # Berg and Lane (2011) 60-80 mg/mL, albumin is 66.5 kDa, pretend all protein is albumin to get uM
  )
{
# Internal function for constructing data object given to JAGS:
  build_mydata <- function(this.data)
  {
    #mg/mL -> g/L is 1:1
    #kDa -> g/mol is *1000
    #g/mol -> M is g/L/MW
    #M <- uM is /1000000
    Test.Nominal.Conc <- unique(this.data$Test.Nominal.Conc) # uM frank parent concentration
    if (length(Test.Nominal.Conc)>1) stop("Multiple test concentrations.")
# Each calibration could be a unique string (such as a date):
    unique.cal <- sort(unique(this.data[,"Calibration"]))
    Num.cal <- length(unique.cal)
# TIME ZERO
    T0.data <- subset(this.data,Sample.Type=="T0")
    T0.df <- unique(T0.data[,"Dilution.Factor"])
    if (length(T0.df)>1) stop("Multiple T0 dilution factors.")
    T0.obs <- T0.data[,"Response"]
# Convert calibrations to sequential integers:
    T0.cal <- sapply(T0.data[,"Calibration"],
                        function(x) which(unique.cal %in% x))
    Num.T0.obs <- length(T0.obs)
# Calibration Curve
    CC.data <- subset(this.data,Sample.Type=="CC")
    CC.df <- unique(CC.data[,"Dilution.Factor"])
    if (length(CC.df)>1) stop("Multiple CC dilution factors.")
    CC.obs <- CC.data[,"Response"]
# Convert calibrations to sequential integers:
    CC.cal <- sapply(CC.data[,"Calibration"],
                        function(x) which(unique.cal %in% x))
    CC.conc <- CC.data[,"Std.Conc"]
    Num.CC.obs <- length(CC.obs)
# PBS
    PBS.data <- subset(this.data,Sample.Type=="PBS")
    PBS.df <- unique(PBS.data[,"Dilution.Factor"])
    if (length(PBS.df)>1) stop("Multiple PBS dilution factors.")
    PBS.obs <-PBS.data[,"Response"]
# Convert calibrations to sequential integers:
    PBS.cal <- sapply(PBS.data[,"Calibration"],
                        function(x) which(unique.cal %in% x))
    Num.PBS.obs <- length(PBS.obs)
# PLASMA
    Plasma.data <- subset(this.data,Sample.Type=="Plasma")
    Plasma.df <- unique(Plasma.data[,"Dilution.Factor"])
    if (length(Plasma.df)>1) stop("Multiple plasma dilution factors.")
    Plasma.obs <- Plasma.data[,"Response"]
# Convert calibrations to sequential integers:
    Plasma.cal <- sapply(Plasma.data[,"Calibration"],
                        function(x) which(unique.cal %in% x))
    Num.Plasma.obs <- length(Plasma.obs)
# Match the PBS and Plasma replicate measurments:
    PBS.rep <- paste0(PBS.data[,"Calibration"],PBS.data[,"Replicate"])
    Plasma.rep <- paste0(Plasma.data[,"Calibration"],Plasma.data[,"Replicate"])
    unique.rep <- sort(unique(c(PBS.rep,Plasma.rep)))
    Num.rep <- length(unique.rep)
# Convert replicates to sequential integers:
    PBS.rep <- sapply(PBS.rep, function(x) which(unique.rep %in% x))
    Plasma.rep <- sapply(Plasma.rep, function(x) which(unique.rep %in% x))
    Assay.Protein.Percent <- Plasma.data[!duplicated(Plasma.data$Replicate),
                              "Percent.Physiologic.Plasma"]
# NO PLASMA BLANK
    NoPlasma.Blank.data <- subset(this.data, Sample.Type=="NoPlasma.Blank")
    NoPlasma.Blank.df <- unique(NoPlasma.Blank.data[,"Dilution.Factor"])
    if (length(NoPlasma.Blank.df)>1) stop("Multiple blank dilution factors.")
    NoPlasma.Blank.obs <- NoPlasma.Blank.data[,"Response"]
# Convert calibrations to sequential integers:
    NoPlasma.Blank.cal <- sapply(NoPlasma.Blank.data[,"Calibration"],
                        function(x) which(unique.cal %in% x))
    Num.NoPlasma.Blank.obs <- length(NoPlasma.Blank.obs)
    if (Num.NoPlasma.Blank.obs == 0) {
      NoPlasma.Blank.df <- 0
      NoPlasma.Blank.obs <- 0
      NoPlasma.Blank.cal <- 0
    }
# PLASMA BLANK
    Plasma.Blank.data <- subset(this.data, Sample.Type=="Plasma.Blank")
    Plasma.Blank.df <- unique(Plasma.Blank.data[,"Dilution.Factor"])
    if (length(Plasma.Blank.df)>1) stop("Multiple blank dilution factors.")
    Plasma.Blank.obs <- Plasma.Blank.data[,"Response"]
# Convert calibrations to sequential integers:
    Plasma.Blank.cal <- sapply(Plasma.Blank.data[,"Calibration"],
                        function(x) which(unique.cal %in% x))
    Num.Plasma.Blank.obs <- length(Plasma.Blank.obs)
    if (!any(is.na(Plasma.Blank.data[,"Replicate"])))
    {
      Plasma.Blank.rep <- paste0(Plasma.Blank.data[,"Calibration"],
                                 Plasma.Blank.data[,"Replicate"])
# Convert replicates to sequential integers:
      Plasma.Blank.rep <- sapply(Plasma.Blank.rep, function(x)
                               which(unique.rep %in% x))
    } else if(length(unique(Plasma.Blank.data[,"Percent.Physiologic.Plasma"]))==1)
    {
      Plasma.Blank.rep <- rep(1, Num.Plasma.Blank.obs)
    } else browser()

    return(list(
# Describe assay:
      'Test.Nominal.Conc' = Test.Nominal.Conc,
      'Num.cal' = Num.cal,
      'Physiological.Protein.Conc' = Physiological.Protein.Conc,
      'Assay.Protein.Percent' = Assay.Protein.Percent,
# Blank data:
      'Num.Plasma.Blank.obs' = Num.Plasma.Blank.obs,
      'Plasma.Blank.obs' = Plasma.Blank.obs,
      'Plasma.Blank.cal' = Plasma.Blank.cal,
      'Plasma.Blank.df' = Plasma.Blank.df,
      'Plasma.Blank.rep' = Plasma.Blank.rep,
      'Num.NoPlasma.Blank.obs' = Num.NoPlasma.Blank.obs,
      'NoPlasma.Blank.obs' = NoPlasma.Blank.obs,
      'NoPlasma.Blank.cal' = NoPlasma.Blank.cal,
      'NoPlasma.Blank.df' = NoPlasma.Blank.df,
## Callibration.curve.data:
      'Num.CC.obs' = Num.CC.obs,
      'CC.conc' = CC.conc,
      'CC.obs' = CC.obs,
      'CC.cal' = CC.cal,
      'CC.df' = CC.df,
## Stability data:
      'Num.T0.obs' = Num.T0.obs,
      'T0.obs' = T0.obs,
      'T0.cal' = T0.cal,
      'T0.df' = T0.df,
#       'Stability.data' = Stability.data[,"ISTDResponseRatio"],
#      'Num.Stability.obs' = Num.Stability.obs,
## Equilibriation data:
#      'EQ1.data' = EQ1.data[,"ISTDResponseRatio"],
#      'Num.EQ1.obs' = Num.EQ1.obs,
#      'EQ2.data' = EQ2.data[,"ISTDResponseRatio"],
#      'Num.EQ2.obs' = Num.EQ2.obs,
# RED data:
      'Num.rep' = Num.rep,
# PBS data:
      'Num.PBS.obs' = Num.PBS.obs,
      'PBS.obs' = PBS.obs,
      'PBS.cal' = PBS.cal,
      'PBS.df' = PBS.df,
      'PBS.rep' = PBS.rep,
# Plasma data:
      'Num.Plasma.obs' = Num.Plasma.obs,
      'Plasma.obs' = Plasma.obs,
      'Plasma.cal' = Plasma.cal,
      'Plasma.df' = Plasma.df,
      'Plasma.rep' = Plasma.rep
    ))
  }

  initfunction <- function(chain)
  {
    seed <- as.numeric(paste(rep(chain,6),sep="",collapse=""))
    set.seed(seed)

    return(list(
# Random number seed:
      .RNG.seed=seed,
      .RNG.name="base::Super-Duper",
# Parameters that may vary between calibrations:
#      log.const.analytic.sd =runif(mydata$Num.cal,-5,-0.5),
#      log.hetero.analytic.slope = runif(mydata$Num.cal,-5,-0.5),
      log.const.analytic.sd = log10(runif(mydata$Num.cal,0,0.1)),
      log.hetero.analytic.slope = log10(runif(mydata$Num.cal,0,0.1)),
      background = rep(0,mydata$Num.cal),
      C.thresh = runif(mydata$Num.cal, 0, 0.1),
      log.calibration = rep(0,mydata$Num.cal),
      log.Plasma.Interference = log10(runif(mydata$Num.cal,0,0.1)),
# Statistics characterizing the measurement:
      log.Kd= runif(1,-8,4),
      C.missing = runif(mydata$Num.rep,0,mydata[["Test.Nominal.Conc"]])
    ))
  }

  if (!is.null(TEMP.DIR))
  {
    current.dir <- getwd()
    setwd(TEMP.DIR)
  }

  MS.data <- read.csv(file=paste(FILENAME,"-fup-RED-Level2.tsv",sep=""),
    sep="\t",header=T)
  MS.data <- subset(MS.data,!is.na(Compound.Name))
  MS.data <- subset(MS.data,!is.na(Response))

# Standardize the column names:
  sample.col <- "Lab.Sample.Name"
  date.col <- "Date"
  compound.col <- "Compound.Name"
  dtxsid.col <- "DTXSID"
  lab.compound.col <- "Lab.Compound.Name"
  type.col <- "Sample.Type"
  dilution.col <- "Dilution.Factor"
  replicate.col <- "Replicate"
  cal.col <- "Calibration"
  istd.name.col <- "ISTD.Name"
  istd.conc.col <- "ISTD.Conc"
  istd.col <- "ISTD.Area"
  std.conc.col <- "Std.Conc"
  Test.Nominal.Conc.col <- "Test.Nominal.Conc"
  plasma.percent.col <- "Percent.Physiologic.Plasma"
  time.col <- "Time"
  area.col <- "Area"
  analysis.method.col <- "Analysis.Method"
  analysis.instrument.col <- "Analysis.Instrument"
  analysis.parameters.col <- "Analysis.Parameters"
  note.col <- "Note"
  level0.file.col <- "Level0.File"
  level0.sheet.col <- "Level0.Sheet"

# For a properly formatted level 2 file we should have all these columns:
  cols <-c(
    sample.col,
    date.col,
    compound.col,
    dtxsid.col,
    lab.compound.col,
    type.col,
    dilution.col,
    replicate.col,
    cal.col,
    istd.name.col,
    istd.conc.col,
    istd.col,
    std.conc.col,
    Test.Nominal.Conc.col,
    plasma.percent.col,
    time.col,
    area.col,
    analysis.method.col,
    analysis.instrument.col,
    analysis.parameters.col,
    note.col,
    level0.file.col,
    level0.sheet.col,
    "Response",
    good.col)
# Throw error if not all columns present with expected names:
  if (!(all(cols %in% colnames(MS.data))))
  {
    warning("Run format_fup_red first (level 1) then curate to (level 2).")
    stop(paste("Missing columns named:",
      paste(cols[!(cols%in%colnames(MS.data))],collapse=", ")))
  }

  # Only include the data types used:
  MS.data <- subset(MS.data,MS.data[,type.col] %in% c(
    "Plasma.Blank","NoPlasma.Blank","PBS","Plasma","T0","Stability","EQ1","EQ2","CC"))

  # Only used verified data:
  unverified.data <- subset(MS.data, MS.data[,good.col] != "Y")
  write.table(unverified.data, file=paste(
    FILENAME,"-fup-RED-Level2-heldout.tsv",sep=""),
    sep="\t",
    row.names=F,
    quote=F)
  MS.data <- subset(MS.data, MS.data[,good.col] == "Y")

  # Clean up data:
  MS.data <- subset(MS.data,!is.na(Response))
  MS.data[MS.data$Response<0,"Response"] <- 0

  # Because of the possibility of crashes we save the results one chemical at a time:
  OUTPUT.FILE <- paste(FILENAME,"-fup-RED-Level4.tsv",sep="")

  # Check to see if we crashed earlier, if so, don't redo something that already is done
  if (!file.exists(OUTPUT.FILE))
  {
    Results <- NULL
  } else {
    Results <- read.table(OUTPUT.FILE,sep="\t",stringsAsFactors=F,header=T)
  }

  # Make a cluster if using multiple cores:
  if (NUM.CORES>1)
  {
    CPU.cluster <- makeCluster(min(NUM.CORES,NUM.CHAINS))
  } else CPU.cluster <-NA

  coda.out <- list()
  ignored.data <- NULL
  for (this.compound in unique(MS.data[,"Compound.Name"]))
    if (!(this.compound %in% Results[,"Compound.Name"]))
    {
      this.subset <- subset(MS.data,MS.data[,"Compound.Name"]==this.compound)
      this.dtxsid <- this.subset$DTXSID[1]
      this.lab.name <- this.subset[1,lab.compound.col]

      # provide running output of where we are in the list:
      print(paste(
        this.compound,
        " (",
        which(unique(MS.data[,compound.col])==this.compound),
        " of ",
        length(unique(MS.data[,compound.col])),
        ")",
        sep=""))

      REQUIRED.DATA.TYPES <- c("Plasma","PBS","Plasma.Blank","NoPlasma.Blank")
      if (all(REQUIRED.DATA.TYPES %in% this.subset[,type.col]))
      {
        mydata <- build_mydata(this.subset)
        if (!is.null(mydata))
        {
          # Use random number seed for reproducibility
          set.seed(RANDOM.SEED)

          # write out arguments to runjags:
          save(this.compound, mydata ,initfunction,
            file=paste(FILENAME,"-fup-RED-PREJAGS.RData",sep=""))

          # Run JAGS:
          coda.out[[this.compound]] <- autorun.jags(
            fup_RED_model,
            n.chains = NUM.CHAINS,
            method="parallel",
            cl=CPU.cluster,
            summarise=T,
            inits = initfunction,
            startburnin = 25000,
            startsample = 25000,
            max.time="5m",
            crash.retry=2,
            adapt=10000,
            psrf.target = 1.1,
            thin.sample=2000,
            data = mydata,
            jags = JAGS.PATH,
            monitor = c(
  # Chemical analysis parameters:
              'const.analytic.sd',
              'hetero.analytic.slope',
              'C.thresh',
              'log.calibration',
              'background',
              'Plasma.Interference',
  # Measurement parameters:
              'C.missing',
              'Kd',
              'Fup'))

          sim.mcmc <- coda.out[[this.compound]]$mcmc[[1]]
          for (i in 2:NUM.CHAINS) sim.mcmc <- rbind(sim.mcmc,coda.out$mcmc[[i]])
          results <- apply(sim.mcmc,2,function(x) quantile(x,c(0.025,0.5,0.975)))

          Fup.point <- signif(
            (mean(mydata$PBS.obs)*mydata$PBS.df -
             mean(mydata$NoPlasma.Blank.obs)*mydata$NoPlasma.Blank.df) /
            (mean(mydata$Plasma.obs)*mydata$Plasma.df -
             mean(mydata$Plasma.Blank.obs)*mydata$Plasma.Blank.df),
             4)

          new.results <- data.frame(Compound.Name=this.compound,
                                    Lab.Compound.Name=this.lab.name,
                                    DTXSID=this.dtxsid,
                                    Fup.point=Fup.point,
                                    stringsAsFactors=F)
          new.results[,c("Fup.Med","Fup.Low","Fup.High")] <-
            sapply(results[c(2,1,3),"Fup"],
            function(x) signif(x,4))

          print(paste("Final results for ",
            this.compound,
            " (",
            which(unique(MS.data[,compound.col])==this.compound),
            " of ",
            length(unique(MS.data[,compound.col])),
            ")",
            sep=""))
          print(results)
          print(new.results)

          Results <- rbind(Results,new.results)

          write.table(Results,
            file=paste(OUTPUT.FILE,sep=""),
            sep="\t",
            row.names=F,
            quote=F)
        }
      } else {
        ignored.data <- rbind(ignored.data, MSdata)
      }
    }

  if (!is.null(TEMP.DIR))
  {
    setwd(current.dir)
  }
  stopCluster(CPU.cluster)

  write.table(ignored.data,
    file=paste(FILENAME,"-fup-RED-Level2-ignoredbayes.tsv",sep=""),
    sep="\t",
    row.names=F,
    quote=F)

  View(Results)
  save(Results,
    file=paste(FILENAME,"-fup-RED-Level4Analysis-",Sys.Date(),".RData",sep=""))

  return(list(Results=Results,coda=coda.out))
}


