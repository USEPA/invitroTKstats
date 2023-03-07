fup_RED_model <- "
model {

# Measurement Model:
  # Mass-spec calibration:
  for (i in 1:Num.cal)
  {
    # Priors:
    log.const.analytic.sd[i] ~ dnorm(-0.1,0.01)
    log.hetero.analytic.slope[i] ~ dnorm(-2,0.01)
    C.thresh[i] ~ dunif(0,Test.Nominal.Conc[1]/10)
    log.calibration[i] ~ dnorm(0,0.01)
    # Scale conversions:
    const.analytic.sd[i] <- 10^log.const.analytic.sd[i]
    hetero.analytic.slope[i] <- 10^log.hetero.analytic.slope[i]
    calibration[i] <- 10^log.calibration[i]
    # Concentrations below this value are not detectable:
    background[i] <- calibration[i]*C.thresh[i]
  }
  
# Likelihood for the blank observations:
  for (i in 1:Num.cal)
  {
    Blank.pred[i] <- background[i]  
    Blank.prec[i] <- (const.analytic.sd[i]+hetero.analytic.slope[i]*(Blank.pred[i]))^(-2)
  }
  for (i in 1:Num.blank.obs) {
    Blank.obs[i] ~ dnorm(Blank.pred[Blank.cal[i]],Blank.prec[Blank.cal[i]])
  }
  
# Likelihood for the T0 observations:  
  for (i in 1:Num.T0.obs)
  {
# Mass spec response as a function of diluted concentration:  
    T0.pred[i] <- 
      calibration[T0.cal[i]]* 
      (Nominal.Test.Conc/T0.df - C.thresh[T0.cal[i]]) *
      step(Nominal.Test.Conc/T0.df - C.thresh[T0.cal[i]]) +
      background[T0.cal[i]] 
# Heteroskedastic precision:
    T0.prec[i] <- (const.analytic.sd[T0.cal[i]] +
      hetero.analytic.slope[T0.cal[i]] * T0.pred[i])^(-2)
# Model for the observation:
    T0.obs[i] ~ dnorm(T0.pred[i],T0.prec[i]])
  }

# Likelihood for the RED plasma protein binding assay::
  log.Fup ~ dunif(-10,0)
  Fup <- 10^log.Fup
# Concentrtion in each replicate:
  for (i in 1:Num.rep) 
  {
# Missing (bound to walls/membrane) chemical:
    C.missing[i] ~ dunif(0, Nominal.Test.Conc)
# Unbound concentration in both wells:
    C.u[i] <- Fup/(Fup+1)*(Nominal.Test.Conc-C.missing[i])
# Bound concentration in plasma well:
    C.b[i] <- (1-Fup)/Fup*C.u[i]
# Total concentration in plasma well:
    C.total[i] <- C.b[i] + C.u[i] 
  }
# Likelihood for the PBS observations:  
  for (i in 1:Num.PBS.obs)
  {
# Mass spec response as a function of diluted concentration:  
    PBS.conc[i] <- C.u[PBS.rep[i]]
    PBS.pred[i] <- 
      calibration[PBS.cal[i]]* 
      (PBS.conc[i]/PBS.df - C.thresh[PBS.cal[i]]) *
      step(PBS.conc[i]/PBS.df - C.thresh[PBS.cal[i]]) +
      background[PBS.cal[i]] 
# Heteroskedastic precision:
    PBS.prec[i] <- (const.analytic.sd[PBS.cal[i]] +
      hetero.analytic.slope[PBS.cal[i]] * PBS.pred[i])^(-2)
# Model for the observation:
    PBS.obs[i] ~ dnorm(PBS.pred[i], PBS.prec[i]])
  }
# Likelihood for the Plasma observations:  
  for (i in 1:Num.Plasma.obs)
  {
# Mass spec response as a function of diluted concentration:  
    Plasma.conc[i] <- C.total[Plasma.rep[i]]
    Plasma.pred[i] <- 
      calibration[Plasma.cal[i]]* 
      (Plasma.conc[i]/Plamsa.df - C.thresh[Plasma.cal[i]]) *
      step(Plasma.conc[i]/Plasma.df - C.thresh[Plasma.cal[i]]) +
      background[Plasma.cal[i]] 
# Heteroskedastic precision:
    Plasma.prec[i] <- (const.analytic.sd[Plasma.cal[i]] +
      hetero.analytic.slope[Plasma.cal[i]] * Plasma.pred[i])^(-2)
# Model for the observation:
    Plasma.obs[i] ~ dnorm(Plasma.pred[i], Plasma.prec[i]])
  }
}
"

#' Calculate fraction unbound in plasma from rapid equilibruim dialysis data
#' 
#' The data frame of observations should be annotated according to
#' of these types:
#' \tabular{rrrrr}{
#'   Sample Blank (no chemical, just plasma) \tab Blank\cr
#'   Time zero chemical and plasma \tab T0\cr
#'   Equilibrium chemical in phosphate-buffered well (no plasma) \tab PBS\cr
#'   Equilibrium chemical in plasma well \tab Plasma\cr
#' }
#'
#' @param MS.data A data frame containing mass-spectrometry peak areas,
#' indication of chemical identiy, and measurment type.
#'
#' @param this.conc The plasma protein concentration relative to physiologic
#' levels (default 100\%)
#'
#' @param FILENAME A string used to identify outputs of the function call.
#' (defaults to "BASE_Model_Results")
#'
#' @param TEMP.DIR An optional directory where file writing may be faster.
#'
#' @param JAGS.PATH The file path to JAGS.
#'
#' @param NUM.CHAINS The number of Markov Chains to use. This allows evaluation
#' of convergence according to Gelman and Rubin diagnostic.
#'
#' @param NUM.CORES The number of processors to use (default 2)
#'
#' @param RANDOM.SEED The seed used by the random number generator 
#' (default 1111) 
#'
#' @param sample.col Which column of MS.data indicates the unique mass 
#' spectrometry (MS) sample name used by the laboratory. (Defaults to 
#' "Lab.Sample.Name")
#' 
#' @param lab.compound.col Which column of MS.data indicates The test compound 
#' name used by the laboratory (Defaults to "Lab.Compound.Name")
#' 
#' @param dtxsid.col Which column of MS.data indicates EPA's DSSTox Structure 
#' ID (\url{http://comptox.epa.gov/dashboard}) (Defaults to "DTXSID")
#' 
#' @param date.col Which column of MS.data indicates the laboratory measurment
#' date (Defaults to "Date")
#' 
#' @param compound.col Which column of MS.data indicates the test compound
#' (Defaults to "Compound.Name")
#' 
#' @param area.col Which column of MS.data indicates the target analyte (that 
#' is, the test compound) MS peak area (Defaults to "Area")
#' 
#' @param series.col Which column of MS.data indicates the "series", that is
#' a simultaneous replicate (Defaults to "Series")
#' 
#' @param type.col Which column of MS.data indicates the sample type (see table
#' above)(Defaults to "Sample.Type")
#' 
#' @param cal.col Which column of MS.data indicates the MS calibration -- for
#' instance different machines on the same day or different days with the same
#' MS analyzer (Defaults to "Cal")
#' 
#' @param dilution.col Which column of MS.data indicates how many times the
#' sample was diluted before MS analysis (Defaults to "Dilution.Factor")
#' 
#' @param istd.col Which column of MS.data indicates the MS peak area for the
#' internal standard (Defaults to "ISTD.Area")
#' 
#' @param istd.name.col Which column of MS.data indicates identity of the 
#' internal standard (Defaults to "ISTD.Name")
#' 
#' @param istd.conc.col Which column of MS.data indicates the concentration of
#' the internal standard (Defaults to "ISTD.Conc")
#' 
#' @param nominal.test.conc.col Which column of MS.data indicates the intended
#' test chemical concentration at time zero (Defaults to "Test.Target.Conc")       
#'
#' @return A data.frame containing quunantiles of the Bayesian posteriors 
#'
#' @references
#' \insertRef{waters2008validation}{invitroTKstats}
#' \insertRef{wambaugh2019assessing}{invitroTKstats}
#'
#' @author John Wambaugh and Chantel Nicolas
#' 
#' @export calc_fup_red_base
calc_fup_red <- function(
  FILENAME, 
  TEMP.DIR = NULL,
  NUM.CHAINS=5, 
  NUM.CORES=2,
  RANDOM.SEED=1111,
  good.col="Verified",
  JAGS.PATH = NA
  )
{
# Internal function for constructing data object given to JAGS:
  build_mydata <- function(this.data)
  {
# Each calibration could be a unique string (such as a date):
    unique.cal <- sort(unique(this.data[,"Calibration"]))
    Num.cal <- length(unique.cal)
#    
    Blank.data <- subset(this.data,Type=="Blank")
    Blank.df <- unique(Blank.data[,"Dilution.Factor"])
    if (length(Blank.df)>1) stop("Multiple blank dilution factors.") 
    Blank.obs <- Blank.data[,"Response"]
# Convert calibrations to sequential integers:
    Blank.cal <- sapply(Blank.data[,"Calibration"],
                        function(x) which(unique.cal %in% x))
    Num.Blank.obs <- length(Blank.obs)
#
    T0.data <- subset(this.data,Type=="T0")
    T0.df <- T0.data[,"Dilution.Factor"]
    if (length(T0.df)>1) stop("Multiple T0 dilution factors.") 
    T0.obs <- T0.data[,"Response"]
# Convert calibrations to sequential integers:
    T0.cal <- sapply(T0.data[,"Calibration"],
                        function(x) which(unique.cal %in% x))
    Num.T0.obs <- legnth(T0.obs)
#
    PBS.data <- subset(this.data,Type=="PBS")
    PBS.df <- PBS.data[,"Dilution.Factor"]
    if (length(PBS.df)>1) stop("Multiple PBS dilution factors.") 
    PBS.obs <-PBS.data[,"Response"]
# Convert calibrations to sequential integers:
    PBS.cal <- sapply(PBS.data[,"Calibration"],
                        function(x) which(unique.cal %in% x))
    Num.PBS.obs <- length(PBS.obs)
#
    Plasma.data <- subset(this.data,Type=="Plasma")
    Plasma.df <- Plasma.data[,"Dilution.Factor"]
    if (length(Plasma.df)>1) stop("Multiple plasma dilution factors.") 
    Plasma.obs <- Plasma.data[,"Response"]
# Convert calibrations to sequential integers:
    Plasma.cal <- sapply(Plasma.data[,"Calibration"],
                        function(x) which(unique.cal %in% x))
    Num.Plasma.obs <- length(Plasma.obs)
# Match the PBS and Plasma replicate measurments:
    PBS.rep <- paste0(PBS.data[,"Calibration"],PBS.data[,"Replicate"])
    Plasma.rep <- paste0(Plasma.data[,"Calibration"],Plasma.data[,"Replicate"])
    unique.rep <- sort(unique(c(PBS.rep,plasma.rep)))
    Num.rep <- length(unique.rep)
# Convert replicates to sequential integers:
    PBS.rep <- sapply(PBS.rep, function(x) which(unique.rep %in% x))
    Plasma.rep <- sapply(Plasma.rep, function(x) which(unique.rep %in% x))
    #mg/mL -> g/L is 1:1
    #kDa -> g/mol is *1000
    #g/mol -> M is g/L/MW
    #M <- uM is /1000000
    PPB100 <- 70/(66.5*1000)*1000000 # Berg and Lane (2011) 60-80 mg/mL, albumin is 66.5 kDa, pretend all protein is albumin to get uM
    Test.Nominal.Conc <- unique(this.data$Nominal.Test.Con) # uM frank parent concentration
    if (length(Test.Nominal.Conc)>1) stop("Multiple test concentrations.")

    return(list(                
# Describe assay:
      'Test.Nominal.Conc' = Test.Nominal.Conc,
      'Num.cal' = Num.cal,
# Blank data:
      'Num.blank.obs' = Num.blank.obs,
      'Blank.obs' = Blank.obs,
      'Blank.cal' = Blank.cal,
      'Blank.df' = Blank.df,'
## Callibration.curve.data:
#      'Num.cc' = Num.cc.obs,
#      'cc.obs.conc' = cc.obs.conc,
#      'cc.obs' = cc.obs,
#      'cc.obs.cal' = cc.obs.cal,
#      'cc.obs.Dilution.Factor' = cc.obs.df,
## Stability data:
#      'T0.data' = T0.data[,"ISTDResponseRatio"],
#      'Num.T0.obs' = Num.T0.obs,
#      'Stability.data' = Stability.data[,"ISTDResponseRatio"],
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
      'Plasma.rep' = Plasma.rep,
    ))
  }

  initfunction <- function(chain)
  {
    return(list(
# Random number seed:
      .RNG.seed=seed,
      .RNG.name="base::Super-Duper",
# Parameters that may vary between calibrations:
      log.const.analytic.sd =runif(mydata$Num.cal,-5,-0.5),
      log.hetero.analytic.slope = runif(mydata$Num.cal,-5,-0.5),
      C.thresh = runif(mydata$Num.cal, 0, 0.1),
      log.calibration = rep(0,mydata$Num.cal),
# Statisticas characterizing the measurment:
      log.Fup = log10(runif(1,0,1)),
      C.missing = runif(1,0,mydata[["Test.Nominal.Conc"]])
    ))
  }
  
  if (!is.null(TEMP.DIR)) 
  {
    current.dir <- getwd()
    setwd(TEMP.DIR)
  }

  MS.data <- read.csv(file=paste(FILENAME,"-PPB-RED-Level2.tsv",sep=""), 
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
  nominal.test.conc.col <- "Nominal.Test.Conc"
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
    nominal.test.conc.col,
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
    "Blank","PBS","Plasma","T0","Stability","EQ1","EQ2","CC"))
  
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
  for (this.compound in unique(MS.data[,compound.col]))
    if (!(this.compound %in% Results[,compound.col]))
    {
      this.subset <- subset(MS.data,MS.data[,compound.col]==this.compound)
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
          jags = findjags(),
          monitor = c(
            'log.const.analytic.sd',
            'hetero.analytic.slope.factor',
            'Fup',
            'C.thresh',
            'background',
            'calibration'))
        
        sim.mcmc <- coda.out[[MSdata[,"CAS"][1]]]$mcmc[[1]]
        for (i in 2:NUM.CHAINS) sim.mcmc <- rbind(sim.mcmc,coda.out$mcmc[[i]])
        results <- apply(sim.mcmc,2,function(x) quantile(x,c(0.025,0.5,0.975)))
    
        Fup.point <- 2/5 *
          (mean(mydata$PBS.data) - mean(mydata$Blank.data)) /
          (mean(mydata$Plasma.data) - mean(mydata$Blank.data))
        
        new.results <- data.frame(Name=MSdata[,"Name"][1],
                                  DTXSID=MSdata[,"DTXSID"][1],
                                  CAS=MSdata[,"CAS"][1],
                                  CompoundName=this.compound,
                                  Fup.point=Fup.point,
                                  stringsAsFactors=F)
        new.results[,c("Fup.Med","Fup.Low","Fup.High")] <- results[c(2,1,3),"Fup"]
    
        print(new.results)
    
        Results <- rbind(Results,new.results)
    
        write.table(Results, 
          file=paste(OUTPUT.FILE,sep=""),
          sep="\t",
          row.names=F,
          quote=F)
      }    
    }
  
  if (!is.null(TEMP.DIR)) 
  {
    setwd(current.dir)
  }
  stopCluster(CPU.cluster)
  
  View(Results)
  save(Results,
    file=paste(FILENAME,"-fup-RED-Level4Analysis-",Sys.Date(),".RData",sep=""))

  return(list(Results=Results,coda=coda.out))
}


