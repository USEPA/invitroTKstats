#' Formatting function for X-axis in log10-scale
#'
#' @param x (Character) String to be formatted. 
#'
#' @return Text with desired expression. Replace any scientific e notation to ten notation, 
#' simplify 10^01 to 10 and 10^0 to 1.
#'
#' @importFrom scales scientific_format
#'
scientific_10 <- function(x) {
  out <- gsub("1e", "10^", scientific_format()(x))
  out <- gsub("\\+","",out)
  out <- gsub("10\\^01","10",out)
  out <- parse(text=gsub("10\\^00","1",out))
}

#' Heaviside
#'
#' @param x (Numeric) A numeric vector.
#'
#' @param threshold (Numeric) A threshold value used to compare to elements in \code{x}. (Defaults to 0.)
#'
#' @return A vector of 1 and 0 in which 1 indicates the element in \code{x} is larger or equal to the threshold.
#'
#'
Heaviside <- function(x, threshold=0)
{
  out <- rep(0,length(x))
  out[x >= threshold] <- 1
  return(out)
}

#' Convert a runjags-class object to a list
#'
#' @param runjagsdata.in (\code{runjags} Object) MCMC results from \code{autorun.jags}.
#'
#' @return A list object containing MCMC results from the provided runjags object.
#'
#'
runjagsdata.to.list <- function(runjagsdata.in)
{
  temp <- strsplit(runjagsdata.in,"\n")[[1]]
  list.out <- list()
  for (i in 1:length(temp))
  {
    temp2 <- strsplit(temp[i]," <- ")[[1]]
    list.out[[gsub("\"","",temp2[1])]] <- eval(parse(text=temp2[2]))
  }
  return(list.out)
}


#' Build Data Object for Intrinsic Hepatic Clearance (Clint) Bayesian Model
#'
#' @param this.cvt (Data Frame) Subset of data containing all "Cvst" sample observations of one test compound.
#' @param this.data (Data Frame) Subset of data containing all observations of one test compound.
#' @param decrease.prob (Numeric) Prior probability that a chemical will decrease in
#' the assay.
#' @param saturate.prob (Numeric) Prior probability that a chemicals rate of metabolism
#' will decrease between 1 and 10 uM.
#' @param degrade.prob (Numeric) Prior probability that a chemical will be unstable
#' (that is, degrade abiotically) in the assay.
#'
#' @return A named list to be passed into the Bayesian model. 
#'
build_mydata_clint <- function(this.cvt, this.data, decrease.prob, saturate.prob, degrade.prob)
{
  #
  # What concentrations were tested (1 and 10 uM typical):
  #
  # Establish a vector of unique nominal test concentrations:
  Test.conc <- sort(unique(unique(this.cvt[,"Clint.Assay.Conc"])))
  Num.conc <- length(Test.conc)
  #
  # How many separate mass-spec calibrations were made:
  #
  # Cal.name is used for matching the observations to the calibrations, but is not
  # passed on to JAGS
  Cal.name <- unique(this.data$Calibration)
  # Identify the number of calibrations:
  Num.cal <- length(Cal.name)
  #
  # CvT Obs:
  #
  # Extract the observations
  this.cvt <- subset(this.data, Sample.Type=="Cvst")
  obs <-  this.cvt[!is.na(this.cvt[,"Time"]), "Response"]
  Num.obs <- length(obs)
  obs.time <- this.cvt[!is.na(this.cvt[,"Time"]), "Time"]
  obs.df <- this.cvt[!is.na(this.cvt[,"Time"]), "Dilution.Factor"]
  obs.conc <- rep(NA, Num.obs)
  for (this.conc in Test.conc)
  {
    obs.conc[this.cvt[
      !is.na(this.cvt[,"Time"]), "Clint.Assay.Conc"] == this.conc] <-
      which(Test.conc == this.conc)
  }
  # Match observations to correct calibration curve:
  obs.cal <- rep(NA, Num.obs)
  for (this.cal in unique(this.cvt$Calibration))
  {
    obs.cal[this.cvt$Calibration == this.cal] <-
      which(Cal.name == this.cal)
  }
  #
  # Blanks (hepatocytes, no chemical):
  #
  # Identify the blanks (observation time should be NA):
  this.blanks <- subset(this.data, Sample.Type=="Blank")
  blank.obs <- this.blanks[,"Response"]
  blank.df <- this.blanks[,"Dilution.Factor"]
  Num.blanks <- length(blank.obs)
  # Create a dummy vector to keep JAGS happy:
  if (Num.blanks == 0) {
    blank.obs <- c(-99,-99)
    blank.df <- c(-99,-99)
  } else if (Num.blanks == 1) {
    blank.obs <- c(blank.obs, -99)
    blank.df <- c(blank.df, -99)
  }
  # Match the blanks to correct calibration curve:
  if (Num.blanks > 0) blank.cal <- rep(NA, Num.blanks)
  else blank.cal <- c(-99,-99)
  for (this.cal in unique(this.blanks$Calibration))
  {
    blank.cal[this.blanks$Calibration == this.cal] <-
      which(Cal.name == this.cal)
  }
  #
  # Inactivated hepatocytes
  #
  # Get the inactive hepatocyte data (if any):
  this.abio <- subset(this.data, Sample.Type=="Inactive" &
                        !is.na(Time))
  Num.abio.obs <- dim(this.abio)[1]
  if (Num.abio.obs > 0)
  {
    abio.obs <-  this.abio[!is.na(this.abio[,"Time"]), "Response"]
    abio.obs.time <- this.abio[!is.na(this.abio[,"Time"]), "Time"]
    abio.obs.df <- this.abio[!is.na(this.abio[,"Time"]), "Dilution.Factor"]
    abio.obs.conc <- rep(NA, Num.abio.obs)
    for (this.conc in Test.conc)
    {
      abio.obs.conc[this.abio[
        !is.na(this.abio[,"Time"]), "Clint.Assay.Conc"] == this.conc] <-
        which(Test.conc == this.conc)
    }
    abio.obs.cal <- rep(NA, Num.abio.obs)
    for (this.cal in unique(this.abio$Calibration))
    {
      abio.obs.cal[this.abio$Calibration == this.cal] <-
        which(Cal.name == this.cal)
    }
  } else {
    abio.obs <- c(-99,-99)
    abio.obs.conc <- c(-99,-99)
    abio.obs.time <- c(-99,-99)
    abio.obs.cal <- c(-99,-99)
    abio.obs.df<- c(-99,-99)
  }
  #
  # Calibration curve measurements
  #
  # Get the calibration curves (if any):
  this.cc <- subset(this.data, Sample.Type=="CC" &
                      !is.na(Test.Compound.Conc))
  Num.cc.obs <- dim(this.cc)[1]
  if (Num.cc.obs > 0)
  {
    cc.obs <- this.cc[, "Response"]
    cc.obs.conc <- this.cc[, "Test.Compound.Conc"]
    cc.obs.df <- this.cc[, "Dilution.Factor"]
    cc.obs.cal <- rep(NA, Num.cc.obs)
    for (this.cal in unique(this.cc[,"Calibration"]))
    {
      cc.obs.cal[this.cc[,"Calibration"] == this.cal] <-
        which(Cal.name == this.cal)
    }
  } else {
    cc.obs <- c(-99,-99)
    cc.obs.conc <- c(-99,-99)
    cc.obs.cal <- c(-99,-99)
    cc.obs.df <- c(-99,-99)
  }
  
  return(mydata <- list('obs' = obs,
                        # Describe assay:
                        'Test.Nominal.Conc' = Test.conc,
                        'Num.cal' = Num.cal,
                        # Cvt data:
                        'Num.obs' = Num.obs,
                        'obs.conc' = obs.conc,
                        'obs.time' = obs.time,
                        'obs.cal' = obs.cal,
                        'obs.Dilution.Factor' = obs.df,
                        # Blank data:
                        'Num.blank.obs' = Num.blanks,
                        'Blank.obs' = blank.obs,
                        'Blank.cal' = blank.cal,
                        'Blank.Dilution.Factor' = blank.df,
                        # Callibration.curve.data:
                        'Num.cc' = Num.cc.obs,
                        'cc.obs.conc' = cc.obs.conc,
                        'cc.obs' = cc.obs,
                        'cc.obs.cal' = cc.obs.cal,
                        'cc.obs.Dilution.Factor' = cc.obs.df,
                        # Abiotic degradation data:
                        'Num.abio.obs' = Num.abio.obs,
                        'abio.obs' = abio.obs,
                        'abio.obs.conc' = abio.obs.conc,
                        'abio.obs.time' = abio.obs.time,
                        'abio.obs.cal' = abio.obs.cal,
                        'abio.obs.Dilution.Factor' = abio.obs.df,
                        # Priors for decrease/saturation/degradation:
                        'DECREASE.PROB' = decrease.prob,
                        'SATURATE.PROB' = saturate.prob,
                        'DEGRADE.PROB' = degrade.prob
  ))
}

#' Set Initial Values for Intrinsic Hepatic Clearance (Clint) Bayesian Model
#' 
#' @param mydata (List) Output of \code{build_mydata_clint}.
#' @param chain (Numeric) The number of Markov Chains to use.
#' 
#' @return A list of initial values.
#' 
initfunction_clint <- function(mydata, chain)
{
  seed <- as.numeric(paste(rep(chain,6),sep="",collapse=""))
  set.seed(seed)
  
  return(list(
    # Random number seed:
    .RNG.seed=seed,
    .RNG.name="base::Super-Duper",
    # Parameters that may vary between calibrations:
    log.const.analytic.sd = log10(runif(mydata$Num.cal,0,0.1)),
    log.hetero.analytic.slope = log10(runif(mydata$Num.cal,0,0.1)),
    C.thresh = runif(mydata$Num.cal, 0, 0.1),
    log.calibration = rep(0,mydata$Num.cal),
    background = rep(0,mydata$Num.cal),
    # Statistics characterizing the measurement:
    decreases = rbinom(1,1,0.5),
    degrades = rbinom(1,1,0.5),
    bio.rate = runif(1,0.05,0.25),
    abio.rate = runif(1,0.05,0.25),
    saturates = rbinom(1,1,0.5),
    saturation = runif(1,0,1)
  ))
}

#' Build Data Object for Fup RED Bayesian Model
#' 
#' @param this.data (Data Frame) Subset of data containing all observations of one test compound.
#' @param Physiological.Protein.Conc (Numeric) The assumed physiological protein concentration 
#' for plasma protein binding calculations. 
#' 
#' @return A named list to be passed into the Bayesian model.
#' 
build_mydata_fup_red <- function(this.data, Physiological.Protein.Conc)
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
  Num.T0.obs <- nrow(T0.data)
  if(Num.T0.obs > 0){
    T0.df <- unique(T0.data[,"Dilution.Factor"])
    if (length(T0.df)>1) stop("Multiple T0 dilution factors.")
    T0.obs <- T0.data[,"Response"]
    # Convert calibrations to sequential integers:
    T0.cal <- sapply(T0.data[,"Calibration"],
                     function(x) which(unique.cal %in% x))
  }else{
    T0.df <- c(-99,-99)
    T0.obs <- c(-99,-99)
    T0.cal <- c(-99,-99)
  }
  
  
  # Calibration Curve
  #
  # Get the calibration curves (if any):
  CC.data <- subset(this.data,Sample.Type=="CC" & !is.na(Test.Compound.Conc)) # may need to also add & !is.na("Test.Compound.Conc")
  Num.CC.obs <- nrow(CC.data)
  if(Num.CC.obs > 0){
    CC.df <- unique(CC.data[,"Dilution.Factor"])
    if (length(CC.df)>1) stop("Multiple CC dilution factors.")
    CC.obs <- CC.data[,"Response"]
    # Convert calibrations to sequential integers:
    CC.cal <- sapply(CC.data[,"Calibration"],
                     function(x) which(unique.cal %in% x))
    CC.conc <- CC.data[,"Test.Compound.Conc"]
  }else{
    CC.obs <- c(-99,-99)
    CC.cal <- c(-99,-99)
    CC.conc <- c(-99,-99)
    CC.df <- c(-99,-99)
  }
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
  PBS.rep <- paste0(PBS.data[,"Calibration"],PBS.data[,"Technical.Replicates"])
  Plasma.rep <- paste0(Plasma.data[,"Calibration"],Plasma.data[,"Technical.Replicates"])
  unique.rep <- sort(unique(c(PBS.rep,Plasma.rep)))
  Num.rep <- length(unique.rep)
  # Convert replicates to sequential integers:
  PBS.rep <- sapply(PBS.rep, function(x) which(unique.rep %in% x))
  Plasma.rep <- sapply(Plasma.rep, function(x) which(unique.rep %in% x))
  Assay.Protein.Percent <- Plasma.data[!duplicated(Plasma.data$Technical.Replicates),
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
    NoPlasma.Blank.obs <- c(-99,-99)
    NoPlasma.Blank.cal <- c(-99,-99)
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
  if (!any(is.na(Plasma.Blank.data[,"Technical.Replicates"])))
  {
    Plasma.Blank.rep <- paste0(Plasma.Blank.data[,"Calibration"],
                               Plasma.Blank.data[,"Technical.Replicates"])
    # Convert replicates to sequential integers:
    Plasma.Blank.rep <- sapply(Plasma.Blank.rep, function(x)
      which(unique.rep %in% x))
  } else if(length(unique(Plasma.Blank.data[,"Percent.Physiologic.Plasma"]))==1)
  {
    Plasma.Blank.rep <- rep(1, Num.Plasma.Blank.obs)
  } else{
    stop("build_mydata_fup_red - `Plasma.Blank.rep` cannot be allocated due to lack of `Technical.Replicates` information and multiple `Percent.Physiologic.Plasma` values for `Plasma.Blank` samples.")
    # browser()
  } 
  
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

#' Set Initial Values for Fup RED Bayesian Model
#' 
#' @param mydata (List) Output of \code{build_mydata_fup_red}.
#' @param chain (Numeric) The number of Markov Chains to use.
#' 
#' @return A list of initial values.
#' 
initfunction_fup_red <- function(mydata, chain)
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

#' Build Data Object for Fup UC Bayesian Model
#'
#' @param MS.data (Data Frame) Subset of data containing all observations of one test compound.
#' @param CC.data (Data Frame) Subset of data containing observations of calibration curves samples.
#' @param T1.data (Data Frame) Subset of data containing observations of Whole Plasma T1h Samples.
#' @param T5.data (Data Frame) Subset of data containing observations of Whole Plasma T5h Samples.
#' @param AF.data (Data Frame) Subset of data containing observations of Aqueous Fraction samples.
#' 
#' @return A named list to be passed into the Bayesian model. 
build_mydata_fup_uc <- function(MS.data, CC.data, T1.data, T5.data, AF.data){
  
    all.cal <- unique(MS.data[,"Calibration"])
    Num.cal <- length(all.cal)        
    #
    #
    #
    #CC.data <- MS.data[MS.data[,type.col]=="CC",]
    Num.cc.obs <- dim(CC.data)[1]
    CC.data$Obs.Conc <- seq(1,Num.cc.obs)
    Conc <- CC.data[,"Test.Compound.Conc"]
    Dilution.Factor <- CC.data[,"Dilution.Factor"]
    #
    #
    #  Each series contains T1, T5, and AF data
    #T1.data <- MS.data[MS.data[,type.col]=="T1",]
    Num.T1.obs <- dim(T1.data)[1]
    #T5.data <- MS.data[MS.data[,type.col]=="T5",]
    Num.T5.obs <- dim(T5.data)[1]
    #AF.data <- MS.data[MS.data[,type.col]=="AF",]
    Num.AF.obs <- dim(AF.data)[1]
    Num.series <- 0
    all.series <- NULL
    Test.Nominal.Conc <- NULL
    for (i in 1:Num.cal)
    {
      these.series <- unique(T5.data[
        T5.data[,"Calibration"]==all.cal[i],
        "Biological.Replicates"])
      Num.series <- Num.series + length(these.series) 
      T1.data[
        T1.data[,"Calibration"]==all.cal[i],
        "Biological.Replicates"] <- paste(all.cal[i],
                             T1.data[                          
                               T1.data[,"Calibration"]==all.cal[i],
                               "Biological.Replicates"],
                             sep="-")
      T5.data[
        T5.data[,"Calibration"]==all.cal[i],
        "Biological.Replicates"] <- paste(all.cal[i],
                             T5.data[                          
                               T5.data[,"Calibration"]==all.cal[i],
                               "Biological.Replicates"],
                             sep="-")
      AF.data[
        AF.data[,"Calibration"]==all.cal[i],
        "Biological.Replicates"] <- paste(all.cal[i],
                             AF.data[
                               AF.data[,"Calibration"]==all.cal[i],
                               "Biological.Replicates"],
                             sep="-")
      all.series <- c(all.series,paste(all.cal[i],these.series,sep="-"))
      Test.Nominal.Conc[i] <- mean(T1.data[
        T1.data[,"Calibration"]==all.cal[i],
        "UC.Assay.T1.Conc"],na.rm=T)
    }
    # There is one initial concentration per series, even if there are
    # multiple observations of that series:
    for (i in 1:Num.series)
    {
      T1.data[T1.data$Biological.Replicates==all.series[i],"Obs.Conc"] <- 
        Num.cc.obs + i
      T5.data[T5.data$Biological.Replicates==all.series[i],"Obs.Conc"] <- 
        Num.cc.obs + 1*Num.series + i
      AF.data[AF.data$Biological.Replicates==all.series[i],"Obs.Conc"] <-   
        Num.cc.obs + 2*Num.series + i
    }
    # There are three total concentrations per series (T1, T5, and AF):
    Conc <- c(Conc,rep(NA,3*Num.series))
    #
    #
    #
    UC.obs <- rbind(CC.data,T1.data,T5.data,AF.data)
    Num.obs <- dim(UC.obs)[1]
    for (i in 1:Num.cal)
    {
      UC.obs[UC.obs[,"Calibration"]==all.cal[i],"Obs.Cal"] <- i
    }
    #
    #
    #
    return(list(   
      'Num.cal' = Num.cal,            
      'Num.obs' = Num.obs,
      "Response.obs" = UC.obs[,"Response"],
      "obs.conc" = UC.obs[,"Obs.Conc"],
      "obs.cal" = UC.obs[,"Obs.Cal"],
      "Conc" = Conc,
      "Num.cc.obs" = Num.cc.obs,
      "Num.series" = Num.series,
      "Dilution.Factor" = UC.obs[,"Dilution.Factor"],
      "Test.Nominal.Conc" = Test.Nominal.Conc
    ))
}


#' Set Initial Values for Fup UC Bayesian Model
#' 
#' @param mydata (List) Output of \code{build_mydata_clint}.
#' @param chain (Numeric) The number of Markov Chains to use.
#' 
#' @return A list of initial values.
#' 
initfunction_fup_uc <- function(mydata, chain)
{
  seed <- as.numeric(paste(rep(chain,6),sep="",collapse=""))
  set.seed(seed)
  cal.coeff <- lm(
    mydata$Response.obs[1:mydata$Num.cc.obs]~
      mydata$Conc[1:mydata$Num.cc.obs])[["coefficients"]]
  slope <- as.numeric(cal.coeff[2])
  intercept <- as.numeric(cal.coeff[1])
  
  # We need a vector with NA's for all the values that are not sampled, but 
  # initial values for the concentrations that are inferred (the T1's):
  init.Conc <- rep(NA,mydata$Num.cc.obs+mydata$Num.series*3)
  # Set initial values for the T1's:
  init.Conc[(mydata$Num.cc.obs+1):
              (mydata$Num.cc.obs+mydata$Num.series)] <- 
    mydata$Test.Nominal.Conc
  
  return(list(
    .RNG.seed=seed,
    .RNG.name="base::Super-Duper",
    # Parameters that may vary between calibrations:
    #      log.const.analytic.sd =runif(mydata$Num.cal,0.5,1),
    #      log.hetero.analytic.slope = runif(mydata$Num.cal,-5,-3),
    log.const.analytic.sd = log10(runif(mydata$Num.cal,0,0.1)),
    log.hetero.analytic.slope = log10(runif(mydata$Num.cal,0,0.1)),
    # Average across all the calibrations (the sampler will vary these):
    C.thresh = rep(
      min(
        max(10^-8,abs(intercept)/slope),
        mydata$Test.Nominal.Conc/10,na.rm=TRUE),
      mydata$Num.cal),
    background = rep(0,mydata$Num.cal),
    log.calibration = rep(max(
      min(-2.95,
          log10(max(0,
                    slope))),
      1.95),
      mydata$Num.cal),
    # There is only one Fup per chemical:
    log.Fup = log10(runif(1,0,1)),
    # There is only one Fstable per chemical:
    log.Floss = runif(1,-4,-2),
    # Set the initial concentrations:
    Conc = init.Conc
  ))
}
