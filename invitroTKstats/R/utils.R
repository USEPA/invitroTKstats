#' Formatting function for X-axis in log10-scale
#'
#' @param x
#'
#' @return text with desired expression
#'
#' @importFrom scales scientific_format
#'
#' @export calc_fup_uc
scientific_10 <- function(x) {
  out <- gsub("1e", "10^", scientific_format()(x))
  out <- gsub("\\+","",out)
  out <- gsub("10\\^01","10",out)
  out <- parse(text=gsub("10\\^00","1",out))
}

#' Heaviside
#'
#' @param x a vector
#'
#' @param threshold defaults to 0
#'
#' @return  a vector of 1 and 0, 1 indicates input element is above threshold
#'
#'
#' @export Heaviside
Heaviside <- function(x, threshold=0)
{
  out <- rep(0,length(x))
  out[x >= threshold] <- 1
  return(out)
}

#' Convert a runjags-class object to a list
#'
#' @param runjagsdata.in MCMC results from autorun.jags(), object of class runjags
#'
#' @return List object containing MCMC results from the provided runjags object.
#'
#'
#' @export runjagsdata.to.list
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
  obs <-  this.cvt[!is.na(this.cvt[,time.col]), "Response"]
  Num.obs <- length(obs)
  obs.time <- this.cvt[!is.na(this.cvt[,time.col]), "Time"]
  obs.df <- this.cvt[!is.na(this.cvt[,time.col]), "Dilution.Factor"]
  obs.conc <- rep(NA, Num.obs)
  for (this.conc in Test.conc)
  {
    obs.conc[this.cvt[
      !is.na(this.cvt[,time.col]), "Clint.Assay.Conc"] == this.conc] <-
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
    abio.obs <-  this.abio[!is.na(this.abio[,time.col]), "Response"]
    abio.obs.time <- this.abio[!is.na(this.abio[,time.col]), "Time"]
    abio.obs.df <- this.abio[!is.na(this.abio[,time.col]), "Dilution.Factor"]
    abio.obs.conc <- rep(NA, Num.abio.obs)
    for (this.conc in Test.conc)
    {
      abio.obs.conc[this.abio[
        !is.na(this.abio[,time.col]), "Clint.Assay.Conc"] == this.conc] <-
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
                      !is.na(Std.Conc))
  Num.cc.obs <- dim(this.cc)[1]
  if (Num.cc.obs > 0)
  {
    cc.obs <- this.cc[, "Response"]
    cc.obs.conc <- this.cc[, "Std.Conc"]
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


initfunction_clint <- function(chain)
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
    # Statistics characterizing the measurment:
    decreases = rbinom(1,1,0.5),
    degrades = rbinom(1,1,0.5),
    bio.rate = runif(1,0.05,0.25),
    abio.rate = runif(1,0.05,0.25),
    saturates = rbinom(1,1,0.5),
    saturation = runif(1,0,1)
  ))
}
