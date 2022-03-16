Clint_model <- "
model {
# Measurement Model:
  # Mass-spec calibration:
  for (i in 1:Num.cal)
  {
    # Priors:
    log.const.analytic.sd[i] ~ dnorm(-0.1,0.01)
    log.hetero.analytic.slope[i] ~ dnorm(-2,0.01)
    C.thresh[i] ~ dunif(0,Test.Nominal.Conc[i]/10)
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

# Likelihood for the calibration curve observations:
  for (i in 1:Num.cc)
  {
# The parameters for calibration curve
    cc.slope[i] <- calibration[cc.obs.cal[i]]
    cc.intercept[i] <- background[cc.obs.cal[i]]
# Mass spec response as a function of diluted concentration:        
    cc.pred[i] <- 
      cc.slope[i] * 
      (cc.obs.conc[i]/cc.obs.Dilution.Factor[i] - C.thresh[cc.obs.cal[i]]) *
      step(cc.obs.conc[i]/cc.obs.Dilution.Factor[i] - C.thresh[cc.obs.cal[i]]) +
      cc.intercept[i] 
# Heteroskedastic precision:
    cc.prec[i] <- (const.analytic.sd[cc.obs.cal[i]] +
      hetero.analytic.slope[cc.obs.cal[i]] * cc.pred[i])^(-2)
# Model for the observation:
    cc.obs[i] ~ dnorm(cc.pred[i],cc.prec[i])
  }

# Clearance model:

# Decreases indicates whether the concentration decreases (1 is yes, 0 is no):
  decreases ~ dbern(0.5) 
# Degrades indicates whether we think abiotic degradation is a factor
# ()1 is yes, 0 is no):
  degrades ~ dbern(0.5) 
# Slope is the clearance rate at the lower concentration:
  bio.rate ~ dunif(0,10)
# In addition to biological elimination, we also check for abiotic elimination
# (for example, degradation) which we can distinguish if we have data from
# inactivated hepatocutes (Num.abio.obs > 0):
  abio.rate ~ dunif(0,10)
# Total elimination rate is a sum of both:
  slope[1] <- decreases * (bio.rate + degrades * abio.rate)
# Actual biological elimination rate:
  bio.slope[1] <- decreases * bio.rate
# Saturates is whether the bio clearance rate decreases (1 is yes, 0 is no):
  saturates ~ dbern(0.5) 
# Saturation is how much the bio clearance rate decreases at the higher conc:
  saturation ~ dunif(0,1)
# Calculate a slope at the highe concentration:
  slope[2] <- decreases * (degrades * abio.rate +
    bio.rate*(1 - saturates*saturation))
# Actual biological elimination rate:
  bio.slope[2] <- decreases * bio.rate*(1 - saturates*saturation)
# The observations are normally distributed (heteroskedastic error):
  for (i in 1:Num.obs)
  {
  # Exponential decay:
    C[i] <- Test.Nominal.Conc[obs.conc[i]] *
      exp(-slope[obs.conc[i]]*obs.time[i])
  # MS prediction:
    obs.pred[i] <- calibration[obs.cal[i]] *
      (C[i]/obs.Dilution.Factor[i] - C.thresh[obs.cal[i]]) *
      step(C[i]/obs.Dilution.Factor[i] - C.thresh[obs.cal[i]]) + 
      background[obs.cal[i]]    
    obs.prec[i] <- (const.analytic.sd[obs.cal[i]] + 
      hetero.analytic.slope[obs.cal[i]]*obs.pred[i])^(-2)
    
    obs[i] ~ dnorm(obs.pred[i], obs.prec[i])
  }
  
# non-biological clearance model:
  abio.slope <- decreases * abio.rate
# The observations are normally distributed (heteroskedastic error):
  for (i in 1:Num.abio.obs)
  {
  # Exponential decay:
    abio.C[i] <- Test.Nominal.Conc[abio.obs.conc[i]] *
      exp(-abio.slope[abio.obs.conc[i]]*abio.obs.time[i])
  # MS prediction:
    abio.obs.pred[i] <- calibration[abio.obs.cal[i]] *
      (abio.C[i]/abio.obs.Dilution.Factor[i] - C.thresh[abio.obs.cal[i]]) *
      step(abio.C[i]/abio.obs.Dilution.Factor[i] - C.thresh[abio.obs.cal[i]]) + 
      background[abio.obs.cal[i]]    
    abio.obs.prec[i] <- (const.analytic.sd[abio.obs.cal[i]] + 
      hetero.analytic.slope[abio.obs.cal[i]]*abio.obs.pred[i])^(-2)
    
    abio.obs[i] ~ dnorm(abio.obs.pred[i], abio.obs.prec[i])
  }
}
"

#' Calculate intrinsic hepatic clearance
#'
#' This function use describing mass spectrometry (MS) peak areas
#' from samples collected as part of in vitro measurement of chemical clearance
#' as characterized by disappearance of parent compound over time when incubated
#' with primary hepatocytes (Shibata et al., 2000)
#'
#' Data are read from a "Level2" text file that should have been formatted and created 
#' by \code{\link{format_fup_red}} (this is the "Level1" file). The Level1 file
#' should have been curated and had a column added with the value "Y" indicating
#' that each row is verified as usable for analysis (that is, the Level2 file).
#' 
#' The data frame of observations should be annotated according to
#' of these types:
#' \tabular{rrrrr}{
#'   Blank \tab Blank\cr
#'   Hepatocyte inciubation concentration \tab Conc\cr
#' }
#'
#' Clint is calculated using \code{\link{lm}} to perform a linear regression of
#' MS response as a function of time.
#'
#' @param FILENAME A string used to identify the input file, whatever the
#' argument given, "-PPB-RED-Level2.tsv" is appended (defaults to "MYDATA")
#' 
#' @param good.col Name of a column indicating which rows have been verified for 
#' analysis, indicated by a "Y" (Defaults to "Verified")
#'
#' @return \item{data.frame}{A data.frame in standardized format} 
#'
#' @author John Wambaugh
#' 
#' @examples
#' 
#' library(invitroTKstats)
#' 
#' clint <- wambaugh2019.clint
#' clint$Date <- "2019"
#' clint$Sample.Type <- "Blank"
#' clint$Time..mins. <- as.numeric(clint$Time..mins.)
#' clint[!is.na(clint$Time..mins.),"Sample.Type"] <- "Cvst"
#' clint$ISTD.Name <- "Bucetin, Propranolol, and Diclofenac"
#' clint$ISTD.Conc <- 1
#' clint$Dilution.Factor <- 1
#' clint[is.na(clint$FileName),"FileName"]<-"Wambaugh2019"
#' clint$Hep.Density <- 0.5
#' clint$Analysis.Method <- "LC or GC" 
#' clint$Analysis.Instrument <- "No Idea"
#' clint$Analysis.Parameters <- "None"
#' 
#' level1 <- format_clint(clint,
#'   FILENAME="Wambaugh2019",
#'   sample.col="Sample.Name",
#'   compound.col="Preferred.Name",
#'   lab.compound.col="Name",
#'   time.col="Time..mins.",
#'   cal.col="FileName")
#'
#' level2 <- level1
#' level2$Verified <- "Y"
#' 
#' # All data (allows test for saturation):
#' write.table(level2,
#'   file="Wambaugh2019-Clint-Level2.tsv",
#'   sep="\t",
#'   row.names=F,
#'   quote=F)
#' 
#' level4 <- calc_clint_point(FILENAME="Wambaugh2019")
#'
#' @export calc_clint
calc_clint <- function(FILENAME, good.col="Verified")
{
# Internal function for constructing data object given to JAGS:
  build_mydata <- function(this.data)
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
    Num.blanks <- length(blank.obs)
    # Create a dummy vector to keep JAGS happy:
    if (Num.blanks == 0) blank.obs <- c(NA, NA)
    if (Num.blanks == 1) blank.obs <- c(blank.obs, NA)
# Match the blanks to correct calibration curve:
    if (Num.blanks > 0) blank.cal <- rep(NA, Num.blanks)
    else blank.cal <- c(NA, NA)
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
      abio.obs <- c(NA,NA)
      abio.obs.conc <- c(NA,NA)
      abio.obs.time <- c(NA,NA)
      abio.obs.cal <- c(NA,NA)
      abio.obs.df<- c(NA,NA)
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
      cc.obs <- c(NA,NA)
      cc.obs.conc <- c(NA,NA)
      cc.obs.cal <- c(NA,NA)
      cc.obs.df <- c(NA,NA)
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
      'abio.obs.Dilution.Factor' = abio.obs.df
      ))
  }
 
  # function to initialize a Markov chain:
  initfunction <- function(chain)
  {
    seed <- as.numeric(paste(rep(chain,6),sep="",collapse=""))
    set.seed(seed)

    return(list(
      .RNG.seed=seed,
      .RNG.name="base::Super-Duper",
# Parameters that may vary between calibrations:
      log.const.analytic.sd =runif(mydata$Num.cal,-5,-0.5),
      log.hetero.analytic.slope = runif(mydata$Num.cal,-5,-0.5),
      C.thresh = runif(mydata$Num.cal, 0, 0.1),
      log.calibration = rep(0,mydata$Num.cal),
      decreases = rbinom(1,1,0.5),
      degrades = rbinom(1,1,0.5),
      bio.rate = runif(1,0,1/15),
      abio.rate = runif(1,0,1/15),
      saturates = rbinom(1,1,0.5),
      saturation = runif(1,0,1)
    ))
  }

  clint.data <- read.csv(file=paste(FILENAME,"-Clint-Level2.tsv",sep=""), 
    sep="\t",header=T)  
  clint.data <- subset(clint.data,!is.na(Compound.Name))
  clint.data <- subset(clint.data,!is.na(Response))

# Standardize the column names:
  sample.col <- "Lab.Sample.Name"
  date.col <- "Date"
  compound.col <- "Compound.Name"
  dtxsid.col <- "DTXSID"
  lab.compound.col <- "Lab.Compound.Name"
  type.col <- "Sample.Type"
  dilution.col <- "Dilution.Factor"
  cal.col <- "Calibration"
  istd.name.col <- "ISTD.Name"
  istd.conc.col <- "ISTD.Conc"
  istd.col <- "ISTD.Area"
  density.col <- "Hep.Density"
  std.conc.col <- "Std.Conc"
  clint.assay.conc.col <- "Clint.Assay.Conc"
  time.col <- "Time"
  area.col <- "Area"

# We need all these columns in clint.data
  cols <-c(
    sample.col,
    date.col,
    compound.col,
    dtxsid.col,
    lab.compound.col,
    type.col,
    dilution.col,
    cal.col,
    istd.name.col,
    istd.conc.col,
    istd.col,
    density.col,
    std.conc.col,
    clint.assay.conc.col,
    time.col,
    area.col)
      
  # Check for missing columns
  if (!(all(cols %in% colnames(clint.data))))
  {
    warning("Run format_clint first (level 1) then curate to (level 2).")
    stop(paste("Missing columns named:",
      paste(cols[!(cols%in%colnames(clint.data))],collapse=", ")))
  }

  # Only include the data types used:
  clint.data <- subset(clint.data,clint.data[,type.col] %in% c(
    "Blank","Cvst","CC","Inactive"))
  
  # Only used verfied data:
  clint.data <- subset(clint.data, clint.data[,good.col] == "Y")

  # Clean up data:
  clint.data <- subset(clint.data,!is.na(Response))
  clint.data[clint.data$Response<0,"Response"] <- 0
   
  # Because of the possibility of crashes we save the results one chemical at a time:
  OUTPUT.FILE <- paste(FILENAME,"-Clint-Level4.tsv",sep="")

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
  for (this.compound in unique(clint.data[,compound.col]))
    if (!(this.compound %in% Results[,compound.col]))
    {
      this.subset <- subset(clint.data,clint.data[,compound.col]==this.compound)
      this.dtxsid <- this.subset$DTXSID[1]
      this.lab.name <- this.subset[1,lab.compound.col]

      this.row <- c(this.subset[1,c(compound.col,dtxsid.col)],
        data.frame(Calibration="All Data",
          Clint=NaN,
          Clint.pValue=NaN))
      this.cvt <- subset(this.subset,Sample.Type=="Cvst")
      this.blank <- subset(this.subset,Sample.Type=="Blank")
      if (length(unique(this.cvt$Dilution.Factor))>1) browser()
      df.cvt <- this.cvt$Dilution.Factor[1]
      if (length(unique(this.cvt$Hep.Density))>1) browser()
      hep.density <- this.cvt$Hep.Density[1]
    
      # provide running output of where we are in the list:
      print(paste(
        this.compound,
        " (",
        which(unique(clint.data[,compound.col])==this.compound),
        " of ",
        length(unique(clint.data[,compound.col])),
        ")",
        sep=""))
  
  
      mydata <- build_mydata(this.subset)
      if (!is.null(mydata))
      {
        # Use random number seed for reproducibility
        set.seed(RANDOM.SEED)
        
        # write out arguments to runjags:
        save(this.compound,mydata,initfunction,
          file=paste(FILENAME,"-Clint-PREJAGS.RData",sep=""))
          
        # Run JAGS:
        coda.out[[this.compound]] <-  autorun.jags(Clint_model, 
                           n.chains = NUM.CHAINS,
                           method="parallel", 
                           cl=CPU.cluster,
                           summarise=T,
                           inits = initfunction,
                           max.time="300s",
                           startsample=4000,
                           adapt=5000,
                           startburnin=20000,
                           psrf.target = 1.1,
                           thin=5,
                           thin.sample=2000,
                           data = mydata,
                           jags = JAGS.PATH,
                           monitor = c(
                             'log.const.analytic.sd',
                             'hetero.analytic.slope',
                             'C.thresh',
                             'calibration',
                             'background'))
        
        coda.out[[this.compound]] <-extend.jags(coda.out[[this.compound]],
                              drop.monitor = c(
                             'log.const.analytic.sd',
                             'hetero.analytic.slope',
                             'C.thresh',
                             'calibration',
                             'background'), 
                              add.monitor = c(
                                'slope',
                                'decreases',
                                'saturates',
                                "degrades"))
 
        sim.mcmc <- coda.out[[this.compound]]$mcmc[[1]]
        for (i in 2:NUM.CHAINS) sim.mcmc <- rbind(sim.mcmc,coda.out[[this.compound]]$mcmc[[i]])
        results <- apply(sim.mcmc,2,function(x) signif(quantile(x,c(0.025,0.5,0.975)),3))
        results <- as.data.frame(results)
        
        # Convert disappareance rates (1/h)to 
        # heaptic clearance (uL/min/10^6 hepatocytes)):
        hep.density <- this.subset[1,density.col]
        
        # Calculate a Clint only for 1 and 10 uM (if tested)
        if (1 %in% mydata$Test.Nominal.Conc)
        {
          index <- which(mydata$Test.Nominal.Conc == 1)
          results[,"Clint.1"] <- 1000 *
            results[,paste("slope[",index,"]",sep="")] / hep.density
        } else results[,"Clint.1"] <- NA
        if (10 %in% mydata$Test.Nominal.Conc)
        {
          index <- which(mydata$Test.Nominal.Conc == 10)
          results[,"Clint.10"] <- 1000 *
            results[,paste("slope[",index,"]",sep="")] / hep.density
        } else results[,"Clint.10"] <- NA
    
        # Round to 3 sig figs:
        for (i in dim(results)[2])
          results[,i] <- signif(results[,i],3)

        # Calculate a Clint "pvalue" from probability that we observed a decrease:
        results[,"Clint.pValue"] <- sum(sim.mcmc[,"decreases"]==0)/dim(sim.mcmc)[1]
         
        # Calculate a "pvalue" for saturation probability that we observed
        # a lower Clint at higher conc:
        results[,"Sat.pValue"] <- sum(sim.mcmc[,"saturates"]==0)/dim(sim.mcmc)[1]

        # Calculate a "pvalue" for abiotic degradation:
        results[,"degrades.pValue"] <- sum(sim.mcmc[,"degrades"]==0)/dim(sim.mcmc)[1]    

        # Create a row of formatted results:
        new.results <- t(data.frame(c(this.compound,this.dtxsid,this.lab.name),stringsAsFactors=F))
        colnames(new.results) <- c(compound.col, dtxsid.col, lab.compound.col)
        for (this.param in c("Clint.1","Clint.10","Clint.pValue","Sat.pValue"))
        {
          new.results <- cbind.data.frame(new.results,
            t(as.data.frame(as.numeric(results[c(2,1,3), this.param]))))
          last.col <- length(new.results)
          colnames(new.results)[(last.col-2):last.col] <- c(
            paste(this.param,".Med",sep=""),
            paste(this.param,".Low",sep=""),
            paste(this.param,".High",sep=""))
        }
        rownames(new.results) <- this.compound
         
        print(results)
    
        Results <- rbind(Results,new.results)
    
        write.table(Results, 
          file=paste(OUTPUT.FILE,sep=""),
          sep="\t",
          row.names=F,
          quote=F)
      }    
    }  
  
  stopCluster(CPU.cluster)
  
  View(Results)
  save(Results,
    file=paste(FILENAME,"-Clint-Level4Analysis-",Sys.Date(),".RData",sep=""))

  return(list(Results=Results,coda=coda.out))  
}


