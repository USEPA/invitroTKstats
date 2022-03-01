Clint_model <- "
model {
# Measurement Model:
  log.const.analytic.sd ~ dnorm(-2,0.1)
  const.analytic.sd <- 10^log.const.analytic.sd
  log.hetero.analytic.slope.factor ~ dnorm(1,.1)
  hetero.analytic.slope.factor <- 10^log.hetero.analytic.slope.factor
  hetero.analytic.slope <- hetero.analytic.slope.factor*const.analytic.sd
  for (i in 1:2)
  {
# Number of pixels in the background signal:
    background[i] ~ dunif(0,1000)         
    log.calibration[j] ~ dnorm(0,.1)
# calibration*concentration = number of pixles:
    calibration[i] <- 10^log.calibration[i]
# Concentrations below this value are not detectable:
    log.C.thresh[j] ~ dnorm(log(C.frank/2),0.25)T(log(C.frank/1000),log(10*C.frank))
    C.thresh[j] <- exp(log.C.thresh[j])
    Blank.pred[i] <- background[i]  
    Blank.prec[i] <- (const.analytic.sd+hetero.analytic.slope*(Blank.pred[i]))^(-2)
  }

# Likelihood for the blank observations:
  for (i in 1:Num.blank.obs) {
    Blank.obs[i] ~ dnorm(Blank.pred[Blank.conc[i]],Blank.prec[Blank.conc[i]])
  }

# Clearance model:

# Specify the nominal conc:
  C0[1] <- 1
  C0[2] <- 10

# Decreases indicates whether or not the concentration decreases (1 is yes, 0 is no):
  decreases ~ dbern(0.5) 
# Slope is the clearance rate at the lower concentration (fastest slope we can 
#identify is assumed to be 99.4% gone in the first 15 minutes):
  rate ~ dunif(0,-5/15)
  slope[1] <- decreases*rate
# Saturates is whether or not the clearance rate decreases (1 is yes, 0 is no):
  saturates ~ dbern(0.5) 
# Saturation is how much the clearance rate decreases at the higher conc:
  saturation ~ dunif(0,1)
  slope[2] <- slope[1]*(1 - saturates*saturation)

# The observations are normally distributed (heteroskedastic error):
  for (i in 1:Num.obs)
  {
    C[i] <- C0[obs.conc[i]]*exp(-slope[obs.conc[i]]*obs.time[i])
    obs.pred[i] <- calibration[obs.conc[i]]*(C[i]-C.thresh[obs.conc[i]])*step(C[i]-C.thresh[obs.conc[i]])+ background[obs.conc[i]]    
    obs.prec[i] <- (const.analytic.sd+hetero.analytic.slope*obs.pred[i])^(-2)
    
    obs[i] ~ dnorm(obs.pred[i],obs.prec[i])
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
# Identify the blanks (observation time should be NA):    
    this.blanks <- subset(this.data, Sample.Type=="Blank")
    blank.obs <- this.blanks[,"Response"]
    Num.blanks <- length(blank.obs)
# If there are no blanks create a single blank with the average blank response:
    if (Num.blanks==0 | all(is.na(blank.obs)))
    {
      all.blanks <- subset(clint.data,is.na(this.data) & !is.na(area.col))
      blank.obs <- median(all.blanks$ISTDResponseRatio,na.rm=T)
      Num.blanks <- 1
      warning("Mean blank across data set used because of missing blank data.")
    }
    blank.conc <- rep(2,Num.blanks)
    blank.conc[this.blanks$Conc ==1] <- 1

# Separate the 1 and 10 uM data so we can order obs by concentration:
    this.cvt <- subset(this.data, Sample.Type=="Cvst")
    this.1data <- subset(this.cvt, Conc==1)
    this.10data <- subset(this.cvt, Conc==10)
# Order the data so that all the 1 uM come first:
    obs <- c(this.1data[!is.na(this.1data[,time.col]), "Response"],
      this.1data[!is.na(this.1data[,time.col]), "Response"])
    Num.obs <- length(obs)
    Num.obs1 <- dim(subset(this.1data, !is.na(Time)))[1]
    obs.time <- c(this.1data[!is.na(this.1data[,time.col]), time.col],
      this.1data[!is.na(this.1data[,time.col]), time.col])
    obs.conc <- c(rep(1,Num.obs1),rep(2,Num.obs-Num.obs1))
      
    return(mydata <- list('obs' = obs,
      'obs.conc' = obs.conc,
      'obs.time' = obs.time,
      'Num.obs' = Num.obs,
      'Blank.obs' = blank.obs,
      'Blank.conc' = blank.conc,
      'Num.blank.obs' = Num.blanks
      ))
  }
 
# Function to organize the data to get initial slope estimate: 
  make.fit.data <- function(mydata)
  {
    fit.data <- cbind(mydata$obs,mydata$obs.time,mydata$obs.conc)
    colnames(fit.data) <- c("obs","time","conc")
    fit.data <- as.data.frame(fit.data)
    fit.data1 <- subset(fit.data,conc==1 & obs>0)
  
    if (dim(fit.data1)[1]<2) fit.data1 <- subset(fit.data,conc==2 & obs>0)
    if (dim(fit.data1)[1]<2) fit.data <- NULL
    else fit.data <- fit.data1
   
    return(fit.data)
  }
  
  # function to initialize a Markov chain:
  initfunction <- function(chain)
  {
    background <- rep(0,2)
    calibration <- rep(1,2)
    log.C.thresh <- rep(NA,2)
    for (this.conc in 1:2)
    {
      Blank.data <- mydata[["Blank.obs"]][data[["Blank.conc"]]==this.conc]
      T0.data <- mydata[["obs"]][
        data[["obs.conc"]]==this.conc & mydata[["obs.time"]]==0]
      background[this.conc] <- runif(1,min(Blank.data),max(Blank.data))
      background[is.na(background)] <- 0
      calibration[this.conc] <- runif(1,min(T0.data),
        max(T0.data))/c(1,10)[this.conc]-background[this.conc]        
      calibration[this.conc] <- min(max(calibration[this.conc],10^-3),10^2)
      log.C.thresh[this.conc] <- log(runif(1,0.01,1))
    }
    calibration[is.na(calibration)] <- 1
    if (!is.null(fit.data))
    { 
      if (dim(fit.data)[1]>1 & any(fit.data$time>0))
      {
        fit.data$obs <- fit.data$obs*runif(length(fit.data$obs),0.9,1.1)
        fit1 <- lm(log(obs/calibration[1])~time,fit.data)
        if (-fit1[["coefficients"]][2] > 0) 
        {
          rate <- -fit1[["coefficients"]][2]
          decreases <- 1
        } else {
          rate <- 0
          decreases <- 0
        }
      } else {
        decreases <- 1
        rate <- 1/15
      }
    } else {
      rate <- 0
      decreases <- 0
    }
    return(list(
      .RNG.seed=as.numeric(paste(rep(chain,6),sep="",collapse="")),
      .RNG.name="base::Super-Duper",
      log.const.analytic.sd =runif(1,0,0.1),
      hetero.analytic.slope.factor =runif(1,0, 1),
      background = background,
      log.calibration = log10(calibration),
      decreases = decreases,
      rate = rate,
      saturates = 0,
      saturation = 0.5,
      log.C.thresh = log.C.thresh
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
  conc.col <- "Conc"
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
    conc.col,
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
    "Blank","Cvst"))
  
  # Only used verfied data:
  clint.data <- subset(clint.data, clint.data[,good.col] == "Y")

  # Clean up data:
  clint.data <- subset(clint.data,!is.na(Response))
  clint.data[clint.data$Response<0,"Response"] <- 0
   
  # Because of the possibility of crashes we save the results one chemical at a time:
  OUTPUT.FILE <- paste(FILENAME,"-Level4.tsv",sep="")

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
      this.dtxsid <- this.subset$dtxsid[1]
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
        fit.data <- make.fit.data(mydata)
        # Use random number seed for reproducibility
        set.seed(RANDOM.SEED)
        
        # write out arguments to runjags:
        save(this.compound,mydata,initfunction,
          file=paste(FILENAME,"-PREJAGS.RData",sep=""))
          
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
                             'hetero.analytic.slope.factor',
                             'C.thresh',
                             'calibration',
                             'background'))
        
        coda.out[[this.compound]] <-extend.jags(coda.out[[this.cas]],
                              drop.monitor = c(
                                'log.const.analytic.sd','
                                hetero.analytic.slope.factor'), 
                              add.monitor = c(
                                'slope',
                                'decreases',
                                'saturates'))
 
        sim.mcmc <- coda.out[[this.compound]]$mcmc[[1]]
        for (i in 2:NUM.CHAINS) sim.mcmc <- rbind(sim.mcmc,coda.out[[this.compound]]$mcmc[[i]])
        results <- apply(sim.mcmc,2,function(x) signif(quantile(x,c(0.025,0.5,0.975)),3))
    
        browser()
        new.results <- t(data.frame(c(this.compound,this.dtxsid,this.lab.name),stringsAsFactors=F))
        colnames(new.results) <- c("Compound","DTXSID","Lab.Compound.Name")
        new.results <- cbind.data.frame(new.results,
          t(as.data.frame(as.numeric(results[c(2,1,3),"slope"]))))
        colnames(new.results)[4:6] <- c(
          "Fup.Med",
          "Fup.Low",
          "Fup.High")

        rownames(new.results) <- this.compound
    
        print(mydata$Num.obs)
        print(mydata$Response.obs)
        print(results)
    
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
    file=paste(FILENAME,"-Clint-Level4Analysis-",Sys.Date(),".RData",sep=""))

  return(list(Results=Results,coda=coda.out))  
}


