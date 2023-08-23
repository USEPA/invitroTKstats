UC_PPB_model <- "
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
  
  # Mass-spec observations:  
  for (i in 1:Num.obs) 
  {
# The parameters for calibration curve
    slope[i] <- calibration[obs.cal[i]]
    intercept[i] <- background[obs.cal[i]]
# Mass spec response as a function of diluted concentration:        
    Response.pred[i] <- 
      slope[i] * 
      ((Conc[obs.conc[i]] - C.thresh[obs.cal[i]]) *
      step(Conc[obs.conc[i]] - C.thresh[obs.cal[i]]) +
      intercept[i])/Dilution.Factor[i] 
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
  # Prior on Fstable:
  log.Floss ~ dunif(-6, 0)
  Fstable <- 1 - 10^log.Floss
    
# The conc's we don't know are for the T1, T5, and AF
  for (i in (Num.cc.obs +1):(Num.cc.obs + Num.series)) 
  {
  # Priors for T1 samples for ultra centrigugation UC):
    Conc[i] ~ dnorm(Test.Nominal.Conc[obs.cal[i]],
      100)
  # The T5 samples after potential breakdown:
    Conc[i + Num.series] <- Fstable * Conc[i]
  # Aqueous fraction concentrations for stable chemical at T5:
    Conc[i + 2*Num.series] <- Fup * Conc[i + Num.series]
  }   
}
"

#' Calculate fraction unbound in plasma from ultracentrifugation data
#'
#' The data frame of observations should be annotated according to
#' of these types:
#' \tabular{rrrrr}{
#'   Calibration Curve \tab CC\cr
#'   Ultracentrifugation Aqueous Fraction \tab UC\cr
#'   Whole Plasma T1h Sample  \tab T1\cr
#'   Whole Plasma T5h Sample \tab T5\cr
#' }
#' We don't currently use the T1 data, but you must have CC, AF, and T5 data.
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
#' @param PPB.data A data frame containing mass-spectrometry peak areas,
#' indication of chemical identity, and measurement type. The data frame should
#' contain columns with names specified by the following arguments:
#' 
#' @param sample.col Which column of PPB.data indicates the unique mass 
#' spectrometry (MS) sample name used by the laboratory. (Defaults to 
#' "Lab.Sample.Name")
#' 
#' @param lab.compound.col Which column of PPB.data indicates The test compound 
#' name used by the laboratory (Defaults to "Lab.Compound.Name")
#' 
#' @param dtxsid.col Which column of PPB.data indicates EPA's DSSTox Structure 
#' ID (\url{http://comptox.epa.gov/dashboard}) (Defaults to "DTXSID")
#' 
#' @param date.col Which column of PPB.data indicates the laboratory measurement
#' date (Defaults to "Date")
#' 
#' @param compound.col Which column of PPB.data indicates the test compound
#' (Defaults to "Compound.Name")
#' 
#' @param area.col Which column of PPB.data indicates the target analyte (that 
#' is, the test compound) MS peak area (Defaults to "Area")
#' 
#' @param series.col Which column of PPB.data indicates the "series", that is
#' a simultaneous replicate (Defaults to "Series")
#' 
#' @param type.col Which column of PPB.data indicates the sample type (see table
#' above)(Defaults to "Sample.Type")
#' 
#' @param cal.col Which column of PPB.data indicates the MS calibration -- for
#' instance different machines on the same day or different days with the same
#' MS analyzer (Defaults to "Cal")
#'
#' #param std.conc.col Which column indicates the intended concentration 
#' of the test chemical for calibration curves (Defaults to "Standard.Conc")
#' 
#' @param dilution.col Which column of PPB.data indicates how many times the
#' sample was diluted before MS analysis (Defaults to "Dilution.Factor")
#' 
#' @param istd.col Which column of PPB.data indicates the MS peak area for the
#' internal standard (Defaults to "ISTD.Area")
#' 
#' @param istd.name.col Which column of PPB.data indicates identity of the 
#' internal standard (Defaults to "ISTD.Name")
#' 
#' @param istd.conc.col Which column of PPB.data indicates the concentration of
#' the internal standard (Defaults to "ISTD.Conc")
#' 
#' @param uc.assay.conc.col Which column of PPB.data indicates the intended
#' test chemical concentration at time zero (Defaults to "UC.Assay.Conc") 
#'
#' @return A data.frame containing quunantiles of the Bayesian posteriors 
#'
#' @author John Wambaugh and Chantel Nicolas
#' 
#' @import parallel 
#' @import runjags
#' 
#' @export calc_fup_uc
calc_fup_uc <- function(PPB.data,
  FILENAME = "UC_Model_Results",
  TEMP.DIR = NULL,
  NUM.CHAINS=5, 
  NUM.CORES=2,
  RANDOM.SEED=1111,
  good.col="Verified",
  JAGS.PATH = NA
  )
{

# local function to give each chain it's own starting values:
  initfunction <- function(chain)
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
        
  if (!is.null(TEMP.DIR)) 
  {
    current.dir <- getwd()
    setwd(TEMP.DIR)
  }
  
  PPB.data <- read.csv(file=paste(FILENAME,"-fup-UC-Level2.tsv",sep=""), 
    sep="\t",header=T)  
  PPB.data <- subset(PPB.data,!is.na(Compound.Name))
  PPB.data <- subset(PPB.data,!is.na(Response))

  # Standardize the column names:
    sample.col <- "Lab.Sample.Name"
    date.col <- "Date"
    compound.col <- "Compound.Name"
    dtxsid.col <- "DTXSID"
    lab.compound.col <- "Lab.Compound.Name"
    type.col <- "Sample.Type"
    dilution.col <- "Dilution.Factor"
    cal.col <- "Calibration"
    std.conc.col <- "Standard.Conc"
    uc.assay.conc.col <- "UC.Assay.T1.Conc"
    istd.name.col <- "ISTD.Name"
    istd.conc.col <- "ISTD.Conc"
    istd.col <- "ISTD.Area"
    series.col <- "Series"
    area.col <- "Area"
    analysis.method.col <- "Analysis.Method"
    analysis.instrument.col <- "Analysis.Instrument"
    analysis.parameters.col <- "Analysis.Parameters" 
    note.col <- "Note"

# For a properly formatted level 2 file we should have all these columns:
  cols <-c(
    sample.col,
    date.col,
    compound.col,
    dtxsid.col,
    lab.compound.col,
    type.col,
    dilution.col,
    cal.col,
    std.conc.col,
    uc.assay.conc.col,
    istd.name.col,
    istd.conc.col,
    istd.col,
    series.col,
    area.col,
    analysis.method.col,
    analysis.instrument.col,
    analysis.parameters.col,
    note.col,
    "Response",
    good.col)
  if (!(all(cols %in% colnames(PPB.data))))
  {
    warning("Run format_fup_uc first (level 1) then curate to level 2.")
    stop(paste("Missing columns named:",
      paste(cols[!(cols%in%colnames(PPB.data))],collapse=", ")))
  }

  # Only include the data types used:
  PPB.data <- subset(PPB.data,PPB.data[,type.col] %in% c(
    "CC",
    "T1",
    "T5",
    "AF"))   
  
  # Only used verified data:
  unverified.data <- subset(PPB.data, PPB.data[,good.col] != "Y")
  write.table(unverified.data, file=paste(
    FILENAME,"-fup-UC-Level2-heldout.tsv",sep=""),
    sep="\t",
    row.names=F,
    quote=F)
  PPB.data <- subset(PPB.data, PPB.data[,good.col] == "Y")
  
  PPB.data <- as.data.frame(PPB.data)
  all.blanks <- subset(PPB.data,!is.na(eval(area.col)))
  
  OUTPUT.FILE <- paste(FILENAME,"-fup-UC-Level4.tsv",sep="")

  set.seed(RANDOM.SEED)
  if (!file.exists(OUTPUT.FILE))
  {
    Results <- NULL
  } else {
    Results <- read.table(OUTPUT.FILE,sep="\t",stringsAsFactors=F,header=T)
  }

  if (NUM.CORES>1)
  {
    CPU.cluster <- makeCluster(min(NUM.CORES,NUM.CHAINS))
  } else CPU.cluster <-NA
  
  coda.out <- list()
  ignored.data <- NULL
  for (this.compound in  unique(PPB.data[,compound.col]))
    if (!(this.compound %in% Results[,"Compound"]))
    {
      this.name <- PPB.data[PPB.data[,compound.col]==this.compound,compound.col][1]
      this.dtxsid <- PPB.data[PPB.data[,compound.col]==this.compound,dtxsid.col][1]
      this.lab.name <- PPB.data[PPB.data[,compound.col]==this.compound,lab.compound.col][1]
      print(paste(
        this.name,
        " (",
        which(unique(PPB.data[,compound.col])==this.compound),
        " of ",
        length(unique(PPB.data[,compound.col])),
        ")",
        sep=""))
      MS.data <- PPB.data[PPB.data[,compound.col]==this.compound,]
    
      for (this.series in unique(MS.data[,series.col]))
        if (!is.na(this.series))
        {
          this.series.subset <- subset(MS.data,MS.data[,series.col]==this.series)
          for (this.cal in unique(this.series.subset[,cal.col]))
            if (!is.na(this.cal))
            {
              this.cal.subset <- subset(this.series.subset,          
                                      this.series.subset[,cal.col]==this.cal)            
              if (!all(c("T1","T5","AF") %in% this.cal.subset[,type.col]))
              {
                # Have to handle the NA series values for CC data:
                series.values <- MS.data[,series.col]
                # Assign a dummy value to the NA's
                series.values[is.na(series.values)]<-"Cat"
                # Identify the bad series from the cal and add to ignored.data:
                ignored.data <- rbind(ignored.data, subset(MS.data,
                                      series.values == this.series & 
                                      MS.data[,cal.col]==this.cal))
                # Remove the bad series:
                MS.data <- subset(MS.data,
                                 series.values != this.series |
                                 MS.data[,cal.col]!=this.cal)
                print(paste("Dropped series",this.series,"from cal",
                           this.cal,"for incomplete data."))
              }
            } 
        }
    
      if (any(MS.data[,type.col]=="CC") &
          any(MS.data[,type.col]=="T1") &
          any(MS.data[,type.col]=="T5") &
          any(MS.data[,type.col]=="AF"))
      {
        all.cal <- unique(MS.data[,cal.col])
        Num.cal <- length(all.cal)        
  #
  #
  #
        CC.data <- MS.data[MS.data[,type.col]=="CC",]
        Num.cc.obs <- dim(CC.data)[1]
        CC.data$Obs.Conc <- seq(1,Num.cc.obs)
        Conc <- CC.data[,std.conc.col]
        Dilution.Factor <- CC.data[,dilution.col]
  #
  #
  #  Each series contains T1, T5, and AF data
        T1.data <- MS.data[MS.data[,type.col]=="T1",]
        Num.T1.obs <- dim(T1.data)[1]
        T5.data <- MS.data[MS.data[,type.col]=="T5",]
        Num.T5.obs <- dim(T5.data)[1]
        AF.data <- MS.data[MS.data[,type.col]=="AF",]
        Num.AF.obs <- dim(AF.data)[1]
        Num.series <- 0
        all.series <- NULL
        Test.Nominal.Conc <- NULL
        for (i in 1:Num.cal)
        {
          these.series <- unique(T5.data[
            T5.data[,cal.col]==all.cal[i],
            series.col])
          Num.series <- Num.series + length(these.series) 
          T1.data[
            T1.data[,cal.col]==all.cal[i],
            series.col] <- paste(all.cal[i],
             T1.data[                          
               T1.data[,cal.col]==all.cal[i],
               series.col],
             sep="-")
          T5.data[
            T5.data[,cal.col]==all.cal[i],
            series.col] <- paste(all.cal[i],
             T5.data[                          
               T5.data[,cal.col]==all.cal[i],
               series.col],
             sep="-")
          AF.data[
            AF.data[,cal.col]==all.cal[i],
            series.col] <- paste(all.cal[i],
             AF.data[
               AF.data[,cal.col]==all.cal[i],
               series.col],
             sep="-")
          all.series <- c(all.series,paste(all.cal[i],these.series,sep="-"))
          Test.Nominal.Conc[i] <- mean(T1.data[
            T1.data[,cal.col]==all.cal[i],
            uc.assay.conc.col],na.rm=T)
        }
        # There is one initial concentration per series, even if there are
        # multiple observations of that series:
        for (i in 1:Num.series)
        {
          T1.data[T1.data$Series==all.series[i],"Obs.Conc"] <- 
            Num.cc.obs + i
          T5.data[T5.data$Series==all.series[i],"Obs.Conc"] <- 
            Num.cc.obs + 1*Num.series + i
          AF.data[AF.data$Series==all.series[i],"Obs.Conc"] <-   
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
          UC.obs[UC.obs[,cal.col]==all.cal[i],"Obs.Cal"] <- i
        }
  #
  #
  #
        mydata <- list(   
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
        )
      
        save(this.compound,mydata,UC_PPB_model,initfunction,
          file=paste(FILENAME,"-Fup-UC-PREJAGS.RData",sep=""))  
        coda.out[[this.compound]] <- autorun.jags(
          UC_PPB_model, 
          n.chains = NUM.CHAINS,
          method="parallel", 
          cl=CPU.cluster,
          summarise=T,
          inits = initfunction,
          startburnin = 50000, 
          startsample = 50000, 
          max.time="1h",
          crash.retry=2,
          adapt=15000,
          psrf.target = 1.05,
          thin.sample=2000,
          data = mydata,
          jags = findjags(look_in=JAGS.PATH),
          monitor = c(
          # Chemical analysis parameters:
            'const.analytic.sd',
            'hetero.analytic.slope',
            'C.thresh',
            'log.calibration',
            'background',
          # Measurement parameters:
            'Fup',
            "Fstable"
            ))

        sim.mcmc <- coda.out[[this.compound]]$mcmc[[1]]
        for (i in 2:NUM.CHAINS) sim.mcmc <- rbind(sim.mcmc,coda.out[[this.compound]]$mcmc[[i]])
        results <- apply(sim.mcmc,2,function(x) signif(quantile(x,c(0.025,0.5,0.975)),3))
    
        new.results <- t(data.frame(c(this.compound,this.dtxsid,this.lab.name),stringsAsFactors=F))
        colnames(new.results) <- c("Compound","DTXSID","Lab.Compound.Name")
         new.results <- cbind.data.frame(new.results,
        t(as.data.frame(as.numeric(results[c(2,1,3),"Fstable"]))))
        colnames(new.results)[4:6] <- c(
          "Fstable.Med",
          "Fstable.Low",
          "Fstable.High")
        new.results <- cbind.data.frame(new.results,
          t(as.data.frame(as.numeric(results[c(2,1,3),"Fup"]))))
        colnames(new.results)[7:9] <- c(
          "Fup.Med",
          "Fup.Low",
          "Fup.High")
        new.results[,"Fup.point"] <- signif(mean(AF.data[,"Response"] *
          AF.data[,"Dilution.Factor"]) / mean(T5.data[,"Response"] *
          T5.data[,"Dilution.Factor"]),3)
        rownames(new.results) <- this.compound
    
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
      } else {
        ignored.data <- rbind(ignored.data, MS.data)
      }   
    }
  
  if (!is.null(TEMP.DIR)) 
  {
    setwd(current.dir)
  }
  stopCluster(CPU.cluster)

  write.table(ignored.data, 
    file=paste(FILENAME,"-fup-UC-Level2-ignoredbayes.tsv",sep=""),
    sep="\t",
    row.names=F,
    quote=F)
  
  View(Results)
  save(Results,
    file=paste(FILENAME,"-fup-UC-Level4Analysis-",Sys.Date(),".RData",sep=""))

  return(list(Results=Results,coda=coda.out))  
}


