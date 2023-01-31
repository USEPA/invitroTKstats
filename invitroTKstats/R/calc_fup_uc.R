UC_PPB_model <- "
model {

# Measurement Model:
  # Mass-spec calibration:
  for (i in 1:Num.cal)
  {
    # Priors:
    log.const.analytic.sd[i] ~ dnorm(-0.1,0.01)
    log.hetero.analytic.slope[i] ~ dnorm(-2,0.01)
    log.C.thresh[i] ~ dnorm(log(Test.Nominal.Conc[i]/10), 0.01)
    log.calibration[i] ~ dnorm(0,0.01)
    # Scale conversions:
    const.analytic.sd[i] <- 10^log.const.analytic.sd[i]
    hetero.analytic.slope[i] <- 10^log.hetero.analytic.slope[i]
    C.thresh[i] <- 10^log.C.thresh[i]
    calibration[i] <- 10^log.calibration[i]
    # Concentrations below this value are not detectable:
    background[i] <- calibration[i]*C.thresh[i]
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

# The conc's we don't know are for the AF and T5 samples
  for (i in (Num.cc.obs +1):(Num.cc.obs + Num.series)) 
  {
  # Priors for whole samples for ultra centrigugation UC):
    Conc[i] ~ dnorm(Test.Nominal.Conc[obs.cal[i]],
      100)
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
#' indication of chemical identiy, and measurment type. The data frame should
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
#' @param date.col Which column of PPB.data indicates the laboratory measurment
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
#' #param std.conc.col Which column indictes the intended concentration 
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
# initial values for the concentrations that are inferred (the T5's):
    init.Conc <- rep(NA,mydata$Num.cc.obs+mydata$Num.series*2)
    init.Conc[(mydata$Num.cc.obs+1):(mydata$Num.cc.obs+mydata$Num.series)] <- 
      mydata$Test.Nominal.Conc[
        mydata$obs.cal[(mydata$Num.cc.obs+1):(mydata$Num.cc.obs+mydata$Num.series)]]
      
    return(list(
      .RNG.seed=seed,
      .RNG.name="base::Super-Duper",
# Parameters that may vary between calibrations:
      log.const.analytic.sd =runif(mydata$Num.cal,-1.5,0.5),
      log.hetero.analytic.slope = runif(mydata$Num.cal,-5,-0.5),
# Average across all the calibrations (the sampler will vary these):
      log.C.thresh = log10(rep(
                     min(
                         max(10^-8,intercept/slope),
                         mydata$Test.Nominal.Conc/10,na.rm=TRUE),
                     mydata$Num.cal)),
      log.calibration = rep(max(min(-2.95,log10(max(0,slope))),1.95),mydata$Num.cal),
# There is only one Fup per chemical:
      log.Fup = log10(runif(1,0,1)),
# Set the initial concentrations:
      Conc = init.Conc
    ))
  }
        
  if (!is.null(TEMP.DIR)) 
  {
    current.dir <- getwd()
    setwd(TEMP.DIR)
  }
  
  PPB.data <- read.csv(file=paste(FILENAME,"-PPB-UC-Level2.tsv",sep=""), 
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
    FILENAME,"-PPB-UC-Level2-heldout.tsv",sep=""),
    sep="\t",
    row.names=F,
    quote=F)
  PPB.data <- subset(PPB.data, PPB.data[,good.col] == "Y")
  
  PPB.data <- as.data.frame(PPB.data)
  all.blanks <- subset(PPB.data,!is.na(eval(area.col)))
  
  OUTPUT.FILE <- paste(FILENAME,"-PPB-UC-Level4.tsv",sep="")

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
      MSdata <- PPB.data[PPB.data[,compound.col]==this.compound,]
    
      for (this.series in unique(MSdata[,series.col]))
        if (!is.na(this.series))
        {
          this.series.subset <- subset(MSdata,MSdata[,series.col]==this.series)
          for (this.cal in unique(this.series.subset[,cal.col]))
            if (!is.na(this.cal))
            {
              this.cal.subset <- subset(this.series.subset,          
                                      this.series.subset[,cal.col]==this.cal)            
              if (!all(c("T1","T5","AF") %in% this.cal.subset[,type.col]))
              {
                # Have to handle the NA series values for CC data:
                series.values <- MSdata[,series.col]
                # Assign a dummy value to the NA's
                series.values[is.na(series.values)]<-"Cat"
                # Identify the bad series from the cal and add to ignored.data:
                ignored.data <- rbind(ignored.data, subset(MSdata,
                                      series.values == this.series & 
                                      MSdata[,cal.col]==this.cal))
                # Remove the bad series:
                MSdata <- subset(MSdata,
                                 series.values != this.series |
                                 MSdata[,cal.col]!=this.cal)
                print(paste("Dropped series",this.series,"from cal",
                           this.cal,"for incomplete data."))
              }
            } 
        }
    
      if (any(MSdata[,type.col]=="CC") &
          any(MSdata[,type.col]=="T1") &
          any(MSdata[,type.col]=="T5") &
          any(MSdata[,type.col]=="AF"))
      {
        all.cal <- unique(MSdata[,cal.col])
        Num.cal <- length(all.cal)        
  #
  #
  #
        CC.data <- MSdata[MSdata[,type.col]=="CC",]
        Num.cc.obs <- dim(CC.data)[1]
        CC.data$Obs.Conc <- seq(1,Num.cc.obs)
        Conc <- CC.data[,std.conc.col]
        Dilution.Factor <- CC.data[,dilution.col]
  #
  #
  #  Each series (currently) contains T5 and AF data
        T1.data <- MSdata[MSdata[,type.col]=="T1",]
        T5.data <- MSdata[MSdata[,type.col]=="T5",]
        AF.data <- MSdata[MSdata[,type.col]=="AF",]
        Num.series <- 0
        all.series <- NULL
        Test.Nominal.Conc <- NULL
        for (i in 1:Num.cal)
        {
          these.series <- unique(T5.data[
            T5.data[,cal.col]==all.cal[i],
            series.col])
          Num.series <- Num.series + length(these.series) 
          T5.data[
            T5.data[,cal.col]==all.cal[i],
            "Series"] <- paste(all.cal[i],
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
        for (i in 1:Num.series)
        {
          T5.data[T5.data$Series==all.series[i],"Obs.Conc"] <- Num.cc.obs + i
          AF.data[AF.data$Series==all.series[i],"Obs.Conc"] <- Num.cc.obs + Num.series + i
        }
        Conc <- c(Conc,rep(NA,2*Num.series))
        Dilution.Factor <- c(Dilution.Factor,
          T5.data[,dilution.col],
          AF.data[,dilution.col])
  #
  #
  #
        UC.obs <- rbind(CC.data,T5.data,AF.data)
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
          "Dilution.Factor" = Dilution.Factor,
          "Test.Nominal.Conc" = Test.Nominal.Conc
        )
      
        save(this.compound,mydata,UC_PPB_model,initfunction,
          file=paste(FILENAME,"-FupUC-PREJAGS.RData",sep=""))  
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
            'Fup',
            'log.C.thresh',
            'log.const.analytic.sd',
            'log.hetero.analytic.slope',
            'log.calibration'
            ))
            #,
            #"Response.pred",
            #"Conc"))

        sim.mcmc <- coda.out[[this.compound]]$mcmc[[1]]
        for (i in 2:NUM.CHAINS) sim.mcmc <- rbind(sim.mcmc,coda.out[[this.compound]]$mcmc[[i]])
        results <- apply(sim.mcmc,2,function(x) signif(quantile(x,c(0.025,0.5,0.975)),3))
    
        new.results <- t(data.frame(c(this.compound,this.dtxsid,this.lab.name),stringsAsFactors=F))
        colnames(new.results) <- c("Compound","DTXSID","Lab.Compound.Name")
        new.results <- cbind.data.frame(new.results,
          t(as.data.frame(as.numeric(results[c(2,1,3),"Fup"]))))
        colnames(new.results)[4:6] <- c(
          "Fup.Med",
          "Fup.Low",
          "Fup.High")
        new.results[,"Fup.point"] <- signif(mean(AF.data[,"Response"] *
          AF.data[,"Dilution.Factor"]) / mean(T5.data[,"Response"] *
          T5.data[,"Dilution.Factor"]),3)
        rownames(new.results) <- this.compound
    
        print(mydata$Num.obs)
        print(mydata$Response.obs)
        print(results)
     
 #       if (results[1,"Fup"]<1e-8) browser()
    
        Results <- rbind(Results,new.results)
    
        write.table(Results, 
          file=paste(OUTPUT.FILE,sep=""),
          sep="\t",
          row.names=F,
          quote=F)
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
    file=paste(FILENAME,"-PPB-UC-Level2-ignoredbayes.tsv",sep=""),
    sep="\t",
    row.names=F,
    quote=F)
  
  View(Results)
  save(Results,
    file=paste(FILENAME,"-UC-Fup-Level4Analysis-",Sys.Date(),".RData",sep=""))

  return(list(Results=Results,coda=coda.out))  
}


