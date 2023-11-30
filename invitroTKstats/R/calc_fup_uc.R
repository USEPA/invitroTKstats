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

#' Calculate the Fraction Unbound in Plasma (Fup) from Ultracentrifugation (UC) Data
#'
#' This function estimates the fraction unbound in plasma (Fup) and credible
#' intervals with a Bayesian modeling approach, via MCMC simulations.
#' Data used in modeling is collected from Ultracentrifugation (UC) Fup assays.
#' Fup and the credible interval are calculated from the MCMC posterior samples
#' and the function returns a summary table along with the full set of
#' MCMC results.
#' 
#' The input to this function should be "Level-2" data. Level-2 data is Level-1,
#' data formatted with the \code{\link{format_fup_uc}} function, and curated
#' with a verification column. "Y" in the verification column indicates the
#' data row is valid for analysis. 
#' 
#' Note: By default, this function writes files to the user's current working
#' directory. Users must specify an alternative path with the `TEMP.DIR`
#' argument if they want the files exported to another path. Exported files 
#' include the summary table (.RData), JAGS model (.RData), and any "unverified" 
#' data excluded from the analysis (.tsv).
#' 
#' The data frame of observations should be annotated according to
#' these types:
#' \tabular{rrrrr}{
#'   Calibration Curve \tab CC\cr
#'   Ultracentrifugation Aqueous Fraction \tab AF\cr
#'   Whole Plasma T1h Sample  \tab T1\cr
#'   Whole Plasma T5h Sample \tab T5\cr
#' }
#' We don't currently use the T1 data, but CC, AF, and T5 data are required.
#'
#' @param FILENAME (Character) A string used to identify the input Level-2 file.
#' "<FILENAME>-fup-UC-Level2.tsv". (Defaults to "UC_Model_Results")
#'
#' @param TEMP.DIR (Character) Alternative directory to save output files. By
#' default, i.e. unspecified, all files will be exported to the user's current
#' working directory. (Defaults to `NULL`.)
#'
#' @param NUM.CHAINS (Numeric) The number of Markov Chains to use. (Defaults to 5.)
#'
#' @param NUM.CORES (Numeric) The number of computer processors to use for
#' parallel computing. (Defaults to 2.)
#'
#' @param RANDOM.SEED (Numeric) The seed used by the random number generator.
#' (Defaults to 1111.)
#' 
#' @param good.col (Character) Column name indicating which rows have been
#' verified, data rows valid for analysis are indicated with a "Y".
#' (Defaults to "Verified".)
#' 
#' @param JAGS.PATH (Character) Computer specific file path to JAGS software.
#' (Defaults to `NA`.)
#' 
#' @return A list of two objects: 
#' \enumerate{
#'    \item{Results: A data frame with Bayesian estimated fraction unbound
#'    in plasma (Fup) and credible intervals for all compounds in the input file.
#'    Column includes:
#'    Compound.Name - compound name,
#'    Lab.Compound.Name - compound name used by the laboratory,
#'    DTXSID - EPA's DSSTox Structure ID,
#'    Fup.point - point estimate of Fup,
#'    Fup.Med - Posterior median,
#'    Fup.Low - 2.5th quantile,
#'    and Fup.High - 97.5th quantile.}
#'    \item{coda: A runjags-class object containing results from JAGS model.}
#' }
#'
#' @author John Wambaugh and Chantel Nicolas
#' 
#' @examples 
#' # Level-2 file
#' write.table(level2,
#'   file="KreutzPFAS-fup-UC-Level2.tsv",
#'   sep="\t",
#'   row.names=F,
#'   quote=F)
#' 
#' # JAGS.PATH should be changed to user's specific computer file path to JAGS software.
#' level4 <- calc_fup_uc(FILENAME="KreutzPFAS",
#'                        NUM.CORES=8,
#'                        JAGS.PATH="C:/Users/jwambaug/AppData/Local/JAGS/JAGS-4.3.0/x64")
#' 
#' @import parallel 
#' @import runjags
#' 
#' @import coda
#'
#' @import Rdpack
#'
#' @export calc_fup_uc
calc_fup_uc <- function(
  FILENAME = "UC_Model_Results",
  TEMP.DIR = NULL,
  NUM.CHAINS=5, 
  NUM.CORES=2,
  RANDOM.SEED=1111,
  good.col="Verified",
  JAGS.PATH = NA
  )
{
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
       
        CC.data <- MS.data[MS.data[,type.col]=="CC",]
        T1.data <- MS.data[MS.data[,type.col]=="T1",]
        T5.data <- MS.data[MS.data[,type.col]=="T5",]
        AF.data <- MS.data[MS.data[,type.col]=="AF",]
        mydata <- build_mydata_fup_uc(MS.data, CC.data, T1.data, T5.data, AF.data)
        
        init_vals <- function(chain) initfunction_fup_uc(mydata=mydata, chain = chain)
        # write out arguments to runjags:
        save(this.compound,mydata,UC_PPB_model,init_vals,
          file=paste(FILENAME,"-Fup-UC-PREJAGS.RData",sep=""))  
        
        coda.out[[this.compound]] <- autorun.jags(
          UC_PPB_model, 
          n.chains = NUM.CHAINS,
          method="parallel", 
          cl=CPU.cluster,
          summarise=T,
          inits = init_vals,
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


