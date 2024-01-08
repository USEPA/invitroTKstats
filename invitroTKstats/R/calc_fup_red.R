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

#' Calculate Fraction Unbound in Plasma (Fup) From Rapid Equilibrium Dialysis (RED) Data
#'
#' This function estimates the fraction unbound in plasma (Fup) with Bayesian
#' modeling on Rapid Equilibrium Dialysis (RED) data. Both Fup and
#' the credible interval are estimated from posterior samples of the MCMC.
#' A summary table along with the full set of MCMC results is returned from
#' the function.
#' 
#' The input data should be "Level-2". Level-2 data is Level-1 data, formatted 
#' with the \code{\link{format_fup_red}} function, and curated with a
#' verification column where "Y" indicates a data row is valid for analysis. 
#' 
#' Note: By default, this function writes files to the user's current working directory 
#' unless another path is specified in the \code{TEMP.DIR} argument. Files saved include 
#' the summary table (.RData), JAGS model (.RData), and any "unverified" data 
#' excluded from the analysis (.tsv).
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
#' @param FILENAME (Character) A string used to identify the input Level-2 file.
#' "<FILENAME>-fup-RED-Level2.tsv".
#' 
#' @param data.in (Data Frame) A Level-2 data frame containing
#' mass-spectrometry peak areas, indication of chemical identity,
#' and measurement type.
#'
#' @param TEMP.DIR (Character) Temporary directory to save intermediate files. 
#' If \code{NULL}, all files will be written to the current working directory. 
#' (Defaults to \code{NULL}.)
#'
#' @param JAGS.PATH (Character) Computer specific file path to JAGS software. (Defaults to \code{NA}.)
#'
#' @param NUM.CHAINS (Numeric) The number of Markov Chains to use. (Defaults to 5.)
#'
#' @param NUM.CORES (Numeric) The number of processors to use for parallel computing. (Defaults to 2.)
#'
#' @param RANDOM.SEED The seed used by the random number generator.
#' (Defaults to 1111.)
#'
#' @param good.col (Character) Column name indicating which rows have been 
#' verified for analysis, valid data rows are indicated with "Y". (Defaults to "Verified".)
#'
#' @param Physiological.Protein.Conc (Numeric) The assumed physiological protein concentration 
#' for plasma protein binding calculations. (Defaults to 70/(66.5*1000)*1000000.
#' According to Berg and Lane (2011): 60-80 mg/mL, albumin is 66.5 kDa,
#' assume all protein is albumin to estimate default in uM.) 
#' 
#' @param output.res (Logical) When set to \code{TRUE}, the result 
#' table (Level-4) will be exported as a .RData file. 
#' (Defaults to \code{TRUE}.)
#' 
#' @param save.MCMC (Logical) When set to \code{TRUE}, will export the MCMC results
#' as an .RData file. (Defaults to \code{FALSE}.)
#' 
#' @param INPUT.DIR (Character) Path to the directory where the input level-2 file exists. 
#' If \code{NULL}, looking for the input level-2 file in the current working
#' directory. (Defaults to \code{NULL}.)
#' 
#' @param OUTPUT.DIR (Character) Path to the directory to save the output file. 
#' If \code{NULL}, the output file will be saved to the current working
#' directory. (Defaults to \code{NULL}.)
#'
#' @return A list of two objects: 
#' \enumerate{
#'    \item{Results: A data frame with the Bayesian estimated fraction unbound in plasma (Fup) 
#'    and credible interval for all compounds in the input file. Column includes:
#'    Compound.Name - compound name, Lab.Compound.Name - compound name used by 
#'    the laboratory, DTXSID - EPA's DSSTox Structure ID, Fup.point - point estimate of Fup,
#'    Fup.Med - posterior median, Fup.Low - 2.5th quantile, and Fup.High - 97.5th quantile}
#'    \item{coda: A runjags-class object containing results from JAGS model.}
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
#' # JAGS.PATH should be changed to user's specific computer file path to JAGS software.
#' level4 <- calc_fup_red(FILENAME="SmeltzPFAS",
#'                        NUM.CORES=8,
#'                        JAGS.PATH="C:/Users/jwambaug/AppData/Local/JAGS/JAGS-4.3.0/x64")
#'
#' @import coda
#'
#' @import Rdpack
#'
#' @export calc_fup_red
calc_fup_red <- function(
  FILENAME,
  data.in,
  TEMP.DIR = NULL,
  NUM.CHAINS=5,
  NUM.CORES=2,
  RANDOM.SEED=1111,
  good.col="Verified",
  JAGS.PATH = NA,
  Physiological.Protein.Conc = 70/(66.5*1000)*1000000, # Berg and Lane (2011) 60-80 mg/mL, albumin is 66.5 kDa, pretend all protein is albumin to get uM
  output.res = TRUE,
  save.MCMC = FALSE,
  INPUT.DIR=NULL, 
  OUTPUT.DIR = NULL
  )
{
  
  if (!missing(data.in)) {
    MS.data <- as.data.frame(data.in)
    } else if (!is.null(INPUT.DIR)) {
      MS.data <- read.csv(file=paste0(INPUT.DIR, "/", FILENAME,"-fup-RED-Level2.tsv"),
                        sep="\t",header=T)
      } else {
        MS.data <- read.csv(file=paste0(FILENAME,"-fup-RED-Level2.tsv"),
                        sep="\t",header=T)
        }
  
  MS.data <- subset(MS.data,!is.na(Compound.Name))
  MS.data <- subset(MS.data,!is.na(Response))
  
  if (!is.null(TEMP.DIR))
  {
    current.dir <- getwd()
    setwd(TEMP.DIR)
  }

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
  write.table(unverified.data, file=paste0(
    FILENAME,"-fup-RED-Level2-heldout.tsv"),
    sep="\t",
    row.names=F,
    quote=F)
  MS.data <- subset(MS.data, MS.data[,good.col] == "Y")

  # Clean up data:
  MS.data <- subset(MS.data,!is.na(Response))
  MS.data[MS.data$Response<0,"Response"] <- 0

  # Because of the possibility of crashes we save the results one chemical at a time:
  OUTPUT.FILE <- paste0(FILENAME,"-fup-RED-Level4.tsv")

  # Check to see if we crashed earlier, if so, don't redo something that already is done
  if (!file.exists(OUTPUT.FILE))
  {
    Results <- NULL
  } else {
    Results <- read.table(OUTPUT.FILE,sep="\t",stringsAsFactors=F,header=T)
  }
  
  # Safety check for parallel computation 
  MAX.CORES <- detectCores(logical = F) - 1
  if (NUM.CORES > MAX.CORES) stop(paste0("Specified NUM.CORES = ", NUM.CORES, " cores for parallel computing exceeds the allowable number of cores, that is ",
                                         MAX.CORES, 
                                         ", and may bog down your machine! (Max cores is based on the total number of available computing cores minus one for overhead.)"))
  if (NUM.CHAINS > 10) warning("Specified number of chains is greater than 10 and may be excessive for computational time.")

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
        mydata <- build_mydata_fup_red(this.subset, Physiological.Protein.Conc)
        if (!is.null(mydata))
        {
          # Use random number seed for reproducibility
          set.seed(RANDOM.SEED)

          init_vals <- function(chain) initfunction_fup_red(mydata=mydata, chain = chain)
          # write out arguments to runjags:
          save(this.compound, mydata ,init_vals,
            file=paste0(FILENAME,"-fup-RED-PREJAGS.RData"))

          # Run JAGS:
          coda.out[[this.compound]] <- autorun.jags(
            fup_RED_model,
            n.chains = NUM.CHAINS,
            method="parallel",
            cl=CPU.cluster,
            summarise=T,
            inits = init_vals,
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
            file=paste0(OUTPUT.FILE),
            sep="\t",
            row.names=F,
            quote=F)
        }
      } else {
        ignored.data <- rbind(ignored.data, this.subset)
      }
    }

  if (!is.null(TEMP.DIR))
  {
    setwd(current.dir)
  }
  stopCluster(CPU.cluster)

  View(Results)
  
  if (output.res) {
    # Write out a "level 4" result table:
    # Determine the path for output
    if (!is.null(OUTPUT.DIR)) {
      file.path <- OUTPUT.DIR
    } else if (!is.null(INPUT.DIR)) {
      file.path <- INPUT.DIR
    } else {
      file.path <- getwd()
    }
    
    save(Results,
      file=paste0(file.path, "/", FILENAME,"-fup-RED-Level4Analysis-",Sys.Date(),".RData"))
    
    cat(paste0("A Level-4 file named ",FILENAME,"-fup-RED-Level4Analysis-",Sys.Date(),".RData", 
                " has been exported to the following directory: ", file.path), "\n")
   
    # Save ignored data if there is any
    if (!is.null(ignored.data)) {
      write.table(ignored.data,
                file=paste0(file.path, "/", FILENAME,"-fup-RED-Level2-ignoredbayes.tsv"),
                sep="\t",
                row.names=F,
                quote=F)
      cat(paste0("A subset of ignored data named ",FILENAME,"-fup-RED-Level2-ignoredbayes.tsv", 
               " has been exported to the following directory: ", file.path), "\n")
    
    }
    
    # Write out the MCMC results separately 
    if (save.MCMC){
      if (length(coda.out) != 0) {
        save(coda.out,
             file=paste0(file.path, "/", FILENAME,"-fup-RED-Level4-MCMC-Results-",Sys.Date(),".RData"))
      } else {
        cat("No MCMC results to be saved.\n")
      }
    }
  }
  
  return(list(Results=Results,coda=coda.out))
}


