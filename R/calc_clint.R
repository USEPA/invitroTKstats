Clint_model <- "
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

# Likelihood for the blank observations:
  # obtain the blank predictions and precisions for each calibration
  for (i in 1:Num.cal)
  {
    Blank.pred[i] <- background[i]/Blank.Dilution.Factor[i]
    Blank.prec[i] <- (const.analytic.sd[i]+hetero.analytic.slope[i]*(Blank.pred[i]))^(-2)
  }
  # obtain the estimated observations for each blank sample available
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
      (cc.obs.conc[i]/cc.obs.Dilution.Factor[i] -
      C.thresh[cc.obs.cal[i]]/cc.obs.Dilution.Factor[i]) *
      step(cc.obs.conc[i]/cc.obs.Dilution.Factor[i] -
      C.thresh[cc.obs.cal[i]]/cc.obs.Dilution.Factor[i]) +
      cc.intercept[i]/cc.obs.Dilution.Factor[i]
# Heteroskedastic precision:
    cc.prec[i] <- (const.analytic.sd[cc.obs.cal[i]] +
      hetero.analytic.slope[cc.obs.cal[i]] * cc.pred[i])^(-2)
# Model for the observation:
    cc.obs[i] ~ dnorm(cc.pred[i],cc.prec[i])
  }

# Clearance model:

# Decreases indicates whether the concentration decreases (1 is yes, 0 is no):
  decreases ~ dbern(DECREASE.PROB)
# Degrades indicates whether we think abiotic degradation is a factor
# (1 is yes, 0 is no), p=0.05 assumes 5 percent of chemicals degrade:
  degrades ~ dbern(DEGRADE.PROB)
# Slope is the clearance rate at the lower concentration:
  bio.rate ~ dunif(0,-log(C.thresh[1]/Test.Nominal.Conc[1]))
# In addition to biological elimination, we also check for abiotic elimination
# (for example, degradation) which we can distinguish if we have data from
# inactivated hepatocutes (Num.abio.obs > 0):
  abio.rate ~ dunif(0,-log(C.thresh[1]/Test.Nominal.Conc[1]))
# Total elimination rate is a sum of both:
  slope[1] <- decreases * bio.rate + degrades * abio.rate
# Actual biological elimination rate:
  bio.slope[1] <- decreases * bio.rate
# Saturates is whether the bio clearance rate decreases (1 is yes, 0 is no)
# p=0.25 assumes 25 percent of chemicals show some metabolic saturation between
# 1 and 10 uM:
  saturates ~ dbern(SATURATE.PROB)
# Saturation is how much the bio clearance rate decreases at the higher conc:
  saturation ~ dunif(0,1)
# Calculate a slope at the higher concentration:
  slope[2] <- degrades * abio.rate +
    decreases * bio.rate * (1 - saturates*saturation)
# Actual biological elimination rate:
  bio.slope[2] <- decreases * bio.rate * (1 - saturates*saturation)
# The observations are normally distributed (heteroskedastic error):
  for (i in 1:Num.obs)
  {
  # Exponential decay:
    C[i] <- Test.Nominal.Conc[obs.conc[i]] *
      exp(-slope[obs.conc[i]]*obs.time[i])
  # MS prediction:
    obs.pred[i] <- calibration[obs.cal[i]] *
      (C[i]/obs.Dilution.Factor[i] -
      C.thresh[obs.cal[i]]/obs.Dilution.Factor[i]) *
      step(C[i]/obs.Dilution.Factor[i] -
      C.thresh[obs.cal[i]]/obs.Dilution.Factor[i]) +
      background[obs.cal[i]]/obs.Dilution.Factor[i]
    obs.prec[i] <- (const.analytic.sd[obs.cal[i]] +
      hetero.analytic.slope[obs.cal[i]]*obs.pred[i])^(-2)

    obs[i] ~ dnorm(obs.pred[i], obs.prec[i])
  }

# non-biological clearance model:
  abio.slope <- degrades * abio.rate
# The observations are normally distributed (heteroskedastic error):
  for (i in 1:Num.abio.obs)
  {
  # Exponential decay:
    abio.C[i] <- Test.Nominal.Conc[abio.obs.conc[i]] *
      exp(-abio.slope*abio.obs.time[i])
  # MS prediction:
    abio.obs.pred[i] <- calibration[abio.obs.cal[i]] *
      (abio.C[i]/abio.obs.Dilution.Factor[i] -
      C.thresh[abio.obs.cal[i]]/abio.obs.Dilution.Factor[i]) *
      step(abio.C[i]/abio.obs.Dilution.Factor[i] -
      C.thresh[abio.obs.cal[i]]/abio.obs.Dilution.Factor[i]) +
      background[abio.obs.cal[i]]/abio.obs.Dilution.Factor[i]
    abio.obs.prec[i] <- (const.analytic.sd[abio.obs.cal[i]] +
      hetero.analytic.slope[abio.obs.cal[i]]*abio.obs.pred[i])^(-2)

    abio.obs[i] ~ dnorm(abio.obs.pred[i], abio.obs.prec[i])
  }
}
"

#' Calculate Intrinsic Hepatic Clearance (Clint) with Bayesian Modeling (Level-4)
#'
#' This function estimates the intrinsic hepatic clearance (Clint) with Bayesian
#' modeling on Hepatocyte Incubation data \insertCite{shibata2002prediction}{invitroTKstats}.
#' Clint and the credible intervals,
#' at both 1 and 10 uM (if tested), are estimated from posterior samples of the MCMC.
#' A summary table (level-4) along with the full set of MCMC results is returned from
#' the function.
#' 
#' The input to this function should be "level-2" data. Level-2 data is level-1,
#' data formatted with the \code{\link{format_clint}} function, and curated
#' with a verification column. "Y" in the verification column indicates the
#' data row is valid for analysis. 
#' 
#' Note: By default, this function writes files to the user's per-session temporary 
#' directory. This temporary directory is a per-session directory whose path can 
#' be found with the following code: \code{tempdir()}. For more details, see
#' \url{https://www.collinberke.com/til/posts/2023-10-24-temp-directories/}.
#' 
#' Users must specify an alternative path with the \code{TEMP.DIR}
#' argument if they want the intermediate files exported to another path. Exported 
#' intermediate files include the summary results table (.tsv), JAGS model (.RData), and any "unverified" 
#' data excluded from the analysis (.tsv). Users must specify an alternative path with the \code{OUTPUT.DIR} argument if they 
#' want the final output file exported to another path. The exported final output
#' file is the summary results table (.RData). 
#' 
#' As a best practice, \code{INPUT.DIR} (when importing a .tsv file) and/or \code{OUTPUT.DIR} 
#' should be specified to simplify the process of importing and exporting files. 
#' This practice ensures that the exported files can easily be found and will 
#' not be exported to a temporary directory.
#'
#' The data frame of observations should be annotated according to
#' these types:
#' \tabular{rl}{
#'   Blank \tab Cell free blank with media\cr
#'   CC \tab Cell free calibration curve \cr
#'   Cvst \tab Hepatocyte incubation concentration vs. time \cr
#'   Inactive \tab Concentration vs. time data with inactivated hepatocytes \cr
#' }
#' We currently require Cvst data. Blank, CC, and Inactive data are optional.
#'
#' Clint is calculated using \code{\link{lm}} to perform a linear regression of
#' MS response as a function of time.
#' 
#' Additional User Notification(s):
#' 
#' \itemize{
#'    \item{runjags::findjags() may not work as \code{JAGS.PATH} argument. Instead,
#'    may need to manually remove the trailing path such that \code{JAGS.PATH} only
#'    contains path information through "/x64" (e.g. \code{JAGS.PATH} = "/Program Files/JAGS/JAGS-4.3.1/x64").}
#' }  
#'
#' @param FILENAME (Character) A string used to identify the input level-2 file,
#' "<FILENAME>-Clint-Level2.tsv", and to name the exported model results. 
#' This argument is required no matter which method of specifying input data is used. 
#' (Defaults to \code{NULL}.)
#'
#' @param data.in (Data Frame) A level-2 data frame generated from the 
#' \code{format_clint} function with a verification column added by 
#' \code{sample_verification}. Complement with manual verification if needed.
#'
#' @param TEMP.DIR (Character) Temporary directory to save intermediate files. 
#' If \code{NULL}, all files will be written to the user's per-session
#' temporary directory. 
#' (Defaults to \code{NULL}.)
#' 
#' @param NUM.CHAINS (Numeric) The number of Markov Chains to use. (Defaults to 5.)
#'
#' @param NUM.CORES (Numeric) The number of processors to use for parallel computing. (Defaults to 2.)
#'
#' @param RANDOM.SEED (Numeric) The seed used by the random number generator.
#' (Defaults to 1111.)
#'
#' @param good.col (Character) Column name indicating which rows have been 
#' verified for analysis, valid data rows are indicated with "Y". (Defaults to "Verified".)
#' 
#' @param JAGS.PATH (Character) Computer specific file path to JAGS software. 
#' (Defaults to \code{NA}.)
#' 
#' @param decrease.prob (Numeric) Prior probability that a chemical will decrease in
#' the assay. (Defaults to 0.5.)
#'
#' @param saturate.prob (Numeric) Prior probability that a chemicals rate of metabolism
#' will decrease between 1 and 10 uM. (Defaults to 0.25.)
#'
#' @param degrade.prob (Numeric) Prior probability that a chemical will be unstable
#' (that is, degrade abiotically) in the assay. (defaults to 0.05.)
#' 
#' @param save.MCMC (Logical) When set to \code{TRUE}, will export the MCMC results 
#' as an .RData file. (Defaults to \code{FALSE}.)
#' 
#' @param sig.figs (Numeric) The number of significant figures to round the exported unverified data (level-2).
#' The exported result table (level-4) is left unrounded for reproducibility.
#' (Note: console print statements are also rounded to specified significant figures.)
#' (Defaults to \code{3}.)
#' 
#' @param INPUT.DIR (Character) Path to the directory where the input level-2 file exists. 
#' If \code{NULL}, looking for the input level-2 file in the current working
#' directory. (Defaults to \code{NULL}.)
#' 
#' @param OUTPUT.DIR (Character) Path to the directory to save the output file. 
#' If \code{NULL}, the output file will be saved to the user's per-session 
#' temporary directory or \code{INPUT.DIR} if specified. 
#' (Defaults to \code{NULL}.)
#'
#' @return A list of two objects: 
#' \enumerate{
#'    \item{Results: A level-4 data frame with the Bayesian estimated intrinsic hepatic clearance (Clint)
#'    for 1 and 10 uM and credible intervals for all compounds in the input file. Column includes:
#'    Compound.Name - compound name, Lab.Compound.Name - compound name used by 
#'    the laboratory, DTXSID - EPA's DSSTox Structure ID, Clint.1.Med/Clint.10.Med - posterior median, 
#'    Clint.1.Low/Clint.10.Low - 2.5th quantile, Clint.1.High/Clint.10.High - 97.5th quantile, 
#'    Clint.pValue, Sat.pValue, degrades.pValue - "p-values" estimated from the probabilities of 
#'    observing decreases, saturations, and abiotic degradations in all posterior samples.}
#'    \item{coda: A runjags-class object containing results from JAGS model.}
#' }
#'
#' @author John Wambaugh
#'
#' @examples
#' ## Example 1: loading level-2 using data.in and export all files to the user's
#' ## temporary directory
#' \dontrun{
#' level2 <- invitroTKstats::clint_L2
#' 
#' # JAGS.PATH should be changed to user's specific computer file path to JAGS software.
#' # findJAGS() from runjags package is a handy function to find JAGS path automatically.
#' # In certain circumstances or cases, one may need to provide the absolute path to JAGS.
#' path.to.JAGS <- runjags::findJAGS()
#' level4 <- calc_clint(FILENAME = "Example1",
#'                      data.in = level2,
#'                      NUM.CORES=2,
#'                      JAGS.PATH=path.to.JAGS)
#' }
#' 
#' ## Example 2: importing level-2 from a .tsv file and export all files to same 
#' ## location as INPUT.DIR 
#' \dontrun{
#' # Refer to sample_verification help file for how to export level-2 data to a directory.
#' # JAGS.PATH should be changed to user's specific computer file path to JAGS software.
#' # findJAGS() from runjags package is a handy function to find JAGS path automatically.
#' # In certain circumstances or cases, one may need to provide the absolute path to JAGS.
#' # Will need to replace FILENAME and INPUT.DIR with name prefix and location of level-2 'tsv'.
#' path.to.JAGS <- runjags::findJAGS()
#' level4 <- calc_clint(# e.g. replace with "Examples" from "Examples-Clint-Level2.tsv"
#'                      FILENAME="<level-2 FILENAME prefix>",
#'                      NUM.CORES=2,
#'                      JAGS.PATH=path.to.JAGS,
#'                      INPUT.DIR = "<level-2 FILE LOCATION>")
#' }
#' 
#' @references
#' \insertRef{shibata2002prediction}{invitroTKstats}
#'
#' @import Rdpack coda runjags parallel
#' @importFrom utils read.csv write.table read.delim
#'
#' @export calc_clint
calc_clint <- function(
  FILENAME,
  data.in,
  TEMP.DIR = NULL,
  NUM.CHAINS=5,
  NUM.CORES=2,
  RANDOM.SEED=1111,
  good.col="Verified",
  JAGS.PATH = NA,
  decrease.prob = 0.5,
  saturate.prob = 0.25,
  degrade.prob = 0.05,
  save.MCMC = FALSE,
  sig.figs = 3, 
  INPUT.DIR=NULL, 
  OUTPUT.DIR = NULL
  )
{
  
  #assigning global variables
  Compound.Name <- Response <- Sample.Type <- NULL
  
  
  if (!missing(data.in)) {
    if (missing(FILENAME)) stop("FILENAME is required to save the model results. Please provide input for this argument.")
    MS.data <- as.data.frame(data.in)
    } else if (!is.null(INPUT.DIR)) {
      MS.data <- read.csv(file=paste0(INPUT.DIR, "/", FILENAME,"-Clint-Level2.tsv"),
                          sep="\t",header=T)
      } else {
        MS.data <- read.csv(file=paste0(FILENAME,"-Clint-Level2.tsv"),
                          sep="\t",header=T)
        }
  
  MS.data <- subset(MS.data,!is.na(Compound.Name))
  MS.data <- subset(MS.data,!is.na(Response))
  
  # save the current working directory 
  current.dir <- getwd()
  
  if (!is.null(TEMP.DIR)) # set working directory to user specified TEMP.DIR
  {
    setwd(TEMP.DIR)
  } else # set working directory to per-session tempdir()
  {
    setwd(tempdir())
  }

  clint.cols <- c(L1.common.cols,
                  time.col = "Time",
                  test.conc.col = "Test.Compound.Conc",
                  test.nominal.conc.col = "Test.Nominal.Conc",
                  density.col = "Hep.Density"
  )
  list2env(as.list(clint.cols), envir = environment())
  cols <- c(unlist(mget(names(clint.cols))), "Response", good.col)
  
  # Check for missing columns
  if (!any(c("Biological.Replicates", "Technical.Replicates") %in% colnames(MS.data)))
    stop("Need at least one column representing replication, i.e. Biological.Replicates or Technical.Replicates. Run format_clint first (level-1) then curate to (level-2).")
  
  if (!(all(cols %in% colnames(MS.data))))
  {
    warning("Run format_clint first (level-1) then curate to (level-2).")
    stop(paste("Missing columns named:",
      paste(cols[!(cols%in%colnames(MS.data))],collapse=", ")))
  }

  # Only include the data types used:
  MS.data <- subset(MS.data,MS.data[,type.col] %in% c(
    "Blank","Cvst","CC","Inactive"))

  # Only used verified data:
  unverified.data <- subset(MS.data, MS.data[,good.col] != "Y")
  # Round unverified data 
  if (!is.null(sig.figs)){
    unverified.data[,"Area"] <- signif(unverified.data[,"Area"], sig.figs)
    unverified.data[,"ISTD.Area"] <- signif(unverified.data[,"ISTD.Area"], sig.figs)
    unverified.data[,"Response"] <- signif(unverified.data[,"Response"], sig.figs)
    cat(paste0("\nHeldout L2 data to export has been rounded to ", sig.figs, " significant figures.\n"))
  }
  write.table(unverified.data, file=paste0(
    FILENAME,"-Clint-Level2-heldout.tsv"),
    sep="\t",
    row.names=F,
    quote=F)
  MS.data <- subset(MS.data, MS.data[,good.col] == "Y")

  # Clean up data:
  MS.data <- subset(MS.data,!is.na(Response))
  MS.data[MS.data$Response<0,"Response"] <- 0

  # Because of the possibility of crashes we save the results one chemical at a time:
  OUTPUT.FILE <- paste(FILENAME,"-Clint-Level4.tsv",sep="")

  # Check to see if we crashed earlier, if so, don't redo something that already is done
  if (!file.exists(OUTPUT.FILE))
  {
    Results <- NULL
  } else {
    Results <- read.delim(OUTPUT.FILE,stringsAsFactors=FALSE)
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
  for (this.compound in unique(MS.data[,compound.col]))
    if (!(this.compound %in% Results[,compound.col]))
    {
      this.subset <- subset(MS.data,MS.data[,compound.col]==this.compound)
      this.dtxsid <- this.subset$DTXSID[1]
      this.lab.name <- this.subset[1,lab.compound.col]

      this.row <- c(this.subset[1,c(compound.col,dtxsid.col)],
        data.frame(Calibration="All Data",
          Clint=NaN,
          Clint.pValue=NaN))
      this.cvt <- subset(this.subset,Sample.Type=="Cvst")
      this.blank <- subset(this.subset,Sample.Type=="Blank")
    #  if (length(unique(this.cvt$Dilution.Factor))>1) browser()
    #  df.cvt <- this.cvt$Dilution.Factor[1]
      if (length(unique(this.cvt$Hep.Density))>1){
        stop("calc_clint - Cvst samples for `",this.compound,"` have more than one `Hep.Density`.")
        # browser()
      } 
      hep.density <- this.cvt$Hep.Density[1]

      # provide running output of where we are in the list:
      print(paste(
        this.compound,
        " (",
        which(unique(MS.data[,compound.col])==this.compound),
        " of ",
        length(unique(MS.data[,compound.col])),
        ")",
        sep=""))

      mydata <- build_mydata_clint(this.cvt, this.subset, decrease.prob, saturate.prob, degrade.prob)
      if (!is.null(mydata))
      {
        # Use random number seed for reproducibility
        set.seed(RANDOM.SEED)
        
        init_vals <- function(chain) initfunction_clint(mydata=mydata, chain = chain)
        # write out arguments to runjags:
        save(this.compound,mydata,init_vals,
        file=paste0(FILENAME,"-Clint-PREJAGS.RData"))

        # Run JAGS:
        coda.out[[this.compound]] <-  autorun.jags(
                           Clint_model,
                           n.chains = NUM.CHAINS,
                           method="parallel",
                           cl=CPU.cluster,
                           summarise=T,
                           inits = init_vals,
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
                           # Chemical analysis parameters:
                            'const.analytic.sd',
                            'hetero.analytic.slope',
                            'C.thresh',
                            'log.calibration',
                            'background'))

  # We don't follow the measurement parameters for convergence because they
  # are discrete and the convergence diagnostics don't work:
        coda.out[[this.compound]] <-extend.jags(coda.out[[this.compound]],
                              drop.monitor = c(
                              # Chemical analysis parameters:
                                'const.analytic.sd',
                                'hetero.analytic.slope',
                                'C.thresh',
                                'log.calibration',
                                'background'),
                              add.monitor = c(
                              # Measurement parameters
                                'bio.slope',
                                'decreases',
                                'saturates',
                                "degrades"))

        sim.mcmc <- coda.out[[this.compound]]$mcmc[[1]]
        for (i in 2:NUM.CHAINS) sim.mcmc <- rbind(sim.mcmc,coda.out[[this.compound]]$mcmc[[i]])
        results <- apply(sim.mcmc,2,function(x) quantile(x,c(0.025,0.5,0.975)))
        results <- as.data.frame(results)

        # Convert disappareance rates (1/h)to
        # heaptic clearance (uL/min/10^6 hepatocytes)):
        hep.density <- this.subset[1,density.col]

        # Calculate a Clint only for 1 and 10 uM (if tested)
        if (1 %in% mydata$Test.Nominal.Conc)
        {
          index <- which(mydata$Test.Nominal.Conc == 1)
          results[,"Clint.1"] <- 1000 *
            results[,paste("bio.slope[",index,"]",sep="")] / hep.density / 60
        } else results[,"Clint.1"] <- NA
        if (10 %in% mydata$Test.Nominal.Conc)
        {
          index <- which(mydata$Test.Nominal.Conc == 10)
          results[,"Clint.10"] <- 1000 *
            results[,paste("bio.slope[",index,"]",sep="")] / hep.density / 60
        } else results[,"Clint.10"] <- NA
        
        # Create a row of formatted results:
        new.results <- t(data.frame(c(this.compound,this.dtxsid,this.lab.name),stringsAsFactors=F))
        colnames(new.results) <- c(compound.col, dtxsid.col, lab.compound.col)
        for (this.param in c("Clint.1","Clint.10"))
        {
          new.results <- cbind.data.frame(new.results,
            t(as.data.frame(as.numeric(results[c(2,1,3), this.param]))))
          last.col <- length(new.results)
          colnames(new.results)[(last.col-2):last.col] <- c(
            paste(this.param,".Med",sep=""),
            paste(this.param,".Low",sep=""),
            paste(this.param,".High",sep=""))
        }

        # Calculate a Clint "pvalue" from probability that we observed a decrease:
        new.results[,"Clint.pValue"] <- 
          sum(sim.mcmc[,"decreases"]==0)/dim(sim.mcmc)[1]

        # Calculate a "pvalue" for saturation probability that we observed
        # a lower Clint at higher conc:
        new.results[,"Sat.pValue"] <- 
          sum(sim.mcmc[,"saturates"]==0)/dim(sim.mcmc)[1]

        # Calculate a "pvalue" for abiotic degradation:
        new.results[,"degrades.pValue"] <- 
          sum(sim.mcmc[,"degrades"]==0)/dim(sim.mcmc)[1]

        rownames(new.results) <- this.compound
        
        # round results and new.results for printing
        rounded.results <- results
        rounded.new.results <- new.results 
        
        if (!is.null(sig.figs)){
          for (this.col in 1:ncol(rounded.results)){
            rounded.results[,this.col] <- signif(rounded.results[,this.col], sig.figs)
          }
          round.cols <- colnames(rounded.new.results)[!colnames(rounded.new.results) %in% c("Compound.Name","DTXSID","Lab.Compound.Name")]
          for (this.col in round.cols){
            rounded.new.results[,this.col] <- signif(rounded.new.results[,this.col], sig.figs)
          }
        }
  
        print(paste("Final results for ",
          this.compound,
          " (",
          which(unique(MS.data[,compound.col])==this.compound),
          " of ",
          length(unique(MS.data[,compound.col])),
          ")",
          sep=""))
        print(rounded.results)
        print(rounded.new.results)

        Results <- rbind(Results,new.results)
        
        write.table(Results,
          file=paste0(OUTPUT.FILE),
          sep="\t",
          row.names=F,
          quote=F)
      }
    }
  
  # set working directory back to original 
  setwd(current.dir)

  stopCluster(CPU.cluster)

  #View(Results)
  
  if (!is.null(OUTPUT.DIR)) { # export output file to OUTPUT.DIR (OUTPUT.DIR specified) 
    file.path <- OUTPUT.DIR
  } else if (!is.null(INPUT.DIR)) { # export output file to INPUT.DIR (OUTPUT.DIR not specified)
    file.path <- INPUT.DIR
  } else { # export output file to tempdir() (OUTPUT.DIR & INPUT.DIR not specified)
    file.path <- tempdir()
  }
  
  save(Results,
       file=paste0(file.path, "/", FILENAME,"-Clint-Level4Analysis-",Sys.Date(),".RData"))
  cat(paste0("A level-4 file named ",FILENAME,"-Clint-Level4Analysis-",Sys.Date(),".RData", 
             " has been exported to the following directory: ", file.path), "\n")
  
  if (save.MCMC){
    if (length(coda.out) != 0) {
      save(coda.out,
           file=paste0(file.path, "/", FILENAME,"-Clint-Level4-MCMC-Results-",Sys.Date(),".RData"))
      } else {
        cat("No MCMC results to be saved.\n")
      }
    }
  
  return(list(Results=Results,coda=coda.out))
}


