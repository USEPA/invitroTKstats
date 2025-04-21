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

#' Calculate Fraction Unbound in Plasma (Fup) from Ultracentrifugation (UC) Data
#' with Bayesian Modeling (Level-4)
#'
#' This function estimates the fraction unbound in plasma (Fup) and credible
#' intervals with a Bayesian modeling approach, via MCMC simulations.
#' Data used in modeling is collected from Ultracentrifugation (UC) Fup assays.
#' Fup and the credible interval are calculated from the MCMC posterior samples
#' and the function returns a summary table (level-4) along with the full set of
#' MCMC results.
#' 
#' The input to this function should be "level-2" data. Level-2 data is level-1,
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
#' We currently require CC, AF, and T5 data. T1 data are optional.
#' 
#' Note: runjags::findjags() may not work as \code{JAGS.PATH} argument. Instead, 
#' may need to manually remove the trailing path such that \code{JAGS.PATH} only 
#' contains path information through "/x64" (e.g. \code{JAGS.PATH} = "/Program Files/JAGS/JAGS-4.3.1/x64").
#'
#' @param FILENAME (Character) A string used to identify the input level-2 file,
#' "<FILENAME>-fup-UC-Level2.tsv", and to name the exported model results. 
#' This argument is required no matter which method of specifying input data is used. 
#' (Defaults to \code{NULL}.)
#' 
#' @param data.in A level-2 data frame generated from the 
#' \code{format_fup_uc} function with a verification column added by 
#' \code{sample_verification}. Complement with manual verification if needed.
#'
#' @param TEMP.DIR (Character) Temporary directory to save intermediate files. If 
#' \code{NULL}, all files will be written to the current working directory.
#' (Defaults to \code{NULL}.)
#'
#' @param NUM.CHAINS (Numeric) The number of Markov Chains to use. (Defaults to 5.)
#'
#' @param NUM.CORES (Numeric) The number of processors to use for
#' parallel computing. (Defaults to 2.)
#'
#' @param RANDOM.SEED (Numeric) The seed used by the random number generator.
#' (Defaults to 1111.)
#' 
#' @param good.col (Character) Column name indicating which rows have been
#' verified for analysis, valid data rows are indicated with "Y".
#' (Defaults to "Verified".)
#' 
#' @param JAGS.PATH (Character) Computer specific file path to JAGS software.
#' (Defaults to `NA`.)
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
#' If \code{NULL}, the output file will be saved to the current working
#' directory. (Defaults to \code{NULL}.)
#' 
#' @return A list of two objects: 
#' \enumerate{
#'    \item{Results: A level-4 data frame with Bayesian estimated fraction unbound
#'    in plasma (Fup) and credible intervals for all compounds in the input file.
#'    Column includes:
#'    Compound.Name - compound name,
#'    Lab.Compound.Name - compound name used by the laboratory,
#'    DTXSID - EPA's DSSTox Structure ID,
#'    Fup.point - point estimate of Fup,
#'    Fup.Med - posterior median,
#'    Fup.Low - 2.5th quantile,
#'    Fup.High - 97.5th quantile,
#'    Fstable.Med - posterior median of stability fraction,
#'    Fstable.Low - 2.5th quantile,
#'    Fstable.High - 97.5th quantile.}
#'    \item{coda: A runjags-class object containing results from JAGS model.}
#' }
#'
#' @author John Wambaugh and Chantel Nicolas
#' 
#' @examples 
#' ## Example 1: loading level-2 using data.in
#' \dontrun{
#' level2 <- invitroTKstats::fup_uc_L2
#' 
#' # JAGS.PATH should be changed to user's specific computer file path to JAGS software.
#' # findJAGS() from runjags package is a handy function to find JAGS path automatically.
#' # In certain circumstances or cases, one may need to provide the absolute path to JAGS.
#' path.to.JAGS <- runjags::findJAGS()
#' level4 <- calc_fup_uc(FILENAME = "Example1",
#'                       data.in = level2,
#'                       NUM.CORES=2,
#'                       JAGS.PATH=path.to.JAGS)
#' }
#' 
#' ## Example 2: importing level-2 from a .tsv file
#' \dontrun{
#' # Refer to sample_verification help file for how to export level-2 data to a directory.
#' # JAGS.PATH should be changed to user's specific computer file path to JAGS software.
#' # findJAGS() from runjags package is a handy function to find JAGS path automatically.
#' # In certain circumstances or cases, one may need to provide the absolute path to JAGS.
#' # Will need to replace FILENAME and INPUT.DIR with name prefix and location of level-2 'tsv'.
#' path.to.JAGS <- runjags::findJAGS()
#' level4 <- calc_fup_uc(# e.g. replace with "Examples" from "Examples-fup-UC-Level2.tsv"
#'                       FILENAME="<level-2 FILENAME prefix>",
#'                       NUM.CORES=2,
#'                       JAGS.PATH=path.to.JAGS,
#'                       INPUT.DIR = "<level-2 FILE LOCATION>")
#' }
#' 
#' @import parallel 
#' @import runjags
#' @import coda
#' @import Rdpack
#' 
#'
#' @export calc_fup_uc
calc_fup_uc <- function(
  FILENAME,
  data.in,
  TEMP.DIR = NULL,
  NUM.CHAINS=5, 
  NUM.CORES=2,
  RANDOM.SEED=1111,
  good.col="Verified",
  JAGS.PATH = NA,
  save.MCMC = FALSE,
  sig.figs = 3, 
  INPUT.DIR=NULL, 
  OUTPUT.DIR = NULL
  )
{
  
  #assigning global variables
  Compound.Name <- Response <- NULL
  
  if (!missing(data.in)) {
    if (missing(FILENAME)) stop("FILENAME is required to save the model results. Please provide input for this argument.")
    PPB.data <- as.data.frame(data.in)
    } else if (!is.null(INPUT.DIR)) {
      PPB.data <- read.csv(file=paste0(INPUT.DIR, "/",FILENAME,"-fup-UC-Level2.tsv"), 
                         sep="\t",header=T)  
      } else {
        PPB.data <- read.csv(file=paste0(FILENAME,"-fup-UC-Level2.tsv"), 
                         sep="\t",header=T)  
        }
  
  PPB.data <- subset(PPB.data,!is.na(Compound.Name))
  PPB.data <- subset(PPB.data,!is.na(Response))
  
  
  if (!is.null(TEMP.DIR)) 
  {
    current.dir <- getwd()
    setwd(TEMP.DIR)
  }
  
  fup.uc.cols <- c(L1.common.cols,
                   test.conc.col = "Test.Compound.Conc",
                   test.nominal.conc.col = "Test.Nominal.Conc"
  )
  list2env(as.list(fup.uc.cols), envir = environment())
  cols <- c(unlist(mget(names(fup.uc.cols))), "Response", good.col)
  
  if (!any(c("Biological.Replicates", "Technical.Replicates") %in% colnames(PPB.data)))
    stop("Need at least one column representing replication, i.e. Biological.Replicates or Technical.Replicates. Run format_fup_uc first (level-1) then curate to (level-2).")
  
  if (!(all(cols %in% colnames(PPB.data))))
  {
    warning("Run format_fup_uc first (level-1) then curate to level-2.")
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
  # Round unverified data 
  if (!is.null(sig.figs)){
    unverified.data[,"Area"] <- signif(unverified.data[,"Area"], sig.figs)
    unverified.data[,"ISTD.Area"] <- signif(unverified.data[,"ISTD.Area"], sig.figs)
    unverified.data[,"Response"] <- signif(unverified.data[,"Response"], sig.figs)
    cat(paste0("\nHeldout L2 data to export has been rounded to ", sig.figs, " significant figures.\n"))
  }
  write.table(unverified.data, file=paste0(
    FILENAME,"-fup-UC-Level2-heldout.tsv"),
    sep="\t",
    row.names=F,
    quote=F)
  PPB.data <- subset(PPB.data, PPB.data[,good.col] == "Y")
  
  PPB.data <- as.data.frame(PPB.data)
  all.blanks <- subset(PPB.data,!is.na(eval(area.col)))
  
  OUTPUT.FILE <- paste0(FILENAME,"-fup-UC-Level4.tsv")

  set.seed(RANDOM.SEED)
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
    
      for (this.series in unique(MS.data[,"Biological.Replicates"]))
        if (!is.na(this.series))
        {
          this.series.subset <- subset(MS.data,MS.data[,"Biological.Replicates"]==this.series)
          for (this.cal in unique(this.series.subset[,cal.col]))
            if (!is.na(this.cal))
            {
              this.cal.subset <- subset(this.series.subset,          
                                      this.series.subset[,cal.col]==this.cal)            
              if (!all(c("T1","T5","AF") %in% this.cal.subset[,type.col]))
              {
                # Have to handle the NA series values for CC data:
                series.values <- MS.data[,"Biological.Replicates"]
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
          file=paste0(FILENAME,"-Fup-UC-PREJAGS.RData"))  
        
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
        results <- apply(sim.mcmc,2,function(x) quantile(x,c(0.025,0.5,0.975)))
    
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
        new.results[,"Fup.point"] <- mean(AF.data[,"Response"] *
          AF.data[,"Dilution.Factor"]) / mean(T5.data[,"Response"] *
          T5.data[,"Dilution.Factor"])
        rownames(new.results) <- this.compound
    
        # round results and new.results for printing
        rounded.results <- results
        rounded.new.results <- new.results 
        
        if (!is.null(sig.figs)){
          for (this.col in 1:ncol(rounded.results)){
            rounded.results[,this.col] <- signif(rounded.results[,this.col], sig.figs)
          }
          round.cols <- colnames(rounded.new.results)[!colnames(rounded.new.results) %in% c("Compound","DTXSID","Lab.Compound.Name")]
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
      } else {
        ignored.data <- rbind(ignored.data, MS.data)
      }   
    }
  
  if (!is.null(TEMP.DIR)) 
  {
    setwd(current.dir)
  }
  stopCluster(CPU.cluster)

  #View(Results)
  
  # Write out a "level-4" result table:
  # Determine the path for output
  if (!is.null(OUTPUT.DIR)) {
    file.path <- OUTPUT.DIR
    } else if (!is.null(INPUT.DIR)) {
      file.path <- INPUT.DIR
      } else {
        file.path <- getwd()
      }

  save(Results,
    file=paste0(file.path, "/", FILENAME,"-fup-UC-Level4Analysis-",Sys.Date(),".RData"))
  cat(paste0("A level-4 file named ",FILENAME,"-fup-UC-Level4Analysis-",Sys.Date(),".RData", 
             " has been exported to the following directory: ", file.path), "\n")
    
  # Save ignored data if there is any
  if (!is.null(ignored.data)) {
    write.table(ignored.data, 
                file=paste0(file.path, "/", FILENAME,"-fup-UC-Level2-ignoredbayes.tsv"),
                sep="\t",
                row.names=F,
                quote=F)
    cat(paste0("A subset of ignored data named ",FILENAME,"-fup-UC-Level2-ignoredbayes.tsv", 
               " has been exported to the following directory: ", file.path), "\n")
    }
    
  # Write out the MCMC results separately 
  if (save.MCMC){
    if (length(coda.out) != 0) {
      save(coda.out,
           file=paste0(file.path, "/", FILENAME,"-fup-UC-Level4-MCMC-Results-",Sys.Date(),".RData"))
      } else {
        cat("No MCMC results to be saved.\n")
      }
    }
  
  return(list(Results=Results,coda=coda.out))  
}


