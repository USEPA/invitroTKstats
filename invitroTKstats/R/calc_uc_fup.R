UC_PPB_model <- "
model {

# Measurement Model:
  # Mass-spec calibration:
  for (i in 1:Num.cal)
  {
    # Priors:
    log.const.analytic.sd[i] ~ dunif(-5, 0)
    log.hetero.analytic.slope.factor[i] ~ dunif(-5, 3)
    background[i] ~ dunif(0, 10000)   
    log.calibration[i] ~ dunif(-2, 2)
    # Concentrations below this value are not detectable:
    log.C.thresh[i] ~ dunif(-10, 1) 
    # Scale conversions:
    const.analytic.sd[i] <- 10^log.const.analytic.sd[i]
    hetero.analytic.slope.factor[i] <- 10^log.hetero.analytic.slope.factor[i]
    hetero.analytic.slope[i] <- hetero.analytic.slope.factor[i]*const.analytic.sd[i]
    calibration[i] <- 10^log.calibration[i]
    C.thresh[i] <- 10^log.C.thresh[i]
  }
  
  # Mass-spec observations:  
  for (i in 1:Num.obs) 
  {
    Response.pred[i] <- 
      calibration[obs.cal[i]]*
      (Conc[obs.conc[i]]/Dilution.Factor[obs.conc[i]] - C.thresh[obs.cal[i]])*
      step(Conc[obs.conc[i]]/Dilution.Factor[obs.conc[i]] - C.thresh[obs.cal[i]]) + 
      background[obs.cal[i]]
    Response.prec[i] <- (const.analytic.sd[obs.cal[i]] +
      hetero.analytic.slope[obs.cal[i]]*
      calibration[obs.cal[i]]*
      Conc[obs.conc[i]]/Dilution.Factor[obs.conc[i]])^(-2)
    Response.obs[i] ~ dnorm(Response.pred[i],Response.prec[i])
  }
  
  
# Binding Model:
  # Prior on Fup: 
  log.Fup ~ dunif(-10, 0)
  # Scale conversion:
  Fup <- 10^log.Fup

  for (i in (Num.cc.obs +1):(Num.cc.obs + Num.series)) 
  {
  # Priors for whole samples for ultra centrigugation UC):
    Conc[i] ~ dunif(0.01,1000)
  # Aqueous fraction concentrations for UC samples:
    Conc[i+Num.series] <- Fup * Conc[i]
  }   
}
"

#' This funcion calculates fraction unbound in plasma
#'
#' The data frame of observations should be annotated according to
#' of these types:
#' \tabular{rrrrr}{
#'   Calibration Curve \tab CC\cr
#'   Ultra-centrifugation Aqueous Fraction \tab UC\cr
#'   Whole Plasma T1h Sample  \tab T1\cr
#'   Whole Plasma T5h Sample \tab T5\cr
#' }
#'
#' @param PPB.data A data frame containing mass-spectrometry peak areas,
#' indication of chemical identiy, and measurment type.
#' @param this.conc The plasma protein concentration relative to physiologic
#' levels (default 100%)
#' @param FILENAME A string used to identify outputs of the function call.
#' (defaults to "BASE_Model_Results")
#' @param TEMP.DIR An optional directory where file writing may be faster.
#' @param JAGS.PATH The file path to JAGS.
#' @param NUM.CHAINS The number of Markov Chains to use. This allows evaluation
#' of convergence according to Gelman and Rubin diagnostic.
#' @param NUM.CORES The number of processors to use (default 2)
#' @param RANDOM.SEED The seed used by the random number generator 
#' (default 1111)
#'
#' @return \item{data.frame}{A data.frame containing quunantiles of the 
#' Bayesian posteriors} 
#'
#' @author John Wambaugh and Chantel Nicolas
#' 
#' @import parallel, runjags
#' 
#' @export calc_fup_base
calc_uc_fup <- function(PPB.data,
  FILENAME = "UC_Model_Results",
  TEMP.DIR = NULL,
  NUM.CHAINS=5, 
  NUM.CORES=2,
  RANDOM.SEED=1111,
  compound.col="Compound.Name",
  response.col="Response",
  type.col="Sample.Type",
  compound.conc.col="Nominal.Conc",
  cal.col="Cal",
  dilution.col="Dilution.Factor"
  )
{
  if (!is.null(TEMP.DIR)) 
  {
    current.dir <- getwd()
    setwd(TEMP.DIR)
  }
  
  all.blanks <- subset(PPB.data,!is.na(eval(response.col)))

  OUTPUT.FILE <- paste(FILENAME,".txt",sep="")

  set.seed(RANDOM.SEED)
  if (!file.exists(OUTPUT.FILE))
  {
    Results <- NULL
  } else {
    Results <- read.table(OUTPUT.FILE,sep=" ",stringsAsFactors=F,header=T)
  }

  if (NUM.CORES>1)
  {
    CPU.cluster <- makeCluster(min(NUM.CORES,NUM.CHAINS))
  } else CPU.cluster <-NA
  
  coda.out <- list()
  for (this.compound in  unique(PPB.data[,compound.col]))
    if (!(this.compound %in% Results[,compound.col]))
    {
      this.name <- PPB.data[PPB.data[,compound.col]==this.compound,compound.col][1]
      
      print(paste(
        this.name,
        " (",
        which(unique(PPB.data[,compound.col])==this.compound),
        " of ",
        length(unique(PPB.data[,compound.col])),
        ")",
        sep=""))
      MSdata <- PPB.data[PPB.data[,compound.col]==this.compound,]
    
      if (any(MSdata[,type.col]=="CC") &
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
        Conc <- CC.data[,compound.conc.col]
        Dilution.Factor <- CC.data[,dilution.col]
  #
  #
  #  Each series (currently) contains T5 and AF data
        T5.data <- MSdata[MSdata[,type.col]=="T5",]
        AF.data <- MSdata[MSdata[,type.col]=="AF",]
        all.series <- unique(T5.data$Series) 
        Num.series <- length(all.series)
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
          "Response.obs" = UC.obs[,response.col],
          "obs.conc" = UC.obs[,"Obs.Conc"],
          "obs.cal" = UC.obs[,"Obs.Cal"],
          "Conc" = Conc,
          "Num.cc.obs" = Num.cc.obs,
          "Num.series" = Num.series,
          "Dilution.Factor" = Dilution.Factor
        )

        initfunction <- function(chain)
        {
          seed <- as.numeric(paste(rep(chain,6),sep="",collapse=""))
          set.seed(seed)
          
          return(list(
            .RNG.seed=seed,
            .RNG.name="base::Super-Duper",
# Parameters that may vary between calibrations:
            log.const.analytic.sd =runif(Num.cal,-2.2,-1.8),
            log.hetero.analytic.slope.factor = runif(Num.cal,-0.2,0.2),
            background = rlnorm(Num.cal,log(10^-3),1),
            log.calibration = runif(Num.cal,-0.1,0.1),
            log.C.thresh = runif(Num.cal,-3,-2),
# There is only one Fup per chemical:
            log.Fup = log10(runif(1,0,1))
          ))
        }
        
        save(this.compound,mydata,UC_PPB_model,initfunction,
          file=paste(FILENAME,"-PREJAGS.RData",sep=""))  
        coda.out[[this.compound]] <- autorun.jags(
          UC_PPB_model, 
          n.chains = NUM.CHAINS,
          method="parallel", 
          cl=CPU.cluster,
          summarise=T,
          inits = initfunction,
          startburnin = 25000, 
          startsample = 25000, 
          max.time="5m",
          crash.retry=2,
          adapt=10000,
          psrf.target = 1.1,
          thin.sample=2000,
          data = mydata,
          jags = findjags(),
          monitor = c(
            'const.analytic.sd',
            'hetero.analytic.slope',
            'Fup',
            'C.thresh',
            'background',
            'calibration',
            "Conc"))

        sim.mcmc <- coda.out[[this.compound]]$mcmc[[1]]
        for (i in 2:NUM.CHAINS) sim.mcmc <- rbind(sim.mcmc,coda.out[[this.compound]]$mcmc[[i]])
        results <- apply(sim.mcmc,2,function(x) quantile(x,c(0.025,0.5,0.975)))
    
        new.results <- data.frame(this.compound,stringsAsFactors=F)
        colnames(new.results) <- compound.col
        new.results[,c(
          "Fup.Med",
          "Fup.Low",
          "Fup.High")] <- results[c(2,1,3),"Fup"]
    
        print(results)
    
        Results <- rbind(Results,new.results)
    
        write.table(Results,file=OUTPUT.FILE,sep=" ",row.names=FALSE)
      }    
    }
  
  if (!is.null(TEMP.DIR)) 
  {
    setwd(current.dir)
  }
  stopCluster(CPU.cluster)
  
  View(Results)
  save(Results,
    file=paste("UC-Fup-Analysis-",Sys.Date(),".RData",sep=""))

  return(Results)  
}


