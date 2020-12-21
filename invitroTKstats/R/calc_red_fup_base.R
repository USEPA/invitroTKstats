PPB_model_base <- "
model {

# Measurement Model:
  log.const.analytic.sd ~ dnorm(-2,0.1)
  const.analytic.sd <- 10^log.const.analytic.sd
  log.hetero.analytic.slope.factor ~ dnorm(1,.1)
  hetero.analytic.slope.factor <- 10^log.hetero.analytic.slope.factor
  hetero.analytic.slope <- hetero.analytic.slope.factor*const.analytic.sd
  background ~ dunif(0,1000)   
  log.calibration ~ dnorm(0,.1)
  calibration <- 10^log.calibration
# Concentrations below this value are not detectable:
  log.C.thresh ~ dnorm(log(C.frank/2),0.25)T(log(C.frank/1000),log(10*C.frank))
  C.thresh <- exp(log.C.thresh)
  Blank.pred <- background  
  Blank.prec <-  (const.analytic.sd+hetero.analytic.slope*(Blank.pred))^(-2)
  for (i in 1:Num.Blank.obs) {
    Blank.data[i] ~ dnorm(Blank.pred,Blank.prec)
  }
  T0.pred <- calibration*(C.frank/5-C.thresh) + background
  T0.prec <- (const.analytic.sd+hetero.analytic.slope*(T0.pred))^(-2)
  for (i in 1:Num.T0.obs) {
# Must include the dilution factor (5 for Plasma)
    T0.data[i] ~ dnorm(T0.pred,T0.prec) 
  }
# Must include the dilution factor (2 for PBS)
  PBS.pred <- calibration*(C.u[1]/2 - C.thresh)*step(C.u[1]/2 - C.thresh) + background
  PBS.prec <- (const.analytic.sd+hetero.analytic.slope*(PBS.pred))^(-2)
  for (i in 1:Num.PBS.obs) {
      PBS.data[i] ~dnorm(PBS.pred,PBS.prec)
  }
# Must include the dilution factor (5 for Plasma)
  Plasma.pred <- calibration*(C.upb[1]/5 - C.thresh)*step(C.upb[1]/5 - C.thresh) + background
  Plasma.prec <- (const.analytic.sd+hetero.analytic.slope*(Plasma.pred))^(-2)
  for (i in 1:Num.Plasma.obs) {
      Plasma.data[i] ~dnorm(Plasma.pred,Plasma.prec)
  }
  
# Binding Model:
  log.Fup ~ dunif(-10,0)
  Fup <- 10^log.Fup
#  for (i in 1:3) {
# Missing (bound to walls/membrane) chemical:
  C.missing[1] ~ dunif(0,C.frank)
# Unbound concentration in both wells:
  C.u[1] <- Fup/(Fup+1)*(C.frank-C.missing)
# Bound concentration in plasma well:
  C.b[1] <- (1-Fup)/Fup*C.u[1]
# Toal concentration in plasma well:
  C.upb[1] <- C.b[1] + C.u[1] 
  
}
"

#' Calculate fraction unbound in plasma from rapid equilibruim dialysis data
#' 
#' The data frame of observations should be annotated according to
#' of these types:
#' \tabular{rrrrr}{
#'   Sample Blank (no chemical, just plasma) \tab Blank\cr
#'   Time zero chemical and plasma \tab T0\cr
#'   Equilibrium chemical in phosphate-buffered well (no plasma) \tab PBS\cr
#'   Equilibrium chemical in plasma well \tab Plasma\cr
#' }
#'
#' @param PPB.data A data frame containing mass-spectrometry peak areas,
#' indication of chemical identiy, and measurment type.
#'
#' @param this.conc The plasma protein concentration relative to physiologic
#' levels (default 100%)
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
#'
#' @return A data.frame containing quunantiles of the Bayesian posteriors 
#'
#' @author John Wambaugh and Chantel Nicolas
#' 
#' @export calc_fup_red_base
calc_fup_red_base <- function(PPB.data,
  this.conc = 100,
  FILENAME = "BASE_Model_Results",
  TEMP.DIR = NULL,
#  JAGS.PATH ="C:/Program Files/JAGS/JAGS-4.3.0/x64/bin/",
  NUM.CHAINS=5, 
  NUM.CORES=2,
  sample.col="Lab.Sample.Name",
  lab.compound.col="Lab.Compound.Name",
  dtxsid.col="DTXSID",
  date.col="Date",
  compound.col="Compound.Name",
  area.col="Area",
  series.col="Series",
  type.col="Sample.Type",
  compound.conc.col="Nominal.Conc",
  cal.col="Cal",
  dilution.col="Dilution.Factor",
  istd.col="ISTD.Area",
  istd.name.col="ISTD.Name",
  istd.conc.col="ISTD.Conc",
  nominal.test.conc.col="Test.Target.Conc" 
  )
{
  if (!is.null(TEMP.DIR)) 
  {
    current.dir <- getwd()
    setwd(TEMP.DIR)
  }
  
  PPB.data <- as.data.frame(PPB.data)
  all.blanks <- subset(PPB.data,!is.na(eval(area.col)))

  # Only include the data types used:
  PPB.data <- subset(PPB.data,PPB.data[,type.col] %in% c(
    "Blank",
    "T0",
    "PBS",
    "Plasma"))
  
  # Organize the columns:
  PPB.data <- PPB.data[,c(
    sample.col,
    date.col,
    compound.col,
    dtxsid.col,
    lab.compound.col,
    type.col,
    dilution.col,
    cal.col,
    compound.conc.col,
    nominal.test.conc.col,
    istd.name.col,
    istd.conc.col,
    istd.col,
    series.col,
    area.col)]
    
  # Rename the columns:
    sample.col <- "Lab.Sample.Name"
    date.col <- "Date"
    compound.col <- "Compound.Name"
    dtxsid.col <- "DTXSID"
    lab.compound.col <- "Lab.Compound.Name"
    type.col <- "Sample.Type"
    dilution.col <- "Dilution.Factor"
    cal.col <- "Calibration"
    compound.conc.col <- "Standard.Conc"
    nominal.test.conc.col <- "Test.Target.Conc"
    istd.name.col <- "ISTD.Name"
    istd.conc.col <- "ISTD.Conc"
    istd.col <- "ISTD.Area"
    series.col <- "Series"
    area.col <- "Area"
  colnames(PPB.data) <- c(
    sample.col,
    date.col,
    compound.col,
    dtxsid.col,
    lab.compound.col,
    type.col,
    dilution.col,
    cal.col,
    compound.conc.col,
    nominal.test.conc.col,
    istd.name.col,
    istd.conc.col,
    istd.col,
    series.col,
    area.col)
  
  # calculate the reponse:
  PPB.data[,"Response"] <- PPB.data[,area.col] /
     PPB.data[,istd.col] *  PPB.data[,istd.conc.col]
  
# Write out a "level 1" file (data organized into a standard format):  
  write.table(, 
    file=paste(FILENAME,"-Level1.tsv",sep=""),
    sep="\t",
    row.names=F,
    quote=F)

  OUTPUT.FILE <- paste(FILENAME,"-Level2.tsv",sep="")

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
  
  for (this.compound in  unique(PPB.data$CompoundName))
    if (!(this.compound %in% Results[,"CompoundName"]))
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
    # Can't use blanks that are NA:
      MSdata <- subset(MSdata,!is.na(Response) | MSdata[,type.col] != "Blank")
    # Delete any concentrations that are NA:
      MSdata <- subset(MSdata,!is.na(Conc))
    
      if (any(MSdata$Type=="T0") &
          any(MSdata$Type=="PBS") &
          any(MSdata$Type=="Plasma") &
          !(all(MSdata$ISTDResponseRatio==0)))
      {
        for (this.type in c("Blank","T0","PBS","Plasma"))
        {
          if (dim(subset(MSdata,Type==this.type & Conc==this.conc))[1]==0)
          {
             new.row <- MSdata[1,]
             if (this.type=="Blank")
             {
               new.row$ISTDResponseRatio <- median(all.blanks$ISTDResponseRatio,na.rm=T)
               new.row$SampleName <- "MedianBlank"
             } else {
               new.row$SampleName <- "DummyObservation"
               new.row$ISTDResponseRatio <- NA
             }
             new.row$Conc <- this.conc
             new.row$Type <- this.type
          #       new.row$Recovered <- NA
             MSdata <- rbind(MSdata,new.row)
          }
        }
    
  #
  #
  #
        Blank.data <- subset(MSdata,Type=="Blank")
        Blank.data <- subset(Blank.data,!is.na(Blank.data[,"ISTDResponseRatio"]))
        Num.Blank.obs <- dim(Blank.data)[1]
    
  #
  #
  #
        T0.data <- subset(MSdata,Type=="T0")
        Num.T0.obs <- dim(T0.data)[1]
    
  #
  #
  #
        PBS.data <- subset(MSdata,Type=="PBS")
        Num.PBS.obs <- dim(PBS.data)[1]
    
  #
  #
  #
        Plasma.data <- subset(MSdata,Type=="Plasma")
        Num.Plasma.obs <- dim(Plasma.data)[1]
    
        #mg/mL -> g/L is 1:1
        #kDa -> g/mol is *1000
        #g/mol -> M is g/L/MW
        #M <- uM is /1000000
        PPB100 <- 70/(66.5*1000)*1000000 # Berg and Lane (2011) 60-80 mg/mL, albumin is 66.5 kDa, pretend all protein is albumin to get uM
        C.frank <- 5 # uM frank parent concentration
    
        mydata <- list(                
          'Num.Blank.obs' = Num.Blank.obs,
          'T0.data' = T0.data[,"ISTDResponseRatio"],
          'Blank.data' = Blank.data[,"ISTDResponseRatio"],
          'Num.T0.obs' = Num.T0.obs,
          'PBS.data' = PBS.data[,"ISTDResponseRatio"],
          'Num.PBS.obs' = Num.PBS.obs,
          'Plasma.data' = Plasma.data[,"ISTDResponseRatio"],
          'Num.Plasma.obs' = Num.Plasma.obs,
          'C.frank' = C.frank
        )

        initfunction <- function(chain)
        {
          BG <- mean(data[["Blank.data"]],na.rm=T)
          BG[BG<0] <- 0
          background <- rlnorm(1,log(BG/2+10^-6),1)
          calibration <- max(
            10^-3.5,
            (mean(mydata$T0.data)-background)*5/mydata$C.frank)    
          Fup <- max(
            min(
              (calibration *
              mean(mydata$PBS.data,na.rm=T) /
              2-background) /
              (calibration *
              mean(mydata$Plasma.data,na.rm=T) /
              5-background),
              1),
            2*10^-5,
            na.rm=T)
          C.missing <- runif(1,0,mydata[["C.frank"]])
    
          return(list(
            .RNG.seed=as.numeric(paste(rep(chain,6),sep="",collapse="")),
            .RNG.name="base::Super-Duper",
            log.const.analytic.sd =runif(1,0.10,1),
            log.hetero.analytic.slope.factor = log10(runif(1,0, 1)),
            background = background,
            log.calibration = log10(calibration),
            log.Fup = log10(Fup),
            C.missing = C.missing
          ))
        }
        
        save.image(CODA.FILE)  
        coda.out[[MSdata[,"CAS"][1]]] <- autorun.jags(
          PPB_model_base, 
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
            'log.const.analytic.sd',
            'hetero.analytic.slope.factor',
            'Fup',
            'C.thresh',
            'background',
            'calibration'))
        
        sim.mcmc <- coda.out[[MSdata[,"CAS"][1]]]$mcmc[[1]]
        for (i in 2:NUM.CHAINS) sim.mcmc <- rbind(sim.mcmc,coda.out$mcmc[[i]])
        results <- apply(sim.mcmc,2,function(x) quantile(x,c(0.025,0.5,0.975)))
    
        Fup.point <- 2/5 *
          (mean(mydata$PBS.data) - mean(mydata$Blank.data)) /
          (mean(mydata$Plasma.data) - mean(mydata$Blank.data))
        
        new.results <- data.frame(Name=MSdata[,"Name"][1],
                                  DTXSID=MSdata[,"DTXSID"][1],
                                  CAS=MSdata[,"CAS"][1],
                                  CompoundName=this.compound,
                                  Fup.point=Fup.point,
                                  stringsAsFactors=F)
        new.results[,c("Fup.Med","Fup.Low","Fup.High")] <- results[c(2,1,3),"Fup"]
    
        print(new.results)
    
        Results <- rbind(Results,new.results)
    
        write.table(Results,file=OUTPUT.FILE,sep=" ",row.names=FALSE)
      }    
    }
  
  if (!is.null(TEMP.DIR)) 
  {
    setwd(current.dir)
  }
  stopCluster(CPU.cluster)
  
  BASE_Model_Results <- Results
  
  View(BASE_Model_Results)
  save(BASE_Model_Results,
    file=paste("PPB-base-analysis-",Sys.Date(),".RData",sep=""))

  return(BASE_Model_Results)  
}


