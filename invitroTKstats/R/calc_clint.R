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
      all.blanks <- subset(Clint.data,is.na(this.data) & !is.na(area.col))
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
    Num.obs1 <- dim(subset(this1data, !is.na(Time)))[1]
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
  initfunction <- function(mydata, fit.data, chain)
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
   
  OUTPUT.FILE <- paste(FILENAME,"-Level4.tsv",sep="")

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
  for (this.chem in unique(clint.data[,compound.col]))
    if (!(this.compound %in% Results[,compound.col]))
    {
      this.subset <- subset(clint.data,clint.data[,compound.col]==this.chem)
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
      
    # Check to make sure there are data for PBS and plasma: 
      if (dim(this.cvt)[1] > 1 & dim(this.blank)[1] > 1)
      {
  
      }
    
    
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
        Conc <- CC.data[,compound.conc.col]
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
            nominal.test.conc.col],na.rm=T)
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
          init.Conc <- rep(NA,Num.cc.obs+Num.series*2)
          init.Conc[(Num.cc.obs+1):(Num.cc.obs+Num.series)] <- 
            mydata$Test.Nominal.Conc[
              mydata$obs.cal[(Num.cc.obs+1):(Num.cc.obs+Num.series)]]
            
          return(list(
            .RNG.seed=seed,
            .RNG.name="base::Super-Duper",
# Parameters that may vary between calibrations:
            log.const.analytic.sd =runif(Num.cal,-5,-0.5),
            log.hetero.analytic.slope = runif(Num.cal,-5,-0.5),
# Average across all the calibrations (the sampler will vary these):
            C.thresh = rep(min(max(0,intercept/slope),Test.Nominal.Conc/10),Num.cal),
            log.calibration = rep(max(min(-2.95,log10(slope)),1.95),Num.cal),
# There is only one Fup per chemical:
            log.Fup = log10(runif(1,0,1)),
# Set the initial concentrations:
            Conc = init.Conc
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
          startsample = 50000, 
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
            'calibration'))
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
    
        Results <- rbind(Results,new.results)
    
        write.table(Results, 
          file=paste(OUTPUT.FILE,sep=""),
          sep="\t",
          row.names=F,
          quote=F)
      }    
    }  
  
  
  
  
  
  
  
  
  
  for (this.chem in unique(clint.data[,compound.col]))
    if (!(this.compound %in% Results[,"CompoundName"]))
    {
      this.subset <- subset(clint.data,clint.data[,compound.col]==this.chem)
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
      
      print(paste(
        this.chem,
        " (",
        which( unique(clint.data[,compound.col]==this.chem)),
        " of ",
        length( unique(clint.data[,compound.col])),
        ")",
        sep=""))

  # Check to make sure there are data for PBS and plasma: 
    if (dim(this.cvt)[1] > 1 & dim(this.blank)[1] > 1)
    {
      this.data <- rbind(this.blank,this.cvt)
      this.data[this.data$Sample.Type=="Blank","Conc"] <- 0
      this.data[this.data$Sample.Type=="Blank","Time"] <- 0
      this.data[this.data$Sample.Type=="Cvst","Response"] <-
        this.data[this.data$Sample.Type=="Cvst","Response"]*df.cvt
      min.response <- sort(unique(this.data$Response))
      min.response <- min.response[min.response!=0]
      min.response <- min.response[1]
      this.data[this.data$Response==0,"Response"] <- min.response/2




      num.chem <- num.chem + 1
      
      this.fit <- try(mle(lldecay,
        start=list(cal=1,k_elim=0.1,sigma=1),
        lower=list(cal=0,k_elim=0,sigma = 0.0001)))
      this.null <- try(mle(lldecay,
        start=list(cal=1,sigma=1),
        lower=list(cal=0,sigma = 0.0001),
        fixed=list(k_elim=0)))      
      
      if (class(this.fit)!="try-error" & class(this.null)!="try-error")
      {
        this.row$Clint <- 1000*coef(this.fit)["k_elim"]/hep.density
        this.row$Clint.pValue <- min(exp(-(AIC(this.null)-AIC(this.fit))),1)
        this.row$Fit <- paste(paste(unique(this.data$Conc),collapse=", "),"uM")
        this.row$AIC <- AIC(this.fit)
        this.row$AIC.Null <- AIC(this.null)
        this.row$Clint.1 <- NA
        this.row$Clint.10 <- NA
        this.row$AIC.Sat <- NA
        this.row$Sat.pValue <- NA
        if (all(c(1,10)%in%unique(this.data$Conc)))
        {
          this.sat.fit <- try(mle(llsatdecay,
            start=list(cal=1,k_elim=0.1,sigma=1,sat=1),
            lower=list(cal=0,k_elim=0,sigma = 0.0001,sat=0),
            upper=list(sat=1)))
          if (class(this.sat.fit)!="try-error")
          {
            this.row$Clint.1 <- 1000*coef(this.sat.fit)["k_elim"]/hep.density
            this.row$Clint.10 <- 1000*coef(this.sat.fit)["k_elim"]*
              coef(this.sat.fit)["sat"]/hep.density
            this.row$AIC.Sat <- AIC(this.sat.fit)
            if (this.row$Clint.pValue==1) test.AIC <- this.row$AIC.Null
            else test.AIC <- this.row$AIC
            this.row$Sat.pValue <- min(exp(-(test.AIC-AIC(this.sat.fit))),1)
          }
        }
        out.table <- rbind(out.table, this.row)
        print(paste(
          this.row$Compound.Name,
          "Cl_int =",
          signif(this.row$Clint,3),
          "uL/min/million hepatocytes, p-Value =",
          signif(this.row$Clint.pValue,3),
          "."
          ))
      }
    }  
  }

  rownames(out.table) <- make.names(out.table$Compound.Name, unique=TRUE)
  out.table <- apply(out.table,2,unlist) 
  out.table[,"Clint"] <- signif(as.numeric(out.table[,"Clint"]),3) 
  out.table[,"Clint.1"] <- signif(as.numeric(out.table[,"Clint.1"]),3) 
  out.table[,"Clint.10"] <- signif(as.numeric(out.table[,"Clint.10"]),3) 
  out.table[,"Clint.pValue"] <- signif(as.numeric(out.table[,"Clint.pValue"]),3) 
  out.table[,"AIC"] <- signif(as.numeric(out.table[,"AIC"]),3)
  out.table[,"AIC.Null"] <- signif(as.numeric(out.table[,"AIC.Null"]),3)
  out.table[,"AIC.Sat"] <- signif(as.numeric(out.table[,"AIC.Sat"]),3)
  out.table[,"Sat.pValue"] <- signif(as.numeric(out.table[,"Sat.pValue"]),3) 
   
  
  out.table <- as.data.frame(out.table)
  out.table$Clint <- as.numeric(out.table$Clint)
  out.table$Clint.pValue <- as.numeric(out.table$Clint.pValue)
    
# Write out a "level 3" file (data organized into a standard format):  
  write.table(out.table, 
    file=paste(FILENAME,"-Clint-Level3.tsv",sep=""),
    sep="\t",
    row.names=F,
    quote=F)
 
  print(paste("Intrinsic clearance (Clint) calculated for",num.chem,"chemicals."))
  print(paste("Intrinsic clearance (Clint) calculated for",num.cal,"measurements."))

  return(out.table)  
  
  if (!is.null(TEMP.DIR)) 
  {
    setwd(current.dir)
  }
  stopCluster(CPU.cluster)

  View(out.table)

  return(out.table)  
}




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


