#0.5 * 10^-6 viable cells/mL
library(runjags)
library(parallel)
library(gdata)
library(coda)
library(ggplot2)

# Get rid of anything in the workspace:
rm(list=ls()) 

#JAGS.PATH <- "C:/Program Files/JAGS/JAGS-4.2.0/x64/bin/"
JAGS.PATH <- "C:/Program Files/JAGS/JAGS-4.3.0/x64/bin/"
RJAGS.MODEL.FILE <- "Z:/Research Projects/CeetoxHumanHTTK/Clint/clearance-model-052318.jags"
#FILENAME <- paste("CLint-",Sys.Date(),sep="")
FILENAME <- "CLint-2018-10-17"
Output.File <- paste(FILENAME,"-Results.txt",sep="")

RANDOM.NUMBER.SEED <- 1111

#setwd("Z:/Research Projects/CeetoxHumanHTTK/Clint")

NUM.CHAINS <- 5
NUM.CORES <- 12


if (!file.exists(Output.File)) 
{
  clearance.table <- NULL
} else {                                        
  clearance.table <- read.table(Output.File,sep=" ",stringsAsFactors=F,header=T)
}
                  
setwd("c:/users/jwambaug/Rwork")
CODA.FILE <- paste(FILENAME,"-CODA.RData",sep="")
if (!file.exists(CODA.FILE)) 
{
  Clintrawdata <- read.xls("L:/Lab/NCCT_ExpoCast/ExpoCast2018/HTTKNewData/Summary/raw-clint-data.xlsx",stringsAsFactors=F)
  length(unique(Clintrawdata$Name))
  
  source("L:/Lab/NCCT_ExpoCast/ExpoCast2018/HTTKNewData/Summary/build-httk-master-list-061318.R")
  

  
  Clint.data <- subset(Clintrawdata,Name!="")
  Clint.data$SampleID <- Clint.data$Name
  Clint.data$Area <- as.numeric(Clint.data$Area)
  Clint.data$ISTD.Area <- as.numeric(Clint.data$ISTD.Area)
  Clint.data$Time <- as.numeric(Clint.data$Time)
  Clint.data$ISTDResponseRatio <- as.numeric(Clint.data$ISTDResponseRatio)                                     
  
  for(this.compound in unique(Clint.data$SampleID))
  {
    if (this.compound %in% master.table$EPA_SAMPLE_ID)
    {
      Clint.data[Clint.data$SampleID==this.compound,"CAS"] <- master.table[master.table$EPA_SAMPLE_ID==this.compound,"CASRN"]
      Clint.data[Clint.data$SampleID==this.compound,"Name"] <- master.table[master.table$EPA_SAMPLE_ID==this.compound,"Preferred_Name"]
      Clint.data[Clint.data$SampleID==this.compound,"DTXSID"] <- master.table[master.table$EPA_SAMPLE_ID==this.compound,"DTXSID"]
#    } else if (this.compound == "Verapamil")
#    {
#      Clint.data[Clint.data$SampleID==this.compound,"CAS"] <- "52-53-9" 
#      Clint.data[Clint.data$SampleID==this.compound,"Name"] <- "Verapamil"
#      Clint.data[Clint.data$SampleID==this.compound,"DTXSID"] <- "DTXSID9041152"
#    } else if (this.compound == "7-OH-Coumarin")
#    {
#      Clint.data[Clint.data$SampleID==this.compound,"CAS"] <- "93-35-6" 
#      Clint.data[Clint.data$SampleID==this.compound,"Name"] <- "7-OH-Coumarin"
#      Clint.data[Clint.data$SampleID==this.compound,"DTXSID"] <- "DTXSID5052626"
#    } else if (this.compound == "Midazolam")
#    {
#      Clint.data[Clint.data$SampleID==this.compound,"CAS"] <- "59467-70-8" 
#      Clint.data[Clint.data$SampleID==this.compound,"Name"] <- "Midazolam"
#      Clint.data[Clint.data$SampleID==this.compound,"DTXSID"] <- "DTXSID5023320"
#    } else if (this.compound == "Warfarin")
#    {
#      Clint.data[Clint.data$SampleID==this.compound,"CAS"] <- "81-81-2" 
#      Clint.data[Clint.data$SampleID==this.compound,"Name"] <- "Warfarin"
#      Clint.data[Clint.data$SampleID==this.compound,"DTXSID"] <- "DTXSID5023742"
#    } else if (this.compound == "7-Hydroxy4trifluoromethyl coumarin")
#    {
#      Clint.data[Clint.data$SampleID==this.compound,"CAS"] <- "81-81-2" 
#      Clint.data[Clint.data$SampleID==this.compound,"Name"] <- "Warfarin"
#      Clint.data[Clint.data$SampleID==this.compound,"DTXSID"] <- "DTXSID5023742"
    } else {
       Clint.data <- subset(Clint.data,SampleID!=this.compound)
       print(paste("Data for",this.compound,"omitted."))
    }
  }
  colnames(Clint.data)[colnames(Clint.data)=="Sample.Name"] <- "SampleName"
#  colnames(Clint.data)[colnames(Clint.data)=="Name"] <- "CompoundName"



  coda.out <- list()
} else {
  load(CODA.FILE)
}

build_mydata <- function(this.sample)
{
  this.1data <- subset(Clint.data,Name==this.sample & Conc==1)
  this.1data$obs <- this.1data$ISTDResponseRatio
  if (all(is.na(this.1data$obs))) this.1data$obs <- this.1data$Area/this.1data$ISTD.Area
  this.10data <- subset(Clint.data,Name==this.sample & Conc==10)
  this.10data$obs <- this.10data$ISTDResponseRatio
  if (all(is.na(this.10data$obs))) this.10data$obs <- this.10data$Area/this.10data$ISTD.Area
  
  
  blank.obs <- c(subset(this.1data,is.na(Time..mins.))$obs,subset(this.10data,is.na(Time..mins.))$obs)
  Num.blanks <- length(blank.obs)
  if (Num.blanks==0 | all(is.na(blank.obs)))
  {
    all.blanks <- subset(Clint.data,is.na(Time)&!is.na(Area))
    blank.obs <- median(all.blanks$ISTDResponseRatio,na.rm=T)
    Num.blanks <- 1
  }
  blank.conc <- rep(2,length(blank.obs))
  if (length(subset(this.1data,is.na(Time..mins.))$obs)>0) blank.conc[1:length(subset(this.1data,is.na(Time..mins.))$obs)] <- 1
  
  obs <- c(subset(this.1data,!is.na(Time..mins.))$obs,subset(this.10data,!is.na(Time..mins.))$obs)
  Num.obs <- length(obs)
  obs.time <- c(subset(this.1data,!is.na(Time..mins.))$Time..mins.,subset(this.10data,!is.na(Time..mins.))$Time..mins.)
  
  if (Num.obs>0 & Num.blanks>0 & any(!is.na(obs)) & length(unique(obs.time))>1)
  {
    obs.conc <- rep(2,length(obs))
    if (length(subset(this.1data,!is.na(Time..mins.))$obs)>0) obs.conc[1:length(subset(this.1data,!is.na(Time..mins.))$obs)] <- 1
    
    return(mydata <- list('obs' = obs,
                   'obs.conc' = obs.conc,
                   'obs.time' = obs.time,
                   'Num.obs' = Num.obs,
                   'Blank.obs' = blank.obs,
                   'Blank.conc' = blank.conc,
                   'Num.blank.obs' = Num.blanks
                   ))
    } else return(NULL)
}

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

if (NUM.CORES>1)
{
  cl2 <- makeCluster(min(NUM.CORES,NUM.CHAINS))
} else cl2 <-NA
#pdf(file=paste(FILENAME,".pdf",sep=""))
for (this.TXCode in sort(unique(Clint.data$SampleID)))
if (!(this.TXCode %in% c("BF00175289","EV0001852"))) # Background higher than observations
{
#  setwd("Z:/Research Projects/CeetoxHumanHTTK/Clint")
  this.sample <- Clint.data[Clint.data$SampleID==this.TXCode,"Name"][1]      
  this.cas <- Clint.data[Clint.data$Name==this.sample,"CAS"][1]
  this.name <- this.sample
#  this.TXCode <- Clint.data[Clint.data$Name==this.sample,"SampleID"][1]
  this.DTXSID <- Clint.data[Clint.data$Name==this.sample,"DTXSID"][1]
  
  print(paste(this.name," (",which(sort(unique(Clint.data$SampleID))==this.TXCode)," of ",length(unique(Clint.data$SampleID)),")",sep=""))
  
  mydata <- build_mydata(this.sample)
  if (!is.null(mydata))
  {
    fit.data <- make.fit.data(mydata)
    set.seed(RANDOM.NUMBER.SEED) 
  
    initfunction <- function(chain)
    {
      background <- rep(0,2)
      calibration <- rep(1,2)
      log.C.thresh <- rep(NA,2)
      for (this.conc in 1:2)
      {
        Blank.data <- data[["Blank.obs"]][data[["Blank.conc"]]==this.conc]
        T0.data <- data[["obs"]][data[["obs.conc"]]==this.conc & data[["obs.time"]]==0]
        background[this.conc] <- runif(1,min(Blank.data),max(Blank.data))
        background[is.na(background)] <- 0
        calibration[this.conc] <- runif(1,min(T0.data),max(T0.data))/c(1,10)[this.conc]-background[this.conc]        
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
  
    #if (this.name %in% c("Chlorpromazine hydrochloride","Amiodarone hydrochloride")) browser()
  
    if (!(this.sample %in% clearance.table[,"Name"]))
    {                      
      save.image(CODA.FILE)  
      coda.out[[this.cas]]  <- autorun.jags(RJAGS.MODEL.FILE, 
                         n.chains = NUM.CHAINS,
                         method="parallel", 
                         cl=cl2,
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
                         monitor = c('log.const.analytic.sd','hetero.analytic.slope.factor','C.thresh','calibration','background'))
      
      coda.out[[this.cas]] <-extend.jags(coda.out[[this.cas]],
                            drop.monitor = c('log.const.analytic.sd','hetero.analytic.slope.factor'), 
                            add.monitor = c('slope','decreases','saturates'))
    } 
    
    sim.mcmc <- coda.out[[this.cas]]$mcmc[[1]]
    for (i in 2:NUM.CHAINS) sim.mcmc <- rbind(sim.mcmc,coda.out[[this.cas]]$mcmc[[i]])
    results <- apply(sim.mcmc,2,function(x) quantile(x,c(0.025,0.5,0.975)))
  
    if (!(this.sample %in% clearance.table[,"Name"]))
    {                      
      new.row <- as.data.frame(this.TXCode,stringsAsFactors=F)
      new.row <- cbind(new.row,as.data.frame(this.DTXSID,stringsAsFactors=F))
      new.row <- cbind(new.row,as.data.frame(this.cas,stringsAsFactors=F))
      new.row <- cbind(new.row,as.data.frame(this.name,stringsAsFactors=F))
      
      new.row <- cbind(new.row,as.data.frame(mean(sim.mcmc[,"decreases"])))
      
      new.row <- cbind(new.row,as.data.frame(mean(sim.mcmc[,"saturates"])))
      new.row <- cbind(new.row,as.data.frame(t(results[2,"slope[1]"])))
      new.row <- cbind(new.row,as.data.frame(t(results[2,"slope[2]"])))
      new.row <- cbind(new.row,as.data.frame(t(results[c(2,1,3),"slope[1]"]))/0.5*1000)
      new.row <- cbind(new.row,as.data.frame(t(results[c(2,1,3),"slope[2]"]))/0.5*1000)

      this.1data <- subset(Clint.data,Name==this.sample & Conc==1)
      this.1data$obs <- this.1data$ISTDResponseRatio
      if (all(is.na(this.1data$obs))) this.1data$obs <- this.1data$Area/this.1data$ISTD.Area
      this.10data <- subset(Clint.data,Name==this.sample & Conc==10)
      this.10data$obs <- this.10data$ISTDResponseRatio
      if (all(is.na(this.10data$obs))) this.10data$obs <- this.10data$Area/this.10data$ISTD.Area
      sub1 <- subset(this.1data,!is.na(Time)&!is.na(ISTDResponseRatio)&is.finite(log(ISTDResponseRatio)))
      if (dim(sub1)[1]>0&length(unique(sub1$Time))>1&any(sub1$ISTDResponseRatio>0))
      {
        fit1 <- lm(log(ISTDResponseRatio)~Time,sub1)
        fit1stats <- summary(fit1)$fstatistic
        fit1p <- pf(fit1stats[1],fit1stats[2],fit1stats[3],lower.tail=F)
        if (!is.na(fit1p))
        {
          if (fit1p<0.05) 
          {
            new.row <- cbind(new.row,-fit1$coefficients["Time"]/0.5*1000)
          } else new.row <- cbind(new.row,0)
        } else new.row <- cbind(new.row,NA)
      } else new.row <- cbind(new.row,NA)
      sub10 <- subset(this.10data,!is.na(Time)&!is.na(ISTDResponseRatio)&is.finite(log(ISTDResponseRatio)))
      if (dim(sub10)[1]>0&length(unique(sub10$Time))>1&any(sub10$ISTDResponseRatio>0))
      {
        fit10 <- lm(log(ISTDResponseRatio)~Time,sub10)
        fit10stats <- summary(fit10)$fstatistic
        fit10p <- pf(fit10stats[1],fit10stats[2],fit10stats[3],lower.tail=F)
        if (!is.na(fit10p))
        {
          if (fit10p<0.05) 
          { 
            new.row <- cbind(new.row,-fit10$coefficients["Time"]/0.5*1000)
          } else new.row <- cbind(new.row,0)
        } else new.row <- cbind(new.row,NA)
      } else new.row <- cbind(new.row,NA)
      colnames(new.row) <- c("TXCode","DTXSID","CAS","Name","Decreases.Prob","Saturates.Prob","Slope.1uM.Median","Slope.10uM.Median","CLint.1uM.Median","CLint.1uM.Low95th","CLint.1uM.High95th","CLint.10uM.Median","CLint.10uM.Low95th","CLint.10uM.High95th","CLint.1uM.Point","CLint.10uM.Point")
#browser()      
      print(new.row)
      clearance.table <- rbind(clearance.table,new.row) 
      write.table(clearance.table,file=Output.File,sep=" ",row.names=FALSE)
    }
  }
}

for (this.cas in unique(clearance.table$CAS))
  if (this.cas %in% master.table$CASRN)
  clearance.table[clearance.table$CAS==this.cas,"DTXSID"] <- unique(master.table[master.table$CASRN==this.cas,"DTXSID"])
  
setwd("Z:/Research Projects/CeetoxHumanHTTK/Clint")
                  
save.image(CODA.FILE)  
if (!is.na(cl2)) stopCluster(cl2)

pdf(file=paste(FILENAME,".pdf",sep=""))
for (this.cas in unique(clearance.table$CAS))
{
  if (dim(subset(Clint.data[Clint.data$CAS==this.cas,],!is.na(ISTDResponseRatio)))[1]>0)
  {
    this.sample <- Clint.data[Clint.data$CAS==this.cas,"Name"][1]
    this.name <- this.sample
    this.TXCode <- Clint.data[Clint.data$Name==this.sample,"CompoundName"][1]
    
    mydata<-build_mydata(this.sample)
    if (!is.null(mydata))
    {
      fit.data <- make.fit.data(mydata)
      
      # Get the calibration estimates from the MCMC:
      sim.mcmc <- coda.out[[this.cas]]$mcmc[[1]]
      for (i in 2:NUM.CHAINS) sim.mcmc <- rbind(sim.mcmc,coda.out[[this.cas]]$mcmc[[i]])
      results <- apply(sim.mcmc,2,function(x) quantile(x,c(0.025,0.5,0.975)))
      int1 <- results[2,"calibration[1]"]*(1-results[2,"C.thresh[1]"])+results[2,"background[1]"]
      int2 <- results[2,"calibration[2]"]*(10-results[2,"C.thresh[2]"])+results[2,"background[2]"]
      slope1 <- -results[2,"slope[1]"]
      slope2 <- -results[2,"slope[2]"]
      
      # Make data frame for ggplot:
      this.data <- as.data.frame(cbind(mydata$obs,mydata$obs.time,mydata$obs.conc))
      
      # Calibrate the data:
      this.data[this.data[,3]==1,1] <- this.data[this.data[,3]==1,1]/int1
      this.data[this.data[,3]==2,1] <- this.data[this.data[,3]==2,1]/int2*10
      this.data[this.data[,3]==2,3] <- 10
      this.data[,3] <- as.factor(this.data[,3])
      colnames(this.data) <- c("Concentration","Time","Test.Conc")
      
      # Ploting parameters:
      maxconc <- max(this.data$Concentration)
      minconc <- min(this.data$Concentration[this.data$Concentration>0])
      maxtime <- max(mydata$obs.time)
      this.data.plot <- ggplot(this.data, aes(Time,Concentration)) +
        geom_point(aes(colour=Test.Conc,shape=Test.Conc),size=3) +
        geom_segment(aes(x = 0, y = 1,xend=maxtime,yend=1*exp(slope1*maxtime))) +
        geom_segment(aes(x = 0, y = 10,xend=maxtime,yend=10*exp(slope2*maxtime))) + 
        scale_y_log10() +
        ylab("Concentration (uM)")+
        xlab("Time (h)")+
        ggtitle(paste(this.name," (",this.cas,")",sep="")) +
        annotate("text",label=paste(paste("Probability Decreases: ",signif(mean(sim.mcmc[,"decreases"],3))*100,"%",sep=""),
          paste("Probability Saturates: ",signif(mean(sim.mcmc[,"saturates"],3))*100,"%",sep=""),
          paste("CLint (1uM): ",signif(results[2,"slope[1]"]/60/0.5*1000,3), "(",signif(results[1,"slope[1]"]/60/0.5*1000,3),"-",signif(results[3,"slope[1]"]/60/0.5*1000,3),")",sep=""),
          paste("CLint (10uM): ",signif(results[2,"slope[2]"]/60/0.5*1000,3), "(",signif(results[1,"slope[2]"]/60/0.5*1000,3),"-",signif(results[3,"slope[2]"]/60/0.5*1000,3),")",sep=""),sep="\n"),
          0, min(minconc,1*exp(slope1*maxtime)), hjust=0, vjust=0)
   #     theme(plot.margin=unit(c(5,1,0.5,0.5),"lines"))
   this.data.plot <- this.data.plot+ annotate("text",x=-Inf,y=-Inf,hjust=0,vjust=0,label="Text annotation")
        print(this.data.plot)
    }
  }
}
dev.off()


clearance.table$Fit <- "Decreasing"
clearance.table[clearance.table$CLint.1uM.Median == 0,"Fit"] <- "Bayesian Zero"
#clearance.table[is.na(clearance.table$CLint.1uM.Point == 0),"Fit"] <- "Point Est. Missing"
for (i in 1:dim(clearance.table)[1])
  if (!is.na(clearance.table[i,"CLint.1uM.Point"]))
  {
    if (clearance.table[i,"CLint.1uM.Point"] == 0) clearance.table[i,"Fit"] <- "Point Est. Zero"
  }

zeroval <- 10^-1
performance.fig <- ggplot(clearance.table, aes(CLint.1uM.Point+zeroval ,CLint.1uM.Median+zeroval,colour=Fit)) +
  geom_point(size=3) +
  geom_segment(aes(x=CLint.1uM.Point+zeroval,xend=CLint.1uM.Point+zeroval,y=CLint.1uM.Low95th+zeroval,yend=CLint.1uM.High95th+zeroval),size=1)+
  scale_y_log10() + 
   scale_x_log10() +
  xlab(expression(paste(CL[int]," (",mu,"L/min/",10^{6}," Hep.) Point Estimate",sep="")))+
  ylab(expression(paste(CL[int]," (",mu,"L/min/",10^{6}," Hep.) Bayesian",sep="")))+
  geom_abline(intercept = 0, slope = 1,linetype="dashed")+
      theme(legend.position="bottom", text  = element_text(size=18))+ 
    scale_colour_discrete(name="Assay Result")

print(performance.fig)

pdf(file=paste(FILENAME,"-summary.pdf",sep=""))
print(performance.fig)
dev.off()



write.csv(clearance.table,row.names=F,file=paste("Processed-Data-",Sys.Date(),".txt",sep=""))
save(clearance.table,file=paste("ClintData-",Sys.Date(),".RData",sep=""))
save.image(file=paste("Finished-Clearance-Analysi-",Sys.Date(),".RData",sep=""))
