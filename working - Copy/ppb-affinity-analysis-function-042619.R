# Script created by Chantel Nicolas and revised by John Wambaugh
library(runjags)
library(parallel)
library(gdata)
library(coda)
library(ggplot2)

# Get rid of anything in the workspace:
rm(list=ls()) 

setwd("Z:/Research Projects/CeetoxHumanHTTK/PPB-JAGS-Analyses/")

# The model file that is used by JAGS:
JAGS.PATH <- "C:/Program Files/JAGS/JAGS-4.3.0/x64/bin/"
RJAGS.MODEL.FILE <- "Z:/Research Projects/CeetoxHumanHTTK/PPB-JAGS-Analyses/PPB-model-affinity-101318.jags"
#FILENAME <- paste("AFFINITY_Model_Results-",Sys.Date(),sep="")
FILENAME <- "AFFINITY_Model_Results-2019-04-26"
OUTPUT.FILE <- paste(FILENAME,".txt",sep="")
# The number of Markov Chains used (should be >2 by hopefully <= number of CPU cores):
NUM.CHAINS <- 5
NUM.CORES <- 12
RANDOM.SEED <- 1111

setwd("c:/users/jwambaug/temp")
CODA.FILE <- paste(FILENAME,"-CODA.RData",sep="")
if (!file.exists(CODA.FILE)) 
{
  
  #I. Get MS Data
  
  PPBrawdata <- read.xls("L:/Lab/NCCT_ExpoCast/ExpoCast2018/HTTKNewData/Summary/raw-ppb-data.xlsx",stringsAsFactors=F)
  length(unique(PPBrawdata$CompoundName))
  
  PPB.data <- PPBrawdata
  colnames(PPB.data)[colnames(PPB.data)=="Protein"] <- "Conc"
  PPB.data$Type<-"Blank"
  PPB.data[(regexpr("plasma",tolower(PPB.data$SampleName))!=-1),"Type"] <- "Plasma"
  PPB.data[(regexpr("_pl_",tolower(PPB.data$SampleName))!=-1),"Type"] <- "Plasma"
  PPB.data[(regexpr("t0",tolower(PPB.data$SampleName))!=-1),"Type"] <- "T0"
  PPB.data[(regexpr("blank",tolower(PPB.data$SampleName))!=-1),"Type"] <- "Blank"
  PPB.data[(regexpr("pbs",tolower(PPB.data$SampleName))!=-1),"Type"] <- "PBS"
  PPB.data[(regexpr("buffer",tolower(PPB.data$SampleName))!=-1),"Type"] <- "PBS"
  PPB.data[regexpr("-100p",tolower(PPB.data$CompoundName))!=-1,"CompoundName"] <- unlist(lapply(strsplit(PPB.data[regexpr("-100p",tolower(PPB.data$CompoundName))!=-1,"CompoundName"],"-"),function(x) x[[1]]))
  PPB.data[regexpr("-10p",tolower(PPB.data$CompoundName))!=-1,"CompoundName"] <- unlist(lapply(strsplit(PPB.data[regexpr("-10p",tolower(PPB.data$CompoundName))!=-1,"CompoundName"],"-"),function(x) x[[1]]))
  PPB.data[regexpr("-30p",tolower(PPB.data$CompoundName))!=-1,"CompoundName"] <- unlist(lapply(strsplit(PPB.data[regexpr("-30p",tolower(PPB.data$CompoundName))!=-1,"CompoundName"],"-"),function(x) x[[1]]))

  source("L:/Lab/NCCT_ExpoCast/ExpoCast2018/HTTKNewData/Summary/build-httk-master-list-061318.R")
  
  PPB.data <- subset(PPB.data,CompoundName!="")

  for(this.compound in unique(PPB.data$CompoundName))
  {
    if (this.compound %in% master.table$EPA_SAMPLE_ID)
    {
      PPB.data[PPB.data$CompoundName==this.compound,"CAS"] <- master.table[master.table$EPA_SAMPLE_ID==this.compound,"CASRN"]
      PPB.data[PPB.data$CompoundName==this.compound,"Name"] <- master.table[master.table$EPA_SAMPLE_ID==this.compound,"Preferred_Name"]
      PPB.data[PPB.data$CompoundName==this.compound,"DTXSID"] <- master.table[master.table$EPA_SAMPLE_ID==this.compound,"DTXSID"]
    } else if (this.compound == "Verapamil")
    {
      PPB.data[PPB.data$CompoundName==this.compound,"CAS"] <- "52-53-9" 
      PPB.data[PPB.data$CompoundName==this.compound,"Name"] <- "Verapamil"
      PPB.data[PPB.data$CompoundName==this.compound,"DTXSID"] <- "DTXSID9041152"
    } else if (this.compound == "Warfarin")
    {
      PPB.data[PPB.data$CompoundName==this.compound,"CAS"] <- "81-81-2" 
      PPB.data[PPB.data$CompoundName==this.compound,"Name"] <- "Warfarin"
      PPB.data[PPB.data$CompoundName==this.compound,"DTXSID"] <- "DTXSID5023742"
    } else if (this.compound == "Propranolol")
    {
      PPB.data[PPB.data$CompoundName==this.compound,"CAS"] <- "525-66-6" 
      PPB.data[PPB.data$CompoundName==this.compound,"Name"] <- "Propranolol"
      PPB.data[PPB.data$CompoundName==this.compound,"DTXSID"] <- "DTXSID6023525"
    }
  }
  PPB.data$ISTDResponseRatio <- as.numeric(PPB.data$ISTDResponseRatio)
  #Set NA's to 0 (consistent with using a step function for signal):
  PPB.data[is.na(PPB.data$ISTDResponseRatio),"ISTDResponseRatio"] <-0
  # Compound not detected with MS:
  PPB.data <- subset(PPB.data,CompoundName!="EV0001328")
  PPB.data[PPB.data$CompoundName=="propranolol","CompoundName"] <- "Propranolol"

# Some protein concentrations were mislabeled:
  mislabeled.chems <- c("BF00175258","BF00175270","EV0000613","EV0000634","EV0000635","EV0000679")
  for (this.compound in mislabeled.chems)
  {
    PPB.data[PPB.data$CompoundName==this.compound&PPB.data$Conc==10,"Conc"] <- 0
    PPB.data[PPB.data$CompoundName==this.compound&PPB.data$Conc==100,"Conc"] <- 10
    PPB.data[PPB.data$CompoundName==this.compound&PPB.data$Conc==0,"Conc"] <- 100
  }  





  # II. WRITE OUT DATA
  write.csv(PPB.data,row.names=F,file=paste("PPB-all-rawdata-",Sys.Date(),".txt",sep=""))
  coda.out <- list()
} else {
  load(CODA.FILE)
}

# III. RUN JAGS MODEL

all.blanks <- subset(PPB.data,!is.na(Area))

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

for (this.compound in  unique(PPB.data$CompoundName))
if (!(this.compound %in% c(Results[,"CompoundName"],"BF00174490")))
{
  this.name <- PPB.data[PPB.data$CompoundName==this.compound,"Name"][1]
  
  print(paste(this.name," (",which(unique(PPB.data$CompoundName)==this.compound)," of ",length(unique(PPB.data$CompoundName)),")",sep=""))
  MSdata <- subset(PPB.data,CompoundName==this.compound)
# Can't use blanks that are NA:
  MSdata <- subset(MSdata,!is.na(ISTDResponseRatio) | Type!="Blank")
# Delete any concentrations that are NA:
  MSdata <- subset(MSdata,!is.na(Conc) | Type!="Blank")

  for (this.conc in c(10,30,100))
  {
    if (all(subset(MSdata,Conc==this.conc)$ISTDResponseRatio==0|is.na(subset(MSdata,Conc==this.conc)$ISTDResponseRatio)))
    {
      MSdata <- subset(MSdata,Conc!=this.conc)
    }
  }

  if (any(MSdata$Type=="T0") &
      any(MSdata$Type=="PBS") &
      any(MSdata$Type=="Plasma"))
#       &
#      length(unique(MSdata$Conc))==3)
  {
    for (this.conc in c(10,30,100))
#    if (dim(subset(MSdata,Conc==this.conc))[1]>0)
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
    }



    Blank.data <- subset(MSdata,Type=="Blank")
    Num.Blank.obs <- c(sum(Blank.data$Conc==10),
                       sum(Blank.data$Conc==30),
                       sum(Blank.data$Conc==100))
        
    
    T0.data <- subset(MSdata,Type=="T0")
    Num.T0.obs <- c(sum(T0.data$Conc==10),
                    sum(T0.data$Conc==30),
                    sum(T0.data$Conc==100))
    
    PBS.data <- subset(MSdata,Type=="PBS")
    Num.PBS.obs <- c(sum(PBS.data$Conc==10),
                       sum(PBS.data$Conc==30),
                       sum(PBS.data$Conc==100))


    Plasma.data <- subset(MSdata,Type=="Plasma")
    Num.Plasma.obs <- c(sum(Plasma.data$Conc==10),
                       sum(Plasma.data$Conc==30),
                       sum(Plasma.data$Conc==100))

    for (this.conc in c(10,30,100))
    {
      if (subset(MSdata,Conc==this.conc&Type=="Blank")$SampleName[1] == "MedianBlank") Num.Blank.obs[this.conc==c(10,30,100)] <- 1
      if (subset(MSdata,Conc==this.conc&Type=="T0")$SampleName[1] == "DummyObservation") Num.T0.obs[this.conc==c(10,30,100)] <- 0
      if (subset(MSdata,Conc==this.conc&Type=="PBS")$SampleName[1] == "DummyObservation") Num.PBS.obs[this.conc==c(10,30,100)] <- 0
      if (subset(MSdata,Conc==this.conc&Type=="Plasma")$SampleName[1] == "DummyObservation") Num.Plasma.obs[this.conc==c(10,30,100)] <- 0
    }  
    
    #mg/mL -> g/L is 1:1
    #kDa -> g/mol is *1000
    #g/mol -> M is g/L/MW
    #M <- uM is /1000000
    PPB100 <- 70/(66.5*1000)*1000000 # Berg and Lane (2011) 60-80 mg/mL, albumin is 66.5 kDa, pretend all protein is albumin to get uM
    C.frank <- 5 # uM frank parent concentration

    mydata <- list(                
      'Num.Blank.obs' = Num.Blank.obs,
      'T0.data' = rbind(subset(T0.data,Conc==10)$ISTDResponseRatio,subset(T0.data,Conc==30)$ISTDResponseRatio,subset(T0.data,Conc==100)$ISTDResponseRatio),
      'Blank.data' = rbind(subset(Blank.data,Conc==10)$ISTDResponseRatio,subset(Blank.data,Conc==30)$ISTDResponseRatio,subset(Blank.data,Conc==100)$ISTDResponseRatio),
      'Num.T0.obs' = Num.T0.obs,
      'PBS.data' = rbind(subset(PBS.data,Conc==10)$ISTDResponseRatio,subset(PBS.data,Conc==30)$ISTDResponseRatio,subset(PBS.data,Conc==100)$ISTDResponseRatio),
      'Num.PBS.obs' = Num.PBS.obs,
      'Plasma.data' = rbind(subset(Plasma.data,Conc==10)$ISTDResponseRatio,subset(Plasma.data,Conc==30)$ISTDResponseRatio,subset(Plasma.data,Conc==100)$ISTDResponseRatio),
      'Num.Plasma.obs' = Num.Plasma.obs,
      'C.frank' = C.frank,
      'C.protein' = c(PPB100*0.1,PPB100*0.3,PPB100)
    )
        
    initfunction <- function(chain)
    {
#      if (dim(mydata[["Blank.data"]])[1]>1)
#      {  
        BG <- apply(mydata[["Blank.data"]],1,function(x) mean(x))
#      } else {
#        BG <- mean(mydata[["Blank.data"]])
#      }
      BG[is.na(BG)] <- 0
      BG[BG<0] <- 0
      background <- rlnorm(length(BG),log(BG/2+10^-6),1)
 #     if (dim(mydata[["T0.data"]])[1]>1)
 #     {
        calibration <- apply(mydata$T0.data-background,1, function(x) max(10^-3.5,(mean(x)*5/mydata$C.frank)))
 #     } else {
 #       calibration <- max(10^-3.5,(mean(mydata$T0.data-background)*5/mydata$C.frank)) 
 #     }
      calibration[is.na(calibration)] <- 1
      calibration[calibration<=0] <- 10^-3.5
      C.missing <- runif(length(BG),0,mydata[["C.frank"]])
      log.C.thresh <- log(C.frank*runif(length(BG),0.01,0.5))

      return(list(
        .RNG.seed=as.numeric(paste(rep(chain,6),sep="",collapse="")),
        .RNG.name="base::Super-Duper",
        log.const.analytic.sd =runif(1,0.10,1),
        log.hetero.analytic.slope.factor = log10(runif(1,0, 1)),
        background = background,
        log.calibration = log10(calibration),
        log.C.thresh = log.C.thresh,
        log.Kd = runif(1,-3,0),
        C.missing = C.missing
      ))
    }
    
    save.image(CODA.FILE)  
    coda.out[[MSdata[,"CAS"][1]]] <- autorun.jags(RJAGS.MODEL.FILE, 
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
                             jags = JAGS.PATH,
                             monitor = c('Fup','Kd','log.const.analytic.sd','hetero.analytic.slope.factor','C.thresh','calibration','background','C.missing'))
    
    sim.mcmc <- coda.out[[MSdata[,"CAS"][1]]]$mcmc[[1]]
    for (i in 2:NUM.CHAINS) sim.mcmc <- rbind(sim.mcmc,coda.out[[MSdata[,"CAS"][1]]$mcmc[[i]])
    results <- apply(sim.mcmc,2,function(x) quantile(x,c(0.025,0.5,0.975)))

    Fup.point <- 2/5*(mean(mydata$PBS.data)-mean(mydata$Blank.data))/(mean(mydata$Plasma.data)-mean(mydata$Blank.data))
    
    new.results <- data.frame(Name=MSdata[,"Name"][1],
                              DTXSID=MSdata[,"DTXSID"][1],
                              CAS=MSdata[,"CAS"][1],
                              CompoundName=this.compound,
                              Fup.point=Fup.point,
                              stringsAsFactors=F)
    new.results[,c("Fup.Med","Fup.Low","Fup.High")] <- results[c(2,1,3),"Fup"]
    new.results[,c("Kd.Med","Kd.Low","Kd.High")] <- results[c(2,1,3),"Kd"]


    print(new.results)

    Results <- rbind(Results,new.results)

    write.table(Results,file=OUTPUT.FILE,sep=" ",row.names=FALSE)
  }    
}
setwd("Z:/Research Projects/CeetoxHumanHTTK/PPB-JAGS-Analyses/")
save.image(CODA.FILE)
stopCluster(CPU.cluster)

AFFINITY_Model_Results <- Results

missing <- subset(AFFINITY_Model_Results,is.na(Name))$CompoundName
for (this.id in missing)
{
  if (this.id %in% master.table$EPA_SAMPLE_ID)
  {
    AFFINITY_Model_Results[AFFINITY_Model_Results$CompoundName==this.id,"Name"] <- master.table[master.table$EPA_SAMPLE_ID==this.id,"Preferred_Name"]
    AFFINITY_Model_Results[AFFINITY_Model_Results$CompoundName==this.id,"DTXSID"] <- master.table[master.table$EPA_SAMPLE_ID==this.id,"DTXSID"]
    AFFINITY_Model_Results[AFFINITY_Model_Results$CompoundName==this.id,"CAS"] <- master.table[master.table$EPA_SAMPLE_ID==this.id,"CASRN"]
  }
}
View(AFFINITY_Model_Results)
subset(AFFINITY_Model_Results,(Fup.High-Fup.Low)>0.95&!(CompoundName%in%mislabeled.chems))

save.image(paste("PPB-affinity-analysis-",Sys.Date(),".RData",sep=""))
savehistory(paste("PPB-affinity-analysis-",Sys.Date(),".Rhistory",sep=""))


Results <- NULL
for (this.id in names(coda.out))
{
    sim.mcmc <- coda.out[[this.id]]$mcmc[[1]]
    for (i in 2:NUM.CHAINS) sim.mcmc <- rbind(sim.mcmc,coda.out[[this.id]]$mcmc[[i]])
    results <- apply(sim.mcmc,2,function(x) quantile(x,c(0.025,0.5,0.975)))

    Fup.point <- 2/5*(mean(mydata$PBS.data)-mean(mydata$Blank.data))/(mean(mydata$Plasma.data)-mean(mydata$Blank.data))
    
    new.results <- data.frame(Name=PPB.data[PPB.data$CAS==this.id,"Name"][1],
                              DTXSID=PPB.data[PPB.data$CAS==this.id,"DTXSID"][1],
                              CAS=this.id,
                              CompoundName=PPB.data[PPB.data$CAS==this.id,"CompoundName"][1],
                              Fup.point=Fup.point,
                              stringsAsFactors=F)
    new.results[,c("Fup.Med","Fup.Low","Fup.High")] <- results[c(2,1,3),"Fup"]
    new.results[,c("Kd.Med","Kd.Low","Kd.High")] <- results[c(2,1,3),"Kd"]
    new.results[,c("Fmissing_10.Med","Fmissing_10.Low","Fmissing_10.High")] <- results[c(2,1,3),"C.missing[1]"]/C.frank
    new.results[,c("Fmissing_30.Med","Fmissing_30.Low","Fmissing_30.High")] <- results[c(2,1,3),"C.missing[2]"]/C.frank
    new.results[,c("Fmissing_100.Med","Fmissing_100.Low","Fmissing_100.High")] <- results[c(2,1,3),"C.missing[3]"]/C.frank

    print(new.results)

    Results <- rbind(Results,new.results)

    write.table(Results,file="cmissing-043019.txt",sep=" ",row.names=FALSE)
}






