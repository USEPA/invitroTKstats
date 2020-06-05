library(ggplot2)
library(scales)

scientific_10 <- function(x) {                                  
  out <- gsub("1e", "10^", scientific_format()(x))              
  out <- gsub("\\+","",out)                                     
  out <- gsub("10\\^01","10",out)                               
  out <- parse(text=gsub("10\\^00","1",out))                    
}  

plot_uc_results <- function(dat,bayes,compound,cal,MW,quad.cal=NULL,cal.col="Cal",name.col="Compound.Name")
{
  sim.mcmc <- bayes[[compound]]$mcmc[[1]]
  for (i in 2:length(bayes[[compound]]$mcmc))
  {
    sim.mcmc <- rbind(sim.mcmc,bayes[[compound]]$mcmc[[i]])
  }
  results <- apply(sim.mcmc,2,function(x) quantile(x,c(0.025,0.5,0.975)))
        
  this.data <- dat[dat[,name.col]==compound,]
  all.cal <- unique(this.data[,cal.col])
  this.cal <- which(all.cal == cal)
  this.data <- dat[dat[,cal.col]==cal,]
                         
  cal.slope <- results["50%",paste("calibration[",this.cal,"]",sep="")]
  cal.int <- results["50%",paste("background[",this.cal,"]",sep="")]
  
  plotcal.curves <- data.frame(Conc=10^seq(-4,1,by=0.1))
  cal.curves$Bayesian <- cal.slope*cal.curves$Conc + cal.int

  if (!is.null(quad.cal))
  {
    a <- quad.cal[["avar"]]
    b <- quad.cal[["bvar"]]
    cvar <- quad.cal[["cvar"]]
    cal.curves$Quadratic <- a*(cal.curves$Conc*305/MW)^2 + b*(cal.curves$Conc*305/MW)+ cvar
  } else cal.curves$Quadratic <-1000
  
  plinear <- ggplot(this.data, aes(Response, Nominal.Conc)) +
    geom_point(size=3) +
    scale_y_continuous(limits=c(0,max(subset(all.data,Cal=="010720")$Nominal.Conc,na.rm=T))) +
    scale_x_continuous(limits=c(0,max(subset(all.data,Cal=="010720")$Response,na.rm=T))) +
    ylab("Nominal Concentration (uM)") +
    xlab("Normalized Response") +
    geom_vline(xintercept = subset(this.data,Sample.Type=="AF")$Response,color="Red") +
    geom_vline(xintercept = subset(this.data,Sample.Type=="T5")$Response,color="Blue") +
  #  geom_abline(slope=1/4.267595, intercept=0.00295*(1-1/4.267595),linetype="dashed") +
    geom_line(data=cal.curves,aes(x=Bayesian,y=Conc),linetype="dashed")+
    geom_line(data=cal.curves,aes(x=Quadratic,y=Conc),linetype="dotted")+
    ggtitle(paste(compound,cal)) +
    theme(text = element_text(size=20))
    
  plog <- ggplot(this.data, aes(Response, Nominal.Conc)) +
    geom_point(size=3) +
    scale_y_log10(label=scientific_10) +
    scale_x_log10(label=scientific_10) +
    ylab("Nominal Concentration (uM)") +
    xlab("Normalized Response") +
    geom_vline(xintercept = subset(this.data,Sample.Type=="AF")$Response,color="Red") +
    geom_vline(xintercept = subset(this.data,Sample.Type=="T5")$Response,color="Blue") +
    geom_line(data=cal.curves,aes(x=Bayesian,y=Conc),linetype="dashed")+
    geom_line(data=cal.curves,aes(x=Quadratic,y=Conc),linetype="dotted")+
    ggtitle(paste(compound,cal)) +
    theme(text = element_text(size=20))

  return(list(plinear = plinear, plog = plog))
}
