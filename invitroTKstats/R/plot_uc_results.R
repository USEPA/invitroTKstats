#' Plot Mass Spectrometry Responses for Fraction Unbound in Plasma Data from
#' Ultracentrifugation (UC) with Calibration Curves
#'
#' This function generates two scatter plots of mass spectrometry (MS) responses 
#' by nominal concentration for one chemical, overlaid with Bayesian calibration curves.
#' The MS responses were collected from measurement of fraction unbound in plasma (Fup)
#' using ultracentrifugation (UC). The two scatter plots are differentiated by the unit of axes. 
#' One plot uses axes in the original unit and one plot uses axes in log-10 scale. 
#' Users can use the returned plots for a reality check on the data collected and the 
#' Bayesian results, see if the data behave as expected.
#'
#' @param dat (Data Frame) Level-2 Ultracentrifugation (UC) Fup assays data. 
#' This should be the same Level-2 data used to generate \code{bayes}.
#'
#' @param bayes (List) MCMC results returned from \code{calc_fup_uc}.
#'
#' @param compound (Character) Name of the compound to be plotted.
#'
#' @param cal (Character) Which calibration the data to be plotted is based on.
#'
#' @param MW (Numeric) Molecular weight.
#'
#' @param quad.cal (List or Named Vector) Quadratic calibration curve results from Mass Spec. 
#' (Defaults to \code{NULL}.)
#'
#' @param cal.col (Character) Column name from \code{dat} containing the 
#' the calibration information. Typically, this uses indices or dates to represent 
#' if the analyses were done on different machines on 
#' the same day or on different days with the same MS analyzer. (Defaults to "Cal".)
#'
#' @param name.col (Character) Column name from \code{dat} containing compound names. 
#' (Defaults to "Compound.Name".)
#'
#' @param uc.dilute (Numeric) Number of times the sample was diluted before MS 
#' analysis. (Defaults to 5.)
#'
#' @param af.dilute (Numeric) Number of times the Aqueous Fraction samples were diluted before MS 
#' analysis. (Defaults to 2.)
#'
#' @return A list of two scatter plots of MS responses versus nominal concentration,
#' overlaid with a Bayesian calibration curve. One plot displays the axes in the original unit and 
#' one displays the axes in log-10 scale.
#'
#' @author John Wambaugh
#'
#' @import ggplot2
#'
#' @export plot_uc_results
plot_uc_results <- function(dat,bayes,compound,cal,MW,quad.cal=NULL,cal.col="Cal",name.col="Compound.Name",uc.dilute=5,af.dilute=2)
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

  cal.slope.low <- results["50%",paste("calibration.low[",this.cal,"]",sep="")]
  cal.slope.high <- results["50%",paste("calibration.high[",this.cal,"]",sep="")]
  cal.int <- results["50%",paste("background[",this.cal,"]",sep="")]
  C.cal.thresh <- results["50%",paste("C.cal.thresh[",this.cal,"]",sep="")]

  fup.025 <- signif(results["2.5%","Fup"],2)
  fup.med <- signif(results["50%","Fup"],2)
  fup.975<- signif(results["97.5%","Fup"],2)
  fup <- paste("Fup = ",fup.med," (",fup.025," - ",fup.975,")",sep="")

  input.data <- runjagsdata.to.list(bayes[[compound]]$data)
  Num.cc.obs <- input.data[["Num.cc.obs"]]
  Num.series <- input.data[["Num.series"]]
  # Dilution.factor <- input.data[["Dilution.Factor"]]

  UC.conc <- results["50%",paste("Conc[",Num.cc.obs+this.cal,"]",sep="")]
  AF.conc <- results["50%",paste("Conc[",Num.cc.obs+this.cal+Num.series,"]",sep="")]

  cal.curves <- data.frame(Conc=10^seq(-4,1,by=0.1))
  cal.curves$Bayesian <- (cal.slope.low +
      cal.slope.high *
      Heaviside(cal.curves$Conc,threshold=C.cal.thresh)) *
      cal.curves$Conc +
    cal.int -
    C.cal.thresh * cal.slope.high *
      Heaviside(cal.curves$Conc,threshold=C.cal.thresh)

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
    ylab("Sample Concentration (uM)") +
    xlab("Normalized Response") +
    geom_vline(xintercept = subset(this.data,Sample.Type=="AF")$Response,color="Red") +
    geom_vline(xintercept = subset(this.data,Sample.Type=="T5")$Response,color="Blue") +
    geom_hline(yintercept = AF.conc/af.dilute,color="Red",linetype="dashed") +
    annotate("text", x = 1.1, y = 0.01+AF.conc/af.dilute, label = paste(
      "AF = ",
      signif(AF.conc,3),
      " (",
      af.dilute,
      "-fold)",
      sep="")) +
    geom_hline(yintercept = UC.conc/uc.dilute ,color="Blue",linetype="dashed") +
    annotate("text", x = 1.1, y = 0.01+UC.conc/uc.dilute, label = paste(
      "T5 =",
      signif(UC.conc,3),
      " (",
      uc.dilute,
      "-fold)",
      sep="")) +
  #  geom_abline(slope=1/4.267595, intercept=0.00295*(1-1/4.267595),linetype="dashed") +
    geom_line(data=cal.curves,aes(x=Bayesian,y=Conc),linetype="dashed")+
    geom_line(data=cal.curves,aes(x=Quadratic,y=Conc),linetype="dotted")+
    ggtitle(paste(compound,cal,fup)) +
    theme(text = element_text(size=18))

  plog <- ggplot(this.data, aes(Response, Nominal.Conc)) +
    geom_point(size=3) +
    scale_y_log10(label=scientific_10) +
    scale_x_log10(label=scientific_10) +
    ylab("Sample Concentration (uM)") +
    xlab("Normalized Response") +
    geom_vline(xintercept = subset(this.data,Sample.Type=="AF")$Response,color="Red") +
    geom_vline(xintercept = subset(this.data,Sample.Type=="T5")$Response,color="Blue") +
    geom_hline(yintercept = (AF.conc/af.dilute),color="Red",linetype="dashed") +
    annotate("text", x = 5, y = 1.5*AF.conc/af.dilute, label = paste("AF =",signif(AF.conc,3))) +
    geom_hline(yintercept = UC.conc/uc.dilute ,color="Blue",linetype="dashed") +
    annotate("text", x = 5, y = 1.5*UC.conc/uc.dilute, label = paste("T5 =",signif(UC.conc,3))) +
    geom_line(data=cal.curves,aes(x=Bayesian,y=Conc),linetype="dashed")+
    geom_line(data=cal.curves,aes(x=Quadratic,y=Conc),linetype="dotted")+
    ggtitle(paste(compound,cal,fup)) +
    theme(text = element_text(size=18))

  return(list(plinear = plinear, plog = plog))
}
