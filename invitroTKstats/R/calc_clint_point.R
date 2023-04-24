#' Calculate a point estimate of intrinsic hepatic clearance
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
#' argument given, "-Clint-Level3.tsv" is appended (defaults to "MYDATA")
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
#' level3 <- calc_clint_point(FILENAME="Wambaugh2019")
#'  
#' # Just 1 uM data:
#' write.table(subset(level2,Conc==1),
#'   file="Wambaugh2019-1-Clint-Level2.tsv",
#'   sep="\t",
#'   row.names=F,
#'   quote=F)
#' 
#' level3.1 <- calc_clint_point(FILENAME="Wambaugh2019-1")
#' 
#' # Just 10 uM data:
#' write.table(subset(level2,Conc==10),
#'   file="Wambaugh2019-10-Clint-Level2.tsv",
#'   sep="\t",
#'   row.names=F,
#'   quote=F) 
#' 
#' level3.10 <- calc_clint_point(FILENAME="Wambaugh2019-10")
#'
#' @references
#' \insertRef{shibata2002prediction}{invitroTKstats}
#'
#' @importFrom stats4 mle coef AIC
#'
#' @export calc_clint_point
calc_clint_point <- function(FILENAME, good.col="Verified")
{
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
  std.conc.col <- "Std.Conc"
  clint.assay.conc.col <- "Clint.Assay.Conc"
  time.col <- "Time"
  area.col <- "Area"
  analysis.method.col <- "Analysis.Method"
  analysis.instrument.col <- "Analysis.Instrument"
  analysis.parameters.col <- "Analysis.Parameters" 
  note.col <- "Note"

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
    std.conc.col,
    clint.assay.conc.col,
    time.col,
    area.col,
    analysis.method.col,
    analysis.instrument.col,
    analysis.parameters.col,
    note.col
    )
      
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
  clint.data[clint.data$Sample.Type=="Blank" & is.na(clint.data$Time),"Time"] <- 0
                  
  out.table <-NULL
  num.chem <- 0
  num.cal <- 0
  decay <- function(time.min,conc,cal,k_elim) cal*conc*exp(-k_elim*time.min)
  lldecay <- function(cal,k_elim,sigma)
  {
    if (sigma < 0.0001) sigma <- 0.0001
    N <- dim(this.data)[1]
    pred <- decay(
      time.min=this.data$Time,
      conc=this.data$Clint.Assay.Conc,
      cal=cal,
      k_elim=k_elim)
    ll <- log(1/sigma/sqrt(2*pi))*N
    res <- pred-this.data$Response
    ll <- ll+sum(-1/2*res^2/sigma^2)
    if (is.na(ll)) browser()
    return(-ll)
  }
  satdecay <- function(time.min,conc,cal,k_elim,sat) cal*conc*exp(-k_elim*ifelse(conc==10,sat,1)*time.min)
  llsatdecay <- function(cal,k_elim,sigma,sat)
  {
    if (sigma < 0.0001) sigma <- 0.0001
    N <- dim(this.data)[1]
    pred <- satdecay(
      time.min=this.data$Time,
      conc=this.data$Clint.Assay.Conc,
      cal=cal,
      k_elim=k_elim,
      sat=sat)
    ll <- log(1/sigma/sqrt(2*pi))*N
    res <- pred-this.data$Response
    ll <- ll+sum(-1/2*res^2/sigma^2)
    if (is.na(ll)) browser()
    return(-ll)
  }

  for (this.chem in unique(clint.data[,compound.col]))
  {     
    this.subset <- subset(clint.data,clint.data[,compound.col]==this.chem)
    this.dtxsid <- this.subset$dtxsid[1]
    this.row <- cbind(this.subset[1,c(compound.col,dtxsid.col,lab.compound.col)],
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
      this.data <- rbind(this.blank,this.cvt)
      this.data[this.data$Sample.Type=="Blank","Clint.Assay.Conc"] <- 0
      this.data[this.data$Sample.Type=="Blank","Time"] <- 0
      this.data[this.data$Sample.Type=="Cvst","Response"] <-
        this.data[this.data$Sample.Type=="Cvst","Response"]*df.cvt
      min.response <- sort(unique(this.data$Response))
      min.response <- min.response[min.response!=0]
      min.response <- min.response[1]
      this.data[this.data$Response==0,"Response"] <- min.response/2

      num.chem <- num.chem + 1
      num.cal <- length(unique(this.data[,"Calibration"]))
      
      this.fit <- try(mle(lldecay,
        start=list(cal=1,k_elim=0.1,sigma=1),
        lower=list(cal=0,k_elim=0,sigma = 0.0001)))
      this.null <- try(mle(lldecay,
        start=list(cal=1,sigma=1),
        lower=list(cal=0,sigma = 0.0001),
        fixed=list(k_elim=0)))      
      
      if (class(this.fit)!="try-error" & class(this.null)!="try-error")
      {
        # k_elim has units 1/h, convert to uL/min/10^6 hepatocytes
        # hep density is 10^6 hepatocytes/mL
        this.row$Clint <- 1000*coef(this.fit)["k_elim"]/hep.density/60
        this.row$Clint.pValue <- min(exp(-(AIC(this.null)-AIC(this.fit))),1)
        this.row$Fit <- paste(paste(unique(this.data$Clint.Assay.Conc),collapse=", "),"uM")
        this.row$AIC <- AIC(this.fit)
        this.row$AIC.Null <- AIC(this.null)
        this.row$Clint.1 <- NA
        this.row$Clint.10 <- NA
        this.row$AIC.Sat <- NA
        this.row$Sat.pValue <- NA
        if (all(c(1,10)%in%unique(this.data$Clint.Assay.Conc)))
        {
          this.sat.fit <- try(mle(llsatdecay,
            start=list(cal=1,k_elim=0.1,sigma=1,sat=1),
            lower=list(cal=0,k_elim=0,sigma = 0.0001,sat=0),
            upper=list(sat=1)))
          if (class(this.sat.fit)!="try-error")
          {
            this.row$Clint.1 <- 1000*coef(this.sat.fit)["k_elim"]/hep.density/60
            this.row$Clint.10 <- 1000*coef(this.sat.fit)["k_elim"]*
              coef(this.sat.fit)["sat"]/hep.density/60
            this.row$AIC.Sat <- AIC(this.sat.fit)
            if (this.row$Clint.pValue==1) test.AIC <- this.row$AIC.Null
            else test.AIC <- this.row$AIC
            this.row$Sat.pValue <- min(exp(-(test.AIC-AIC(this.sat.fit))),1)
          } else browser()
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
      } else browser()
    }  
  }

  out.table <- as.data.frame(out.table)
  rownames(out.table) <- make.names(out.table$Compound.Name, unique=TRUE)
  #out.table <- apply(out.table,2,unlist) 
  out.table[,"Clint"] <- signif(as.numeric(out.table[,"Clint"]),3) 
  out.table[,"Clint.1"] <- signif(as.numeric(out.table[,"Clint.1"]),3) 
  out.table[,"Clint.10"] <- signif(as.numeric(out.table[,"Clint.10"]),3) 
  out.table[,"Clint.pValue"] <- signif(as.numeric(out.table[,"Clint.pValue"]),3) 
  out.table[,"AIC"] <- signif(as.numeric(out.table[,"AIC"]),3)
  out.table[,"AIC.Null"] <- signif(as.numeric(out.table[,"AIC.Null"]),3)
  out.table[,"AIC.Sat"] <- signif(as.numeric(out.table[,"AIC.Sat"]),3)
  out.table[,"Sat.pValue"] <- signif(as.numeric(out.table[,"Sat.pValue"]),3) 
    
# Write out a "level 3" file (data organized into a standard format):  
  write.table(out.table, 
    file=paste(FILENAME,"-Clint-Level3.tsv",sep=""),
    sep="\t",
    row.names=F,
    quote=F)
 
  print(paste("Intrinsic clearance (Clint) calculated for",num.chem,"chemicals."))
  print(paste("Intrinsic clearance (Clint) calculated for",num.cal,"measurements."))

  return(out.table)  
}


