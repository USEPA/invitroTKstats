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
#' argument given, "-PPB-RED-Level2.tsv" is appended (defaults to "MYDATA")
#' 
#' @param good.col Name of a column indicating which rows have been verified for 
#' analysis, indicated by a "Y" (Defaults to "Verified")
#'
#' @return \item{data.frame}{A data.frame in standardized format} 
#'
#' @author John Wambaugh
#' 
#' @references
#' Shibata, Yoshihiro, Hiroyuki Takahashi, and Yasuyuki Ishii. "A convenient in 
#' vitro screening method for predicting in vivo drug metabolic clearance using 
#' isolated hepatocytes suspended in serum." Drug metabolism and disposition 
#' 28.12 (2000): 1518-1523.
#'
#' @importFrom stats4 mle coef AIC
#'
#' @export calc_clint_point
calc_clint_point <- function(FILENAME, good.col="Verified")
{
  clint.data <- read.csv(file=paste(FILENAME,"-PPB-RED-Level2.tsv",sep=""), 
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

  out.table <-NULL
  num.chem <- 0
  num.cal <- 0
  decay <- function(time.min,conc,cal,k_elim) cal*conc*exp(-k_elim*time.min)
  lldecay <- function(cal,k_elim,sigma)
  {
    N <- dim(this.data)[1]
    pred <- decay(
      time.min=this.data$Time,
      conc=this.data$Conc,
      cal=cal,
      k_elim=k_elim)
    ll <- log(1/sigma/sqrt(2*pi))*N
    res <- pred-this.data$Response
    ll <- ll+sum(-1/2*res^2/sigma^2)
    return(-ll)
  }

  for (this.chem in unique(clint.data[,compound.col]))
  {     
    this.subset <- subset(clint.data,clint.data[,compound.col]==this.chem)
    this.dtxsid <- this.subset$dtxsid[1]
    this.row <- c(this.subset[1,c(compound.col,dtxsid.col)],
      data.frame(Calibration="All Data",
        Fup=NaN))
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
        lower=list(cal=0,k_elim=0,sigma = 0.01)))
      this.null <- try(mle(lldecay,
        start=list(cal=1,sigma=1),
        lower=list(cal=0,sigma = 0.01),
        fixed=list(k_elim=0)))      
      
      if (class(this.fit)!="try-error" & class(this.null)!="try-error")
      {
        this.row$Clint <- 1000*coef(this.fit)["k_elim"]/hep.density
        this.row$Clint.pValue <- exp(-(AIC(this.null)-AIC(this.fit)))
        this.row$Fit <- "All Concs"
        this.row$AIC <- AIC(this.fit)
        this.row$AIC.Null <- AIC(this.null)
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
  out.table[,"Clint.pValue"] <- signif(as.numeric(out.table[,"Clint.pValue"]),3) 
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
}


