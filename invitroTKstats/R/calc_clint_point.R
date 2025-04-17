#' Calculate a Point Estimate of Intrinsic Hepatic Clearance (Clint) (Level-3)
#'
#' This function calculates a point estimate of intrinsic hepatic clearance (Clint) 
#' using mass spectrometry (MS) peak area data collected as part of in vitro measurement 
#' of chemical clearance, as characterized by the disappearance of parent compound over 
#' time when incubated with primary hepatocytes \insertCite{shibata2002prediction}{invitroTKstats}.
#'
#' The input to this function should be "Level-2" data. Level-2 data is Level-1,
#' data formatted with the \code{\link{format_clint}} function, and curated
#' with a verification column. "Y" in the verification column indicates the
#' data row is valid for analysis. 
#' 
#' The data frame of observations should be annotated according to
#' of these types:
#' \tabular{rrrrr}{
#'   Blank \tab Blank\cr
#'   Hepatocyte incubation concentration vs. time \tab Cvst\cr
#' }
#'
#' Clint is calculated using \code{\link{lm}} to perform a linear regression of
#' MS response as a function of time.
#'
#' @param FILENAME A string used to identify the input Level-2 file.
#' "<FILENAME>-Clint-Level2.tsv".
#' 
#' @param data.in (Data Frame) A Level-2 data frame generated from the 
#' \code{format_clint} function with a verification column added by 
#' \code{sample_verification}. Complement with manual verification if needed. 
#'
#' @param good.col (Character) Column name indicating which rows have been
#' verified, data rows valid for analysis are indicated with a "Y".
#' (Defaults to "Verified".)
#' 
#' @param output.res (Logical) When set to \code{TRUE}, the result 
#' table (Level-3) will be exported the current directory as a .tsv file. 
#' (Defaults to \code{TRUE}.)
#' 
#' @param sig.figs (Numeric) The number of significant figures to round the exported result table (Level-3). 
#' (Note: console print statements are also rounded to specified significant figures.)
#' (Defaults to \code{3}.)
#' 
#' @param INPUT.DIR (Character) Path to the directory where the input level-2 file exists. 
#' If \code{NULL}, looking for the input level-2 file in the current working
#' directory. (Defaults to \code{NULL}.)
#' 
#' @param OUTPUT.DIR (Character) Path to the directory to save the output file. 
#' If \code{NULL}, the output file will be saved to the current working
#' directory or \code{INPUT.DIR} if specified. (Defaults to \code{NULL}.)
#'
#' @return A Level-3 data frame with one row per chemical, containing a point estimate of intrinsic 
#' clearance (Clint), estimates of Clint of assays performed at 1 and 10 uM (if tested), 
#' the p-value and the Akaike Information Criterion (AIC) of the linear 
#' regression fit for all chemicals in the input data frame. 
#' 
#' @author John Wambaugh
#'
#' @examples
#' ## Load example level-2 data
#' level2 <- invitroTKstats::clint_L2
#' 
#' ## scenario 1: 
#' ## input level-2 data from the R session and do not export the result table
#' level3 <- calc_clint_point(data.in = level2, output.res = FALSE)
#' 
#' ## scenario 2: 
#' ## import level-2 data from a 'tsv' file and export the result table
#' \dontrun{
#' ## Refer to sample_verification help file for how to export level-2 data to a directory.
#' ## Unless a different path is specified in OUTPUT.DIR,
#' ## the result table will be saved to the directory specified in INPUT.DIR.
#' ## Will need to replace FILENAME and INPUT.DIR with name prefix and location of level-2 'tsv'.
#' level3 <- calc_clint_point(# e.g. replace with "Examples" from "Examples-Clint-Level2.tsv"
#'                            FILENAME="<level-2 FILENAME prefix>",
#'                            INPUT.DIR = "<level-2 FILE LOCATION>")
#' }
#'
#' @references
#' \insertRef{shibata2002prediction}{invitroTKstats}
#'
#' @importFrom stats4 mle coef AIC
#' @importFrom utils read.csv write.table
#'
#' @import Rdpack
#'
#' @export calc_clint_point
calc_clint_point <- function(
    FILENAME, 
    data.in,
    good.col="Verified", 
    output.res=TRUE, 
    sig.figs = 3,
    INPUT.DIR=NULL, 
    OUTPUT.DIR = NULL)
{
  
  #assigning global variables
  Compound.Name <- Response <- Sample.Type <- Time <- NULL
  
  
  if (!missing(data.in)) {
    clint.data <- as.data.frame(data.in)
  } else if (!is.null(INPUT.DIR)) {
      clint.data <- read.csv(file=paste0(INPUT.DIR, "/", FILENAME,"-Clint-Level2.tsv"),
                             sep="\t",header=T)
    } else {
      clint.data <- read.csv(file=paste0(FILENAME,"-Clint-Level2.tsv"),
                             sep="\t",header=T)
    }
  clint.data <- subset(clint.data,!is.na(Compound.Name))
  clint.data <- subset(clint.data,!is.na(Response))
  
  clint.cols <- c(L1.common.cols,
                  time.col = "Time",
                  test.conc.col = "Test.Compound.Conc",
                  test.nominal.conc.col = "Test.Nominal.Conc",
                  density.col = "Hep.Density"
  )
  list2env(as.list(clint.cols), envir = environment())
  cols <- c(unlist(mget(names(clint.cols))), "Response", good.col)
  
  if (!any(c("Biological.Replicates", "Technical.Replicates") %in% colnames(clint.data)))
    stop("Need at least one column representing replication, i.e. Biological.Replicates or Technical.Replicates. Run format_clint first (level 1) then curate to (level 2).")
  
  if (!(all(cols %in% colnames(clint.data))))
  {
    warning("Run format_clint first (level 1) then curate to (level 2).")
    stop(paste("Missing columns named:",
      paste(cols[!(cols%in%colnames(clint.data))],collapse=", ")))
  }
  

  # Only include the data types used:
  clint.data <- subset(clint.data,clint.data[,type.col] %in% c(
    "Blank","Cvst"))

  # Only used verified data:
  clint.data <- subset(clint.data, clint.data[,good.col] == "Y")

  # Clean up data:
  clint.data <- subset(clint.data,!is.na(Response))
  clint.data[clint.data$Response<0,"Response"] <- 0
  clint.data[clint.data$Sample.Type=="Blank" & is.na(clint.data$Time),"Time"] <- 0

  out.table <-NULL
  num.chem <- 0
  num.cal <- 0
  
  # decay - This function calculates the test compound concentration at time t, using a model of
  # exponential decay with time C(t) = C_0*e^{-mt}, where C_0 is the 
  # test compound concentration at time 0, m is a rate constant (argument k_elim),
  # and t is the incubation time in hour.
  decay <- function(time.hours,conc,cal,k_elim) cal*conc*exp(-k_elim*time.hours)
  
  # Negative log-likelihood of the linear regression fit
  lldecay <- function(cal,k_elim,sigma)
  {
    if (sigma < 0.0001) sigma <- 0.0001
    N <- dim(this.data)[1]
    pred <- decay(
      time.hours=this.data$Time,
      conc=this.data$Test.Nominal.Conc,
      cal=cal,
      k_elim=k_elim)
    ll <- log(1/sigma/sqrt(2*pi))*N
    res <- pred-this.data$Response
    ll <- ll+sum(-1/2*res^2/sigma^2)
    if (is.na(ll)){
      stop("lldecay - Estimated Log-likelihood is `NA`.")
      # browser()
    } 
    return(-ll)
  }
  
  # satdecay - This function calculates the test compound concentration at time t, 
  # using a model of exponential decay with time while considering a saturation probability. 
  # C(t) = C_0*e^{-m*sat*t}. C_0 is the test compound concentration at time 0, m is the elimination 
  # rate constant (argument k_elim), and t is the incubation time in hours. 
  # sat is the probability of saturation, defined as observing a lower clearance at a higher 
  # concentration. At 1 uM, sat is 1, meaning saturation is unlikely at the current concentration 
  # and is going to be observed at a higher concentration. At 10 uM, sat is between 0 and 1, 
  # meaning full saturation may or may not have been reached.
  satdecay <- function(time.hours,conc,cal,k_elim,sat) cal*conc*exp(-k_elim*ifelse(conc==10,sat,1)*time.hours)
  
  # Negative log-likelihood of the linear regression fit
  llsatdecay <- function(cal,k_elim,sigma,sat)
  {
    if (sigma < 0.0001) sigma <- 0.0001
    N <- dim(this.data)[1]
    pred <- satdecay(
      time.hours=this.data$Time,
      conc=this.data$Test.Nominal.Conc,
      cal=cal,
      k_elim=k_elim,
      sat=sat)
    ll <- log(1/sigma/sqrt(2*pi))*N
    res <- pred-this.data$Response
    ll <- ll+sum(-1/2*res^2/sigma^2)
    if (is.na(ll)){
      stop("llsatdecay - Estimated Log-likelihood is `NA`.")
      # browser()
    } 
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
    if (length(unique(this.cvt$Dilution.Factor))>1){
      stop("calc_clint_point - Cvst samples for `",this.chem,"` have more than one `Dilution.Factor`.")
      # browser()
    }
    df.cvt <- this.cvt$Dilution.Factor[1]
    if (length(unique(this.cvt$Hep.Density))>1){
      stop("calc_clint_point - Cvst samples for `",this.chem,"` have more than one `Hep.Density`.")
      # browser()
    } 
    hep.density <- this.cvt$Hep.Density[1]

    if (dim(this.cvt)[1] > 1)
    {
      this.data <- rbind(this.blank,this.cvt)
      this.data[this.data$Sample.Type=="Blank","Test.Nominal.Conc"] <- 0
      this.data[this.data$Sample.Type=="Blank","Time"] <- 0
      this.data[this.data$Sample.Type=="Cvst","Response"] <-
        this.data[this.data$Sample.Type=="Cvst","Response"]*df.cvt
      min.response <- sort(unique(this.data$Response))
      min.response <- min.response[min.response!=0]
      min.response <- min.response[1]
      this.data[this.data$Response==0,"Response"] <- min.response/2

      num.chem <- num.chem + 1
      num.cal <- num.cal + length(unique(this.data[,"Calibration"]))

      this.data$Response <- this.data$Response /
        mean(subset(this.data,Time==0)$Response)
      this.fit <- try(mle(lldecay,
        start=list(cal=1, k_elim=0.1, sigma=0.1),
        lower=list(cal=0,k_elim=0,sigma = 0.0001)))
      this.null <- try(mle(lldecay, 
        start=list(cal=1, sigma=0.1),
        lower=list(cal=0, sigma = 0.0001),
        fixed=list(k_elim=0)))

      if (!inherits(this.fit, "try-error") & !inherits(this.null, "try-error"))
      {
        # k_elim has units 1/h, convert to uL/min/10^6 hepatocytes
        # hep density is 10^6 hepatocytes/mL
        this.row$Clint <- 1000*coef(this.fit)["k_elim"]/hep.density/60
        this.row$Clint.pValue <- min(exp(-(AIC(this.null)-AIC(this.fit))),1)
        this.row$Fit <- paste(paste(unique(this.data$Test.Nominal.Conc),collapse=", "),"uM")
        this.row$AIC <- AIC(this.fit)
        this.row$AIC.Null <- AIC(this.null)
        this.row$Clint.1 <- NA
        this.row$Clint.10 <- NA
        this.row$AIC.Sat <- NA
        this.row$Sat.pValue <- NA
        if (all(c(1,10)%in%unique(this.data$Test.Nominal.Conc)))
        {
          this.sat.fit <- try(mle(llsatdecay,
            start=list(cal=1, k_elim=0.1, sigma=0.1, sat=0.5),
            lower=list(cal=0, k_elim=0, sigma = 0.0001, sat=0),
            upper=list(sat=1)))
          if (!inherits(this.sat.fit, "try-error"))
          {
            this.row$Clint.1 <- 1000*coef(this.sat.fit)["k_elim"]/hep.density/60
            this.row$Clint.10 <- 1000*coef(this.sat.fit)["k_elim"]*
              coef(this.sat.fit)["sat"]/hep.density/60
            this.row$AIC.Sat <- AIC(this.sat.fit)
            if (this.row$Clint.pValue==1) test.AIC <- this.row$AIC.Null
            else test.AIC <- this.row$AIC
            this.row$Sat.pValue <- min(exp(-(test.AIC-AIC(this.sat.fit))),1)
          } else{
            # warning message to users to indicate failed fitting
            warning("calc_clint_point - Saturation decay fit resulted in an error when fitting. Returning `NA` for the following ourputs:\n\t",
                    paste0(c("Clint.1","Clint.10","AIC.Sat","AIC","Sat.pValue"),collapse = ", "))
            # assign NA's to outputs since saturation decay fit failed
            this.row$Clint.1    <- NA
            this.row$Clint.10   <- NA
            this.row$AIC.Sat    <- NA
            this.row$Sat.pValue <- NA
            # browser()
          } 
        }
        if (!is.null(sig.figs)){
          # Print results to desired sig.figs or default of 3  
          print(paste(
            this.row$Compound.Name,
            "Cl_int =",
            signif(this.row$Clint,sig.figs),
            "uL/min/million hepatocytes, p-Value =",
            signif(this.row$Clint.pValue,sig.figs),
            "."
          ))
        } else {
          # If sig.figs = NULL
          print(paste(
            this.row$Compound.Name,
            "Cl_int =",
            this.row$Clint,
            "uL/min/million hepatocytes, p-Value =",
            this.row$Clint.pValue,
            "."
          ))
        }
      } else {
        for (col in c("Fit","AIC","AIC.Null","Clint.1","Clint.10","AIC.Sat","Sat.pValue"))
          this.row[,col] <- NA
        this.row$Clint <- "Linear Regression Failed"
        cat("Linear regression failed for:",this.chem,".\n")
        plot(this.data$Time, this.data$Response)
        # browser()
      }
      out.table <- rbind(out.table, this.row)
    }
  }

  out.table <- as.data.frame(out.table)
  rownames(out.table) <- make.names(out.table$Compound.Name, unique=TRUE)
  #out.table <- apply(out.table,2,unlist)
  out.table[!(out.table[,"Clint"]%in%"Linear Regression Failed"),"Clint"] <-
    as.numeric(out.table[
    !(out.table[,"Clint"]%in%"Linear Regression Failed"),"Clint"])
  out.table[,"Clint.1"] <- as.numeric(out.table[,"Clint.1"])
  out.table[,"Clint.10"] <- as.numeric(out.table[,"Clint.10"])
  out.table[,"Clint.pValue"] <- as.numeric(out.table[,"Clint.pValue"])
  out.table[,"AIC"] <- as.numeric(out.table[,"AIC"])
  out.table[,"AIC.Null"] <- as.numeric(out.table[,"AIC.Null"])
  out.table[,"AIC.Sat"] <- as.numeric(out.table[,"AIC.Sat"])
  out.table[,"Sat.pValue"] <- as.numeric(out.table[,"Sat.pValue"])

  if (output.res) {
    # Determine the path for output
    if (!is.null(OUTPUT.DIR)) {
      file.path <- OUTPUT.DIR
    } else if (!is.null(INPUT.DIR)) {
      file.path <- INPUT.DIR
    } else {
      file.path <- getwd()
    }
    
    rounded.out.table <- out.table
    
    # Round results to desired number of sig figs
    if (!is.null(sig.figs)){
      rounded.out.table[!(rounded.out.table[,"Clint"]%in%"Linear Regression Failed"),"Clint"] <-
        signif(rounded.out.table[
          !(rounded.out.table[,"Clint"]%in%"Linear Regression Failed"),"Clint"],sig.figs)
      rounded.out.table[,"Clint.1"] <- signif(rounded.out.table[,"Clint.1"],sig.figs)
      rounded.out.table[,"Clint.10"] <- signif(rounded.out.table[,"Clint.10"],sig.figs)
      rounded.out.table[,"Clint.pValue"] <- signif(rounded.out.table[,"Clint.pValue"],sig.figs)
      rounded.out.table[,"AIC"] <- signif(rounded.out.table[,"AIC"],sig.figs)
      rounded.out.table[,"AIC.Null"] <- signif(rounded.out.table[,"AIC.Null"],sig.figs)
      rounded.out.table[,"AIC.Sat"] <- signif(rounded.out.table[,"AIC.Sat"],sig.figs)
      rounded.out.table[,"Sat.pValue"] <- signif(rounded.out.table[,"Sat.pValue"],sig.figs)
      cat(paste0("\nData to export has been rounded to ", sig.figs, " significant figures.\n"))
    }
    
    # Write out a "level 3" file:
    write.table(rounded.out.table,
                file=paste0(file.path, "/", FILENAME,"-Clint-Level3.tsv"),
                sep="\t",
                row.names=F,
                quote=F)
    
    # Print notification message stating where the file was output to
    cat(paste0("A Level-3 file named ",FILENAME,"-Clint-Level3.tsv", 
                " has been exported to the following directory: ", file.path), "\n")
  }

  print(paste("Intrinsic clearance (Clint) calculated for",num.chem,"chemicals."))
  print(paste("Intrinsic clearance (Clint) calculated for",num.cal,"measurements."))

  return(out.table)
}


