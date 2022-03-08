  # Here's the new R package for analyzing these data:
  
  library(invitroTKstats)
  
  # There are multiple packages for loading Excel files, but I've been using this
  # one lately:
  library(readxl)
  
  # Change to the directory that has your Excel file 
  # (this is where it is on my computer):
  setwd("c:/users/jwambaug/git/invitroTKstats/working/")
  
  # Assumed dilution factors:
  CC.DILUTE <- 1
  BLANK.DILUTE <- 1
  AF.DILUTE <- 2*16
  T5.DILUTE <- 5*16
  T1.DILUTE <- 5*16
  
  assayinfo <- read_excel(
             "20220201_PFAS-LC_FractionUnbound_MGS.xlsx",
             sheet=1)
  assayinfo <- as.data.frame(assayinfo)
  
  cheminfo <- read_excel(
             "20220201_PFAS-LC_FractionUnbound_MGS.xlsx",
             sheet=2)
  cheminfo <- as.data.frame(cheminfo)
  
  
  
  # Fill in the date column for all rows:
  this.date <- assayinfo[1,"LCMS Analysis Date"]
  this.row <- 1
  while (this.row <= dim(assayinfo)[1])
  {
    if (is.na(assayinfo[this.row,"LCMS Analysis Date"]))
    {
      assayinfo[this.row,"LCMS Analysis Date"] <- this.date 
    } else {
      this.date <- assayinfo[this.row,"LCMS Analysis Date"] 
    }
    this.row <- this.row + 1
  }
  # Get rid of the UTC:
  assayinfo[,1] <- sapply(assayinfo[,1],function(x) gsub("  UTC","",x))
  
  # Define the info for control chemicals:
  controls <- data.frame(
    Compound="n-butylparaben",
    DTXSID="DTXSID3020209",
    ISTD="13C6-n-butylparaben")
                                       
  # All the data (for now) are in different tabs (by date) of a single file 
  # (still going with a list called files but each "file") is now a tab:
  files <- list()
  # This list stores the chemical ID's in each file
  compounds <- list()           
  DTXSIDs <- list()           
  # This list stores the internal standards for each chemical in each file:
  ISTDs <- list()
  
  for (this.date in unique(assayinfo[,1]))
  {
    this.sheet <- gsub("-","",this.date)
    this.index <- which(unique(assayinfo[,1])==this.date)
        
    files[[this.index]] <- as.data.frame(read_excel(
             "20220201_PFAS-LC_FractionUnbound_MGS.xlsx",
             sheet=this.sheet))
    files[[this.index]]$Date <- this.date        
    this.assayinfo <- subset(assayinfo,regexpr(this.date,assayinfo[,1])!=-1)
    
    compounds[[this.index]] <- this.assayinfo[,"Analyte"]
    DTXSIDs[[this.index]] <- this.assayinfo[,"DTXSID"]
    ISTDs[[this.index]] <- this.assayinfo[,"Int Std"]
  }
           
  
  
  # Now loop over the files and annotate the chemicals:
  UC.data <- NULL
  for (this.file in 1:length(files))
  {
    total.rows <- dim(files[[this.file]])[1]
    this.compound <- 1
    
    files[[this.file]]$DTXSID <- DTXSIDs[[this.file]][this.compound] 
    files[[this.file]]$ISTD.Name <- ISTDs[[this.file]][this.compound] 
    for (this.row in 1:total.rows)
    {
      if (!is.na(as.character(files[[this.file]][this.row,1])))
        if (regexpr("Compound",as.character(files[[this.file]][this.row,1]))!=-1)
        {
          this.id <- NA
          # Try to get the id using two spaces after colon:
          this.id <- 
            strsplit(as.character(files[[this.file]][this.row,1]),":  ")[[1]][2]
          # If is.na, try with only one:
          if (is.na(this.id)) this.id <- 
            strsplit(as.character(files[[this.file]][this.row,1]),": ")[[1]][2]
          
          this.chem.index <- regexpr(paste("\\(",this.id,"\\)",sep=""), 
            cheminfo[,2])!=-1
          # check to see if this is one of the test chemicals:
          if (any(this.chem.index))
          {
            this.dtxsid <- cheminfo[this.chem.index,"DTXSID"]
            if (length(this.dtxsid)>1)
            {
              match.strings <- unlist(lapply(
                strsplit(cheminfo[this.chem.index,2],this.id),
                function(x) x[[1]]))
              match.strings <- gsub(" \\(","",match.strings)
              for (this.string in match.strings)
                if (any(regexpr(tolower(this.string),
                  tolower(compounds[[this.file]]))!=-1))
                {
                  this.assay.index <- which(regexpr(tolower(this.string),
                  tolower(compounds[[this.file]]))!=-1)
                  files[[this.file]][this.row:total.rows,"Compound.Name"] <- 
                    compounds[[this.file]][this.assay.index]
                  files[[this.file]][this.row:total.rows,"DTXSID"] <- 
                    DTXSIDs[[this.file]][this.assay.index]
                  this.dtxsid <- DTXSIDs[[this.file]][this.assay.index]
                  files[[this.file]][this.row:total.rows,"ISTD.Name"] <- 
                    ISTDs[[this.file]][this.assay.index]
                }
            }
            this.assay.index <- which(DTXSIDs[[this.file]] == this.dtxsid)
            if (any(this.assay.index))
            {
              files[[this.file]][this.row:total.rows,"Compound.Name"] <- 
                compounds[[this.file]][this.assay.index]
              files[[this.file]][this.row:total.rows,"DTXSID"] <- 
                this.dtxsid
              files[[this.file]][this.row:total.rows,"ISTD.Name"] <- 
                ISTDs[[this.file]][this.assay.index]
            }
          # Check to see if we can find the chemical in the assay table:
          } else if (any(
            regexpr(tolower(this.id), tolower(assayinfo[,"Analyte"]))!=-1 &
            assayinfo[,"LCMS Analysis Date"] == files[[this.file]][this.row,"Date"]))
          {
            this.assay.index <- which(any(
              regexpr(tolower(this.id), tolower(assayinfo[,"Analyte"]))!=-1 &
              assayinfo[,"LCMS Analysis Date"] == files[[this.file]][this.row,"Date"]))
            files[[this.file]][this.row:total.rows,"Compound.Name"] <- 
              compounds[[this.file]][this.assay.index]
            files[[this.file]][this.row:total.rows,"DTXSID"] <- 
              DTXSIDs[[this.file]][this.assay.index]
            files[[this.file]][this.row:total.rows,"ISTD.Name"] <- 
              ISTDs[[this.file]][this.assay.index]          
          } else {
          # Maybe it's a control:
            this.index <- regexpr(tolower(this.id), tolower(controls[,"Compound"]))!=-1
            if (any(this.index))
            {
              this.compound <- this.compound + 1
              files[[this.file]][this.row:total.rows,"Compound.Name"] <- 
                controls[this.index,"Compound"]
              files[[this.file]][this.row:total.rows,"DTXSID"]  <- 
                controls[this.index,"DTXSID"]
              files[[this.file]][this.row:total.rows,"ISTD.Name"]  <- 
                controls[this.index,"ISTD"]
            } else {
            # Assume it's an internal standard:
              files[[this.file]][this.row:total.rows,"Compound.Name"] <- this.id
              files[[this.file]][this.row:total.rows,"DTXSID"] <- "ISTD"
              files[[this.file]][this.row:total.rows,"ISTD.Name"] <- "Eponymous"
            }
          } 
        }
    }
    UC.data <- rbind(UC.data,files[[this.file]])
  }
  dim(UC.data)
  misses <- unique(subset(UC.data,DTXSID=="ISTD")$Compound.Name)
  misses <- misses[!(misses %in% unique(assayinfo[,"Int Std"]))]
  
  miss.table <- NULL
  for (i in 1:length(misses))
  {
    this.miss <- misses[i]
    dates <- unique(subset(UC.data,Compound.Name==this.miss)$Date)
    for (this.date in dates)
    {
      this.row <- data.frame(Compound = this.miss, Date = this.date)
      miss.table <- rbind(miss.table,this.row)
    }
  }
  miss.table[order(miss.table$Date),]
  
  write.csv(miss.table[order(miss.table$Date),],file="chems-not-identified.csv",row.names=FALSE)
  
subset(UC.data, DTXSID=="DTXSID00379925")

 
  
  UC.data <- as.data.frame(UC.data)
  colnames(UC.data)[3:16] <- UC.data[6,3:16]
  
  
  UC.data <- subset(UC.data,!is.na(UC.data[,"Sample Text"]))
  dim(UC.data)
  
  
  # Extract the sample type from column Sample Text:
  UC.data[regexpr("AF",unlist(UC.data[,"Sample Text"]))!=-1,"Sample.Type"] <- "AF"
  UC.data[regexpr("UF",unlist(UC.data[,"Sample Text"]))!=-1,"Sample.Type"] <- "AF"
  UC.data[regexpr("T1",unlist(UC.data[,"Sample Text"]))!=-1,"Sample.Type"] <- "T1"
  UC.data[regexpr("T5",unlist(UC.data[,"Sample Text"]))!=-1,"Sample.Type"] <- "T5"
  UC.data[UC.data$Type == "Standard","Sample.Type"] <- "CC"
  UC.data[UC.data$Type == "Blank","Sample.Type"] <- "Blank"
  
  # Get rid of unused samples (QC samples):
  UC.data <- subset(UC.data,!is.na(Sample.Type))
  dim(UC.data)
 subset(UC.data, DTXSID=="DTXSID00379925"&Date=="2021-06-03") 
  # Extract the concentration of the standard from column Sample Text:
#  UC.data[,"Std.Conc"] <- unlist(lapply(
#    strsplit(gsub(" pg/uL","",unlist(UC.data[,"Sample Text"]))," - "),
#    function(x) ifelse(length(x)==2,as.numeric(x[2]),NA)))
  UC.data[,"Std.Conc"] <- UC.data[,"Std. Conc"]
  UC.data[,"Std.Units"] <- "uM" 
  
  # No replicates:
  UC.data$Series <- 1
#  # Extract the date from the time stamp
#  UC.data$Date <- unlist(lapply(strsplit(unlist(UC.data[,"Name"]),"_PFAS"),
#    function(x) x[1]))
  # Assume all samples analyzed on same date were the same calibration
  UC.data$Cal <- UC.data$Date
  # We'll need to go back and set this per sample type:
  UC.data[UC.data[,"Sample.Type"]=="AF","Dilution.Factor"] <- AF.DILUTE
  UC.data[UC.data[,"Sample.Type"]=="T1","Dilution.Factor"] <- T1.DILUTE
  UC.data[UC.data[,"Sample.Type"]=="T5","Dilution.Factor"] <- T5.DILUTE
  UC.data[UC.data[,"Sample.Type"]=="CC","Dilution.Factor"] <- CC.DILUTE
  UC.data[UC.data[,"Sample.Type"]=="Blank","Dilution.Factor"] <- BLANK.DILUTE
 subset(UC.data, DTXSID=="DTXSID00379925"&Date=="2021-06-03") 
    
  # Treat the blanks as calibration data with concentration 0:
  UC.data[UC.data[,"Sample.Type"]=="Blank","Std.Conc"] <- 0
  
  UC.data[UC.data[,"Sample.Type"]=="Blank","Sample.Type"] <- "CC"
  dim(UC.data)
  
  # Convert to uM:
  #UC.data[,"Std.Conc"] <- as.numeric(unlist(UC.data[,"Std. Conc (nM)"]))/1000
  # Get rid of CC samples that don't have a concentration:
  UC.data <- subset(UC.data,Sample.Type=!"CC" | !is.na(Std.Conc))
  
  dim(UC.data)
  
  
  
  # Should update this with input from Marci:
  UC.data$Analysis.Method <- "UPLC-MS/MS"
  UC.data$Analysis.Instrument <- "Waters Xevo TQ-S micro (QEB0036)"
  # Give Retention time:
  UC.data$Analysis.Parameters <- UC.data$RT
  
  # ISTD cocn 1ppm:
  UC.data$ISTD.Conc <- 1 #ppm
  
  # Test chemical concentration:
  UC.data$Test.Target.Conc <- 10 # uM
  
  
  # Set replicate id's:
  UC.data[,"Replicate"] <- ""
  UC.data[regexpr("_A",UC.data[,"Sample Text"])!=-1,"Replicate"] <- "A"
  UC.data[regexpr("_B",UC.data[,"Sample Text"])!=-1,"Replicate"] <- "B"
  UC.data[regexpr("_C",UC.data[,"Sample Text"])!=-1,"Replicate"] <- "C"
  
  
  # Make the numeric values numeric:
  UC.data[,"Area"] <- as.numeric(unlist(UC.data[,"Area"]))
  UC.data[,"IS Area"] <- as.numeric(unlist(UC.data[,"IS Area"]))
  
  write.table(UC.data,file="SmeltzPFAS/SmeltzPFAS-PPB-UC-Level0.tsv",sep="\t",row.names=FALSE)
  
  level1 <- format_fup_uc(subset(UC.data,DTXSID!="ISTD"),
    FILENAME="SmeltzPFAS/SmeltzPFAS",
    sample.col="Name",
    compound.col="Compound.Name",
    compound.conc.col="Std.Conc", 
    lab.compound.col="Compound.Name", 
    type.col="Sample.Type", 
    istd.col="IS Area",
    note.col="Replicate"
    )
   
  level2 <- level1
  # Taking all data in spreadsheet as human verfied for starters
  level2$Verified <- "Y"
  # From Marci:
  # 07/23/19:
  level2[level2$DTXSID %in% c(
    "DTXSID1032646",
    "DTXSID1067629",
    "DTXSID20375106",
    "DTXSID3037709",
    "DTXSID40380257",
    "DTXSID4059916",
    "DTXSID50375114",
    "DTXSID8031865") &
    level2$Date == "20190723", "Verified"] <- "Exclude"
  level2[level2$DTXSID %in% c(
    "DTXSID5030030",
    "DTXSID8037706") &
    level2$Date == "20190723" &
    level2$Note == "A", "Verified"] <- "Exclude"
  # 01/02/20:
  level2[level2$DTXSID %in% c(
    "DTXSID0060985",
    "DTXSID30170109",
    "DTXSID30382104",
    "DTXSID4059833",
    "DTXSID3031864") &
    level2$Date == "20200102", "Verified"] <- "Exclude"
  # 01/03/20:
  level2[level2$DTXSID %in% c(
    "DDTXSID504693204") &
    level2$Date == "20200103", "Verified"] <- "Exclude"
  level2[level2$DTXSID %in% c(
    "DTXSID3031864",
    "DTXSID3038939",
    "DTXSID30891564",
    "DTXSID5030030") &
    level2$Date == "20200103" &
    level2$Note == "A", "Verified"] <- "Exclude"
  # 01/06/20:
  level2[level2$DTXSID %in% c(
    "DTXSID1037303",
    "DTXSID3031860",
    "DTXSID3059921",
    "DTXSID90868151") &
    level2$Date == "20200106", "Verified"] <- "Exclude"
  level2[level2$DTXSID %in% c(
    "DTXSID8047553") &
    level2$Date == "20200106" &
    level2$Note == "A", "Verified"] <- "Exclude"  
  # 01/07/20:
  level2[level2$DTXSID %in% c(
    "DTXSID20375106",
    "DTXSID3037707",
    "DTXSID3037709",
    "DTXSID30382063",
    "DTXSID50375114") &
    level2$Date == "20200107", "Verified"] <- "Exclude"
  level2[level2$DTXSID %in% c(
    "DTXSID60500450") &
    level2$Date == "20200107" &
    level2$Note == "B", "Verified"] <- "Exclude"    
  # 01/08/20:
  level2[level2$DTXSID %in% c(
    "DDTXSID60380390",
    "DTXSID7027831",
    "DTXSID90315130") &
    level2$Date == "20200108", "Verified"] <- "Exclude"
  level2[level2$DTXSID %in% c(
    "DTXSID1032646",
    "DTXSID20179883",
    "DTXSID30627108") &
    level2$Date == "20200108" &
    level2$Note == "A", "Verified"] <- "Exclude"     
  # 11/23/20:
  level2[level2$DTXSID %in% c(
    "DTXSID4059833") &
    level2$Date == "20201123", "Verified"] <- "Exclude"
  # 11/24/20:
  level2[level2$DTXSID %in% c(
    "DTXSID0060985",
    "DTXSID20375106",
    "DTXSID30382104",
    "DTXSID90315130") &
    level2$Date == "20201124", "Verified"] <- "Exclude"
  # 02/25/21:
  level2[level2$DTXSID %in% c(
    "DTXSID8059926",
    "DTXSID8059928") &
    level2$Date == "20210225", "Verified"] <- "Exclude"
  level2[level2$DTXSID %in% c(
    "DTXSID40880025",
    "DTXSID50226894",
    "DTXSID50379814") &
    level2$Date == "20210225" &
    level2$Note == "A", "Verified"] <- "Exclude"   
  # 03/01/21:
  level2[level2$DTXSID %in% c(
    "DTXSID20861913") &
    level2$Date == "20210301", "Verified"] <- "Exclude"
  level2[level2$DTXSID %in% c(
    "DTXSID30382063") &
    level2$Date == "20210301" &
    level2$Note == "A", "Verified"] <- "Exclude"        
  level2[level2$DTXSID %in% c(
    "DTXSID00880026") &
    level2$Date == "20210301" &
    level2$Note == "B", "Verified"] <- "Exclude"  
  level2[level2$DTXSID %in% c(
    "DTXSID00880026",
    "DTXSID30382063") &
    level2$Date == "20210301" &
    level2$Note == "C", "Verified"] <- "Exclude"  
  # 03/08/21:
  level2[level2$DTXSID %in% c(
    "DTXSID5027140") &
    level2$Date == "20210308", "Verified"] <- "Exclude"
  # 03/11/21:
  level2[level2$DTXSID %in% c(
    "DTXSID00379925",
    "DTXSID40187142",
    "DTXSID70165670",
    "DTXSID70379917",
    "DTXSID80380837",
    "DTXSID80382154") &
    level2$Date == "20210311", "Verified"] <- "Exclude"
  level2[level2$DTXSID %in% c(
    "DTXSID70565479",
    "DTXSID9059915") &
    level2$Date == "20210311" &
    level2$Note == "A", "Verified"] <- "Exclude"   
  # 04/02/21:
  level2[level2$DTXSID %in% c(
    "DTXSID90315130") &
    level2$Date == "20210402", "Verified"] <- "Exclude"
  level2[level2$DTXSID %in% c(
    "DTXSID30382063",
    "DTXSID40880025",
    "DTXSID50226894",
    "DTXSID50379814") &
    level2$Date == "20210402" &
    level2$Note == "A", "Verified"] <- "Exclude"     
  # 06/03/21:
  level2[level2$DTXSID %in% c(
    "DTXSID00880026",
    "DTXSID70565479") &
    level2$Date == "20210603" &
    level2$Note == "C", "Verified"] <- "Exclude"      
        
  write.table(level2,
    file="SmeltzPFAS/SmeltzPFAS-PPB-UC-Level2.tsv",
    sep="\t",
    row.names=F,
    quote=F)
  
  level3 <- calc_fup_uc_point(FILENAME="SmeltzPFAS/SmeltzPFAS") 
  
  library(invitroTKstats)
  setwd("c:/users/jwambaug/git/invitroTKstats/working/")
  level4 <- calc_fup_uc(FILENAME="SmeltzPFAS/SmeltzPFAS") 
