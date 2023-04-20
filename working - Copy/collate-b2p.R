library(readxl)

# Change to the directory that has your Excel file 
# (this is where it is on my computer):
setwd("c:/users/jwambaug/git/invitroTKstats/working/")

PATH <- "CyprotexB2P"

TO1b2p<- NULL
for (this.file in dir(PATH))
  if (this.file !="Problem")
  {
    sheets <- excel_sheets(paste(PATH,"/",this.file,sep=""))
    for (this.sheet in sheets)
    {
      new.data <- suppressMessages(read_excel(paste(PATH,"/",this.file,sep=""),
        skip=2,
        sheet=which(sheets==this.sheet)))
      if (colnames(new.data)[1]!="SampleName")
      {
        if (regexpr("Summary",this.sheet)==-1) print(
          paste(this.file,":",this.sheet,
            "first column name is",
            colnames(new.data)[1],
            "changing skip to 0"))
        suppressMessages(new.data <- read_excel(paste(PATH,"/",this.file,sep=""),
          skip=0,
          sheet=which(sheets==this.sheet)))
      }
      good <- FALSE
      if (regexpr("10uM",this.sheet)!=-1 |
        regexpr("10 uM",this.sheet)!=-1)
      {
        good <- TRUE
        new.data$Test.Conc <- 10
      } else if (regexpr("5uM",this.sheet)!=-1 |
        regexpr("5 uM",this.sheet)!=-1)
      {
        good <- TRUE
        new.data$Test.Conc <- 5
      } else if (this.sheet %in% c("Raw Data","Data", "Raw data") |
        regexpr("Methazolamide",this.sheet)!=-1 |
        regexpr("Verapamil",this.sheet)!=-1)
      {
        good <- TRUE
        new.data$Test.Conc <-NaN
      } 
      if (good)
      {
        colnames(new.data)[colnames(new.data)=="ISTD.Area"] <- "ISTD Area"
        colnames(new.data)[colnames(new.data)=="Filename"] <- "SampleName"
        colnames(new.data)[colnames(new.data)=="Sample Name"] <- "CompoundName"
        colnames(new.data)[colnames(new.data)=="Area Ratio"] <- "ISTDResponseRatio"
        colnames(new.data)[colnames(new.data)=="mass"] <- "Feature"
        colnames(new.data)[colnames(new.data)=="Transition"] <- "Feature"
        new.data <- subset(new.data,!is.na(SampleName) & !is.na(Area))
        new.data <- subset(new.data,SampleName!="SampleName")
        new.data <- subset(new.data,SampleName!="0")
        if (!("Feature" %in% colnames(new.data))) new.data$Feature <- ""
        if (!("Test.Conc" %in% colnames(new.data))) new.data$Test.Conc <- NA
        if (!("CompoundName" %in% colnames(new.data)))
        {
          sample.name <- new.data$SampleName[
            regexpr("DTXSID",new.data$SampleName)!=-1]
          sample.name <- sample.name[1]
          new.data$CompoundName <- strsplit(sample.name,"_")[[1]][1]
        }
        if (any(regexpr("BLANK",new.data$CompoundName)!=-1))
        {
          this.id <- new.data$CompoundName[regexpr("BLANK",new.data$CompoundName)==-1][1]
          new.data[regexpr("BLANK",new.data$CompoundName)!=-1,"CompoundName"] <-
            this.id
        }
        if (any(new.data$CompoundName=="CompoundName")) browser()
        new.data <- new.data[,c(
          "SampleName",
          "CompoundName",
          "Feature",
          "Area",
          "ISTD Area",
          "ISTDResponseRatio",
          "Test.Conc")]
        new.data$TO <- 1
        new.data$FileName <- this.file
        new.data$SheetName <- this.sheet
        if (!is.null(TO1b2p)) new.data <- new.data[,colnames(TO1b2p)] 
        TO1b2p <- rbind(TO1b2p, new.data)
      } else {
        print(paste("Skipped",this.file,":",this.sheet))
      }
    }
  }

TO1b2p[regexpr("_5uM",TO1b2p$CompoundName)!=-1,"Test.Conc"] <- 5
TO1b2p[regexpr("_10uM",TO1b2p$CompoundName)!=-1,"Test.Conc"] <- 10
TO1b2p$CompoundName <- gsub("_5uM","",TO1b2p$CompoundName)
TO1b2p$CompoundName <- gsub("_10uM","",TO1b2p$CompoundName)
TO1b2p$CompoundName <- gsub("_amm","",TO1b2p$CompoundName)
TO1b2p$CompoundName <- gsub("_Ref_Plasma","",TO1b2p$CompoundName)
TO1b2p$CompoundName <- gsub("_Plasma","",TO1b2p$CompoundName)
TO1b2p$CompoundName <- gsub("_1","",TO1b2p$CompoundName)
TO1b2p$CompoundName <- gsub("_2","",TO1b2p$CompoundName)
TO1b2p$CompoundName <- gsub("_3","",TO1b2p$CompoundName)
TO1b2p$CompoundName <- gsub("_4","",TO1b2p$CompoundName)
TO1b2p$CompoundName <- gsub("_","",TO1b2p$CompoundName)

TO1b2p[regexpr("Ref_Plasma",TO1b2p$SampleName)!=-1,"Type"] <- "Blood"
TO1b2p[regexpr("Ref Plasma",TO1b2p$SampleName)!=-1,"Type"] <- "Blood"
TO1b2p[regexpr("RefP",TO1b2p$SampleName)!=-1,"Type"] <- "Blood"
TO1b2p[regexpr("blank",tolower(TO1b2p$SampleName))!=-1,"Type"] <- "Blank"
TO1b2p[regexpr("Ref_Plasma",TO1b2p$SampleName)==-1 &
  regexpr("Ref Plasma",TO1b2p$SampleName)==-1 &
  regexpr("Plasma",TO1b2p$SampleName)!=-1,"Type"] <- "Plasma"


 
length(unique(TO1b2p$CompoundName)) 

write.table(TO1b2p,file="HTTK2TO1-b2p-all.txt",row.names=F,sep="\t")
