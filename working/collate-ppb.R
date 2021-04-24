library(readxl)

# Change to the directory that has your Excel file 
# (this is where it is on my computer):
setwd("c:/users/jwambaug/git/invitroTKstats/working/")

PATH <- "CyprotexFup"

TO1ppb <- NULL
for (this.file in dir(PATH))
  if (this.file !="Problem")
  {
    sheets <- excel_sheets(paste(PATH,"/",this.file,sep=""))
    for (this.sheet in sheets)
    {
      new.data <- read_excel(paste(PATH,"/",this.file,sep=""),
        skip=0,
        sheet=which(sheets==this.sheet))
      good <- FALSE
      if (this.sheet %in% c("Data","Data-All Comps","Control","Data (2)",
        "Data-Raw","Control Data","Raw data Control","Control data"))
      {
        good <- TRUE
      } else if (this.sheet %in% c("DTXSID9047205-100","Data 100%",
        "Raw data 100","100% plasma"))
      {
        good <- TRUE
        new.data$Protein.Conc <- 100
      } else if (this.sheet %in% c("Raw data 30","Data 30%","30% plasma"))
      {
        good <- TRUE
        new.data$Protein.Conc <- 30
      } else if (this.sheet %in% c("Raw Data -10","Data 10%","10% plasma"))
      {
        good <- TRUE
        new.data$Protein.Conc <- 10
      }
      if (good)
      {
        new.data <- read_excel(paste(PATH,"/",this.file,sep=""),skip=0,sheet=2)
        new.data <- subset(new.data,new.data[,2]!="")
        new.data[,1][is.na(new.data[,1])] <- ""
        if (any(new.data[1]=="Client ID"))
        {
          first.row <- which(new.data[,1]=="Client ID")
          colnames(new.data) <- new.data[first.row,]
          new.data <- new.data[(first.row+1):dim(new.data)[1],]
        }
        colnames(new.data)[colnames(new.data)=="ISTD.Area"] <- "ISTD Area"
        colnames(new.data)[colnames(new.data)=="Filename"] <- "SampleName"
        colnames(new.data)[colnames(new.data)=="Sample Name"] <- "CompoundName"
        colnames(new.data)[colnames(new.data)=="Area Ratio"] <- "ISTDResponseRatio"
        colnames(new.data)[colnames(new.data)=="mass"] <- "Feature"
        colnames(new.data)[colnames(new.data)=="Transition"] <- "Feature"
        if (!("Feature" %in% colnames(new.data))) new.data$Feature <- ""
        if (!("Protein.Conc" %in% colnames(new.data))) new.data$Protein.Conc <- NA
        new.data <- new.data[,c(
          "SampleName",
          "CompoundName",
          "Feature",
          "Area",
          "ISTD Area",
          "ISTDResponseRatio")]
        new.data$TO <- 1
        new.data$FileName <- this.file
        new.data$Protein.Conc <- 100
        new.data[regexpr("30%_Plasma",new.data$SampleName)!=-1,"Protein.Conc"] <- 30
        new.data[regexpr("10%_Plasma",new.data$SampleName)!=-1,"Protein.Conc"] <- 10
      if (!is.null(TO1ppb)) new.data <- new.data[,colnames(TO1ppb)] 
      TO1ppb <- rbind(TO1ppb, new.data)
      }  else {
        print(paste(this.file,":",this.sheet))
      }
    }
  }
 
TO1ppb[regexpr("-100",TO1ppb$CompoundName)!=-1,"Protein.Conc"] <- 100
TO1ppb[regexpr("-30",TO1ppb$CompoundName)!=-1,"Protein.Conc"] <- 30
TO1ppb[regexpr("-10",TO1ppb$CompoundName)!=-1,"Protein.Conc"] <- 10
TO1ppb[regexpr("-100",TO1ppb$SampleName)!=-1,"Protein.Conc"] <- 100
TO1ppb[regexpr("-30",TO1ppb$SampleName)!=-1,"Protein.Conc"] <- 30
TO1ppb[regexpr("-10",TO1ppb$SampleName)!=-1,"Protein.Conc"] <- 10
TO1ppb$CompoundName <- gsub("-100","",TO1ppb$CompoundName)
TO1ppb$CompoundName <- gsub("-30","",TO1ppb$CompoundName)
TO1ppb$CompoundName <- gsub("-10","",TO1ppb$CompoundName)
TO1ppb$CompoundName <- gsub("-1","",TO1ppb$CompoundName)
TO1ppb$CompoundName <- gsub("-2","",TO1ppb$CompoundName)

 
length(unique(TO1ppb$CompoundName)) 

write.csv(TO1ppb,file="HTTK2TO1-PPB-all.txt")

  