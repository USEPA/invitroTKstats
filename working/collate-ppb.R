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
      new.data <- suppressMessages(read_excel(paste(PATH,"/",this.file,sep=""),
        skip=0,
        sheet=which(sheets==this.sheet)))
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
        new.data <- subset(new.data,new.data[,2]!="")
        new.data[,1][is.na(new.data[,1])] <- ""
        if (any(new.data[1]=="Client ID"))
        {
          first.row <- which(new.data[,1]=="Client ID")
          colnames(new.data) <- new.data[first.row,]
          new.data <- new.data[(first.row+1):dim(new.data)[1],]
        } else if (any(new.data[1]=="Filename"))
        {
          first.row <- which(new.data[,1]=="Filename")
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
        if (!("Protein.Conc" %in% colnames(new.data))) new.data$Protein.Conc <- NaN
        if (any(regexpr("BLANK",toupper(new.data$CompoundName)!=-1)))
        {
          this.id <- new.data$CompoundName[regexpr("BLANK",
            toupper(new.data$CompoundName))==-1][1]
          new.data[regexpr("BLANK",toupper(new.data$CompoundName))!=-1,
            "CompoundName"] <- this.id
        }
        new.data <- new.data[,c(
          "SampleName",
          "CompoundName",
          "Feature",
          "Area",
          "ISTD Area",
          "ISTDResponseRatio",
          "Protein.Conc")]
        new.data[regexpr("100%_Plasma",new.data$SampleName)!=-1,"Protein.Conc"] <- 100
        new.data[regexpr("30%_Plasma",new.data$SampleName)!=-1,"Protein.Conc"] <- 30
        new.data[regexpr("10%_Plasma",new.data$SampleName)!=-1,"Protein.Conc"] <- 10
        # Sometimes no protein concentration is given in an experiment (yes there's no
        # protein for the specific measurment, but we need to know which experiment it
        # corresponds to:)
        if (any(is.nan(new.data$Protein.Conc)))
        {
          this.row <- 1
          this.conc <- NaN
          while (this.row < dim(new.data)[1])
          {
            if (!is.nan(unlist(new.data[this.row,"Protein.Conc"])))
            {
              this.conc <- new.data[this.row,"Protein.Conc"]
            } else {
              new.data[this.row,"Protein.Conc"] <- this.conc 
            }
            this.row <- this.row + 1
          }
        }
        new.data$TO <- 1
        new.data$FileName <- this.file
        new.data$SheetName <- this.sheet

        if (!is.null(TO1ppb)) new.data <- new.data[,colnames(TO1ppb)] 
        TO1ppb <- rbind(TO1ppb, new.data)
      }  else {
        print("Skipped",paste(this.file,":",this.sheet))
      }
    }
  }
 
TO1ppb[regexpr("-10",TO1ppb$CompoundName)!=-1,"Protein.Conc"] <- 10
TO1ppb[regexpr("-100",TO1ppb$CompoundName)!=-1,"Protein.Conc"] <- 100
TO1ppb[regexpr("-30",TO1ppb$CompoundName)!=-1,"Protein.Conc"] <- 30
TO1ppb[regexpr("-10",TO1ppb$SampleName)!=-1,"Protein.Conc"] <- 10
TO1ppb[regexpr("-100",TO1ppb$SampleName)!=-1,"Protein.Conc"] <- 100
TO1ppb[regexpr("-30",TO1ppb$SampleName)!=-1,"Protein.Conc"] <- 30
TO1ppb$CompoundName <- gsub("-100","",TO1ppb$CompoundName)
TO1ppb$CompoundName <- gsub("-30","",TO1ppb$CompoundName)
TO1ppb$CompoundName <- gsub("-10","",TO1ppb$CompoundName)
TO1ppb$CompoundName <- gsub("-1","",TO1ppb$CompoundName)
TO1ppb$CompoundName <- gsub("-2","",TO1ppb$CompoundName)

TO1ppb[regexpr("_Plasma",TO1ppb$CompoundName)!=-1,"Type"] <- "Plasma"
TO1ppb[regexpr("_TO",TO1ppb$CompoundName)!=-1,"Type"] <- "T0"
TO1ppb[regexpr("_PBS",TO1ppb$CompoundName)!=-1,"Type"] <- "PBS"
TO1ppb[regexpr("Plasma",TO1ppb$SampleName)!=-1,"Type"] <- "Plasma"
TO1ppb[regexpr("TO",TO1ppb$SampleName)!=-1,"Type"] <- "T0"
TO1ppb[regexpr("PBS",TO1ppb$SampleName)!=-1,"Type"] <- "PBS"
TO1ppb[regexpr("Blank",TO1ppb$SampleName)!=-1,"Type"] <- "PBS"
TO1ppb$CompoundName <- gsub("_Plasma","",TO1ppb$CompoundName)
TO1ppb$CompoundName <- gsub("_T0","",TO1ppb$CompoundName)
TO1ppb$CompoundName <- gsub("_PBS","",TO1ppb$CompoundName)
TO1ppb$CompoundName <- gsub("_T4","",TO1ppb$CompoundName)
TO1ppb$CompoundName <- gsub("_1","",TO1ppb$CompoundName)
TO1ppb$CompoundName <- gsub("_2","",TO1ppb$CompoundName)
TO1ppb$CompoundName <- gsub("_","",TO1ppb$CompoundName)





 
length(unique(TO1ppb$CompoundName)) 

write.csv(TO1ppb,file="HTTK2TO1-PPB-all.txt")

  