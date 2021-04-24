library(readxl)

# Change to the directory that has your Excel file 
# (this is where it is on my computer):
setwd("c:/users/jwambaug/git/invitroTKstats/working/")

PATH <- "CyprotexClint"

TO1clint<- NULL
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
      if (this.sheet %in% c("10 uM","10uM","10uM_a","10uM_b",
         "Data 10uM","10 uM raw data","10uM Data","10 uM active",
         "6500 10uM Active","Xevo 10uM Active","5500 10uM Active",
         "Data 10uM control","10 uM control","6500 10uM Control",
         "10uM Control Group 2","Xevo 10uM Control","5500 10uM Control",
         "Xevo-1 10uM Active","Xevo-1 10uM Control",
         "Data - 10uM","Data - 10uM controls","Xevo 10 uM Active",
         "6500 10 uM Active","Data 10uM 5500","Data 10uM Xevo"))
      {
        good <- TRUE
        new.data$Test.Conc <- 10
        new.data$Heat.Control <- 0
      } else if (this.sheet %in% c("1 uM","1uM","1uM_a","1uM_b",
         "Data 1uM","1 uM raw data","1uM Data","1 uM active",
         "6500 1uM Active","Xevo 1uM Active","5500 1uM Active",
         "Data 1uM Control","1 uM control","6500 1uM Control",
         "1uM Control Group 2","Xevo 1uM Control","5500 1uM Control",
         "Xevo-1 1uM Active","Xevo-1 1uM Control",
         "Data - 1uM","Data - 1uM controls","Xevo 1 uM Active",
         "6500 1 uM Active","Data 1uM 5500","Data 1uM Xevo"))
      {
        good <- TRUE
        new.data$Test.Conc <- 1
        new.data$Heat.Control <- 0
      } else if (this.sheet %in% c("10uM_a Inactive","10uM_b Inactive",
        "10uM Inactive","10uM Data - Inactive","Xevo 10 uM Inactive",
        "6500 10 uM Inactive","10 uM HI"," 10 uM HI","10 uM Raw data"))
      {
        good <- TRUE
        new.data$Test.Conc <- 10
        new.data$Heat.Control <- 1
      } else if (this.sheet %in% c("1uM_a Inactive","1uM_b Inactive",
        "1uM Inactive","1uM Data - Inactive","Xevo 1 uM Inactive",
        "6500 1 uM Inactive","1 uM HI"," 1 uM HI","1 uM Raw data"))
      {
        good <- TRUE
        new.data$Test.Conc <- 1
        new.data$Heat.Control <- 1
      } else if (this.sheet %in% c("Data","Data2","Control",
        "DTXSID6025272 10uM Control","Data-Plate 2","Xevo1.PRO}"))
      {
        good <- TRUE
        new.data$Test.Conc <- NA
        new.data$Heat.Control <- NA
      } else if (this.sheet %in% c("Heat inactivated data","Control Inactive"))
      {
        good <- TRUE
        new.data$Test.Conc <- NA
        new.data$Heat.Control <- 1
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
        }
        colnames(new.data)[colnames(new.data)=="ISTD.Area"] <- "ISTD Area"
        colnames(new.data)[colnames(new.data)=="Filename"] <- "SampleName"
        colnames(new.data)[colnames(new.data)=="Sample Name"] <- "CompoundName"
        colnames(new.data)[colnames(new.data)=="Area Ratio"] <- "ISTDResponseRatio"
        colnames(new.data)[colnames(new.data)=="mass"] <- "Feature"
        colnames(new.data)[colnames(new.data)=="Transition"] <- "Feature"
        colnames(new.data)[colnames(new.data)=="Time (mins"] <- "Time"
        if (!("Feature" %in% colnames(new.data))) new.data$Feature <- ""
        if (!("Time" %in% colnames(new.data))) new.data$Time <- NA
        if (!("Test.Conc" %in% colnames(new.data))) new.data$Test.Conc <- NA
        if (!("Heat.Control" %in% colnames(new.data))) new.data$Heat.Control <- NA
        new.data <- new.data[,c(
          "SampleName",
          "CompoundName",
          "Feature",
          "Area",
          "ISTD Area",
          "ISTDResponseRatio",
          "Time",
          "Test.Conc",
          "Heat.Control")]
        new.data$TO <- 1
        new.data$FileName <- this.file
        if (!is.null(TO1clint)) new.data <- new.data[,colnames(TO1clint)] 
        TO1clint <- rbind(TO1clint, new.data)
      } else {
        print(paste(this.file,":",this.sheet))
      }
    }
  }

TO1clint$CompoundName <- gsub("_Human","",TO1clint$CompoundName)
TO1clint[regexpr("_1uM",TO1clint$CompoundName)!=-1,"Test.Conc"] <- 1
TO1clint[regexpr("_10uM",TO1clint$CompoundName)!=-1,"Test.Conc"] <- 10
TO1clint[regexpr("-1uM",TO1clint$CompoundName)!=-1,"Test.Conc"] <- 1
TO1clint[regexpr("-10uM",TO1clint$CompoundName)!=-1,"Test.Conc"] <- 10
TO1clint[regexpr("_1uM",TO1clint$SampleName)!=-1,"Test.Conc"] <- 1
TO1clint[regexpr("_10uM",TO1clint$SampleName)!=-1,"Test.Conc"] <- 10
TO1clint$CompoundName <- gsub("_1uM","",TO1clint$CompoundName)
TO1clint$CompoundName <- gsub("_10uM","",TO1clint$CompoundName)
TO1clint$CompoundName <- gsub("-1uM","",TO1clint$CompoundName)
TO1clint$CompoundName <- gsub("-10uM","",TO1clint$CompoundName)
TO1clint[regexpr("_HI",TO1clint$CompoundName)!=-1,"Heat.Control"] <- 1
TO1clint[regexpr("_HI",TO1clint$SampleName)!=-1,"Heat.Control"] <- 1
TO1clint$CompoundName <- gsub("_HI","",TO1clint$CompoundName)

TO1clint$Time <- -999
TO1clint[regexpr("_120",TO1clint$CompoundName)!=-1,"Time"] <- 120
TO1clint[regexpr("_120",TO1clint$SampleName)!=-1,"Time"] <- 120
TO1clint$CompoundName <- gsub("_120","",TO1clint$CompoundName)
TO1clint[regexpr("_60",TO1clint$CompoundName)!=-1,"Time"] <- 60
TO1clint[regexpr("_60",TO1clint$SampleName)!=-1,"Time"] <- 60
TO1clint$CompoundName <- gsub("_60","",TO1clint$CompoundName)
TO1clint[regexpr("_30",TO1clint$CompoundName)!=-1,"Time"] <- 30
TO1clint[regexpr("_30",TO1clint$SampleName)!=-1,"Time"] <- 30
TO1clint$CompoundName <- gsub("_30","",TO1clint$CompoundName)
TO1clint[regexpr("_15",TO1clint$CompoundName)!=-1,"Time"] <- 15
TO1clint[regexpr("_15",TO1clint$SampleName)!=-1,"Time"] <- 15
TO1clint$CompoundName <- gsub("_15","",TO1clint$CompoundName)
TO1clint[regexpr("_0",TO1clint$CompoundName)!=-1,"Time"] <- 0
TO1clint[regexpr("_0",TO1clint$SampleName)!=-1,"Time"] <- 0
TO1clint[regexpr("Blank",TO1clint$SampleName)!=-1,"Time"] <- NA
TO1clint$CompoundName <- gsub("_0","",TO1clint$CompoundName)
TO1clint$CompoundName <- gsub("_1","",TO1clint$CompoundName)
TO1clint$CompoundName <- gsub("_2","",TO1clint$CompoundName)
TO1clint$CompoundName <- gsub("_3","",TO1clint$CompoundName)
TO1clint[regexpr("_HI",TO1clint$SampleName)!=-1,"Time"] <- 120
TO1clint[TO1clint$CompoundName=="57.1","CompoundName"] <- "DTXSID5020605"   

TO1clint <- subset(TO1clint,!is.na(Area))
                           
 
length(unique(TO1clint$CompoundName)) 

write.csv(TO1clint,file="HTTK2TO1-Clint-all.txt")
