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
        "6500 10 uM Inactive"))
      {
        good <- TRUE
        new.data$Test.Conc <- 10
        new.data$Heat.Control <- 1
      } else if (this.sheet %in% c("1uM_a Inactive","1uM_b Inactive",
        "1uM Inactive","1uM Data - Inactive","Xevo 1 uM Inactive",
        "6500 1 uM Inactive"))
      {
        good <- TRUE
        new.data$Test.Conc <- 1
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
        colnames(new.data)[colnames(new.data)=="mass"] <- "Feature"
        colnames(new.data)[colnames(new.data)=="Transition"] <- "Feature"
        new.data <- new.data[,c(
          "SampleName",
          "CompoundName",
          "Feature",
          "Area",
          "ISTD Area",
          "ISTDResponseRatio")]
        new.data$TO <- 1
        new.data$FileName <- this.file
        if (!is.null(TO1clint)) new.data <- new.data[,colnames(TO1clint)] 
        TO1caco2 <- rbind(TO1clint, new.data)
      } else {
        print(this.sheet)
      }
    }
  }
 
length(unique(TO1clint$CompoundName)) 

write.csv(TO1clint,file="HTTK2TO1-Clint-all.txt")
