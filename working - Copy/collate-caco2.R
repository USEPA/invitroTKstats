library(readxl)

# Change to the directory that has your Excel file 
# (this is where it is on my computer):
setwd("c:/users/jwambaug/git/invitroTKstats/working/")

PATH <- "CyprotexCaco2"

TO1caco2 <- NULL
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
      if (this.sheet %in% c("Data","Data (2)","Raw data 1","Raw data 2",
        "Raw Data","Control Data"))
      {
        good <- TRUE
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
        if (!("Feature" %in% colnames(new.data))) new.data$Feature <- ""
        new.data <- new.data[,c(
          "SampleName",
          "CompoundName",
          "Feature",
          "Area",
          "ISTD Area",
          "ISTDResponseRatio")]
        new.data$TO <- 1
        new.data$FileName <- this.file
        new.data$SheetName <- this.sheet
        if (!is.null(TO1caco2)) new.data <- new.data[,colnames(TO1caco2)] 
        TO1caco2 <- rbind(TO1caco2, new.data)
      }  else {
        print("Skipped",paste(this.file,":",this.sheet))
      }
    }
  }
 
length(unique(TO1caco2$CompoundName)) 

write.table(TO1caco2,file="HTTK2TO1-Caco2-all.txt",row.names=F,sep="\t")
