library(readxl)

# Change to the directory that has your Excel file 
# (this is where it is on my computer):
setwd("c:/users/jwambaug/git/invitroTKstats/working/")

PATH <- "CyprotexFup"

TO1ppb <- NULL
for (this.file in dir(PATH))
  if (this.file !="Problem")
  {
    if (length(excel_sheets(paste(PATH,"/",this.file,sep="")))==4)
    {
      for (this.conc in 1:3)
      {
        new.data <- read_excel(paste(PATH,"/",this.file,sep=""),skip=0,sheet=(1+this.conc))
        new.data <- subset(new.data,new.data[,2]!="")
        new.data[,1][is.na(new.data[,1])] <- ""
        if (any(new.data[1]=="Client ID"))
        {
          first.row <- which(new.data[,1]=="Client ID")
          colnames(new.data) <- new.data[first.row,]
          new.data <- new.data[(first.row+1):dim(new.data)[1],]
        }
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
        new.data$Protein.Conc <- c("10","30","100")[this.conc]
      }
    } else {
      new.data <- read_excel(paste(PATH,"/",this.file,sep=""),skip=0,sheet=2)
      new.data <- subset(new.data,new.data[,2]!="")
      new.data[,1][is.na(new.data[,1])] <- ""
      if (any(new.data[1]=="Client ID"))
      {
        first.row <- which(new.data[,1]=="Client ID")
        colnames(new.data) <- new.data[first.row,]
        new.data <- new.data[(first.row+1):dim(new.data)[1],]
      }
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
      new.data$Protein.Conc <- "100"
      new.data[regexpr("30%_Plasma",new.data$SampleName)!=-1,"Protein.Conc"] <- "30"
      new.data[regexpr("10%_Plasma",new.data$SampleName)!=-1,"Protein.Conc"] <- "10"
    }
    if (!is.null(TO1ppb)) new.data <- new.data[,colnames(TO1ppb)] 
    TO1ppb <- rbind(TO1ppb, new.data)
  }
 
TO1ppb[regexpr("-30",TO1ppb$CompoundName)!=-1,"Protein.Conc"] <- "30"
TO1ppb[regexpr("-10",TO1ppb$CompoundName)!=-1,"Protein.Conc"] <- "10"
TO1ppb$CompoundName <- gsub("-100","",TO1ppb$CompoundName)
TO1ppb$CompoundName <- gsub("-30","",TO1ppb$CompoundName)
TO1ppb$CompoundName <- gsub("-10","",TO1ppb$CompoundName)
TO1ppb$CompoundName <- gsub("-1","",TO1ppb$CompoundName)
TO1ppb$CompoundName <- gsub("-2","",TO1ppb$CompoundName)

 
length(unique(TO1ppb$CompoundName)) 

write.csv(TO1ppb,file="HTTK2TO1-PPB-all.txt")

  