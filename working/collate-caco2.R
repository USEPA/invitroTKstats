library(readxl)

# Change to the directory that has your Excel file 
# (this is where it is on my computer):
setwd("c:/users/jwambaug/git/invitroTKstats/working/")

PATH <- "CyprotexCaco2"

TO1caco2 <- NULL
for (this.file in dir(PATH))
  if (this.file !="Problem")
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
#  colnames(new.data)[1] <- "Compound"
#  colnames(new.data)[10] <- "Test.Article"
#  colnames(new.data)[11] <- "Test.Conc"
  if (!is.null(TO1caco2)) new.data <- new.data[,colnames(TO1caco2)] 
  TO1caco2 <- rbind(TO1caco2, new.data)
}
 
length(unique(TO1caco2$CompoundName)) 

write.csv(TO1caco2,file="HTTK2TO1-Caco2-all.txt")
