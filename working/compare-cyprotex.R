all.chems <- sort(unique(c(
  TO1clint$CompoundName,
  TO1ppb$CompoundName,
  TO1caco2$CompoundName,
  TO1b2p$CompoundName)))

all.chems <- all.chems[regexpr("DTXSID",all.chems)!=-1]

all.data <- NULL
for (this.chem in all.chems)
{
  this.row <- data.frame(DTXSID=this.chem)
  
  clint.data <- subset(TO1clint,CompoundName==this.chem)
  this.row$Clint <- 0
  if (dim(clint.data)[1]>0) this.row$Clint <- 1
  this.row$Clint.1 <- 0
  if (dim(subset(clint.data,Test.Conc==1))[1]>0) this.row$Clint.1 <- 1
  this.row$Clint.10 <- 0
  if (dim(subset(clint.data,Test.Conc==10))[1]>0) this.row$Clint.10 <- 1
  this.row$Clint.HI <- 0
  if (dim(subset(clint.data,Heat.Control==1))[1]>0) this.row$Clint.HI <- 1

  fup.data <- subset(TO1ppb,CompoundName==this.chem)
  this.row$PPB <- 0
  if (dim(fup.data)[1]>0) this.row$PPB <- 1
  this.row$PPB.100 <- 0
  if (dim(subset(fup.data,Protein.Conc==100))[1]>0) this.row$PPB.100 <- 1
  this.row$PPB.30 <- 0
  if (dim(subset(fup.data,Protein.Conc==30))[1]>0) this.row$PPB.30 <- 1
  this.row$PPB.10 <- 0
  if (dim(subset(fup.data,Protein.Conc==10))[1]>0) this.row$PPB.10 <- 1

  caco2.data <- subset(TO1caco2,CompoundName==this.chem)
  this.row$Caco2 <- 0
  if (dim(caco2.data)[1]>0) this.row$Caco2 <- 1

  b2p.data <- subset(TO1b2p,CompoundName==this.chem)
  this.row$B2P <- 0
  if (dim(b2p.data)[1]>0) this.row$B2P <- 1
  
  all.data <- rbind(all.data,this.row)
}

write.csv(TO1ppb,file="HTTK2TO1-datastatus.txt")

