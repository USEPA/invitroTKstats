library(invitroTKstats)

level0 <- library(invitroTKstats)

level0 <- wambaugh2019.clint
level0$Date <- "2019"
level0$Sample.Type <- "Blank"
level0$Time..mins. <- as.numeric(level0$Time..mins.)
level0[!is.na(level0$Time..mins.),"Sample.Type"] <- "Cvst"
level0$ISTD.Name <- "Bucetin, Propranolol, and Diclofenac"
level0$ISTD.Conc <- 1
level0$Dilution.Factor <- 1
level0[is.na(level0$FileName),"FileName"]<-"Wambaugh2019"
level0$Hep.Density <- 0.5
level0$Analysis.Method <- "LCMS or GCMS"
level0$Analysis.Instrument <- "Agilent QQQ or Water.s.Xevo or AB Sciex Qtrap or Agilent GCMS or GCTOF"
level0$Analysis.Parameters <- "Unknown"
level0[is.na(level0$FileName),"FileName"] <- "Unknown"
level0[is.na(level0$TaskOrder),"TaskOrder"] <- "Unknown"
this.cal <- 1
this.row <- 1
while (this.row < dim(level0)[1])
{
  while (level0[this.row,"Preferred.Name"]=="Umbelliferone")
  {
    level0[this.row,"Calibration"] <- this.cal
    this.row <- this.row + 1
    if (is.na(level0[this.row,"Preferred.Name"])) break()
  }
  while (level0[this.row,"Preferred.Name"]!="Umbelliferone")
  {
    level0[this.row,"Calibration"] <- this.cal
    this.row <- this.row + 1
    if (is.na(level0[this.row,"Preferred.Name"])) break()
  }
  this.cal <- this.cal + 1
}

for (this.id in unique(level0$DTXSID))
  if (this.id %in% unique(wambaugh2019.methods$DTXSID))
  {
    if (!is.na(wambaugh2019.methods[wambaugh2019.methods$DTXSID==this.id,"LC"]))
       if (wambaugh2019.methods[wambaugh2019.methods$DTXSID==this.id,"LC"]=="Y")
    {
       level0[level0$DTXSID==this.id,"Analysis.Method"] <- "LCMS"
    }
    else if (!is.na(wambaugh2019.methods[wambaugh2019.methods$DTXSID==this.id,"GC"]))
       if (wambaugh2019.methods[wambaugh2019.methods$DTXSID==this.id,"GC"]=="Y")
    {
       level0[level0$DTXSID==this.id,"Analysis.Method"] <- "GCMS"
    }
    if (!is.na(wambaugh2019.methods[wambaugh2019.methods$DTXSID==this.id,"Agilent.QQQ"]))
       if (wambaugh2019.methods[wambaugh2019.methods$DTXSID==this.id,"Agilent.QQQ"]
      %in% c("Positive (Electrospray)", "Negative (Electrospray)"))
      { 
        level0[level0$DTXSID==this.id,"Analysis.Instrument"] <- "Agilent.QQQ"
        level0[level0$DTXSID==this.id,"Analysis.Parameters"] <-
          wambaugh2019.methods[wambaugh2019.methods$DTXSID==this.id,"Agilent.QQQ"] 
      }
    if (!is.na(wambaugh2019.methods[wambaugh2019.methods$DTXSID==this.id,"Water.s.Xevo"]))
       if (wambaugh2019.methods[wambaugh2019.methods$DTXSID==this.id,"Water.s.Xevo"]
      %in% c("Positive (Electrospray)", "Negative (Electrospray)"))
      { 
        level0[level0$DTXSID==this.id,"Analysis.Instrument"] <- "Water.s.Xevo"
        level0[level0$DTXSID==this.id,"Analysis.Parameters"] <-
          wambaugh2019.methods[wambaugh2019.methods$DTXSID==this.id,"Water.s.Xevo"] 
      }
     if (!is.na(wambaugh2019.methods[wambaugh2019.methods$DTXSID==this.id,"AB.Sciex.Qtrap"]))
       if (wambaugh2019.methods[wambaugh2019.methods$DTXSID==this.id,"AB.Sciex.Qtrap"]
       == "Passed") 
       level0[level0$DTXSID==this.id,"Analysis.Instrument"] <- "AB.Sciex.Qtrap"
     if (!is.na(wambaugh2019.methods[wambaugh2019.methods$DTXSID==this.id,"Agilent.GCMS"]))
       if (wambaugh2019.methods[wambaugh2019.methods$DTXSID==this.id,"Agilent.GCMS"]
       %in% c("Positive (El)","Y"))
       { 
         level0[level0$DTXSID==this.id,"Analysis.Instrument"] <- "Agilent.GCMS"
         level0[level0$DTXSID==this.id,"Analysis.Parameters"] <-
           wambaugh2019.methods[wambaugh2019.methods$DTXSID==this.id,"Agilent.GCMS"]
       }
     if (!is.na(wambaugh2019.methods[wambaugh2019.methods$DTXSID==this.id,"GCTOF"]))
       if (wambaugh2019.methods[wambaugh2019.methods$DTXSID==this.id,"GCTOF"]
       == "Ionizable") 
       level0[level0$DTXSID==this.id,"Analysis.Instrument"] <- "GCTOF"
}
  
level1 <- format_clint(level0,
  FILENAME="Wambaugh2019",
  sample.col="Sample.Name",
  compound.col="Preferred.Name",
  lab.compound.col="Name",
  time.col="Time..mins.",
  cal.col="Calibration"
)

   
level2 <- level1
level2$Verified <- "Y"
  
 write.table(subset(level2, DTXSID %in% unique(level2$DTXSID)[4:13]),
   file="Wambaugh2019-Clint-Level2.tsv",
   sep="\t",
   row.names=F,
   quote=F)
    
#level3 <- calc_clint_point(FILENAME="Wambaugh2019")
library(invitroTKstats)

level4 <- calc_clint(FILENAME="Wambaugh2019")

