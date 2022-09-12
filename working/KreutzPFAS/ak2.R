library(invitroTKstats)

# There are multiple packages for loading Excel files, but I've been using this
# one lately:
library(readxl)

# Change to the directory that has your Excel file 
# (this is where it is on my computer):
#setwd("c:/users/jwambaug/git/invitroTKstats/working/")

# load the data one sheet at a time from the Excel file, we skip the first row
# so that we get more of the column names loaded:
sheet1 <- read_excel(
  "021220_915_965_476_267_906_273_913_899_900_analysis_090920.xlsx",
  sheet=6,
  skip=1)

# 899-225 is on a different sheet:
sheet2<- read_excel(
  "021220_915_965_476_267_906_273_913_899_900_analysis_090920.xlsx",
  sheet=9,
  skip=1)

# 900 has no internal standard recorded
  
# Do each test chemical:
# 915 is given by "Area...14"
sheet1.915 <- extract_level1_fup_uc(
  sheet1,
  chem.name="915",
  area.col=14,
  ISTD.name="BET")
  
# 965 is given by "Area...22", we also need to change the ISTD name
sheet1.965 <- extract_level1_fup_uc(
  sheet1,
  chem.name="965",
  area.col=22,
  ISTD.name="HET")

# 476 is given by "Area...30", we also need to change the ISTD name
sheet1.476 <- extract_level1_fup_uc(
  sheet1,
  chem.name="476",
  area.col=30,
  ISTD.name="MFHET")

# 267 is given by "Area...38", we also need to change the ISTD name
sheet1.267 <- extract_level1_fup_uc(
  sheet1,
  chem.name="267",
  area.col=38,
  ISTD.name="DET")

# 3125 is given by "Area...46", we also need to change the ISTD name
sheet1.3125 <- extract_level1_fup_uc(
  sheet1,
  chem.name="3125",
  area.col=46,
  ISTD.name="DET")

# 906 is given by "Area...54", we also need to change the ISTD name
sheet1.906 <- extract_level1_fup_uc(
  sheet1,
  chem.name="906",
  area.col=54,
  ISTD.name="DET")

# 273 is given by "Area...62", we also need to change the ISTD name
sheet1.273 <- extract_level1_fup_uc(
  sheet1,
  chem.name="273",
  area.col=62,
  ISTD.name="DET")

# 913 is given by "Area...70", we also need to change the ISTD name
sheet1.913 <- extract_level1_fup_uc(
  sheet1,
  chem.name="913",
  area.col=70,
  ISTD.name="DET")

# 4NT is given by "Area...78", we also need to change the ISTD name
sheet1.4NT <- extract_level1_fup_uc(
  sheet1,
  chem.name="4NT",
  area.col=78,
  ISTD.name="4NT13C6")  


# 899 is given by "Area...78", we also need to change the ISTD name
sheet2.899 <- extract_level1_fup_uc(
  sheet2,
  chem.name="899",
  area.col=16,
  ISTD.offset=-4,
  ISTD.name="MFBET")  

AK.UCdata <- rbind(
  sheet1.915,
  sheet1.965,
  sheet1.476,
  sheet1.267,
  sheet1.3125,
  sheet1.906,
  sheet1.273,
  sheet1.913,
  sheet2.899,
  sheet1.4NT
  )

# No replicates:
AK.UCdata$Series <- 1
# Extract the date from the time stamp
AK.UCdata$Date <- unlist(lapply(strsplit(as.character(as.data.frame(AK.UCdata[,3])[,1])," "),function(x) x[1]))
# Assume all samples analyzed on same date were the same calibration
AK.UCdata$Cal <- AK.UCdata$Date
# I need to look up what the ID's for each lab sample were
AK.UCdata$DTXSID <- NA
# We'll need to go back and set this per sample type:
AK.UCdata$Dilution.Factor <- 1
# I don't think we need this:
AK.UCdata$Test.Target.Conc <- NA
# I'm assiming all ISTD were 1 PPM
AK.UCdata$ISTD.Conc <- 1

AK.UCdata$Note <- ""



