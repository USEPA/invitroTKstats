smeltz2023.red <- read.csv("SmeltzPFAS/SmeltzPFAS-fup-RED-Level2.tsv",sep="\t")
smeltz2023.uc <- read.csv("SmeltzPFAS/SmeltzPFAS-PPB-UC-Level2.tsv",sep="\t")
smeltz2023.clint <- read.csv("SmeltzPFAS/SmeltzPFAS-Clint-Level2.tsv",sep="\t")
kreutz2023.uc <- read.csv("KreutzPFAS/KreutzPFAS-fup-UC-Level2.tsv",sep="\t")
kreutz2023.clint <- read.csv("KreutzPFAS/KreutzPFAS-Clint-Level2.tsv",sep="\t")

save(smeltz2023.red,smeltz2023.uc,smeltz2023.clint,
     file="Smeltz2023.RData",
     version=2)

save(kreutz2023.uc,kreutz2023.clint,
     file="Kreutz2023.RData",
     version=2)
