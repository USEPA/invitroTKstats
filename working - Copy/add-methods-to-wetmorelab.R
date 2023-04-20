setwd("c:/users/jwambaug/git/invitroTKstats/working")

load("../invitroTKstats/data/smeltz2020.RData")
smeltz2020$Analysis.Method <- "GC" 
smeltz2020$Analysis.Instrument <- "Whatever the Wetmore lab GC/MS is"
smeltz2020$Analysis.Parameters <- "Retention time or the like"
save(smeltz2020,file="smeltz2020.RData")

load("../invitroTKstats/data/kreutz2020.RData")
kreutz2020$Analysis.Method <- "GC" 
kreutz2020$Analysis.Instrument <- "Whatever the Wetmore lab GC/MS is"
kreutz2020$Analysis.Parameters <- "Retention time or the like"
save(kreutz2020,file="kreutz2020.RData")
