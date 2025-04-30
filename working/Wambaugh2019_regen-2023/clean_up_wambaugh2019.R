# Clean up wambaugh2019
# There are several variables that are not relevant as an example dataset in the 
# invitroTKstats package and need to be removed. (Jira ticket IVTKS-49)

# The previous version of dataset is saved under data-raw as wambaugh2019_old.RData
# load("~/invitrotkstats/invitroTKstats/Data/wambaugh2019.RData")

# If there is a need to re-produce the steps, you will need to run the following code instead 
# to load the previous version of wambaugh2019 dataset.
# load("~/invitrotkstats/invitroTKstats/Data/wambaugh2019_old.RData")

library(dplyr)
# remove irrelevant variables and save the updated datasets
wambaugh2019.clint <- wambaugh2019.clint %>% select(-X, -Transition, -`ln...Remaining`, -TaskOrder)
wambaugh2019.red <- wambaugh2019.red %>% select(-Transition, -Task.Order, -CyprotexEPASetNumber)

save(wambaugh2019.clint, wambaugh2019.red, wambaugh2019.methods, 
     file = "~/invitrotkstats/invitroTKstats/Data/wambaugh2019.RData")
