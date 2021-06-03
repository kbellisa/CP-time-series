#2015 Conner Prairie Time-Series
#by Jack VanSchaik and Kristen Bellisario

#install.packages("lmtest")
#install.packages("fpp")

library(lmtest)
library(fpp) #stationarity, Augmented-Dickey-Fuller Test

source ("/Volumes/depot/bpijanow/data/CGS_Projects/CP_SoundEvent/Scripts/cbindNA.R")
source ("/Volumes/depot/bpijanow/data/CGS_Projects/CP_SoundEvent/Scripts/rbindNA.R")

setwd("/Volumes/depot/bpijanow/data/CGS_Projects/CP_SoundEvent/Results/2015_Results")

###############################
#    2015 COMPARISON          #
###############################

##############################
#   PRE COMMUNITY FUNCTION   #
##############################
PrairiePre <- read.csv("/Volumes/depot/bpijanow/data/CGS_Projects/CP_SoundEvent/Results/2015_Results/Indices/Pre_15116.csv", header=TRUE)
RiparianPre <- read.csv("/Volumes/depot/bpijanow/data/CGS_Projects/CP_SoundEvent/Results/2015_Results/Indices/Pre_9219.csv", header=TRUE)
AgPre <- read.csv("/Volumes/depot/bpijanow/data/CGS_Projects/CP_SoundEvent/2015_Results/Results/Indices/Pre_15106.csv", header=TRUE)

setwd("/Volumes/depot/bpijanow/data/CGS_Projects/CP_SoundEvent/Results/2015_Results/Pre")

#Test all pre-event indices / plot
PreTest <- sapply(11:ncol(PrairiePre), function(x) {
  n <- names(PrairiePre)
  tsvar <- ts(PrairiePre[,x], frequency=24)
  png(paste(n[x], "PrairiePreTSvar.png"))
  plot(tsvar)
  dev.off()
  png(paste(n[x], "PrairiePreDecompose.png"))
  plot(decompose(tsvar))
  dev.off()
  png(paste(n[x], "PrairiePreACF.png"))
  plot(acf(tsvar, lag.max=72)) 
  dev.off()
  adftest <- adf.test(tsvar, alternative="stationary")
  return(adftest)
})

write.csv(PreTest, file="/Volumes/depot/bpijanow/data/CGS_Projects/CP_SoundEvent/Results/2015_Results/ADFtests/PreTestPrairieADF.csv")

PreTest2 <- sapply(11:ncol(RiparianPre), function(x) {
  n <- names(RiparianPre)
  tsvar <- ts(RiparianPre[,x], frequency=24)
  png(paste(n[x], "RiparianPreTSvar.png"))
  plot(tsvar)
  dev.off()
  png(paste(n[x], "RiparianPreDecompose.png"))
  plot(decompose(tsvar))
  dev.off()
  png(paste(n[x], "RiparianPreACF.png"))
  plot(acf(tsvar, lag.max=72)) 
  dev.off()
  adftest2 <- adf.test(tsvar, alternative="stationary")
  return(adftest2)
})

write.csv(PreTest2, file="/Volumes/depot/bpijanow/data/CGS_Projects/CP_SoundEvent/Results/2015_Results/ADFtests/PreTestRiparianADF.csv")

PreTest4 <- sapply(11:ncol(AgPre), function(x) {
  n <- names(AgPre)
  tsvar <- ts(AgPre[,x], frequency=24)
  png(paste(n[x], "AgPreTSvar.png"))
  plot(tsvar)
  dev.off()
  png(paste(n[x],"AgPreDecompose.png"))
  plot(decompose(tsvar))
  dev.off()
  png(paste(n[x],"AgPreACF.png"))
  plot(acf(tsvar, lag.max=72)) 
  dev.off()
  adftest4 <- adf.test(tsvar, alternative="stationary")
  return(adftest4)
})
write.csv(PreTest4, file="/Volumes/depot/bpijanow/data/CGS_Projects/CP_SoundEvent/Results/2015_Results/ADFtests/PreTestAgADF.csv")

##############################
# EVENT COMMUNITY FUNCTION   #
##############################

#TIME SERIES: EVENT
setwd("/Volumes/depot/bpijanow/data/CGS_Projects/CP_SoundEvent/Results/2015_Results/Event")

event9219 <- read.csv("/Volumes/depot/bpijanow/data/CGS_Projects/CP_SoundEvent/Results/2015_Results/Indices/Event_9219.csv", header=TRUE)
event15106 <- read.csv("/Volumes/depot/bpijanow/data/CGS_Projects/CP_SoundEvent/Results/2015_Results/Indices/Event_15106.csv", header=TRUE)
event15116 <- read.csv("/Volumes/depot/bpijanow/data/CGS_Projects/CP_SoundEvent/Results/2015_Results/Indices/Event_15116.csv", header=TRUE)

# Riparian Event:  Test all indices
EventTest9219 <- sapply(2:ncol(event9219), function(x) {
  n <- names(event9219)
  tsvar <- ts(event9219[,x], frequency=24)
  png(paste(n[x], "RipEventTSvar.png"))
  plot(tsvar)
  dev.off()
  png(paste(n[x],"RipEventDecompose.png"))
  plot(decompose(tsvar))
  dev.off()
  png(paste(n[x],"RipEventACF.png"))
  plot(acf(tsvar, lag.max=72)) 
  dev.off()
  adftestE9219 <- adf.test(tsvar, alternative="stationary")
  return(adftestE9219)
})

write.csv(EventTest9219, file="/Volumes/depot/bpijanow/data/CGS_Projects/CP_SoundEvent/Results/2015_Results/Results/ADFtests/EventTest9219ADF.csv")

EventTest15116<- sapply(2:ncol(event15116), function(x) {
  n <- names(event15116)
  tsvar <- ts(event15116[,x], frequency=24)
  png(paste(n[x], "AgEventTSvar.png"))
  plot(tsvar)
  dev.off()
  png(paste(n[x],"AgEventDecompose.png"))
  plot(decompose(tsvar))
  dev.off()
  png(paste(n[x],"AgEventACF.png"))
  plot(acf(tsvar, lag.max=72)) 
  dev.off()
  adftestE15116 <- adf.test(tsvar, alternative="stationary")
  return(adftestE15116)
})

write.csv(EventTest15116, file="/Volumes/depot/bpijanow/data/CGS_Projects/CP_SoundEvent/Results/2015_Results/ADFtests/EventTest15116ADF.csv")

################################
#   POST COMMUNITY FUNCTION    #
################################

setwd("/Volumes/depot/bpijanow/data/CGS_Projects/CP_SoundEvent/Results/2015_Results/Post")

Prairie <- read.csv("/Volumes/depot/bpijanow/data/CGS_Projects/CP_SoundEvent/Results/2015_Results/Indices/Post_Prairie.csv", header=TRUE)
Riparian <- read.csv("/Volumes/depot/bpijanow/data/CGS_Projects/CP_SoundEvent/Results/2015_Results/Indices/Post_Riparian.csv", header=TRUE)
Ag <- read.csv("/Volumes/depot/bpijanow/data/CGS_Projects/CP_SoundEvent/Results/2015_Results/Indices/Post_AgPrairie.csv", header=TRUE)

PostTestPrairie<- sapply(2:ncol(Prairie), function(x) {
  n <- names(Prairie)
  tsvar <- ts(Prairie[,x], frequency=24)
  png(paste(n[x], "PrairiePostTSvar.png"))
  plot(tsvar)
  dev.off()
  png(paste(n[x],"PrairiePostDecompose.png"))
  plot(decompose(tsvar))
  dev.off()
  png(paste(n[x],"PrairiePostACF.png"))
  plot(acf(tsvar, lag.max=72)) 
  dev.off()
  adftestPrPo <- adf.test(tsvar, alternative="stationary")
  return(adftestPrPo)
})

write.csv(PostTestPrairie, file="/Volumes/depot/bpijanow/data/CGS_Projects/CP_SoundEvent/Results/2015_Results/ADFtests/PostTestPrairieADF.csv")

PostTestRiparian<- sapply(2:ncol(Riparian), function(x) {
  n <- names(Riparian)
  tsvar <- ts(Riparian[,x], frequency=24)
  png(paste(n[x], "RiparianPostTSvar.png"))
  plot(tsvar)
  dev.off()
  png(paste(n[x],"RiparianPostDecompose.png"))
  plot(decompose(tsvar))
  dev.off()
  png(paste(n[x],"RiparianPostACF.png"))
  plot(acf(tsvar, lag.max=72)) 
  dev.off()
  adftestRipPo <- adf.test(tsvar, alternative="stationary")
  return(adftestRipPo)
})

write.csv(PostTestRiparian, file="/Volumes/depot/bpijanow/data/CGS_Projects/CP_SoundEvent/Results/2015_Results/ADFtests/PostTestRiparianADF.csv")


PostTestAg<- sapply(2:ncol(Ag), function(x) {
  n <- names(Ag)
  tsvar <- ts(Ag[,x], frequency=24)
  png(paste(n[x], "AgPostTSvar.png"))
  plot(tsvar)
  dev.off()
  png(paste(n[x],"AgPostDecompose.png"))
  plot(decompose(tsvar))
  dev.off()
  png(paste(n[x],"AgPostACF.png"))
  plot(acf(tsvar, lag.max=72)) 
  dev.off()
  adftestAgPo <- adf.test(tsvar, alternative="stationary")
  return(adftestAgPo)
})

write.csv(PostTestAg, file="/Volumes/depot/bpijanow/data/CGS_Projects/CP_SoundEvent/Results/2015_Results/ADFtests/PostTestAgADF.csv")

############################
#   COMBINED TIME SERIES   #
#     PRE, EVENT, POST     #
############################

setwd("/Volumes/depot/bpijanow/data/CGS_Projects/CP_SoundEvent/Results/2015_Results/FullTimeSeriesPlots")

#PRAIRIE COMBINED

#adjust columns for meterological data (check - applicable for both 2015 and 2016?)
PrairiePre2 <- PrairiePre[,10:21]

PrairieTS.comb<- sapply(2:ncol(PrairiePre), function(x) {
  comb.prairie <- (cbind.na(PrairiePre[,x], event15106[,x], Prairie[,x]))
  comb.prairie.3 <- stack(as.data.frame(comb.prairie))
  comb.prairie.3 <- na.omit(comb.prairie.3)
  n <- names(PrairiePre[x])
  tsvar <- ts(comb.prairie.3[,1], frequency=24)
  png(paste(n, "CombPrairieTSvar.png"))
  plot(tsvar)
  dev.off()
  png(paste(n,"CombPrairieDecompose.png"))
  plot(decompose(tsvar))
  dev.off()
  png(paste(n,"CombPrairieACF.png"))
  plot(acf(tsvar, lag.max=240)) 
  dev.off()
  adftestCombPrairie <- adf.test(tsvar, alternative="stationary")
  return(adftestCombPrairie)
})

write.csv(PrairieTS.comb, file="/Volumes/depot/bpijanow/data/CGS_Projects/CP_SoundEvent/Results/2015_Results/ADFtests/PrairieCombinedADF.csv")

#adjust columns for meterological data
RiparianPre2 <- RiparianPre[,10:21]

RiparianTS.comb<- sapply(2:ncol(RiparianPre), function(x) {
  comb.riparian <- (cbind.na(RiparianPre[,x], event9219[,x], Riparian[,x]))
  comb.riparian.3 <- stack(as.data.frame(comb.riparian))
  comb.riparian.3 <- na.omit(comb.riparian.3)
  n <- names(RiparianPre[x])
  tsvar <- ts(comb.riparian.3[,1], frequency=24)
  png(paste(n, "CombRiparianTSvar.png"))
  plot(tsvar)
  dev.off()
  png(paste(n,"CombRiparianDecompose.png"))
  plot(decompose(tsvar))
  dev.off()
  png(paste(n,"CombRiparianACF.png"))
  plot(acf(tsvar, lag.max=240)) 
  dev.off()
  adftestCombRiparian <- adf.test(tsvar, alternative="stationary")
  return(adftestCombRiparian)
})

write.csv(RiparianTS.comb, file="/Volumes/depot/bpijanow/data/CGS_Projects/CP_SoundEvent/Results/2015_Results/ADFtests/RiparianCombinedADF.csv")

#adjust columns for meterological data
AgPre2 <- AgPre[,10:21]

AgTS.comb<- sapply(2:ncol(AgPre), function(x) {
  comb.Ag <- (cbind.na(AgPre[,x], event15116[,x], Ag[,x]))
  comb.Ag.3 <- stack(as.data.frame(comb.Ag))
  comb.Ag.3 <- na.omit(comb.Ag.3)
  n <- names(AgPre[x])
  tsvar <- ts(comb.Ag.3[,1], frequency=24)
  png(paste(n, "CombAgTSvar.png"))
  plot(tsvar)
  dev.off()
  png(paste(n,"CombAgDecompose.png"))
  plot(decompose(tsvar))
  dev.off()
  png(paste(n,"CombAgACF.png"))
  plot(acf(tsvar, lag.max=240)) 
  dev.off()
  adftestCombAg <- adf.test(tsvar, alternative="stationary")
  return(adftestCombAg)
})

write.csv(AgTS.comb, file="/Volumes/depot/bpijanow/data/CGS_Projects/CP_SoundEvent/Results/2015_Results/ADFtests/AgCombinedADF.csv")

