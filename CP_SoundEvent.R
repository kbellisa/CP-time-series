### SCRIPT FOR CIVIL WAR EXPERIMENT
setwd("/Volumes/depot/bpijanow/data/CGS_Projects/CP_SoundEvent_2015")

source ("/Volumes/depot/bpijanow/data/CGS_Projects/CP_SoundEvent_2015/CPSound_kb.R")

install.packages('tuneR')
install.packages('seewave')
install.packages('soundecology')
install.packages("stringr")

library(tuneR)
library(seewave)
library(soundecology)
library(stringr)
library(parallel)
library(tools)
library(dplyr)

#IF FILE IS CORRUPT, CONTINUE PROCESSING
readWave.failsafe <- failwith(default=NA, readWave)

#LOAD FILES FOR EVENT
E <- "Event_9219" #complete
G <- "Event_15099" #complete
H <- "Event_15106" 
I <- "Event_15116"

#LOAD FILES FOR PRE EVENT COMMUNITY
J <- "Pre_comm9219"
K <- "Pre_comm_15099"
L <- "Pre_comm_15106"
M <- "Pre_comm_15116"

#LOAD FILES FOR POST EVENT COMMUNITY
A <- "AgPrairie_CivilWar"
B <- "BermCivilWar"
C <- "Riparian_CivilWar"
D <- "Prairie_CivilWar"

#replace E with actual event
filesE <- paste(H, list.files(H), sep="/")
df5 <- data.frame(filepath=filesE)

#READ WAVE FUNCTION IN A LOOP, RUN INDEX, CHECK FOR CORRUPT FILE, THEN OUTPUT
try.this <- mclapply(1:nrow(df5), function(k) {
  print(k)
  W <- readWave.failsafe(filesE[[k]], from=1, to=60, units="seconds")
  if(is.na(W) == TRUE) {
    indices <- rep(NA, 36)
  }
  else {
    indices <- CPSound(W)
  }
  rm(W)
  return(indices)
}, mc.cores=getOption("mc.cores", 12))


DF <- as.data.frame(matrix(as.numeric(unlist(try.this)), nrow=nrow(df5), byrow=TRUE))
names(DF) <- c("AcouOccupancyR", "Bioac", "Hf", "Ht", "H", "ACI", "AEI", "M", "NDSI", "ADI", "NPIC", "ASA", "BN", "SNR", "AA", "CAE", "ADAE", "RMS")
DF$wav <- substring(df5$filepath, 11)
names(DF) <- gsub(".", "", names(DF), fixed=TRUE)

write.csv(DF, file="/Volumes/depot/bpijanow/data/CGS_Projects/CP_SoundEvent_2015/15106_Event.csv")


