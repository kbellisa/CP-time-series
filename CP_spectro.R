#July 3, 2017 - Kristen Bellisario
#1 Min Spectrograms

library(seewave)
library(parallel)

### GENERATE SPECTROGRAMS - 10 s.
#Pre X Event Post
setwd("/Volumes/bpijanow/data/CGS_Projects/CP_SoundEvent_2015/2015/Pre/")

#15166_PreEvent X/ 9219_PreEventNewX / 15099_PreEvent X / 15106_PreEvent (uh oh)
#Event_15116 X / Event_9219 X / Event_15106 X
#Prairie_CivilWar / Riparian_CivilWar / AgPrairie_CivilWar

batch1 <- "15106_PreEvent"
filesBatch1 <- paste(batch1, list.files(batch1), sep="/")
df <- data.frame(filepath=filesBatch1)
df$filepath <- as.character(df$filepath)

wavs <- mclapply(1:nrow(df), function(j) readWave(df$filepath[j], from=1, to=30, units="seconds"), mc.cores=getOption("mc.cores", 4))
wav.names <- sapply(basename(df$filepath), function(x) strsplit(x, "/"))

#PreEvent X / Event X/ PostEvent
setwd("/Users/kristenbellisario/Documents/R_Dir/CP_SoundEvent_2015/spectrograms/PreEvent")
spec <- lapply(1:length(wavs), function(k) {
  jpg(paste(wav.names[k], ".png", sep=""), height=400, width=600) 
  spectro(wavs[[k]], f=44100, wl=512, flim=c(.5, 8), plot=T)
  dev.off()
})


