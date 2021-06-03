## WEATHER CORRELATION
## Civil War Reenactment - LSE 
## Compute Daily Averages and Pull In Acoustic Indices
## Kristen Bellisario

#########################################################################
# INTERVAL DATA - needs non-parametric test
library(dplyr)
weather2015 <- read.csv("data/2015_Weather.csv", header=T)

#2015 - Indices
dir <- ("/data")
Pre_9219 <- read.csv(paste(dir,"Pre_9219.csv", sep=""))
Pre_15106 <- read.csv(paste(dir,"Pre_15106.csv", sep=""))
Pre_15116 <- read.csv(paste(dir,"Pre_15116.csv", sep=""))

Event_9219 <- read.csv(paste(dir,"Event_9219.csv", sep=""))
Event_15106 <- read.csv(paste(dir,"Event_15106.csv", sep=""))
Event_15116 <- read.csv(paste(dir,"Event_15116.csv", sep=""))

Post_9219 <- read.csv(paste(dir,"Post_AgPrairie.csv", sep=""))
Post_15106 <- read.csv(paste(dir,"Post_Riparian.csv", sep=""))
Post_15116 <- read.csv(paste(dir,"Post_Prairie.csv", sep=""))

# 2015: Compile list of means of 24 hours for Hf for each day, labeled according to experimental period
L9219_15 <- list (
  Pre_9219_Hf <- data.frame(Pre_9219[,5], rep(1:7, each=24)),
  names(Pre_9219_Hf) <- c("Hf", "Day"),
  res1 <- Pre_9219_Hf %>% group_by(Day) %>%
    summarise(value=mean(Hf)),
  Event_9219_Hf <- data.frame(Event_9219[,4], rep(1:2, each=24)),
  names(Event_9219_Hf) <- c("Hf", "Day"),
  res2 <- Event_9219_Hf %>% group_by(Day) %>%
    summarise(value=mean(Hf)),
  Post_9219_Hf <- data.frame(Post_9219[1:192,4], rep(1:8, each=24)),
  names(Post_9219_Hf) <- c("Hf", "Day"),
  res3 <- Post_9219_Hf %>% group_by(Day) %>%
    summarise(value=mean(Hf))
)


#2015 Hf list for sensor 15106
L15106_15 <- list(
  Pre_15106_Hf <- data.frame(Pre_15106[,5], rep(1:7, each=24)),
  names(Pre_15106_Hf) <- c("Hf", "Day"),
  res4 <- Pre_15106_Hf %>% group_by(Day) %>%
    summarise(value=mean(Hf)),
  Event_15106_Hf <- data.frame(Event_15106[,4], rep(1:2, each=24)),
  names(Event_15106_Hf) <- c("Hf", "Day"),
  res5 <- Event_15106_Hf %>% group_by(Day) %>%
    summarise(value=mean(Hf)),
  Post_15106_Hf <- data.frame(Post_15106[1:192,4], rep(1:8, each=24)),
  names( Post_15106_Hf) <- c("Hf", "Day"),
  res6 <- Post_15106_Hf %>% group_by(Day) %>%
    summarise(value=mean(Hf))
)

#2015 Hf list for 15116
L15116_15 <- list(  
  Pre_15116_Hf <- data.frame(Pre_15116[,5], rep(1:7, each=24)),
  names(Pre_15116_Hf) <- c("Hf", "Day"),
  res7 <- Pre_15116_Hf %>% group_by(Day) %>%
    summarise(value=mean(Hf)),
  Event_15116_Hf <- data.frame(Event_15116[,4], rep(1:2, each=24)),
  names( Event_15116_Hf) <- c("Hf", "Day"),
  res8 <- Event_15116_Hf %>% group_by(Day) %>%
    summarise(value=mean(Hf)),
  Post_15116_Hf <- data.frame(Post_15116[1:192,4], rep(1:8, each=24)),
  names(Post_15116_Hf) <- c("Hf", "Day"),
  res9 <- Post_15116_Hf %>% group_by(Day) %>%
    summarise(value=mean(Hf))
)

#combine 9219 / 15106 / 15116 --- include only 2015 res
Hf_total_9219 <- data.frame(rbind(res1, res2, res3, res1.6, res2.6, res3.6))
Hf_total_15106 <- data.frame(rbind(res4, res5, res6, res4.6, res5.6, res6.6))
Hf_total_15116 <- data.frame(rbind(res7, res8, res9, res7.6, res8.6, res9.6))
weather2015[,1] <- as.Date(weather2015[,1], format="%M/%D/%Y")
w2 <- seq.Date(weather2015[,1])
day <- seq(1:17)
day2 <- seq(1:17)

#generalized linear model
#Hf 2015 --> 
install.packages("lmerTest")
library(lme4)
library(lmerTest)
Hf9219 <- lmer(Hf_total_9219[,2] ~ weather_total[,3]*weather_total[,4] + (1|weather_total[,1]), data=weather_total, REML=T)
Hf15106 <- lmer(Hf_total_15106[,2] ~ weather_total[,3]*weather_total[,4] + (1|weather_total[,1]))
Hf15116 <- lmer(Hf_total_15116[,2] ~ weather_total[,3]*weather_total[,4] + (1|weather_total[,1]))
#residuals

# Compare Acoustic Density measurements with interval data from Raw Listening -- proposed bird
ADTest <- read.csv("/data/AcousticDensityTest.csv")
boxplot(ADTest[,3]*10~ADTest[,4])
boxplot(ADTest[,3]*10~ADTest[,5])
plot(ADTest[,3]*10~ADTest[,4])

Cum_9219 <- rbind(Pre_9219[,c(1,5)], Event_9219[,c(1,4)], Post_9219[,c(1,4)])
Cum_11516 <- rbind(Pre_15116[,c(1,5)], Event_15116[,c(1,4)], Post_15116[,c(1,4)])
Cum_15106 <- rbind(Pre_15106[,c(1,5)], Event_15106[,c(1,4)], Post_15106[,c(1,4)])

lo <- weather2015[c(1:16),c(1,4)]
avg <- weather2015[c(1:16),c(1,3)]
hi <- weather2015[c(1:16),c(1,2)]
precip <- weather2015[c(1:16),c(1,16)]
time2 <- seq(1:16)
time <- seq(1:432)


# Weather-Index Plots -----------------------------------------------------
#for saving directly
printdir <- ("UserDefined")

# 9219_weather ------------------------------------------------------------

## 9219 LO

par(mar=c(5, 4, 4, 6) + 0.1)   ## add extra space to right margin of plot within frame
plot(time, Cum_9219[,2], pch=16, axes=FALSE, ylab="Hf Index",xlab="", 
     type="p",col="black", main="Riparian")     ## Plot first set of data and draw its axis
axis(2,ylim=c(0,1),col="black",las=1)  ## las=1 makes horizontal labels
box()
par(new=TRUE)   ## Allow a second plot on the same graph
plot(time2, lo[,2], pch=15,  xlab="", ylab="", ylim=NULL, axes=FALSE, type="l", col="red")   ## Plot the second plot and put axis scale on right
mtext("Low Temperature",side=4,col="red",line=4) ## a little farther out (line=4) to make room for labels
axis(4, ylim=c(0,90), col="red",col.axis="red",las=1)  ## Draw the time axis
axis(1,pretty(range(time2), n=8))
mtext("Days",side=1,col="black",line=2.5, cex=1)  
legend("bottomleft",legend=c("Hf","Lo"),
       text.col=c("black","red"),pch=c(16,15),col=c("black","red"))   ## Add Legend
#png(filename=paste(printdir, "9219-lo.png", sep="")) ## need to export all layers - how

#9219 - AVG
plot(time, Cum_9219[,2], pch=16, axes=FALSE, ylab="Hf Index",xlab="", 
     type="p",col="black", main="Riparian")
axis(2,ylim=c(0,1),col="black",las=1)  ## las=1 makes horizontal labels
box()
par(new=TRUE)
plot(time2, avg[,2], pch=15,  xlab="", ylab="", ylim=NULL, axes=FALSE, type="l", col="red")
mtext("Avg Temperature",side=4,col="red",line=4) 
axis(4, ylim=c(0,90), col="red",col.axis="red",las=1)
axis(1,pretty(range(time2), n=8))
mtext("Days",side=1,col="black",line=2.5, cex=.6)  
legend("bottomleft",legend=c("Hf","Avg"),
       text.col=c("black","red"),pch=c(16,15),col=c("black","red"))

#9219 - HI
plot(time, Cum_9219[,2], pch=16, axes=FALSE, ylab="Hf Index",xlab="", 
     type="p",col="black", main="Riparian")
axis(2,ylim=c(0,1),col="black",las=1)  ## las=1 makes horizontal labels
box()
par(new=TRUE)
plot(time2, hi[,2], pch=15,  xlab="", ylab="", ylim=NULL, axes=FALSE, type="l", col="red")
mtext("High Temperature",side=4,col="red",line=4) 
axis(4, ylim=c(0,90), col="red",col.axis="red",las=1)
axis(1,pretty(range(time2), n=8))
mtext("Days",side=1,col="black",line=2.5, cex=.6)  
legend("bottomleft",legend=c("Hf","High"),
       text.col=c("black","red"),pch=c(16,15),col=c("black","red"))

#9219 - PRECIP
plot(time, Cum_9219[,2], pch=16, axes=FALSE, ylab="Hf Index",xlab="", 
     type="p",col="black", main="Riparian")
axis(2,ylim=c(0,1),col="black",las=1)  ## las=1 makes horizontal labels
box()
par(new=TRUE)
plot(time2, precip[,2], pch=15,  xlab="", ylab="", ylim=NULL, axes=FALSE, type="l", col="red")
mtext("Precipitation",side=4,col="red",line=4) 
axis(4, ylim=c(0,90), col="red",col.axis="red",las=1)
axis(1,pretty(range(time2), n=8))
mtext("Days",side=1,col="black",line=2.5, cex=.6)  
legend("bottomleft",legend=c("Hf","Precip"),
       text.col=c("black","red"),pch=c(16,15),col=c("black","red"))

# 11516_weather -----------------------------------------------------------

#11516-Lo
plot(time, Cum_11516[,2], pch=16, axes=FALSE, ylab="Hf Index",xlab="", 
     type="p",col="black", main="Prairie")
axis(2,ylim=c(0,1),col="black",las=1)  ## las=1 makes horizontal labels
box()
par(new=TRUE)
plot(time2, lo[,2], pch=15,  xlab="", ylab="", ylim=NULL, axes=FALSE, type="l", col="red")
mtext("Low Temperature",side=4,col="red",line=4) 
axis(4, ylim=c(0,90), col="red",col.axis="red",las=1)
axis(1,pretty(range(time2), n=8))
mtext("Days",side=1,col="black",line=2.5, cex=.6)  
legend("bottomleft",legend=c("Hf","Lo"),
       text.col=c("black","red"),pch=c(16,15),col=c("black","red"))

#11516 - AVG
plot(time, Cum_11516[,2], pch=16, axes=FALSE, ylab="Hf Index",xlab="", 
     type="p",col="black", main="Prairie")
axis(2,ylim=c(0,1),col="black",las=1)  ## las=1 makes horizontal labels
box()
par(new=TRUE)
plot(time2, avg[,2], pch=15,  xlab="", ylab="", ylim=NULL, axes=FALSE, type="l", col="red")
mtext("Avg Temperature",side=4,col="red",line=4) 
axis(4, ylim=c(0,90), col="red",col.axis="red",las=1)
axis(1,pretty(range(time2), n=8))
mtext("Days",side=1,col="black",line=2.5, cex=.6)  
legend("bottomleft",legend=c("Hf","Avg"),
       text.col=c("black","red"),pch=c(16,15),col=c("black","red"))

#11516 - HI
plot(time, Cum_11516[,2], pch=16, axes=FALSE, ylab="Hf Index",xlab="", 
     type="p",col="black", main="Prairie")
axis(2,ylim=c(0,1),col="black",las=1)  ## las=1 makes horizontal labels
box()
par(new=TRUE)
plot(time2, hi[,2], pch=15,  xlab="", ylab="", ylim=NULL, axes=FALSE, type="l", col="red")
mtext("High Temperature",side=4,col="red",line=4) 
axis(4, ylim=c(0,90), col="red",col.axis="red",las=1)
axis(1,pretty(range(time2), n=8))
mtext("Days",side=1,col="black",line=2.5, cex=.6)  
legend("bottomleft",legend=c("Hf","High"),
       text.col=c("black","red"),pch=c(16,15),col=c("black","red"))

#11516 - PRECIP
plot(time, Cum_11516[,2], pch=16, axes=FALSE, ylab="Hf Index",xlab="", 
     type="p",col="black", main="Prairie")
axis(2,ylim=c(0,1),col="black",las=1)  ## las=1 makes horizontal labels
box()
par(new=TRUE)
plot(time2, precip[,2], pch=15,  xlab="", ylab="", ylim=NULL, axes=FALSE, type="l", col="red")
mtext("Precipitation",side=4,col="red",line=4) 
axis(4, ylim=c(0,90), col="red",col.axis="red",las=1)
axis(1,pretty(range(time2), n=8))
mtext("Days",side=1,col="black",line=2.5, cex=.6)  
legend("bottomleft",legend=c("Hf","Precip"),
       text.col=c("black","red"),pch=c(16,15),col=c("black","red"))


# 15106_weather -----------------------------------------------------------

#15106-Lo
plot(time, Cum_15106[,2], pch=16, axes=FALSE, ylab="Hf Index",xlab="", 
     type="p",col="black", main="Ag")
axis(2,ylim=c(0,1),col="black",las=1)  ## las=1 makes horizontal labels
box()
par(new=TRUE)
plot(time2, lo[,2], pch=15,  xlab="", ylab="", ylim=NULL, axes=FALSE, type="l", col="red")
mtext("Low Temperature",side=4,col="red",line=4) 
axis(4, ylim=c(0,90), col="red",col.axis="red",las=1)
axis(1,pretty(range(time2), n=8))
mtext("Days",side=1,col="black",line=2.5, cex=.6)  
legend("bottomleft",legend=c("Hf","Lo"),
       text.col=c("black","red"),pch=c(16,15),col=c("black","red"))

#15106 - AVG
plot(time, Cum_15106[,2], pch=16, axes=FALSE, ylab="Hf Index",xlab="", 
     type="p",col="black", main="Ag")
axis(2,ylim=c(0,1),col="black",las=1)  ## las=1 makes horizontal labels
box()
par(new=TRUE)
plot(time2, avg[,2], pch=15,  xlab="", ylab="", ylim=NULL, axes=FALSE, type="l", col="red")
mtext("Avg Temperature",side=4,col="red",line=4) 
axis(4, ylim=c(0,90), col="red",col.axis="red",las=1)
axis(1,pretty(range(time2), n=8))
mtext("Days",side=1,col="black",line=2.5, cex=.6)  
legend("bottomleft",legend=c("Hf","Avg"),
       text.col=c("black","red"),pch=c(16,15),col=c("black","red"))

#15106 - HI
plot(time, Cum_15106[,2], pch=16, axes=FALSE, ylab="Hf Index",xlab="", 
     type="p",col="black", main="Ag")
axis(2,ylim=c(0,1),col="black",las=1)  ## las=1 makes horizontal labels
box()
par(new=TRUE)
plot(time2, hi[,2], pch=15,  xlab="", ylab="", ylim=NULL, axes=FALSE, type="l", col="red")
mtext("High Temperature",side=4,col="red",line=4) 
axis(4, ylim=c(0,90), col="red",col.axis="red",las=1)
axis(1,pretty(range(time2), n=8))
mtext("Days",side=1,col="black",line=2.5, cex=.6)  
legend("bottomleft",legend=c("Hf","High"),
       text.col=c("black","red"),pch=c(16,15),col=c("black","red"))

#15106 - PRECIP
plot(time, Cum_15106[,2], pch=16, axes=FALSE, ylab="Hf Index",xlab="", 
     type="p",col="black", main="Ag")
axis(2,ylim=c(0,1),col="black",las=1)  ## las=1 makes horizontal labels
box()
par(new=TRUE)
plot(time2, precip[,2], pch=15,  xlab="", ylab="", ylim=NULL, axes=FALSE, type="l", col="red")
mtext("Precipitation",side=4,col="red",line=4) 
axis(4, ylim=c(0,90), col="red",col.axis="red",las=1)
axis(1,pretty(range(time2), n=8))
mtext("Days",side=1,col="black",line=2.5, cex=.6)  
legend("bottomleft",legend=c("Hf","Precip"),
       text.col=c("black","red"),pch=c(16,15),col=c("black","red"))
