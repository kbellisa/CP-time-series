## Species Accumulation 
## CP-Chapter for Dissertation
## Kristen Bellisario
## Baseline code for evaluation - not optimized

library(vegan)
library(scales)

# Data Input --------------------------------------------------------------

morning <- read.csv("data/MORNING-Table-1-F.csv")
event <- read.csv("data/EVENT-Table-1.csv")
morning[is.na(morning)] <- 0
morn.tr <- morning[,2]
grassland.morning <- as.data.frame(cbind(morn.tr, morning$BOBO, morning$GRSP, morning$HESP, morning$HOLA, morning$SAVS))
names(grassland.morning) <- c( "TREATMENT", "BOBO", "GRSP", "HESP", "HOLA", "SAVS")
generalist.morning <- morning[,c(-1, -8,-23, -24, -25,-34)]

event[is.na(event)] <- 0
event.tr <- event[,2]
grassland.event <- data.frame(cbind(event.tr, event$BOBO, event$GRSP, event$HESP,  event$HOLA, event$SAVS))
names(grassland.event) <- c( "TREATMENT", "BOBO", "GRSP", "HESP", "HOLA", "SAVS")
generalist.event <- event[,c(-1, -8,-23, -24, -25,-34)]

'SPECIES POOL MODELING OPTION
morning.mat <- as.matrix(morning[,3:42])
specpool(morning.mat) #get gamma diversity
Species     chao  chao.se    jack1 jack1.se    jack2     boot  boot.se   n
All      22 23.59444 2.155261 25.98611 1.993056 25.01036 24.27635 1.337967 288

event.mat <- as.matrix(event[,3:42])
specpool(event.mat)
Species     chao  chao.se   jack1 jack1.se    jack2     boot  boot.se  n
All      32 51.82812 19.77728 40.8125 3.253804 47.56161 35.73388 1.793866 48'

# morning sites -----------------------------------------------------------

pre.rip <- morning[which(morning$TREATMENT=="9219_PRE"),]
pre.pra <- morning[which(morning$TREATMENT=="15116_PRE"),]
pre.ag <- morning[which(morning$TREATMENT=="15106_PRE"),]
event.rip <- morning[which(morning$TREATMENT=="9219_EVENT"),]
event.pra <- morning[which(morning$TREATMENT=="15116_EVENT"),]
event.ag <- morning[which(morning$TREATMENT=="15106_EVENT"),]
post.rip <- morning[which(morning$TREATMENT=="9219_POST"),]
post.pra <- morning[which(morning$TREATMENT=="15116_POST"),]
post.ag <- morning[which(morning$TREATMENT=="15106_POST"),]

# repeat for obligate versus non-obligate for each period
pre.rip.g <- grassland.morning[which(grassland.morning$TREATMENT=="9"),]
pre.pra.g <- grassland.morning[which(grassland.morning$TREATMENT=="6"),]
pre.ag.g <- grassland.morning[which(grassland.morning$TREATMENT=="3"),]

pre.rip.ge <- generalist.morning[which(generalist.morning$TREATMENT=="9219_PRE"),]
pre.pra.ge <- generalist.morning[which(generalist.morning$TREATMENT=="15116_PRE"),]
pre.ag.ge <- generalist.morning[which(generalist.morning$TREATMENT=="15106_PRE"),]

post.rip.g <- grassland.morning[which(grassland.morning$TREATMENT=="8"),]
post.pra.g <- grassland.morning[which(grassland.morning$TREATMENT=="5"),]
post.ag.g <- grassland.morning[which(grassland.morning$TREATMENT=="2"),]

post.rip.ge <- generalist.morning[which(generalist.morning$TREATMENT=="9219_POST"),]
post.pra.ge <- generalist.morning[which(generalist.morning$TREATMENT=="15116_POST"),]
post.ag.ge <- generalist.morning[which(generalist.morning$TREATMENT=="15106_POST"),]

# repeat for obligate versus non-obligate for each period
pre.rip.g.evt <- grassland.event[which(grassland.event$TREATMENT=="9"),]
pre.pra.g.evt <- grassland.event[which(grassland.event$TREATMENT=="6"),]
pre.ag.g.evt <- grassland.event[which(grassland.event$TREATMENT=="3"),]

pre.rip.ge.evt <- generalist.event[which(generalist.event$TREATMENT=="9219-PRE"),]
pre.pra.ge.evt <- generalist.event[which(generalist.event$TREATMENT=="15116-PRE"),]
pre.ag.ge.evt <- generalist.event[which(generalist.event$TREATMENT=="15106-PRE"),]

post.rip.g.evt <- grassland.event[which(grassland.event$TREATMENT=="8"),]
post.pra.g.evt <- grassland.event[which(grassland.event$TREATMENT=="5"),]
post.ag.g.evt <- grassland.event[which(grassland.event$TREATMENT=="2"),]

post.rip.ge.evt <- generalist.event[which(generalist.event$TREATMENT=="9219-POST"),]
post.pra.ge.evt <- generalist.event[which(generalist.event$TREATMENT=="15116-POST"),]
post.ag.ge.evt <- generalist.event[which(generalist.event$TREATMENT=="15106-POST"),]

### TEST IF SIGNIFICANT DIFFERENCE IN PRE-POST FOR GRASSLAND MORNING
pre.ripsum.g <- colSums(pre.rip.g[,-1])
pre.prasum.g <- colSums(pre.pra.g[,-1])
pre.agsum.g <- colSums(pre.ag.g[,-1])

morning.g.exp <- rbind(pre.ripsum.g, pre.prasum.g, pre.agsum.g)

post.ripsum.g <- colSums(post.rip.g[,-1])
post.prasum.g <- colSums(post.pra.g[,-1])
post.agsum.g <- colSums(post.ag.g[,-1])

morning.g.obs <- rbind(post.ripsum.g, post.prasum.g, post.agsum.g)
observed <- colSums(morning.g.obs)
expected <- colSums(morning.g.exp)

chisq.test(matrix(c(observed, expected), nrow=2, ncol=5), simulate.p.value=T)

## GRASSLAND AFTERNOON
pre.ripsum.g.evt <- colSums(pre.rip.g.evt[,-1])
pre.prasum.g.evt <- colSums(pre.pra.g.evt[,-1])
pre.agsum.g.evt<- colSums(pre.ag.g.evt[,-1])

event.g.exp <- rbind(pre.ripsum.g.evt, pre.prasum.g.evt, pre.agsum.g.evt)

post.ripsum.g.evt <- colSums(post.rip.g.evt[,-1])
post.prasum.g.evt <- colSums(post.pra.g.evt[,-1])
post.agsum.g.evt <- colSums(post.ag.g.evt[,-1])

event.g.obs <- rbind(post.ripsum.g.evt, post.prasum.g.evt, post.agsum.g.evt)

observed <- colSums(event.g.obs)
expected <- colSums(event.g.exp)

chisq.test(matrix(c(observed, expected), nrow=2, ncol=5), simulate.p.value=F)
chisq.test(matrix(c(observed, expected), nrow=2, ncol=5), simulate.p.value=T)

## GENERALIST MORNING
pre.ripsum.ge <- colSums(pre.rip.ge[,-1])
pre.prasum.ge <- colSums(pre.pra.ge[,-1])
pre.agsum.ge <- colSums(pre.ag.ge[,-1])

morning.ge.exp <- rbind(pre.ripsum.ge, pre.prasum.ge, pre.agsum.ge)

post.ripsum.ge <- colSums(post.rip.ge[,-1])
post.prasum.ge <- colSums(post.pra.ge[,-1])
post.agsum.ge <- colSums(post.ag.ge[,-1])

morning.ge.obs <- rbind(post.ripsum.ge, post.prasum.ge, post.agsum.ge)

observed <- colSums(morning.ge.obs)
expected <- colSums(morning.ge.exp)

observed <- observed[c(1,6,7,9,11,14,15,18,19,21,23,25,26,27,29,30,31,32),]
expected <- expected[c(1,6,7,9,11,14,15,18,19,21,23,25,26,27,29,30,31,32),]

chisq.test(observed, expected, simulate.p.value=F)
chisq.test(observed, expected, simulate.p.value=T, B=1000)

## GENERALIST AFTERNOON
pre.ripsum.ge.evt <- colSums(pre.rip.ge.evt[,-1])
pre.prasum.ge.evt <- colSums(pre.pra.ge.evt[,-1])
pre.agsum.ge.evt<- colSums(pre.ag.ge.evt[,-1])

event.ge.exp <- rbind(pre.ripsum.ge.evt, pre.prasum.ge.evt, pre.agsum.ge.evt)

post.ripsum.ge.evt <- colSums(post.rip.ge.evt[,-1])
post.prasum.ge.evt <- colSums(post.pra.ge.evt[,-1])
post.agsum.ge.evt <- colSums(post.ag.ge.evt[,-1])

event.ge.obs <- rbind(post.ripsum.ge.evt, post.prasum.ge.evt, post.agsum.ge.evt)

observed <- colSums(event.ge.obs)
expected <- colSums(event.ge.exp)

chisq.test(observed, expected, simulate.p.value=F)
chisq.test(observed, expected, simulate.p.value=T, B=1000)

#sum(((observed-expected)^2)/expected)
# this calculation
#TO COMPARE GROUPS - assign the relative vocal indicator species to each sample period
# Then test normailty with Shapiro wilks test / DONT NEED NORMALITY TEST FOR POISSON DATA
# then chi-square test between pre / post group of species (not per site)

morning.rvi <- data.frame(morning$TREATMENT, morning$FISP, morning$GRSP, morning$HESP, morning$KILL, morning$NOCA, morning$RWBL, morning$SOSP, morning$GHOW, morning$AMRO, morning$EASO)
names(morning.rvi) <- c("TREATMENT", "FISP", "GRSP", "HESP", "KILL", "NOCA", "RWBL", "SOSP", "GHOW", "AMRO", "EASO")
pre.rip.r <- morning.rvi[which(morning.rvi$TREATMENT=="9219_PRE"),]
pre.pra.r <- morning.rvi[which(morning.rvi$TREATMENT=="15116_PRE"),]
pre.ag.r <- morning.rvi[which(morning.rvi$TREATMENT=="15106_PRE"),]

pre.ripsum <- colSums(pre.rip.r[,-1])
pre.prasum <- colSums(pre.pra.r[,-1])
pre.agsum <- colSums(pre.ag.r[,-1])

morning.rvi.exp <- rbind(pre.ripsum, pre.prasum, pre.agsum)

post.rip.r <- morning.rvi[which(morning.rvi$TREATMENT=="9219_POST"),]
post.pra.r <- morning.rvi[which(morning.rvi$TREATMENT=="15116_POST"),]
post.ag.r <- morning.rvi[which(morning.rvi$TREATMENT=="15106_POST"),]

post.ripsum <- colSums(post.rip.r[,-1])
post.prasum <- colSums(post.pra.r[,-1])
post.agsum <- colSums(post.ag.r[,-1])

morning.rvi.obs <- rbind(post.ripsum, post.prasum, post.agsum)
observed <- colSums(morning.rvi.obs)
expected <- colSums(morning.rvi.exp)

chisq.test(matrix(c(observed, expected), nrow=2, ncol=10), simulate.p.value=T)
#sum(((observed-expected)^2)/expected)
# this calculation

event.rvi <- data.frame(event$TREATMENT, event$BLJA, event$BWWA, event$CACH, event$CARW, event$COYE, event$EAME, event$EAWP, event$FISP, event$GCFL, event$GCKI, event$INBU, event$GRSP, event$KILL, event$NOCA, event$RBWO, event$REVI, event$RWBL, event$SAVS, event$ETTI, event$WBNU, event$YTVI)
names(event.rvi) <- c("TREATMENT", "BLJA", "BWWA", "CACH", "CARW", "COYE", "EAME", "EAWP", "FISP", "GCFL", "GCKI", "INBU", "GRSP", "KILL", "NOCA", "RBWO", "REVI", "RWBL", "SAVS", "ETTI", "WBNU", "YTVI")

pre.ripe.r <- event.rvi[which(event.rvi$TREATMENT=="9219-PRE"),]
pre.prae.r <- event.rvi[which(event.rvi$TREATMENT=="15116-PRE"),]
pre.age.r <- event.rvi[which(event.rvi$TREATMENT=="15106-PRE"),]

pre.ripsume <- colSums(pre.ripe.r[,-1])
pre.prasume <- colSums(pre.prae.r[,-1])
pre.agsume <- colSums(pre.age.r[,-1])

event.rvi.exp1 <- rbind(pre.ripsume, pre.prasume, pre.agsume)

post.ripe.r <- event.rvi[which(event.rvi$TREATMENT=="9219-POST"),]
post.prae.r <- event.rvi[which(event.rvi$TREATMENT=="15116-POST"),]
post.age.r <- event.rvi[which(event.rvi$TREATMENT=="15106-POST"),]

post.ripsume <- colSums(post.ripe.r[,-1])
post.prasume <- colSums(post.prae.r[,-1])
post.agsume <- colSums(post.age.r[,-1])

event.rvi.obs1 <- rbind(post.ripsume, post.prasume, post.agsume)

observede <- colSums(event.rvi.obs1)
expectede <- colSums(event.rvi.exp1)

chisq.test(matrix(c(observede, expectede), nrow=2, ncol=21), simulate.p.value=T)

# NOW TEST TO SEE IF THERE IS A DIFFERENCE BETWEEN GRASSLAND & GENERALISTS (in terms of count)
# then chi-square test between pre / post group of species (not per site)

sp1 <- specaccum(pre.rip[,3:length(pre.rip)], method="rarefaction", xvar="effort")
sp2 <- specaccum(event.rip[,3:length(event.rip)], method="rarefaction", xvar="effort")
sp3 <- specaccum(post.rip[1:12,3:length(post.rip)], method="rarefaction", xvar="effort")

sp4 <- specaccum(pre.pra[,3:length(pre.pra)], method="rarefaction", xvar="effort")
sp5 <- specaccum(event.pra[,3:length(event.pra)], method="rarefaction", xvar="effort")
sp6 <- specaccum(post.pra[1:12,3:length(post.pra)], method="rarefaction", xvar="effort")

sp7 <- specaccum(pre.ag[,3:length(pre.ag)], method="rarefaction", xvar="effort")
sp8 <- specaccum(event.ag[,3:length(event.ag)], method="rarefaction", xvar="effort")
sp9 <- specaccum(post.ag[1:12,3:length(post.ag)], method="rarefaction", xvar="effort")

# THIS IS FROM SPECIES POOL MODEL - EXTRACT
'specaccum.list <- list(max(sp1$sites), max(sp1$individuals), max(sp1$richness), max(sp1$sd),
    max(sp2$sites), max(sp2$individuals), max(sp2$richness), max(sp2$sd),
    max(sp3$sites), max(sp3$individuals), max(sp3$richness), max(sp3$sd),
    max(sp4$sites), max(sp4$individuals), max(sp4$richness), max(sp4$sd),
    max(sp5$sites), max(sp5$individuals), max(sp5$richness), max(sp5$sd),
    max(sp6$sites), max(sp6$individuals), max(sp6$richness), max(sp6$sd),
    max(sp7$sites), max(sp7$individuals), max(sp7$richness), max(sp7$sd),
    max(sp8$sites), max(sp8$individuals), max(sp8$richness), max(sp8$sd),
    max(sp9$sites), max(sp9$individuals), max(sp9$richness), max(sp9$sd))

M <- as.data.frame(matrix(data=specaccum.list, ncol=4, byrow=T))
names(M) <- c("sites", "individuals", "richness", "sd")'

# CALCULATE PERCENTAGE PER SPECIES IN PRE AND POST
rip.species <- rbind(pre.rip, post.rip)
rip.sp.tot <- colSums(rip.species[,3:42])
pre.rip.tot <- colSums(pre.rip[,3:42])
post.rip.tot <- colSums(post.rip[,3:42])
pre.rip.perc <- data.frame(t(na.omit(data.frame(pre.rip.tot/rip.sp.tot))))
post.rip.perc <- data.frame(t(na.omit(data.frame(post.rip.tot/rip.sp.tot))))
#pre / post == 12 AMRO, BAOW, CAGO, EASO, FISP, GHOW, RBWO, RWBL, SOSP, THRUSH, TRSW

rip.pre.mat <- as.matrix(pre.rip[,3:42])
rip.pre.mat.spec <- specpool(rip.pre.mat)
'Species     chao  chao.se    jack1 jack1.se    jack2     boot  boot.se  n
All       8 9.465116 2.245601 10.93023 1.691771 10.99834 9.527181 1.223714 43

rip.event.mat <- as.matrix(event.rip[,3:42])
specpool(rip.event.mat)
Species chao  chao.se jack1 jack1.se    jack2    boot boot.se  n
All       7 23.5 21.68525  12.5 3.246793 16.74242 9.22413 1.71584 12'

rip.post.mat <- as.matrix(post.rip[,3:42])
rip.post.mat.spec <- specpool(rip.post.mat)
'Species     chao   chao.se   jack1  jack1.se    jack2   boot   boot.se  n
All       6 6.121951 0.4302995 6.97561 0.9756098 4.217073 6.8829 0.9164939 41'

# TOTAL RIPARIAN SITE
rip.mat <- as.matrix(rip.species)
rip.mat.tot.spec <- specpool(rip.mat)
'Species     chao chao.se    jack1 jack1.se    jack2     boot  boot.se  n
All      14 14.88929 1.44979 16.96429 1.711431 15.07071 15.82198 1.311827 84
pool <- poolaccum(rip.mat)
plot(pool, main="Species Pool Models - Riparian")'

pra.species <- rbind(pre.pra, post.pra)
pra.sp.tot <- colSums(pra.species[,3:42])
pre.pra.tot<-colSums(pre.pra[,3:42])
post.pra.tot <- colSums(post.pra[,3:42])
pre.pra.perc <- data.frame(t(na.omit(data.frame(pre.pra.tot/pra.sp.tot))))
post.pra.perc <- data.frame(t(na.omit(data.frame(post.pra.tot/pra.sp.tot))))
#pre / post pr == 12 AMRO, CONI, EAME, FISP, GHOW, GEESE, GRSP, HESP, KILL, RWBL, SASP, SOSP

pra.mat <- as.matrix(pra.species)
pra.mat.tot.spec <- specpool(pra.mat)
boxplot(specpool(pra.mat))
'Species     chao  chao.se    jack1 jack1.se    jack2     boot  boot.se  n
All      14 17.95238 5.233006 17.95238  1.97619 19.92828 15.83844 1.152321 84'

pra.pre.mat <- as.matrix(pre.pra[,3:42])
pra.pre.mat.spec <- specpool(pra.pre.mat)
'Species     chao chao.se    jack1 jack1.se    jack2     boot   boot.se  n
All       9 16.80952 11.3924 12.90476 1.952381 15.78513 10.59915 0.9950267 42'

pra.post.mat <- as.matrix(post.pra[,3:42])
pra.post.mat.spec <- specpool(pra.post.mat)
'Species     chao chao.se    jack1 jack1.se    jack2     boot  boot.se  n
All       8 15.80952 11.3924 11.90476 1.952381 14.78513 9.627151 1.112923 42'

ag.species <- rbind(pre.ag, post.ag)
ag.sp.tot <- colSums(ag.species[,3:42])
pre.ag.tot<-colSums(pre.ag[,3:42])
post.ag.tot <- colSums(post.ag[,3:42])
pre.ag.perc <- data.frame(t(na.omit(data.frame(pre.ag.tot/ag.sp.tot))))
post.ag.perc <- data.frame(t(na.omit(data.frame(post.ag.tot/ag.sp.tot))))
# pre/post ag == 13 CONI, EAME, FISP, GHOW, GEESE, GRSP, HESP, HOLA, KILL, NOCA, RWBL, SASP, SOSP

ag.pre.mat <- as.matrix(pre.ag[,3:42])
ag.pre.mat.spec <- specpool(ag.pre.mat)
'Species     chao  chao.se    jack1 jack1.se    jack2     boot  boot.se  n
All      11 16.85714 6.995257 14.90476 2.400869 18.71429 12.51967 1.205021 42'

ag.post.mat <- as.matrix(post.ag[,3:42])
ag.post.mat.spec <- specpool(ag.post.mat)
'Species    chao  chao.se    jack1 jack1.se    jack2    boot boot.se  n
All       8 8.97619 1.832637 9.952381 1.380542 9.998839 9.00441 1.03419 42'

ag.mat <- as.matrix(ag.species)
ag.mat.tot.spc <- specpool(ag.mat)
boxplot(specpool(ag.mat))
'Species     chao  chao.se    jack1 jack1.se    jack2     boot  boot.se  n
All      15 29.82143 13.46916 20.92857 2.798961 26.78571 17.24386 1.363296 84'

par(mfrow=c(3,3))
plot(sp1, ci.type="polygon", col="black", lwd=2, ci.lty=0, xlab = "observations", ylab = "number of species", ylim = c(0,12), ci=2, ci.col=alpha("blue", 0.1), main="Pre-Riparian")
plot(sp2, ci.type="polygon", col="black", lwd=2, ci.lty=0, xlab = "observations", ylab = "number of species",ylim = c(0,12), ci=2, ci.col=alpha("blue", 0.1), main="Event-Riparian")
plot(sp3, ci.type="polygon", col="black", lwd=2, ci.lty=0, xlab = "observations", ylab = "number of species",ylim = c(0,12), ci=2, ci.col=alpha("blue", 0.1), main="Post-Riparian")
# boxplot(sp1, col="yellow", add=TRUE, pch="+") if use RANDOM instead of SITE ORDER

plot(sp4, ci.type="polygon", col="black", lwd=2, ci.lty=0, xlab = "observations", ylab = "number of species",ylim = c(0,12), ci=2, ci.col=alpha("blue", 0.1), main="Pre-Prairie")
plot(sp5, ci.type="polygon", col="black", lwd=2, ci.lty=0, xlab = "observations", ylab = "number of species",ylim = c(0,12), ci=2, ci.col=alpha("blue", 0.1), main="Event-Prairie")
plot(sp6, ci.type="polygon", col="black", lwd=2, ci.lty=0, xlab = "observations", ylab = "number of species",ylim = c(0,12), ci=2, ci.col=alpha("blue", 0.1), main="Post-Prairie")

plot(sp7, ci.type="polygon", col="black", lwd=2, ci.lty=0, xlab = "observations", ylab = "number of species",ylim = c(0,12), ci=2, ci.col=alpha("blue", 0.1), main="Pre-Ag")
plot(sp8, ci.type="polygon", col="black", lwd=2, ci.lty=0, xlab = "observations", ylab = "number of species",ylim = c(0,12), ci=2, ci.col=alpha("blue", 0.1), main="Event-Ag")
plot(sp9, ci.type="polygon", col="black", lwd=2, ci.lty=0, xlab = "observations", ylab = "number of species",ylim = c(0,12), ci=2, ci.col=alpha("blue", 0.1), main="Post-Ag")

####
speclist <- data.frame(rbind(sp1$freq, sp2$freq, sp3$freq))

#REPEAT FOR EVENT SPEC ACCUM
pre.ripE <- event[which(event$TREATMENT=="9219-PRE"),]
pre.praE <- event[which(event$TREATMENT=="15116-PRE"),]
pre.agE <- event[which(event$TREATMENT=="15106-PRE"),]
event.ripE <- event[which(event$TREATMENT=="9219-EVENT"),]
event.praE <- event[which(event$TREATMENT=="15116-EVENT"),]
event.agE <- event[which(event$TREATMENT=="15106-EVENT"),]
post.ripE <- event[which(event$TREATMENT=="9219-POST"),]
post.praE <- event[which(event$TREATMENT=="15116-POST"),]
post.agE <- event[which(event$TREATMENT=="15106-POST"),]

# CALCULATE PERCENTAGE PER SPECIES IN PRE AND POST
rip.speciesE <- rbind(pre.ripE, post.ripE)
rip.sp.totE <- colSums(rip.speciesE[,3:51])
pre.rip.totE <- colSums(pre.ripE[,3:51])
post.rip.totE <- colSums(post.ripE[,3:51])
pre.rip.percE <- data.frame(t(na.omit(data.frame(pre.rip.totE/rip.sp.totE))))
post.rip.percE <- data.frame(t(na.omit(data.frame(post.rip.totE/rip.sp.totE))))
#

ripE.pre.mat <- as.matrix(pre.ripE[,3:42])
ripE.pre.mat.spec <- specpool(ripE.pre.mat)

'Species chao  chao.se jack1 jack1.se    jack2     boot  boot.se n
All      13   19 8.831761    16 1.936492 17.66667 14.34375 1.085685 4'

ripE.post.mat <- as.matrix(post.ripE[,3:42])
ripE.post.mat.spec <- specpool(ripE.post.mat)

'Species   chao  chao.se jack1 jack1.se    jack2     boot  boot.se n
All      15 16.875 2.321772 18.75 2.704163 19.58333 16.91016 1.911043 4'

pra.speciesE <- rbind(pre.praE, post.praE)
pra.sp.totE <- colSums(pra.speciesE[,3:42])
pre.pra.totE<-colSums(pre.praE[,3:42])
post.pra.totE <- colSums(post.praE[,3:42])
pre.pra.percE <- data.frame(t(na.omit(data.frame(pre.pra.totE/pra.sp.totE))))
post.pra.percE <- data.frame(t(na.omit(data.frame(post.pra.totE/pra.sp.totE))))

praE.pre.mat <- as.matrix(pre.praE[,3:42])
praE.pre.mat.spec <- specpool(praE.pre.mat)
'Species  chao  chao.se jack1 jack1.se    jack2     boot  boot.se n
All      12 18.75 7.739913  16.5 2.806243 18.83333 14.03125 1.598033 4'

praE.post.mat <- as.matrix(post.praE[,3:42])
praE.post.mat.spec <- specpool(praE.post.mat)
'Species chao   chao.se jack1 jack1.se    jack2     boot  boot.se n
All      11 11.3 0.7056912  12.5 1.620185 11.83333 11.94922 1.588951 4'

ag.speciesE <- rbind(pre.agE, post.agE)
ag.sp.totE <- colSums(ag.speciesE[,3:42])
pre.ag.totE<-colSums(pre.agE[,3:42])
post.ag.totE <- colSums(post.agE[,3:42])
pre.ag.percE <- data.frame(t(na.omit(data.frame(pre.ag.totE/ag.sp.totE))))
post.ag.percE <- data.frame(t(na.omit(data.frame(post.ag.totE/ag.sp.totE))))

agE.pre.mat <- as.matrix(pre.agE[,3:42])
agE.pre.mat.spec <- specpool(agE.pre.mat)
'Species   chao  chao.se jack1 jack1.se jack2     boot  boot.se n
All      12 18.125 6.329906 17.25 3.799671 19.75 14.40625 2.168444 4'

agE.post.mat <- as.matrix(post.agE[,3:42])
agE.post.mat.spec <- specpool(agE.post.mat)
'Species chao  chao.se jack1 jack1.se    jack2     boot  boot.se n
All      14 15.5 2.076656    17 1.936492 17.66667 15.53125 1.214737 4'

sp1E <- specaccum(pre.ripE[,3:length(pre.ripE)], method="rarefaction", xvar="effort")
sp2E <- specaccum(event.ripE[,3:length(event.ripE)], method="rarefaction", xvar="effort")
sp3E <- specaccum(post.ripE[,3:length(post.ripE)], method="rarefaction", xvar="effort")

sp4E <- specaccum(pre.praE[,3:length(pre.praE)], method="rarefaction", xvar="effort")
sp5E <- specaccum(event.praE[,3:length(event.praE)], method="rarefaction", xvar="effort")
sp6E <- specaccum(post.praE[,3:length(post.praE)], method="rarefaction", xvar="effort")

sp7E <- specaccum(pre.agE[,3:length(pre.agE)], method="rarefaction", xvar="effort")
sp8E <- specaccum(event.agE[,3:length(event.agE)], method="rarefaction", xvar="effort")
sp9E <- specaccum(post.agE[,3:length(post.agE)], method="rarefaction", xvar="effort")

plot(sp1E, ci.type="polygon", col="black", lwd=2, ci.lty=0, xlab = "observations", ylab = "number of species",ylim = c(0,20), ci=2, ci.col=alpha("blue", 0.1), main="Pre-Riparian E")
plot(sp2E, ci.type="polygon", col="black", lwd=2, ci.lty=0, xlab = "observations", ylab = "number of species",ylim = c(0,20), ci=2, ci.col=alpha("blue", 0.1), main="Event-Riparian E")
plot(sp3E, ci.type="polygon", col="black", lwd=2, ci.lty=0, xlab = "observations", ylab = "number of species",ylim = c(0,20), ci=2, ci.col=alpha("blue", 0.1), main="Post-Riparian E")
# boxplot(sp1, col="yellow", add=TRUE, pch="+") if use RANDOM instead of SITE ORDER

plot(sp4E, ci.type="polygon", col="black", lwd=2, ci.lty=0, xlab = "observations", ylab = "number of species",ylim = c(0,20), ci=2, ci.col=alpha("blue", 0.1), main="Pre-Prairie E")
plot(sp5E, ci.type="polygon", col="black", lwd=2, ci.lty=0, xlab = "observations", ylab = "number of species",ylim = c(0,20), ci=2, ci.col=alpha("blue", 0.1), main="Event-Prairie E")
plot(sp6E, ci.type="polygon", col="black", lwd=2, ci.lty=0, xlab = "observations", ylab = "number of species",ylim = c(0,20), ci=2, ci.col=alpha("blue", 0.1), main="Post-Prairie E")

plot(sp7E, ci.type="polygon", col="black", lwd=2, ci.lty=0, xlab = "observations", ylab = "number of species",ylim = c(0,20), ci=2, ci.col=alpha("blue", 0.1), main="Pre-AgE")
plot(sp8E, ci.type="polygon", col="black", lwd=2, ci.lty=0, xlab = "observations", ylab = "number of species",ylim = c(0,20), ci=2, ci.col=alpha("blue", 0.1), main="Event-AgE")
plot(sp9E, ci.type="polygon", col="black", lwd=2, ci.lty=0, xlab = "observations", ylab = "number of species",ylim = c(0,20), ci=2, ci.col=alpha("blue", 0.1), main="Post-AgE")

# BARNY SUGGESTS PRE / POST --> predicted Hf compared to alpha diversity / relative vocal indicator'
## PRE CUMULATIVE ALPHA / CUMULATIVE RVI
x3 <- c(.84, .68, .76)
y3 <- c(.813, 0.688, 0.167)
wilcox.test(x3, y3)

'data:  x3 and y3
W = 7, p-value = 0.4
alternative hypothesis: true location shift is not equal to 0'

## POST
x4 <- c(.9, .583, .667)
y4 <- c(.813, 0.688, 0.167)
wilcox.test(x4, y4)

'data:  x4 and y4
W = 5, p-value = 1
alternative hypothesis: true location shift is not equal to 0'

###### SE FOR EACH PRE / POST PER SITE -- repeat for each / ALPHA
# NIGHT SAMPLE - notice 42 rows (or 39 species)
std_err(as.matrix(colSums(post.pra[,3:42])))
mean(as.matrix(colSums(post.pra[,3:42])))

#AFTERNOON SAMPLE - notice 51 row (or 48 species)
std_err(as.matrix(colSums(post.praE[,3:51])))
mean(as.matrix(colSums(post.pra[,3:51])))

#### SE FOR EACH PRE/ POST PER SITE --- repeat for each  / RVI
#9219 RVI N column 3,30 pre.rip, post.rip
pre - se 1 mean 4
post - se 1 mean 1
#9219 RVI A column 11,12,14,18,23,24,25,27,33,34,37,38,40,44,51 preE.rip, postE.rip
pre - se 0.37, mean 2.8
post - se 0.36, mean 2.13
#15106 RVI N 20,23,24,28,30,33,36
pre se 5.85, mean 14.86
post se 5.80, mean 11
#15106 RVI A 14,18,20,25,30,41,44
pre se 0.42 mean 2.29
post se 0.18 mean 3.29
#15116 RVI N 20, 23, 24, 33
pre se 8.62 mean 22.5
post se 9.15, mean 17.5
#15116 RVI A 18,20, 24, 25, 30, 37, 41, 42, 44
pre se 0.41,  mean 2
post se 0.42, mean 2.89


