
rm(list=ls())

library(lme4)
library(yarg)
library(roquefort)
library(influence.ME)

# load functions
setwd("R:/ecocon_d/clg32/GitHub/PREDICTS_WDPA")
source("compare_randoms.R")
source("model_select.R")

source("prep_PA_11_14_for_analysis.R")




### model abundance

# check polynomials

fF <- c("Within_PA") 
fT <- list("ag_suit" = "3", "log_slope" = "3", "log_elevation" = "3")
keepVars <- list()
fI <- character(0)
RS <-  c("Within_PA")

abundance.best.random <- compare_randoms(PA_11_14, "log_abundance",
				fixedFactors=fF,
                        fixedTerms=fT,
			     	keepVars = keepVars,
                       	fixedInteractions=fI,
                        otherRandoms=character(0),
				fixed_RandomSlopes = RS,
                        fitInteractions=FALSE,
				verbose=TRUE)

abundance.best.random$best.random #"(1+Within_PA|SS)+ (1|SSB)"


ab.model <- model_select(all.data  = PA_11_14, 
			     responseVar = "log_abundance", 
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
			     keepVars = keepVars,
                       randomStruct = "(Within_PA|SS) + (1|SSB)",
			     otherRandoms=character(0),
                       verbose=TRUE)

ab.model$final.call
#log_abundance~Within_PA+(Within_PA|SS)+(1|SSB)

m1a <- lmer(log_abundance ~ Within_PA 
	+ (Within_PA|SS)+ (1|SSB), 
	 data =  PA_11_14)
m2a <- lmer(log_abundance ~ 1
	+ (Within_PA|SS)+ (1|SSB), 
	 data =  PA_11_14)
anova(m1a, m2a)

fixef(m1a)
fixef(m1a_)
summary(m1a)
exp(fixef(m1a)[2]) # 1.147

#no need to redo as no non-linear terms
fixef(m1a)

# plot
tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/01_15/simple model abundance.tif",
	width = 10, height = 10, units = "cm", pointsize = 12, res = 300)


labels <- c("Unprotected", "Protected")
y <- as.numeric(fixef(m1a)[2])
se <- as.numeric(se.fixef(m1a)[2])
yplus <- y + se*1.96
yminus <- y - se*1.96
y <-(exp(y)*100)
yplus<-(exp(yplus)*100)
yminus<-(exp(yminus)*100)

points <- c(100, y)
CI <- c(yplus, yminus)

plot(points ~ c(1,2), ylim = c(80,150), xlim = c(0.5,2.5), 
	bty = "l", pch = 16, col = c(1,3), cex = 1.5,
	yaxt = "n", xaxt = "n",
	ylab = "Abundance difference (% ± 95%CI)",
	xlab = "")
axis(1, c(1,2), labels)
axis(2, c(80,100,120,140), c(80,100,120,140))
arrows(2,CI[1],2,CI[2], code = 3, length = 0.03, angle = 90)
abline(h = 100, lty = 2)
points(points ~ c(1,2), pch = 16, col = c(1,3), cex = 1.5)

data <- PA_11_14[,c("Within_PA", "SSS", "log_abundance")]
data <- na.omit(data)
text(2,80, paste("n =", length(data$SSS[which(data$Within_PA == "yes")]), sep = " "))
text(1,80, paste("n =", length(data$SSS[which(data$Within_PA == "no")]), sep = " "))

#get details for master plot
ab.plot1 <- data.frame(label = c("unprotected", "all protected"), est = points, 
		upper = c(100, CI[1]), lower = c(100,CI[2]),
		n.site = c(length(data$SSS[which(data$Within_PA == "no")]), length(data$SSS[which(data$Within_PA == "yes")])))


dev.off()





# model abundance and IUCN_category
fF <- c("IUCN_CAT") 
fT <- list("ag_suit" = "1", "log_slope" = "1", "log_elevation" = "1")
keepVars <- list()
fI <- character(0)
RS <-  c("IUCN_CAT")
# cant converge with nonlinear terms

abundance.best.random.IUCN <- compare_randoms(PA_11_14, "log_abundance",
				fixedFactors=fF,
                        fixedTerms=fT,
			     	keepVars = keepVars,
                       	fixedInteractions=fI,
                        otherRandoms=character(0),
				fixed_RandomSlopes = RS,
                        fitInteractions=FALSE,
				verbose=TRUE)

abundance.best.random.IUCN$best.random #"(1+IUCN_CAT|SS)+ (1|SSB)"

ab.model.IUCN <- model_select(all.data  = PA_11_14, 
			     responseVar = "log_abundance", 
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
			     keepVars = keepVars,
                       randomStruct = "(IUCN_CAT|SS) + (1|SSB)",
			     otherRandoms=character(0),
                       verbose=TRUE)

ab.model.IUCN$final.call
#"log_abundance~1+(IUCN_CAT|SS)+(1|SSB)"


#PA_11_14$IUCN_CAT <- relevel(PA_11_14$IUCN_CAT, "4.5")

m3ai <- lmer(log_abundance ~ 1 
	+ (IUCN_CAT|SS)+ (1|SSB), 
	 data = PA_11_14,
	control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
m4ai <- lmer(log_abundance ~ IUCN_CAT 
	+ (IUCN_CAT|SS)+ (1|SSB), 
	 data = PA_11_14,
	control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))


#with ordinal - gives exactly the same results when compared ord vs not ordinal with other responses 
#but coefficients v different
m3ai <- lmer(log_abundance ~ 1 +log_slope + log_elevation + ag_suit
	+ (IUCN_CAT|SS)+ (1|SSB), 
	 data = PA_11_14_ord)
m4ai <- lmer(log_abundance ~ IUCN_CAT +log_slope + log_elevation + ag_suit
	+ (IUCN_CAT|SS)+ (1|SSB), 
	 data = PA_11_14_ord)


anova(m3ai, m4ai)
#3.6585      3     0.3008
summary(m4ai)


# get p values for individual comparisons

studies <- unique(PA_11_14$SS[which(PA_11_14$IUCN_CAT == "1.5")])
o.pos <- which(PA_11_14$Within_PA == "no" & PA_11_14$SS %in% studies)
i.pos <- which(PA_11_14$IUCN_CAT == "1.5")
data1 <- droplevels(PA_11_14[c(o.pos, i.pos),])

m3ai <- lmer(log_abundance ~ 1 
	+ (IUCN_CAT|SS)+ (1|SSB), 
	 data = data1,
	control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
m4ai <- lmer(log_abundance ~ IUCN_CAT 
	+ (IUCN_CAT|SS)+ (1|SSB), 
	 data = data1,
	control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
anova(m3ai, m4ai)

studies <- unique(PA_11_14$SS[which(PA_11_14$IUCN_CAT == "4.5")])
o.pos <- which(PA_11_14$Within_PA == "no" & PA_11_14$SS %in% studies)
i.pos <- which(PA_11_14$IUCN_CAT == "4.5")
data1 <- droplevels(PA_11_14[c(o.pos, i.pos),])

m3ai <- lmer(log_abundance ~ 1 
	+ (IUCN_CAT|SS)+ (1|SSB), 
	 data = data1,
	control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
m4ai <- lmer(log_abundance ~ IUCN_CAT 
	+ (IUCN_CAT|SS)+ (1|SSB), 
	 data = data1,
	control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
anova(m3ai, m4ai)


studies <- unique(PA_11_14$SS[which(PA_11_14$IUCN_CAT == "7")])
o.pos <- which(PA_11_14$Within_PA == "no" & PA_11_14$SS %in% studies)
i.pos <- which(PA_11_14$IUCN_CAT == "7")
data1 <- droplevels(PA_11_14[c(o.pos, i.pos),])

m3ai <- lmer(log_abundance ~ 1 
	+ (IUCN_CAT|SS)+ (1|SSB), 
	 data = data1,
	control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
m4ai <- lmer(log_abundance ~ IUCN_CAT 
	+ (IUCN_CAT|SS)+ (1|SSB), 
	 data = data1,
	control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
anova(m3ai, m4ai)



# plot 

tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/01_15/simple model abundance IUCN.tif",
	width = 15, height = 12, units = "cm", pointsize = 12, res = 300)


labels <- c("Unprotected", "IUCN III  - VI", "unknown", "IUCN I & II")

PA_11_14$IUCN_CAT <- relevel(PA_11_14$IUCN_CAT, "0")

m4ai <- lmer(log_abundance ~ IUCN_CAT 
	+ (IUCN_CAT|SS) + (1|SSB), 
	 data = PA_11_14)
summary(m4ai)


pos <- c(grep("4.5", names(fixef(m4ai))),grep("7", names(fixef(m4ai))),grep("1.5", names(fixef(m4ai))))
y <- as.numeric(fixef(m4ai)[pos])
se <- as.numeric(se.fixef(m4ai)[pos])
yplus <- y + se*1.96
yminus <- y - se*1.96
y <-(exp(y)*100)
yplus<-(exp(yplus)*100)
yminus<-(exp(yminus)*100)

points <- c(100, y)
CI <- cbind(yplus, yminus)

plot(points ~ c(1,2,3,4), ylim = c(80,180), xlim = c(0.5,4.5),
	bty = "l", pch = 16, col = c(1,3,3,3), cex = 1.5,
	yaxt = "n", xaxt = "n",
	ylab = "Abundance difference (% ± 95%CI)",
	xlab = "")
axis(1,seq(1,length(points),1), labels)
axis(2, c(80,100,120,140,160), c(80,100,120,140,160))
arrows(seq(2,length(points),1),CI[,1],
	seq(2,length(points),1),CI[,2], code = 3, length = 0.03, angle = 90)
abline(h = 100, lty = 2)
points(points ~ c(1,2,3,4), pch = 16, col = c(1,3,3,3), cex = 1.5)

data <- PA_11_14[,c("IUCN_CAT", "SSS", "log_abundance")]
data <- na.omit(data)
text(1, 80, paste("n =", length(data$SSS[which(data$IUCN_CAT == "0")]), sep = " "))
text(2, 80, paste("n =", length(data$SSS[which(data$IUCN_CAT == "4.5")]), sep = " "))
text(3, 80, paste("n =", length(data$SSS[which(data$IUCN_CAT == "7")]), sep = " "))
text(4, 80, paste("n =", length(data$SSS[which(data$IUCN_CAT == "1.5")]), sep = " "))

#add to master plot details

IUCN.plot <- data.frame(label = labels[2:4], est = points[2:4], 
		upper = CI[,1], lower = CI[,2],
		n.site = c(length(data$SSS[which(data$IUCN_CAT == "4.5")]), 
			length(data$SSS[which(data$IUCN_CAT == "7")]),
			length(data$SSS[which(data$IUCN_CAT == "1.5")])))
ab.plot2 <- rbind(ab.plot1, IUCN.plot)


dev.off()



# abundance and zone

tropical <- subset(PA_11_14, Zone == "Tropical")
temperate <- subset(PA_11_14, Zone == "Temperate")

# check polynomials for confounding variables
fF <- c("Within_PA" ) 
fT <- list("ag_suit" = "3", "log_slope" = "3", "log_elevation" = "3")
keepVars <- list()
fI <- character(0)
RS <-  c("Within_PA")


ab.best.random.trop <- compare_randoms(tropical, "log_abundance",
				fixedFactors=fF,
                        fixedTerms=fT,
			     	keepVars = keepVars,
                       	fixedInteractions=fI,
                        otherRandoms=character(0),
				fixed_RandomSlopes = RS,
                        fitInteractions=FALSE,
				verbose=TRUE)

ab.best.random.trop$best.random #"(1+Within_PA|SS)+ (1|SSB)"

ab.best.random.temp <- compare_randoms(temperate, "log_abundance",
				fixedFactors=fF,
                        fixedTerms=fT,
			     	keepVars = keepVars,
                       	fixedInteractions=fI,
                        otherRandoms=character(0),
				fixed_RandomSlopes = RS,
                        fitInteractions=FALSE,
				verbose=TRUE)

ab.best.random.temp$best.random #"(1+Within_PA|SS)+ (1|SSB)"

# get polynomial relationships
ab.model.trop <- model_select(all.data  = tropical, 
			     responseVar = "log_abundance", 
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
			     keepVars = keepVars,
                       randomStruct =ab.best.random.trop$best.random,
			     otherRandoms=character(0),
                       verbose=TRUE)
#log_abundance~poly(ag_suit,1)+poly(log_elevation,2)+Within_PA+(1+Within_PA|SS)+(1|SSB)"

ab.model.temp <- model_select(all.data  = temperate, 
			     responseVar = "log_abundance", 
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
			     keepVars = keepVars,
                       randomStruct = ab.best.random.temp$best.random,
			     otherRandoms=character(0),
                       verbose=TRUE)
#log_abundance~poly(log_slope,2)+(1+Within_PA|SS)+(1|SSB)"


# run models
data.trop <- tropical[,c("Within_PA", "ag_suit", "log_elevation", "log_slope", "SS", "SSB", "log_abundance")]
data.trop <- na.omit(data.trop)
m1aztr <- lmer(log_abundance ~ Within_PA + poly(log_elevation,2)+ ag_suit
	+ (Within_PA|SS)+ (1|SSB), 
	 data = tropical)
m2aztr <- lmer(log_abundance ~ 1 + poly(log_elevation,2) + ag_suit
	+ (Within_PA|SS)+ (1|SSB), 
	 data = tropical)
anova(m1aztr, m2aztr)
#5.5455      1    0.01853

data.temp <- temperate[,c("Within_PA", "ag_suit", "log_elevation", "log_slope", "SS", "SSB", "log_abundance")]
data.temp <- na.omit(data.temp)
m1azte <- lmer(log_abundance ~ Within_PA +poly(log_slope,2) 
	+ (Within_PA|SS)+ (1|SSB), 
	 data = temperate)
m2azte <- lmer(log_abundance ~ 1 +poly(log_slope,2) 
	+ (Within_PA|SS)+ (1|SSB), 
	 data = temperate)
anova(m1azte, m2azte)
#0.2912      1     0.5895


#add results to master plot
aztr.est <- exp(fixef(m1aztr)[2])*100
aztr.upper <- exp(fixef(m1aztr)[2] + 1.96* se.fixef(m1aztr)[2])*100
aztr.lower <- exp(fixef(m1aztr)[2] - 1.96* se.fixef(m1aztr)[2])*100

azte.est <- exp(fixef(m1azte)[2])*100
azte.upper <- exp(fixef(m1azte)[2] + 1.96* se.fixef(m1azte)[2])*100
azte.lower <- exp(fixef(m1azte)[2] - 1.96* se.fixef(m1azte)[2])*100

a.zone <- data.frame(label = c("Tropical", "Temperate"),
				est = c(aztr.est, azte.est), 
				upper = c(aztr.upper, azte.upper), 
				lower = c(aztr.lower, azte.lower), 
				n.site = c(nrow(data.trop[which(data.trop$Within_PA == "yes"),]),
					 nrow(data.temp[which(data.temp$Within_PA == "yes"),])))
ab.plot3 <- rbind(ab.plot2, a.zone)





# abundance and taxon

plants <- subset(PA_11_14, taxon_of_interest == "Plants")
inverts <- subset(PA_11_14, taxon_of_interest == "Invertebrates")
verts <- subset(PA_11_14, taxon_of_interest == "Vertebrates")
nrow(plants)
nrow(inverts)
nrow(verts)

# check polynomials for confounding variables
fF <- c("Within_PA" ) 
fT <- list("ag_suit" = "3", "log_slope" = "3", "log_elevation" = "3")
keepVars <- list()
fI <- character(0)
RS <-  c("Within_PA")


ab.best.random.p <- compare_randoms(plants, "log_abundance",
				fixedFactors=fF,
                        fixedTerms=fT,
			     	keepVars = keepVars,
                       	fixedInteractions=fI,
                        otherRandoms=character(0),
				fixed_RandomSlopes = RS,
                        fitInteractions=FALSE,
				verbose=TRUE)
ab.best.random.p$best.random # "(1+Within_PA|SS)+ (1|SSB)"

ab.best.random.i <- compare_randoms(inverts, "log_abundance",
				fixedFactors=fF,
                        fixedTerms=fT,
			     	keepVars = keepVars,
                       	fixedInteractions=fI,
                        otherRandoms=character(0),
				fixed_RandomSlopes = RS,
                        fitInteractions=FALSE,
				verbose=TRUE)
ab.best.random.i$best.random # "(1+Within_PA|SS)+ (1|SSB)"

ab.best.random.v <- compare_randoms(verts, "log_abundance",
				fixedFactors=fF,
                        fixedTerms=fT,
			     	keepVars = keepVars,
                       	fixedInteractions=fI,
                        otherRandoms=character(0),
				fixed_RandomSlopes = RS,
                        fitInteractions=FALSE,
				verbose=TRUE)
ab.best.random.v$best.random #"(1+Within_PA|SS)+ (1|SSBS)+ (1|SSB)"





# get polynomial relationships
ab.model.p <- model_select(all.data  = plants, 
			     responseVar = "log_abundance", 
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
			     keepVars = keepVars,
                       randomStruct =ab.best.random.p$best.random,
			     otherRandoms=character(0),
                       verbose=TRUE)
#log_abundance~poly(ag_suit,2)+poly(log_elevation,1)+poly(log_slope,3)+(1+Within_PA|SS)

ab.model.i <- model_select(all.data  = inverts, 
			     responseVar = "log_abundance", 
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
			     keepVars = keepVars,
                       randomStruct =ab.best.random.i$best.random,
			     otherRandoms=character(0),
                       verbose=TRUE)
#"log_abundance~poly(log_slope,1)+(1+Within_PA|SS)+(1|SSB)"


ab.model.v <- model_select(all.data  = verts, 
			     responseVar = "log_abundance", 
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
			     keepVars = keepVars,
                       randomStruct =ab.best.random.v$best.random,
			     otherRandoms=character(0),
                       verbose=TRUE)
#"log_abundance~poly(ag_suit,3)+poly(log_slope,1)+(1+Within_PA|SS)"


# run models
data.p <- plants[,c("Within_PA", "ag_suit", "log_elevation", "log_slope", "SS", "SSB", "log_abundance")]
data.p <- na.omit(data.p)
m1txp <- lmer(log_abundance ~ Within_PA +poly(ag_suit,2)+poly(log_elevation,1)+poly(log_slope,3)
	+ (Within_PA|SS)+ (1|SSB), 
	 data = data.p)
m2txp <- lmer(log_abundance ~ 1+poly(ag_suit,2)+poly(log_elevation,1)+poly(log_slope,3)
	+ (Within_PA|SS)+ (1|SSB), 
	 data = data.p)
anova(m1txp , m2txp)
# 0.0539      1,13     0.8164

data.i <- inverts[,c("Within_PA", "ag_suit", "log_elevation", "log_slope", "SS", "SSB", "log_abundance")]
data.i <- na.omit(data.i)
m1txi <- lmer(log_abundance ~ Within_PA +log_slope 
	+ (Within_PA|SS)+ (1|SSB), 
	 data = data.i)
m2txi<- lmer(log_abundance ~ 1 +log_slope
	+ (Within_PA|SS)+ (1|SSB), 
	 data = data.i)
anova(m1txi, m2txi)
#2.1433      1     0.1432

data.v <- verts[,c("Within_PA", "ag_suit", "log_elevation", "log_slope", "SS", "SSB", "log_abundance")]
data.v <- na.omit(data.v)
m1txv <- lmer(log_abundance ~ Within_PA +  poly(ag_suit,3)+poly(log_slope,1)
	+ (Within_PA|SS)+ (1|SSB), 
	 data = data.v)
m2txv <- lmer(log_abundance ~ 1 +  poly(ag_suit,3)+poly(log_slope,1)
	+ (Within_PA|SS)+ (1|SSB), 
	 data = data.v)
anova(m1txv, m2txv)
#5.1244      1    0.02359

#add results to master plot
txp.est <- exp(fixef(m1txp)[2])*100
txp.upper <- exp(fixef(m1txp)[2] + 1.96* se.fixef(m1txp)[2])*100
txp.lower <- exp(fixef(m1txp)[2] - 1.96* se.fixef(m1txp)[2])*100

txi.est <- exp(fixef(m1txi)[2])*100
txi.upper <- exp(fixef(m1txi)[2] + 1.96* se.fixef(m1txi)[2])*100
txi.lower <- exp(fixef(m1txi)[2] - 1.96* se.fixef(m1txi)[2])*100

txv.est <- exp(fixef(m1txv)[2])*100
txv.upper <- exp(fixef(m1txv)[2] + 1.96* se.fixef(m1txv)[2])*100
txv.lower <- exp(fixef(m1txv)[2] - 1.96* se.fixef(m1txv)[2])*100



a.tax <- data.frame(label = c("Plants", "Inverts", "Verts"),
				est = c(txp.est, txi.est, txv.est), 
				upper = c(txp.upper, txi.upper, txv.upper), 
				lower = c(txp.lower, txi.lower, txv.lower), 
				n.site = c(nrow(data.p[which(data.p$Within_PA == "yes"),]), 
					nrow(data.i[which(data.i$Within_PA == "yes"),]), 
					nrow(data.v[which(data.v$Within_PA == "yes"),])))
ab.plot <- rbind(ab.plot3, a.tax)


# master plot

load("\\\\smbhome.uscs.susx.ac.uk\\clg32\\Documents\\PREDICTS\\WDPA analysis\\RData files\\8 landuses\\simple models - abundance.RData")


tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/02_15/simple models abundance.tif",
	width = 23, height = 16, units = "cm", pointsize = 12, res = 300)

trop.col <- rgb(0.9,0,0)
temp.col <- rgb(0,0.1,0.7)
p.col <- rgb(0.2,0.7,0.2)
i.col <- rgb(0,0.7,0.9)
v.col <- rgb(0.9,0.5,0)

par(mar = c(9,6,4,1))
plot(1,1, 
	ylim = c(65,190), xlim = c(0.5,nrow(ab.plot)+1),
	bty = "l", 
	axes = F,
	ylab = "Abundance difference (%)",
	cex.lab = 1.5,
	xlab = "")
arrows(1:nrow(ab.plot),ab.plot$upper,
	col = c(1,1,rep(rgb(0.5, 0.5, 0.5), 3), c(trop.col, temp.col, p.col, i.col, v.col)),
	lwd = 2,
	1:nrow(ab.plot),ab.plot$lower, code = 3, length = 0, angle = 90)
points(ab.plot$est ~ c(1:nrow(ab.plot)),
	pch = c(21, rep(16,4), rep(15,2),rep(17,3)), 
	lwd = 2,
	col = c(1,1,rep(rgb(0.5, 0.5, 0.5), 3), c(trop.col, temp.col, p.col, i.col, v.col)),
	bg = c("white"), 
	cex = 1.5)
abline(v = c(2.5,5.5,7.5), col = 8)
abline(h= 100, lty = 2)
text(1:nrow(ab.plot),65, ab.plot$n.site, srt = 180)
#axis(1, c(1:nrow(ab.plot)), ab.plot$label, cex.axis = 1.5, las = 2)
axis(2, c(80,100,120,140,160,180), c(80,100,120,140,160,180))


dev.off()



### effectiveness estimates


b.a <- exp(fixef(m1a)[2]) -1
b.a.max <- exp(fixef(m1a)[2] + 1.96*se.fixef(m1a)[2])-1
b.a.min <- exp(fixef(m1a)[2] - 1.96*se.fixef(m1a)[2])-1


# if IUCN cat 1 or 2
pos <- which(names(fixef(m4ai))== "IUCN_CAT1.5")
b.aIUCN1 <- exp(fixef(m4ai)[pos])-1
b.aIUCN1.max <- exp(fixef(m4ai)[pos] + 1.96*se.fixef(m4ai)[2])-1
b.aIUCN1.min <- exp(fixef(m4ai)[pos] - 1.96*se.fixef(m4ai)[2])-1

# if IUCN cat  3 to 6
pos <- which(names(fixef(m4ai))== "IUCN_CAT4.5")
b.aIUCN2 <- exp(fixef(m4ai)[pos])-1
b.aIUCN2.max <- exp(fixef(m4ai)[pos] + 1.96*se.fixef(m4ai)[2])-1
b.aIUCN2.min <- exp(fixef(m4ai)[pos] - 1.96*se.fixef(m4ai)[2])-1


b.vals <- c(b.a, b.a.max, b.a.min,
		b.aIUCN1, b.aIUCN1.max, b.aIUCN1.min,
		b.aIUCN2, b.aIUCN2.max, b.aIUCN2.min)

vals <- data.frame(name =  c("b.a", "b.a.max", "b.a.min",
				"b.aIUCN1", "b.aIUCN1.max", "b.aIUCN1.min",
				"b.aIUCN2", "b.aIUCN2.max", "b.aIUCN2.min"),
			b.vals = b.vals, 
			metric = rep("abundance", 9), 
			NPA.abs = rep(NA,length(b.vals)),
			PA.abs = rep(NA,length(b.vals)),
			est = rep(NA,length(b.vals)))


#By 2005, we estimate that human impacts had reduced local richness by an average of 13.6% (95% CI: 9.1 – 17.8%) and 
#total abundance by 10.7% (95% CI: 3.8% gain – 23.7% reduction) compared with pre-impact times. 
#Approximately 60% of the decline in richness was independent of effects on abundance: 
#average rarefied richness has fallen by 8.1% (95% CI: 3.5 – 12.9%).
			
b <- b.vals[1]

for(b in vals$b.vals){

benefit <- as.numeric(b)		# percentage increase in metric in PAs
PA.pct <- 14.7 				# percentage of total land area in PAs
# global loss of biodiversity (from Newbold et al)
global.loss <- NULL
if(vals$metric[which(vals$b.vals == b)] == "sp.rich"){
	global.loss <- 0.136			
	}else if (vals$metric[which(vals$b.vals == b)] == "abundance"){
	global.loss <- 0.107
	}else if (vals$metric[which(vals$b.vals == b)] == "rar.rich"){
	global.loss <- 0.081
	}

NPA.rel <- 1-benefit			# relative biodiversity in unprotected sites
PA.rel <- 1					# biodiversity in protected sites
NPA.pct <- 100-PA.pct			# land area unprotected 
global.int <- 1-global.loss		# overall status of biodiversity relative to pristine


# we want NPA.abs and PA.abs - where these are the biodiv metrics in unprotected and protected relative to pristine
# simultaneous equations are

# global.int = NPA.pct/100*NPA.abs + PA.pct/100*PA.abs #overall loss is loss in PAs and nonPAs relative to pristine
# NPA.abs = PA.abs*NPA.rel			     


# so
#global.int = (NPA.pct/100)*PA.abs*NPA.rel + (PA.pct/100)*PA.abs
# which is the same as
#global.int = PA.abs*(NPA.pct/100*(NPA.rel) + PA.pct/100)

# then
PA.abs <- global.int/(PA.pct/100 + (NPA.pct/100 * (NPA.rel)))
NPA.abs <- PA.abs*NPA.rel
vals$NPA.abs[which(vals$b.vals == b)] <- NPA.abs
vals$PA.abs[which(vals$b.vals == b)] <- PA.abs

# if pristine is 1, then where between 1 and NPA.abs does PA.abs fall?
# get difference between PA.abs and NPA.abs as a percentage of NPA.abs
#((1-NPA.abs)-(1-PA.abs))/(1-NPA.abs)

est <- 1-(1-PA.abs)/(1-NPA.abs)
vals$est[which(vals$b.vals == b)] <- est

# ie
#est <- ((1-NPA.abs)-(1-PA.abs))/(1-NPA.abs)
}

vals


# need to work out, if all PAs were IUCN cat 1 and still 14.7%, instead of what we have
# what would global loss be
# then how much area do we need to get this, given current benefit of PAs

# NPA abs remains the same as in normal benefit scenario
NPA.abs <- vals$NPA.abs[which(vals$name == "b.a")]
# total area protected remains the same
PA.pct <- 14.7
# NPA.rel becomes 1-benefit of protection in IUCN category
pos <- which(names(fixef(m4ai))== "IUCN_CAT1.5")
NPA.rel <- 1 - (exp(fixef(m4ai)[pos])-1)

#given
#global.int = (1 - PA.pct/100)*NPA.abs + (PA.pct/100)*PA.abs
global.int = (1 - PA.pct/100)*NPA.abs + (PA.pct/100)*(NPA.abs/NPA.rel)

#now, given the normal benefit, how much additional PA area needed to get this
# need area in terms of global.int and NPA abs. 
#global.int = (1 - PA.pct/100)*NPA.abs + (PA.pct/100)*NPA.abs/NPA.rel
#global.int == NPA.abs - (PA.pct/100)*NPA.abs  + (PA.pct/100)*NPA.abs/NPA.rel
#global.int = NPA.abs - (PA.pct/100)*(NPA.abs - NPA.abs/NPA.rel)
#(PA.pct/100) = (NPA.abs - global.int)/(NPA.abs - NPA.abs/NPA.rel)

NPA.rel <- 1 - (exp(fixef(m1a)[2]) -1)
100*(NPA.abs - global.int)/(NPA.abs - NPA.abs/NPA.rel)







