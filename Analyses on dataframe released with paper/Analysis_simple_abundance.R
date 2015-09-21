
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

### load saved model objects
#load("R:\\ecocon_d\\clg32\\PREDICTS\\WDPA analysis\\RData files\\8 landuses\\simple models - abundance.RData")


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
                       randomStruct = abundance.best.random$best.random,
			     otherRandoms=character(0),
                       verbose=TRUE)
ab.model$stats
ab.model$final.call
#log_abundance~Within_PA+(Within_PA|SS)+(1|SSB)


data <- ab.model$data
m1a <- lmer(log_abundance ~ Within_PA + (Within_PA|SS)+ (1|SSB), data =  data)
m2a <- lmer(log_abundance ~ 1+ (Within_PA|SS)+ (1|SSB), data =  data)
anova(m1a, m2a)

# comparison I had previously
m3a <- lmer(log_abundance ~ Within_PA + (Within_PA|SS)+ (1|SSB), data =  PA_11_14)
m4a <- lmer(log_abundance ~ 1+ (Within_PA|SS)+ (1|SSB), data =  PA_11_14)
anova(m3a, m4a)
### omitting incomplete values does make a difference


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
                       randomStruct = abundance.best.random.IUCN$best.random,
			     otherRandoms=character(0),
                       verbose=TRUE)
ab.model.IUCN$stats
ab.model.IUCN$final.call
ab.model.IUCN$warnings
#"log_abundance~1+(IUCN_CAT|SS)+(1|SSB)"





# plot 

tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/01_15/simple model abundance IUCN.tif",
	width = 15, height = 12, units = "cm", pointsize = 12, res = 300)


labels <- c("Unprotected", "IUCN III  - VI", "unknown", "IUCN I & II")

data <-  ab.model.IUCN$data
data$IUCN_CAT <- relevel(data$IUCN_CAT, "0")

m4ai <- lmer(log_abundance ~ IUCN_CAT + (IUCN_CAT|SS) + (1|SSB),  data)
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
ab.model.trop$stats
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
ab.model.temp$stats
#log_abundance~poly(log_slope,2)+(1+Within_PA|SS)+(1|SSB)"


# run models

data.trop <- ab.model.trop$data
data.temp <- ab.model.temp$data

m1aztr <- lmer(log_abundance ~ Within_PA + poly(log_elevation,2)+ ag_suit
	+ (Within_PA|SS)+ (1|SSB), 
	 data = data.trop)

m1azte <- lmer(log_abundance ~ Within_PA +poly(log_slope,2) 
	+ (Within_PA|SS)+ (1|SSB), 
	 data = data.temp)


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
ab.model.p$stats
ab.model.p$final.call
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
ab.model.i$stats
ab.model.i$final.call
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
ab.model.v$stats
ab.model.v$final.call
#log_abundance~poly(ag_suit,3)+poly(log_slope,2)+Within_PA+(1+Within_PA|SS)+(1|SSB)"

data.p <- ab.model.p$data
data.i <- ab.model.i$data
data.v <- ab.model.v$data

# run models
m1txp <- lmer(log_abundance ~ Within_PA +poly(ag_suit,2)+poly(log_elevation,1)+poly(log_slope,3)
	+ (Within_PA|SS)+ (1|SSB), 
	 data = data.p)

m1txi <- lmer(log_abundance ~ Within_PA +log_slope 
	+ (Within_PA|SS)+ (1|SSB), 
	 data = data.i)

m1txv <- lmer(log_abundance ~ Within_PA +  poly(ag_suit,3)+poly(log_slope,2)
	+ (Within_PA|SS)+ (1|SSB), 
	 data = data.v)


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

ab.plot$truncated <- ab.plot$upper
ab.plot$truncated[which(ab.plot$upper > 145)] <- 135
ab.plot$dotlines <- rep(0,nrow(ab.plot))
ab.plot$dotlines_end <- rep(0,nrow(ab.plot))
ab.plot$dotlines[which(ab.plot$upper > 145)] <- 135
ab.plot$dotlines_end[which(ab.plot$upper > 145)] <- 140

# master plot

tiff( "R:/ecocon_d/clg32/PREDICTS/WDPA analysis/plots/02_15/simple models abundance.tif",
	width = 10, height = 15, units = "cm", pointsize = 12, res = 300)

trop.col <- rgb(0.9,0,0)
temp.col <- rgb(0,0.1,0.7)
p.col <- rgb(0.2,0.7,0.2)
i.col <- rgb(0,0.7,0.9)
v.col <- rgb(0.9,0.5,0)

par(mar = c(9,6,4,1))
plot(1,1, 
	ylim = c(70,145), xlim = c(0.5,nrow(ab.plot)),
	bty = "l", 
	axes = F,
	ylab = "Abundance difference (%)",
	cex.lab = 1.5,
	xlab = "")
arrows(1:nrow(ab.plot),ab.plot$truncated,1:nrow(ab.plot),ab.plot$lower,
	col = c(1,1,rep(rgb(0.5, 0.5, 0.5), 3), c(trop.col, temp.col, p.col, i.col, v.col)),
	lwd = 2,
	code = 3, length = 0, angle = 90)
arrows(1:nrow(ab.plot),ab.plot$dotlines, 1:nrow(ab.plot),ab.plot$dotlines_end,
	col = c(1,1,rep(rgb(0.5, 0.5, 0.5), 3), c(trop.col, temp.col, p.col, i.col, v.col)),
	lwd = 2, lty = 3,
	code = 3, length = 0, angle = 90)
points(ab.plot$est ~ c(1:nrow(ab.plot)),
	pch = c(21, rep(16,4), rep(15,2),rep(17,3)), 
	lwd = 2, lty = 3,
	col = c(1,1,rep(rgb(0.5, 0.5, 0.5), 3), c(trop.col, temp.col, p.col, i.col, v.col)),
	bg = c("white"), 
	cex = 1.5)
abline(v = c(2.5,5.5,7.5), col = 8)
abline(h= 100, lty = 2)
text(1:nrow(ab.plot),72, ab.plot$n.site, srt = 90)
text(which(ab.plot$upper > 145),143, round(ab.plot$upper[which(ab.plot$upper > 145)]), srt = 90)
#axis(1, c(1:nrow(ab.plot)), ab.plot$label, cex.axis = 1.5, las = 2)
axis(2, c(80,100,120,140), c(80,100,120,140))


dev.off()


#new plot after text editing

ab.plot2 <- ab.plot[1:5,]

tiff( "R:/ecocon_d/clg32/PREDICTS/WDPA analysis/plots/06_15/simple models abundance.tif",
	width = 8, height = 15, units = "cm", pointsize = 12, res = 300)

par(mar = c(9,6,4,2))
plot(1,1, 
	ylim = c(80,145), xlim = c(0.5,nrow(ab.plot2)+0.1),
	bty = "l", 
	axes = F,
	ylab = "Abundance difference (%)",
	cex.lab = 1.5,
	xlab = "")
abline(v = 2.5, lty = 2)
abline(h= 100, lty = 1, col = 8)
arrows(1:nrow(ab.plot2),ab.plot2$truncated,1:nrow(ab.plot2),ab.plot2$lower,
	col = c(1,1,rep(rgb(0.5, 0.5, 0.5), 3)),
	lwd = 2,
	code = 3, length = 0, angle = 90)
arrows(1:nrow(ab.plot2),ab.plot2$dotlines, 1:nrow(ab.plot2),ab.plot2$dotlines_end,
	col = c(1,1,rep(rgb(0.5, 0.5, 0.5), 3)),
	lwd = 2, lty = 3,
	code = 3, length = 0, angle = 90)
points(ab.plot2$est ~ c(1:nrow(ab.plot2)),
	pch = c(21, rep(16,4), rep(15,2),rep(17,3)), 
	lwd = 2,
	col = c(1,1,rep(rgb(0.5, 0.5, 0.5), 3)),
	bg = "white", 
	cex = 1.5)
text(1:nrow(ab.plot2),82, ab.plot2$n.site, srt = 90, cex = 1)
text(which(ab.plot$upper > 145),143, round(ab.plot$upper[which(ab.plot$upper > 145)]), srt = 90, col = rgb(0.5, 0.5, 0.5))
axis(1, c(1:nrow(ab.plot2)), ab.plot2$label, cex.axis = 1.1, las = 2, tick = 0)
axis(2, c(90,100,110,120,130,140), c(-10,0,10,20,30,40))

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
			
b <- b.vals[1]

for(b in vals$b.vals){

benefit <- as.numeric(b)		# percentage increase in metric in PAs
PA.pct <- 15.4				# percentage of total land area in PAs
# global loss of biodiversity (from Newbold et al)
global.loss <- NULL
if(vals$metric[which(vals$b.vals == b)] == "sp.rich"){
	global.loss <- 0.129			
	}else if (vals$metric[which(vals$b.vals == b)] == "abundance"){
	global.loss <- 0.126
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


### working out area needed to achieve similar level of global biodiversity to increasing all existing to cat I-II
# method:
# need to work out, if all PAs were IUCN cat 1 and still same terrestrial %, instead of what we have
# what would global loss be
# then how much area do we need to get this, given current benefit of PAs

# NPA abs remains the same as in normal benefit scenario
NPA.abs <- vals$NPA.abs[which(vals$name == "b.a")]

# total area protected remains the same
PA.pct <- 15.4
# NPA.rel becomes 1-benefit of protection in IUCN category
#just use mean effect, as CI are so huge as to be uninformative 
pos <- which(names(fixef(m4ai))== "IUCN_CAT1.5")
NPA.rel <- 1 - (exp(fixef(m4ai)[pos])-1)

#given
#global.int = (1 - PA.pct/100)*NPA.abs + (PA.pct/100)*PA.abs
global.int = (1 - PA.pct/100)*NPA.abs + (PA.pct/100)*(NPA.abs/NPA.rel)

#now, given the normal benefit, how much additional PA area needed to get this
# need area in terms of global.int and NPA abs. 
#global.int = (1 - PA.pct/100)*NPA.abs + (PA.pct/100)*NPA.abs/NPA.rel
#global.int = NPA.abs - (PA.pct/100)*NPA.abs  + (PA.pct/100)*NPA.abs/NPA.rel
#global.int = NPA.abs - (PA.pct/100)*(NPA.abs - NPA.abs/NPA.rel)
#(PA.pct/100) = (NPA.abs - global.int)/(NPA.abs - NPA.abs/NPA.rel)

NPA.rel <- 1 - (exp(fixef(m1a)[2]) -1)
100*(NPA.abs - global.int)/(NPA.abs - NPA.abs/NPA.rel)
print("max")
NPA.rel <- 1 - (exp(fixef(m1a)[2]+2*se.fixef(m1a)[2]) -1)
100*(NPA.abs - global.int)/(NPA.abs - NPA.abs/NPA.rel)
print("min")
NPA.rel <- 1 - (exp(fixef(m1a)[2]-2*se.fixef(m1a)[2]) -1)
100*(NPA.abs - global.int)/(NPA.abs - NPA.abs/NPA.rel)

# what if all PAs are lower management restriction (i.e. new PAs only 3 - 6 and old PAs degraded)
# ie NPA.rel is only that for III to VI

pos <- which(names(fixef(m4ai))== "IUCN_CAT4.5")
NPA.rel <- 1 - (exp(fixef(m4ai)[pos])-1)
 100*(NPA.abs - global.int)/(NPA.abs - NPA.abs/NPA.rel)











### What if current/IUCN high/IUCN low level of protection and 17% of total area protected? 
# what global state of biodiversity do we get if increase to 17%
# taking absolute biodiversity within PAs as current 


#given
PA.abs <- vals[which(vals$name == "b.a"),"PA.abs"]
NPA.rel <- 1 - (exp(fixef(m1a)[2])-1)
PA.pct <- 17
#global.int = (1 - PA.pct/100)*NPA.abs + (PA.pct/100)*PA.abs
global.int.17current = (1 - PA.pct/100)*PA.abs*NPA.rel + (PA.pct/100)*PA.abs
global.int.17current



# to get intactness where IUCN 1 or 2 and 17%
# change equ from before to be 17% instead of 14.6%
# NPA abs remains the same as in normal benefit scenario
NPA.abs <- vals$NPA.abs[which(vals$name == "b.a")]
# total area protected remains the same
PA.pct <- 17
# NPA.rel becomes 1-benefit of protection in IUCN category
pos <- which(names(fixef(m4ai))== "IUCN_CAT1.5")
NPA.rel <- 1 - (exp(fixef(m4ai)[pos])-1)
global.int = (1 - PA.pct/100)*NPA.abs + (PA.pct/100)*(NPA.abs/NPA.rel)
global.int

# to get intactness where IUCN 3 to 6 and 17%
# NPA abs remains the same as in normal benefit scenario
NPA.abs <- vals$NPA.abs[which(vals$name == "b.a")]
# total area protected remains the same
PA.pct <- 17
# NPA.rel becomes 1-benefit of protection in IUCN category
pos <- which(names(fixef(m4ai))== "IUCN_CAT4.5")
NPA.rel <- 1 - (exp(fixef(m4ai)[pos])-1)
global.int = (1 - PA.pct/100)*NPA.abs + (PA.pct/100)*(NPA.abs/NPA.rel)
global.int

# to get intactness where as current and 17% (check matches first method above)
# NPA abs remains the same as in normal benefit scenario
NPA.abs <- vals$NPA.abs[which(vals$name == "b.a")]
# total area protected remains the same
PA.pct <- 17
# NPA.rel becomes 1-benefit of protection in IUCN category
NPA.rel <- 1 - (exp(fixef(m1a)[2])-1)
global.int = (1 - PA.pct/100)*NPA.abs + (PA.pct/100)*(NPA.abs/NPA.rel)
global.int


# how does this compare to the output of the nature paper models that we started with
global.loss <- 0.107
global.int <- 1-global.loss
global.int



