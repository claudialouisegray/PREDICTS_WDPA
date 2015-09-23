rm(list=ls()) 

library(yarg)
library(roquefort)
library(gamm4)


setwd("R:/ecocon_d/clg32/GitHub/PREDICTS_WDPA")
source("compare_randoms.R")
source("model_select.R")
setwd("R:/ecocon_d/clg32/PREDICTS/WDPA analysis")
PA_11_14 <- read.csv("PREDICTS_WDPA.csv")


size <- aggregate(SSS ~ SS, PA_11_14, length)
hist(size$SSS, breaks = 50, col = 8) # break after 5
length(which(size$SSS <= 5)) # 19 out of 167
length(which(size$SSS <= 10)) # 41 out of 167
length(which(size$SSS <= 25)) # 89 out of 167

to.keep <- size$SS[which(size$SSS > 10)]
no.small <- subset(PA_11_14, SS %in% to.keep)
nrow(no.small)


### model abundance

# check polynomials

fF <- c("Within_PA") 
fT <- list("ag_suit" = "3", "log_slope" = "3", "log_elevation" = "3")
keepVars <- list()
fI <- character(0)
RS <-  c("Within_PA")
#log_abundance~Within_PA+(Within_PA|SS)+(1|SSB)

abundance.best.random <- compare_randoms(no.small, "log_abundance",
				fixedFactors=fF,
        fixedTerms=fT,
        keepVars = keepVars,
        fixedInteractions=fI,
        otherRandoms=character(0),
				fixed_RandomSlopes = RS,
        fitInteractions=FALSE,
				verbose=TRUE)

abundance.best.random$best.random #"(1+Within_PA|SS)+ (1|SSB)"

ab.model <- model_select(all.data  = no.small, 
			     responseVar = "log_abundance", 
			     alpha = 0.05,
           fixedFactors= fF,
           fixedTerms= fT,
			     keepVars = keepVars,
           randomStruct = abundance.best.random$best.random,
			     otherRandoms=character(0),
           verbose=TRUE)
ab.model$stats
#"log_abundance~Within_PA+(1+Within_PA|SS)+(1|SSB)"

data <- ab.model$data

m1a <- lmer(log_abundance ~ Within_PA + (Within_PA|SS)+ (1|SSB), data =  data)

# plot

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

data <- no.small[,c("Within_PA", "SSS", "log_abundance", "log_slope", "log_elevation", "ag_suit")]
data <- na.omit(data)
text(2,80, paste("n =", length(data$SSS[which(data$Within_PA == "yes")]), sep = " "))
text(1,80, paste("n =", length(data$SSS[which(data$Within_PA == "no")]), sep = " "))

#get details for master plot
ab.plot1 <- data.frame(label = c("unprotected", "all protected"), est = points, 
		upper = c(100, CI[1]), lower = c(100,CI[2]),
		n.site = c(length(data$SSS[which(data$Within_PA == "no")]), 
			length(data$SSS[which(data$Within_PA == "yes")])))


# model abundance and IUCN_category
fF <- c("IUCN_CAT") 
fT <- list("ag_suit" = "1", "log_slope" = "1", "log_elevation" = "1")
keepVars <- list()
fI <- character(0)
RS <-  c("IUCN_CAT")


abundance.best.random.IUCN <- compare_randoms(no.small, "log_abundance",
				fixedFactors=fF,
        fixedTerms=fT,
        keepVars = keepVars,
        fixedInteractions=fI,
        otherRandoms=character(0),
				fixed_RandomSlopes = RS,
        fitInteractions=FALSE,
				verbose=TRUE)

abundance.best.random.IUCN$best.random #"(1+IUCN_CAT|SS)+ (1|SSB)"

ab.model.IUCN <- model_select(all.data  = no.small, 
			  responseVar = "log_abundance", 
        alpha = 0.05,
        fixedFactors= fF,
        fixedTerms= fT,
        keepVars = keepVars,
        randomStruct = abundance.best.random.IUCN$best.random,
        otherRandoms=character(0),
        verbose=TRUE)
# convergence issues with nonlinear terms
# convergence issues with starting with linear terms

#test without confounding variables
# both these converge
data <- ab.model.IUCN$data
m3ai <- lmer(log_abundance ~ 1 
	+ (IUCN_CAT|SS)+ (1|SSB), 
	 data = data,
	control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
m4ai <- lmer(log_abundance ~ IUCN_CAT 
	+ (IUCN_CAT|SS)+ (1|SSB), 
	 data = data,
	control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

anova(m3ai, m4ai)
#4.9139      3     0.1782

ab.model.IUCN$stats # similar, so can use stats table
ab.model.IUCN$final.call

# plot 

labels <- c("Unprotected", "IUCN III  - VI", "unknown", "IUCN I & II")

data <- ab.model.IUCN$data
data$IUCN_CAT <- relevel(data$IUCN_CAT, "0")

m4ai <- lmer(log_abundance ~ IUCN_CAT + (IUCN_CAT|SS) + (1|SSB), 
	 data = data)
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

data <- no.small[,c("IUCN_CAT", "SSS", "log_abundance", "log_slope", "log_elevation", "ag_suit")]
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






# abundance and zone

tropical <- subset(no.small, Zone == "Tropical")
temperate <- subset(no.small, Zone == "Temperate")

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

ab.best.random.trop$best.random 

ab.best.random.temp <- compare_randoms(temperate, "log_abundance",
				fixedFactors=fF,
        fixedTerms=fT,
        keepVars = keepVars,
        fixedInteractions=fI,
        otherRandoms=character(0),
				fixed_RandomSlopes = RS,
        fitInteractions=FALSE,
				verbose=TRUE)

ab.best.random.temp$best.random 

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
# "log_abundance~poly(ag_suit,1)+poly(log_elevation,2)+poly(log_slope,2)+Within_PA+(1+Within_PA|SS)+(1|SSB)"


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
#"log_abundance~poly(log_slope,2)+(1+Within_PA|SS)+(1|SSB)"


# run models for plot
data.trop <- ab.model.trop$data
data.temp <- ab.model.temp$data

m1aztr <- lmer(log_abundance ~ Within_PA + poly(log_slope,2) + poly(log_elevation,2)+ ag_suit
	+ (Within_PA|SS)+ (1|SSB), 
	 data = data.trop)

m1azte <- lmer(log_abundance ~ Within_PA +poly(log_slope,2) 
	+ (Within_PA|SS)+ (1|SSB), 
	 data = temperate)


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

plants <- subset(no.small, taxon_of_interest == "Plants")
inverts <- subset(no.small, taxon_of_interest == "Invertebrates")
verts <- subset(no.small, taxon_of_interest == "Vertebrates")
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
ab.best.random.p$best.random # 

ab.best.random.i <- compare_randoms(inverts, "log_abundance",
				fixedFactors=fF,
        fixedTerms=fT,
        keepVars = keepVars,
        fixedInteractions=fI,
        otherRandoms=character(0),
				fixed_RandomSlopes = RS,
        fitInteractions=FALSE,
				verbose=TRUE)
ab.best.random.i$best.random # 

ab.best.random.v <- compare_randoms(verts, "log_abundance",
				fixedFactors=fF,
        fixedTerms=fT,
        keepVars = keepVars,
        fixedInteractions=fI,
        otherRandoms=character(0),
				fixed_RandomSlopes = RS,
        fitInteractions=FALSE,
				verbose=TRUE)
ab.best.random.v$best.random #


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
# "log_abundance~poly(log_elevation,1)+poly(log_slope,3)+(1+Within_PA|SS)+(1|SSB)"

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
#"log_abundance~poly(ag_suit,3)+poly(log_slope,2)+Within_PA+(1+Within_PA|SS)+(1|SSB)"


# run models for plot
data.p <- ab.model.p$data
data.i <- ab.model.i$data
data.v <- ab.model.v$data

m1txp <- lmer(log_abundance ~ Within_PA +poly(log_elevation,1)+poly(log_slope,3)
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


# master plot


tiff( "simple models abundance NoSmall.tif",
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



