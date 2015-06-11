
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


# test effect of protection status

# check polynomials for confounding variables
fF <- c("Within_PA") 
fT <- list("ag_suit" = "3", "log_slope" = "3", "log_elevation" = "3")
keepVars <- list()
fI <- character(0)
RS <-  c("Within_PA")

r.sp.best.random <- compare_randoms(multiple.taxa.PA_11_14, "Richness_rarefied",
				fitFamily = "poisson",
				siteRandom = TRUE,
				fixedFactors=fF,
                        fixedTerms=fT,
			     	keepVars = keepVars,
                       	fixedInteractions=fI,
                        otherRandoms=character(0),
				fixed_RandomSlopes = RS,
                        fitInteractions=FALSE,
				verbose=TRUE)

r.sp.best.random$best.random #  

r.sp.model <- model_select(all.data  = multiple.taxa.PA_11_14, 
			     responseVar = "Richness_rarefied", 
			     fitFamily = "poisson", 
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
			     keepVars = keepVars,
                       randomStruct = r.sp.best.random$best.random,
			     otherRandoms=character(0),
                       verbose=TRUE)
r.sp.model$stats
r.sp.model$final.call  
#Richness_rarefied~poly(ag_suit,3)+poly(log_elevation,3)+(Within_PA|SS)+(1|SSB)+(1|SSBS)

#model for plot
m1 <- glmer(Richness_rarefied ~ Within_PA + poly(ag_suit,3)+poly(log_elevation,3)
	+ (Within_PA|SS)+ (1|SSB) + (1|SSBS), 
	family = "poisson", data = r.sp.model$data)

tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/01_15/simple model rar rich.tif",
	width = 10, height = 10, units = "cm", pointsize = 12, res = 300)

labels <- c("Unprotected", "Protected")
y <- as.numeric(fixef(m1)[2])
se <- as.numeric(se.fixef(m1)[2])
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
	ylab = "Rarefied richness difference (% ± 95%CI)",
	xlab = "")
#data <- multiple.taxa.PA_11_14[,c("Within_PA", "SSS", "Richness_rarefied")]
#data <- na.omit(data)
text(2,80, paste("n =", length(data$SS[which(data$Within_PA == "yes")]), sep = " "))
text(1,80, paste("n =", length(data$SS[which(data$Within_PA == "no")]), sep = " "))
axis(1, c(1,2), labels)
axis(2, c(80,100,120,140), c(80,100,120,140))
arrows(2,CI[1],2,CI[2], code = 3, length = 0.03, angle = 90)
abline(h = 100, lty = 2)
points(points ~ c(1,2), pch = 16, col = c(1,3), cex = 1.5)

dev.off()

#keep points for master plot
r.sp.plot1 <- data.frame(label = c("unprotected", "all protected"), est = points, 
		upper = c(100, CI[1]), lower = c(100,CI[2]),
		n.site = c(length(data$SS[which(data$Within_PA == "no")]), length(data$SS[which(data$Within_PA == "yes")])))


# rarefied richness with IUCN cat

# convergence issues at several stages in model select with all polynomials
fF <- c("IUCN_CAT") 
fT <- list("ag_suit" = "1", "log_slope" = "1", "log_elevation" = "1")
keepVars <- list()
fI <- character(0)
RS <-  c("IUCN_CAT")

r.sp.best.random.IUCN <- compare_randoms(multiple.taxa.PA_11_14, "Richness_rarefied",
				fitFamily = "poisson",
				siteRandom = TRUE,
				fixedFactors=fF,
                        fixedTerms=fT,
			     	keepVars = keepVars,
                       	fixedInteractions=fI,
                        otherRandoms=character(0),
				fixed_RandomSlopes = RS,
                        fitInteractions=FALSE,
				verbose=TRUE)

r.sp.best.random.IUCN$best.random #


r.sp.model.IUCN <- model_select(all.data  = multiple.taxa.PA_11_14, 
			     responseVar = "Richness_rarefied", 
			     fitFamily = "poisson", 
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
			     keepVars = keepVars,
                       randomStruct = r.sp.best.random.IUCN$best.random,
			     otherRandoms=character(0),
                       verbose=TRUE)
r.sp.model.IUCN$stats
r.sp.model.IUCN$final.call
r.sp.model.IUCN$warnings
# still too many convergence warnings

# compare with and without manually 
# test as linear
# convergence warnings for null model
m5i <- glmer(Richness_rarefied ~ 1 + log_slope + ag_suit + log_elevation
	+ (IUCN_CAT|SS) + (1|SSB)+ (1|SSBS), 
	family = "poisson", data = multiple.taxa.PA_11_14,
	control= glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
m6i <- glmer(Richness_rarefied ~ IUCN_CAT + log_slope + ag_suit + log_elevation 
	+ (IUCN_CAT|SS) + (1|SSB)+ (1|SSBS), 
	control= glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),
	family = "poisson", data = multiple.taxa.PA_11_14)
anova(m5i, m6i)
#2.1784      3     0.5362

# test without confounding variables
# lack of confounding variables only a concern if result was significant, as could be due to those other effects
# so ok to drop here as a comparison
m3i <- glmer(Richness_rarefied ~ 1 + 
	+ (IUCN_CAT|SS) + (1|SSB)+ (1|SSBS), 
	family = "poisson", data = multiple.taxa.PA_11_14,
	control= glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
m4i <- glmer(Richness_rarefied ~ IUCN_CAT 
	+ (IUCN_CAT|SS) + (1|SSB)+ (1|SSBS), 
	control= glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),
	family = "poisson", data = multiple.taxa.PA_11_14)
anova(m3i, m4i) 
# no convergence warnings - gives similar result

#compare model coefficients
fixef(m4i)
fixef(m6i)
#similar - therefore trust stats from model select and use those


# plot 

tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/01_15/simple model rar rich IUCN.tif",
	width = 15, height = 12, units = "cm", pointsize = 12, res = 300)


labels <- c("Unprotected", "IUCN III  - VI", "unknown", "IUCN I & II")

levels.IUCN <- levels(multiple.taxa.PA_11_14$IUCN_CAT)

multiple.taxa.PA_11_14$IUCN_CAT <- relevel(multiple.taxa.PA_11_14$IUCN_CAT, "0")

m6i <- glmer(Richness_rarefied~IUCN_CAT+poly(ag_suit,1)+poly(log_elevation,1)+(IUCN_CAT|SS)+(1|SSB)+(1|SSBS), 
	family = "poisson", data = r.sp.model.IUCN$data,
	control= glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

pos <- c(grep("4.5", names(fixef(m6i))),grep("7", names(fixef(m6i))),grep("1.5", names(fixef(m6i))))
y <- as.numeric(fixef(m6i)[pos])
se <- as.numeric(se.fixef(m6i)[pos])
yplus <- y + se*1.96
yminus <- y - se*1.96
y <-(exp(y)*100)
yplus<-(exp(yplus)*100)
yminus<-(exp(yminus)*100)

points <- c(100, y)
CI <- cbind(yplus, yminus)

plot(points ~ c(1,2,3,4), ylim = c(80,150), xlim = c(0.5,4.5),
	bty = "l", pch = 16, col = c(1,3,3,3), cex = 1.5,
	yaxt = "n", xaxt = "n",
	ylab = "Rarefied species richness difference (% ± 95%CI)",
	xlab = "")
axis(1,seq(1,length(points),1), labels)
axis(2, c(80,100,120,140), c(80,100,120,140))

data <- multiple.taxa.PA_11_14[,c("ag_suit", "log_elevation", "log_slope", "Within_PA", "SS", "SSB", "SSBS", "SSS",
		 "Richness_rarefied", "Source_ID", "IUCN_CAT", "Zone", "taxon_of_interest")]
data <- na.omit(data)
text(1, 80, paste("n =", length(data$SSS[which(data$IUCN_CAT == "0")]), sep = " "))
text(2, 80, paste("n =", length(data$SSS[which(data$IUCN_CAT == "4.5")]), sep = " "))
text(3, 80, paste("n =", length(data$SSS[which(data$IUCN_CAT == "7")]), sep = " "))
text(4, 80, paste("n =", length(data$SSS[which(data$IUCN_CAT == "1.5")]), sep = " "))

arrows(seq(2,length(points),1),CI[,1],
	seq(2,length(points),1),CI[,2], code = 3, length = 0.03, angle = 90)
abline(h = 100, lty = 2)
points(points ~ c(1,2,3,4), pch = 16, col = c(1,3,3,3), cex = 1.5)

dev.off()



IUCN.plot <- data.frame(label = labels[2:4], est = points[2:4], 
		upper = CI[,1], lower = CI[,2],
		n.site = c(length(data$SSS[which(data$IUCN_CAT == "4.5")]), 
			length(data$SSS[which(data$IUCN_CAT == "7")]),
			length(data$SSS[which(data$IUCN_CAT == "1.5")])))
r.sp.plot2 <- rbind(r.sp.plot1, IUCN.plot)




# simple species richness for Zone data


sp.tropical <- subset(multiple.taxa.PA_11_14, Zone == "Tropical")
sp.temperate <- subset(multiple.taxa.PA_11_14, Zone == "Temperate")

# check polynomials for confounding variables
fF <- c("Within_PA" ) 
fT <- list("ag_suit" = "3", "log_slope" = "3", "log_elevation" = "3")
keepVars <- list()
fI <- character(0)
RS <-  c("Within_PA")

r.sp.best.random.trop <- compare_randoms(sp.tropical, "Richness_rarefied",
				fitFamily = "poisson",
				siteRandom = TRUE,
				fixedFactors=fF,
                        fixedTerms=fT,
			     	keepVars = keepVars,
                       	fixedInteractions=fI,
                        otherRandoms=character(0),
				fixed_RandomSlopes = RS,
                        fitInteractions=FALSE,
				verbose=TRUE)

r.sp.best.random.trop$best.random #

r.sp.best.random.temp <- compare_randoms(sp.temperate, "Richness_rarefied",
				fitFamily = "poisson",
				siteRandom = TRUE,
				fixedFactors=fF,
                        fixedTerms=fT,
			     	keepVars = keepVars,
                       	fixedInteractions=fI,
                        otherRandoms=character(0),
				fixed_RandomSlopes = RS,
                        fitInteractions=FALSE,
				verbose=TRUE)

r.sp.best.random.temp$best.random #

# get polynomial relationships
sp.model.trop <- model_select(all.data  = sp.tropical, 
			     responseVar = "Richness_rarefied", 
			     fitFamily = "poisson", 
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
			     keepVars = keepVars,
                       randomStruct = r.sp.best.random.trop$best.random,
			     otherRandoms=character(0),
                       verbose=TRUE)
sp.model.trop$stats
sp.model.trop$final.call
#Richness_rarefied~poly(ag_suit,1)+poly(log_slope,3)+(1+Within_PA|SS)+(1|SSBS)"

sp.model.temp <- model_select(all.data  = sp.temperate, 
			     responseVar = "Richness_rarefied", 
			     fitFamily = "poisson", 
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
			     keepVars = keepVars,
                       randomStruct = r.sp.best.random.temp$best.random,
			     otherRandoms=character(0),
                       verbose=TRUE)
sp.model.temp$stats
sp.model.temp$final.call
#"Richness_rarefied~poly(ag_suit,3)+(1|SS)+(1|SSB)"


# run models for plot
# final call cannot be used as Within_PA not necessarily in it

m1ztr <- glmer(Richness_rarefied ~ Within_PA + poly(log_slope,3) + poly(ag_suit,1)
	+ (1+Within_PA|SS)+ (1|SSBS), family = "poisson",
	 data = sp.model.trop$trop)

m1zte <- glmer(Richness_rarefied ~ Within_PA  + poly(ag_suit,3)
	+ (1|SS)+ (1|SSB), family = "poisson", 
	 data =sp.model.temp$temp)



#add results to master plot
ztr.est <- exp(fixef(m1ztr)[2])*100
ztr.upper <- exp(fixef(m1ztr)[2] + 1.96* se.fixef(m1ztr)[2])*100
ztr.lower <- exp(fixef(m1ztr)[2] - 1.96* se.fixef(m1ztr)[2])*100

zte.est <- exp(fixef(m1zte)[2])*100
zte.upper <- exp(fixef(m1zte)[2] + 1.96* se.fixef(m1zte)[2])*100
zte.lower <- exp(fixef(m1zte)[2] - 1.96* se.fixef(m1zte)[2])*100

zone <- data.frame(label = c("Tropical", "Temperate"),
				est = c(ztr.est, zte.est), 
				upper = c(ztr.upper, zte.upper), 
				lower = c(ztr.lower, zte.lower), 
				n.site = c(nrow(data.trop[which(data.trop$Within_PA == "yes"),]), 
				nrow(data.temp[which(data.temp$Within_PA == "yes"),])))
r.sp.plot3 <- rbind(r.sp.plot2, zone)







# species richness and taxon

plants <- subset(multiple.taxa.PA_11_14, taxon_of_interest == "Plants")
inverts <- subset(multiple.taxa.PA_11_14, taxon_of_interest == "Invertebrates")
verts <- subset(multiple.taxa.PA_11_14, taxon_of_interest == "Vertebrates")
nrow(plants)
nrow(inverts)
nrow(verts)

# check polynomials for confounding variables
fF <- c("Within_PA" ) 
fT <- list("ag_suit" = "3", "log_slope" = "3", "log_elevation" = "3")
keepVars <- list()
fI <- character(0)
RS <-  c("Within_PA")


best.random.p <- compare_randoms(plants, "Richness_rarefied",
				fitFamily = "poisson",
				siteRandom = TRUE,
				fixedFactors=fF,
                        fixedTerms=fT,
			     	keepVars = keepVars,
                       	fixedInteractions=fI,
                        otherRandoms=character(0),
				fixed_RandomSlopes = RS,
                        fitInteractions=FALSE,
				verbose=TRUE)
best.random.p$best.random # 

best.random.i <- compare_randoms(inverts, "Richness_rarefied",
				fitFamily = "poisson",
				siteRandom = TRUE,
				fixedFactors=fF,
                        fixedTerms=fT,
			     	keepVars = keepVars,
                       	fixedInteractions=fI,
                        otherRandoms=character(0),
				fixed_RandomSlopes = RS,
                        fitInteractions=FALSE,
				verbose=TRUE)
best.random.i$best.random # 

best.random.v <- compare_randoms(verts, "Richness_rarefied",
				fitFamily = "poisson",
				siteRandom = TRUE,
				fixedFactors=fF,
                        fixedTerms=fT,
			     	keepVars = keepVars,
                       	fixedInteractions=fI,
                        otherRandoms=character(0),
				fixed_RandomSlopes = RS,
                        fitInteractions=FALSE,
				verbose=TRUE)
best.random.v$best.random #





# get polynomial relationships
model.p <- model_select(all.data  = plants, 
				fitFamily = "poisson",
			     responseVar = "Richness_rarefied", 
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
			     keepVars = keepVars,
                       randomStruct =best.random.p$best.random,
			     otherRandoms=character(0),
                       verbose=TRUE)
model.p$stats
model.p$final.call
#"Richness_rarefied~poly(ag_suit,1)+poly(log_elevation,3)+(1+Within_PA|SS)+(1|SSBS)"

model.i <- model_select(all.data  = inverts, 
			     responseVar = "Richness_rarefied", 
				fitFamily = "poisson",
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
			     keepVars = keepVars,
                       randomStruct =best.random.i$best.random,
			     otherRandoms=character(0),
                       verbose=TRUE)
model.i$stats
model.i$final.call
#"Richness_rarefied~1+(1+Within_PA|SS)+(1|SSBS)+(1|SSB)"


model.v <- model_select(all.data  = verts, 
			     responseVar = "Richness_rarefied", 
				fitFamily = "poisson",
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
			     keepVars = keepVars,
                       randomStruct =best.random.v$best.random,
			     otherRandoms=character(0),
                       verbose=TRUE)
model.v$stats
model.v$final.call
#"Richness_rarefied~1+(1+Within_PA|SS)+(1|SSBS)"


# run models for plots

m1txp <- glmer(Richness_rarefied ~ Within_PA + poly(ag_suit,1)+poly(log_elevation,3)
	+ (Within_PA|SS) + (1|SSBS), family = "poisson", 
	 data = model.p$data)

m1txi <- glmer(Richness_rarefied ~ Within_PA 
	+ (Within_PA|SS)+ (1|SSB)+ (1|SSBS), family = "poisson", 
	 data = model.i$data)

m1txv <- glmer(Richness_rarefied ~ Within_PA 
	+(1+Within_PA|SS)+(1|SSBS), family = "poisson", 
	 data = model.v$data)


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



tax <- data.frame(label = c("Plants", "Inverts", "Verts"),
				est = c(txp.est, txi.est, txv.est), 
				upper = c(txp.upper, txi.upper, txv.upper), 
				lower = c(txp.lower, txi.lower, txv.lower), 
				n.site = c(nrow(data.p[which(data.p$Within_PA == "yes"),]), 
					nrow(data.i[which(data.i$Within_PA == "yes"),]), 
					nrow(data.v[which(data.v$Within_PA == "yes"),])))
r.sp.plot <- rbind(r.sp.plot3, tax)

r.sp.plot$truncated <- r.sp.plot$upper
r.sp.plot$truncated[which(r.sp.plot$upper > 145)] <- 135
r.sp.plot$dotlines <- rep(0,nrow(r.sp.plot))
r.sp.plot$dotlines_end <- rep(0,nrow(r.sp.plot))
r.sp.plot$dotlines[which(r.sp.plot$upper > 145)] <- 135
r.sp.plot$dotlines_end[which(r.sp.plot$upper > 145)] <- 140


# master plot

load("\\\\smbhome.uscs.susx.ac.uk\\clg32\\Documents\\PREDICTS\\WDPA analysis\\RData files\\8 landuses\\simple models - rar rich.RData")

tiff( "R:/ecocon_d/clg32/PREDICTS/WDPA analysis/plots/02_15/simple models rar rich.tif",
	width = 10, height = 15, units = "cm", pointsize = 12, res = 300)

trop.col <- rgb(0.9,0,0)
temp.col <- rgb(0,0.1,0.7)
p.col <- rgb(0.2,0.7,0.2)
i.col <- rgb(0,0.7,0.9)
v.col <- rgb(0.9,0.5,0)

par(mar = c(9,6,4,1))
plot(1,1, 
	ylim = c(70,145), xlim = c(0.5,nrow(r.sp.plot)),
	bty = "l", 
	axes = F,
	ylab = "Rarefied richness difference (%)",
	cex.lab = 1.5,
	xlab = "")
arrows(1:nrow(r.sp.plot),r.sp.plot$truncated,
	col = c(1,1,rep(rgb(0.5, 0.5, 0.5), 3), c(trop.col, temp.col, p.col, i.col, v.col)),
	lwd = 2,
	1:nrow(r.sp.plot),r.sp.plot$lower, code = 3, length = 0, angle = 90)
arrows(1:nrow(r.sp.plot),r.sp.plot$dotlines, 1:nrow(r.sp.plot),r.sp.plot$dotlines_end,
	col = c(1,1,rep(rgb(0.5, 0.5, 0.5), 3), c(trop.col, temp.col, p.col, i.col, v.col)),
	lwd = 2, lty = 3,
	code = 3, length = 0, angle = 90)
points(r.sp.plot$est ~ c(1:nrow(r.sp.plot)),
	pch = c(21, rep(16,4), rep(15,2),rep(17,3)), 
	lwd = 2,
	col = c(1,1,rep(rgb(0.5, 0.5, 0.5), 3), c(trop.col, temp.col, p.col, i.col, v.col)),
	bg = "white", 
	cex = 1.5)
abline(v = c(2.5,5.5,7.5), col = 8)
abline(h= 100, lty = 2)
text(1:nrow(r.sp.plot),72, r.sp.plot$n.site, srt = 180)
text(which(r.sp.plot$upper > 145),143, round(r.sp.plot$upper[which(r.sp.plot$upper > 145)]), srt = 90)
#axis(1, c(1:nrow(r.sp.plot)), r.sp.plot$label, cex.axis = 1.5, las = 2)
axis(2, c(80,100,120,140), c(80,100,120,140))


dev.off()



save.image("N:\\Documents\\PREDICTS\\WDPA analysis\\RData files\\8 landuses\\simple models - rar rich.RData")


