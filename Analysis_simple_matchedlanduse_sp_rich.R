rm(list=ls())

library(lme4)
library(yarg)
library(roquefort)

# load functions
setwd("R:/ecocon_d/clg32/GitHub/PREDICTS_WDPA")
source("compare_randoms.R")
source("model_select.R")


# Load dataset 

#load data
source("prep_matched.landuse_for_analysis.R")



#make ordinal datasets

matched.landuse_ord <- matched.landuse
matched.landuse_ord$IUCN_CAT <- factor(matched.landuse_ord$IUCN_CAT, ordered = T, 
	levels = c("0", "4.5", "7", "1.5"))

multiple.taxa.matched.landuse_ord <- multiple.taxa.matched.landuse
multiple.taxa.matched.landuse_ord$IUCN_CAT <- factor(multiple.taxa.matched.landuse_ord$IUCN_CAT, ordered = T, 
	levels = c("0", "4.5", "7", "1.5"))



### model species richness

# check polynomials
fF <- c("Within_PA") 
fT <- list("ag_suit" = "3", "log_slope" = "3", "log_elevation" = "3")
keepVars <- list()
fI <- character(0)
RS <-  c("Within_PA")
#Species_richness~poly(ag_suit,3)+poly(log_elevation,2)+(Within_PA|SS)+(1|SSB)+(1|SSBS)

Species_richness.best.random <- compare_randoms(multiple.taxa.matched.landuse, "Species_richness",
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

Species_richness.best.random$best.random #"(1+Within_PA|SS)+ (1|SSBS)+ (1|SSB)"


s.model <- model_select(all.data  = multiple.taxa.matched.landuse, 
			     responseVar = "Species_richness", 
			     fitFamily = "poisson", 
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
			     keepVars = keepVars,
                       randomStruct = Species_richness.best.random$best.random,
			     otherRandoms=character(0),
                       verbose=TRUE)


data <- multiple.taxa.matched.landuse[,c("ag_suit", "log_elevation", "log_slope", "Within_PA", "SS", "SSB", "SSBS", "Species_richness")]
data <- na.omit(data)
m1 <- glmer(Species_richness ~ Within_PA + log_slope +poly(ag_suit,3)+poly(log_elevation,2)
	+ (Within_PA|SS) + (1|SSB) + (1|SSBS), 
	family = "poisson", data = data)
m2 <- glmer(Species_richness ~ 1 + log_slope + poly(ag_suit,3)+poly(log_elevation,2)
	+ (Within_PA|SS) + (1|SSB) + (1|SSBS), 
	family = "poisson", data = data)
anova(m1, m2)

summary(m1)
# significant difference - now to plot the difference

# PREDICTS plots for landuse - the y axis is % difference
# this is obtained by
# getting estimates and lower and upper CI values
y <- fixef(m1)[2]
se <- se.fixef(m1)[2]  
yplus <- y + se*1.96
yminus <- y - se*1.96

# then backtransform and turn into percentage in one step
y <-(exp(y)*100)-100 
yplus<-(exp(yplus)*100)-100
yminus<-(exp(yminus)*100)-100



### need to understand this step ###
# this approach indicates that the exp backtransformed coefficient for inside PA 
# is actually a proportion of the outside PA value
#(when the to the intercept is outside PA) 
# rather than a number of species units that should be added on to the backtransformed value for the intercept
# I thought that where outside PA is the reference value,
# the coefficient for inside PA is added on to get the absolute value for inside PA
# i.e. the absolute value of species richness outside a PA would be
intercept <-fixef(m1)[1]
exp(intercept)
# in a PA it would be 
y <- fixef(m1)[2]
exp(intercept+y)

# Then, the percentage difference would be
exp(intercept+y)/exp(intercept)*100 - 100

#upper and lower would be
exp(fixef(m1)[1] + fixef(m1)[2]+1.96*se.fixef(m1)[2])/exp(fixef(m1)[1])*100
exp(fixef(m1)[1] + fixef(m1)[2]-1.96*se.fixef(m1)[2])/exp(fixef(m1)[1])*100

# ! its the same as the approach above. Excellent.

# but then, when the binomial family is used, the intercept + y as a percentage of intercept has to be used:
# this is from plot_lu_effects
# intercept<-fixef(model)['(Intercept)']
#   y<-(((1/(1+exp(-(intercept+y))))/(1/(1+exp(-(intercept)))))*100)-100
#   yplus<-(((1/(1+exp(-(intercept+yplus))))/(1/(1+exp(-(intercept)))))*100)-100
#   yminus<-(((1/(1+exp(-(intercept+yminus))))/(1/(1+exp(-(intercept)))))*100)-100



# what I had before dec 14 - following Andys instructions in Cambridge oct 2014
#x <- exp(fixef(m1)[2]) - 1

# plot

tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/01_15/simple model matched.landuse sp rich.tif",
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
	ylab = "Species richness difference (% ± 95%CI)",
	xlab = "")

data <- multiple.taxa.matched.landuse[,c("Within_PA", "SSS", "Species_richness")]
data <- na.omit(data)
text(2,80, paste("n =", length(data$SSS[which(data$Within_PA == "yes")]), sep = " "))
text(1,80, paste("n =", length(data$SSS[which(data$Within_PA == "no")]), sep = " "))


axis(1, c(1,2), labels)
axis(2, c(80,100,120,140), c(80,100,120,140))
arrows(2,CI[1],2,CI[2], code = 3, length = 0.03, angle = 90)
abline(h = 100, lty = 2)
points(points ~ c(1,2), pch = 16, col = c(1,3), cex = 1.5)

dev.off()


#keep points for master plot
sp.plot1 <- data.frame(label = c("unprotected", "all protected"), est = points, 
		upper = c(100, CI[1]), lower = c(100,CI[2]),
		n.site = c(length(data$SSS[which(data$Within_PA == "no")]), length(data$SSS[which(data$Within_PA == "yes")])))





#simple species richness with IUCN cat 


fF <- c("IUCN_CAT") 
fT <- list("ag_suit" = "3", "log_slope" = "3", "log_elevation" = "3")
keepVars <- list()
fI <- character(0)
RS <-  c("IUCN_CAT")
#cant converge with non-linear confounding variables


Species_richness.best.random.IUCN <- compare_randoms(multiple.taxa.matched.landuse, "Species_richness",
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

Species_richness.best.random.IUCN$best.random # "(1+IUCN_CAT|SS)+ (1|SSBS)+ (1|SSB)"

s.model.IUCN <- model_select(all.data  = multiple.taxa.matched.landuse, 
			     responseVar = "Species_richness", 
			     fitFamily = "poisson", 
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
			     keepVars = keepVars,
                       randomStruct = Species_richness.best.random.IUCN$best.random,
			     otherRandoms=character(0),
                       verbose=TRUE)


data <- multiple.taxa.matched.landuse[,c("ag_suit", "log_elevation", "log_slope", "IUCN_CAT", "SS", "SSB", "SSBS", "Species_richness")]
data <- na.omit(data)


m2i <- glmer(Species_richness ~ 1  + log_slope + log_elevation + ag_suit
	+ (IUCN_CAT|SS) + (1|SSB) + (1|SSBS), 
	family = "poisson", data = multiple.taxa.matched.landuse,
	control= glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
m3i <- glmer(Species_richness ~ IUCN_CAT + log_slope + log_elevation + ag_suit
	+ (IUCN_CAT|SS) + (1|SSB) + (1|SSBS), 
	family = "poisson", data = multiple.taxa.matched.landuse,
	control= glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
#tried this m3i with Nelder_Mead - also convergence warnings

#doesnt converge with ordinal data either
m4i <- glmer(Species_richness ~ 1  + log_slope + log_elevation + ag_suit
	+ (IUCN_CAT|SS) + (1|SSB) + (1|SSBS), 
	family = "poisson", data = multiple.taxa.matched.landuse_ord,
	control= glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
m5i <- glmer(Species_richness ~ IUCN_CAT + log_slope + log_elevation + ag_suit
	+ (IUCN_CAT|SS) + (1|SSB) + (1|SSBS), 
	family = "poisson", data = multiple.taxa.matched.landuse_ord,
	control= glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))


anova(m1i, m0i)
anova(m2i, m3i)
summary(m1i)

m5i <- glmer(Species_richness ~ 1  + log_slope + log_elevation + ag_suit
	+ (Within_PA|SS) + (1|SSB) + (1|SSBS), 
	family = "poisson", data = multiple.taxa.matched.landuse,
	control= glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
m6i <- glmer(Species_richness ~ IUCN_CAT + log_slope + log_elevation + ag_suit
	+ (Within_PA|SS) + (1|SSB) + (1|SSBS), 
	family = "poisson", data = multiple.taxa.matched.landuse,
	control= glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
anova(m5i, m6i)

m7i <- glmer(Species_richness ~ 1  + log_slope + log_elevation + ag_suit
	+ (1|SS) + (1|SSB) + (1|SSBS), 
	family = "poisson", data = multiple.taxa.matched.landuse,
	control= glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
m8i <- glmer(Species_richness ~ IUCN_CAT + log_slope + log_elevation + ag_suit
	+ (1|SS) + (1|SSB) + (1|SSBS), 
	family = "poisson", data = multiple.taxa.matched.landuse,
	control= glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
anova(m7i, m8i)


# plot 

tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/01_15/simple model matched.landuse sp rich IUCN.tif",
	width = 15, height = 12, units = "cm", pointsize = 12, res = 300)


labels <- c("Unprotected", "IUCN III  - VI", "unknown", "IUCN I & II")

levels(multiple.taxa.matched.landuse$IUCN_CAT)
multiple.taxa.matched.landuse$IUCN_CAT <- relevel(multiple.taxa.matched.landuse$IUCN_CAT, "0")

m4i <- glmer(Species_richness ~ IUCN_CAT + log_slope + log_elevation + ag_suit
	+ (Within_PA|SS) + (1|SSB) + (1|SSBS), 
	family = "poisson", data = multiple.taxa.matched.landuse,
	control= glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

pos <- c(grep("4.5", names(fixef(m4i))),grep("7", names(fixef(m4i))),grep("1.5", names(fixef(m4i))))
y <- as.numeric(fixef(m4i)[pos])
se <- as.numeric(se.fixef(m4i)[pos])
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
	ylab = "Species richness difference (% ± 95%CI)",
	xlab = "")
axis(1,seq(1,length(points),1), labels)
axis(2, c(80,100,120,140), c(80,100,120,140))


data <- multiple.taxa.matched.landuse[,c("IUCN_CAT", "SSS", "Species_richness")]
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
sp.plot2 <- rbind(sp.plot1, IUCN.plot)



# simple species richness for Zone data


sp.tropical <- subset(multiple.taxa.matched.landuse, Zone == "Tropical")
sp.temperate <- subset(multiple.taxa.matched.landuse, Zone == "Temperate")

# check polynomials for confounding variables
fF <- c("Within_PA" ) 
fT <- list("ag_suit" = "3", "log_slope" = "3", "log_elevation" = "3")
keepVars <- list()
fI <- character(0)
RS <-  c("Within_PA")

Sp.best.random.trop <- compare_randoms(sp.tropical, "Species_richness",
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

Sp.best.random.trop$best.random #"(1+Within_PA|SS)+ (1|SSBS)+ (1|SSB)"

Sp.best.random.temp <- compare_randoms(sp.temperate, "Species_richness",
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

Sp.best.random.temp$best.random #"(1+Within_PA|SS)+ (1|SSBS)+ (1|SSB)"

# get polynomial relationships
sp.model.trop <- model_select(all.data  = sp.tropical, 
			     responseVar = "Species_richness", 
			     fitFamily = "poisson", 
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
			     keepVars = keepVars,
                       randomStruct = Sp.best.random.trop$best.random,
			     otherRandoms=character(0),
                       verbose=TRUE)
#"Species_richness~poly(ag_suit,1)+poly(log_elevation,3)+Within_PA+(1+Within_PA|SS)+(1|SSBS)+(1|SSB)"



sp.model.temp <- model_select(all.data  = sp.temperate, 
			     responseVar = "Species_richness", 
			     fitFamily = "poisson", 
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
			     keepVars = keepVars,
                       randomStruct = Sp.best.random.temp$best.random,
			     otherRandoms=character(0),
                       verbose=TRUE)
#"Species_richness~poly(log_elevation,1)+poly(log_slope,2)+(1+Within_PA|SS)+(1|SSBS)+(1|SSB)"




# run models
data.trop <- sp.tropical[,c("Within_PA", "ag_suit", "log_elevation", "log_slope", "SS", "SSB", "SSBS", "Species_richness")] 
data.trop <- na.omit(data.trop)
m1ztr <- glmer(Species_richness ~ Within_PA + log_slope + poly(ag_suit,1)+poly(log_elevation,3)
	+ (1+Within_PA|SS)+ (1|SSBS)+ (1|SSB), family = "poisson",
	 data = data.trop)
m2ztr <- glmer(Species_richness ~ 1 + log_slope +poly(ag_suit,1)+poly(log_elevation,3)
	+(1+Within_PA|SS)+ (1|SSBS)+ (1|SSB), family = "poisson", 
	 data = data.trop)
anova(m1ztr, m2ztr)
#0.1603      1     0.6889

data.temp <- sp.temperate[,c("Within_PA", "ag_suit", "log_elevation", "log_slope", "SS", "SSB","SSBS", "Species_richness")] 
data.temp <- na.omit(data.temp)
m1zte <- glmer(Species_richness ~ Within_PA +poly(log_slope,2) + poly(log_elevation,1) + ag_suit
	+ (1+Within_PA|SS)+ (1|SSBS)+ (1|SSB), family = "poisson", 
	 data = data.temp)
m2zte <- glmer(Species_richness ~ 1 + poly(log_slope,2) + poly(log_elevation,1) + ag_suit
	+ (1+Within_PA|SS)+ (1|SSBS)+ (1|SSB), family = "poisson", 
	 data = data.temp)
anova(m1zte, m2zte)
# 1.4382      1     0.2304


#add results to master plot
ztr.est <- exp(fixef(m1ztr)[2])*100
ztr.upper <- exp(fixef(m1ztr)[2] + 1.96* se.fixef(m1ztr)[2])*100
ztr.lower <- exp(fixef(m1ztr)[2] - 1.96* se.fixef(m1ztr)[2])*100

zte.est <- exp(fixef(m1zte)[2])*100
zte.upper <- exp(fixef(m1zte)[2] + 1.96* se.fixef(m1zte)[2])*100
zte.lower <- exp(fixef(m1zte)[2] - 1.96* se.fixef(m1zte)[2])*100

a.zone <- data.frame(label = c("Tropical", "Temperate"),
				est = c(ztr.est, zte.est), 
				upper = c(ztr.upper, zte.upper), 
				lower = c(ztr.lower, zte.lower), 
				n.site = c(nrow(data.trop), nrow(data.temp)))
sp.plot3 <- rbind(sp.plot2, a.zone)







# species richness and taxon

plants <- subset(multiple.taxa.matched.landuse, taxon_of_interest == "Plants")
inverts <- subset(multiple.taxa.matched.landuse, taxon_of_interest == "Invertebrates")
verts <- subset(multiple.taxa.matched.landuse, taxon_of_interest == "Vertebrates")
nrow(plants)
nrow(inverts)
nrow(verts)

# check polynomials for confounding variables
fF <- c("Within_PA" ) 
fT <- list("ag_suit" = "3", "log_slope" = "3", "log_elevation" = "3")
keepVars <- list()
fI <- character(0)
RS <-  c("Within_PA")


best.random.p <- compare_randoms(plants, "Species_richness",
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
best.random.p$best.random # "(1+Within_PA|SS)+ (1|SSBS)"

best.random.i <- compare_randoms(inverts, "Species_richness",
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
best.random.i$best.random # "(1+Within_PA|SS)+ (1|SSBS)+ (1|SSB)"

best.random.v <- compare_randoms(verts, "Species_richness",
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
best.random.v$best.random #"(1+Within_PA|SS)+ (1|SSBS)+ (1|SSB)"





# get polynomial relationships
model.p <- model_select(all.data  = plants, 
				fitFamily = "poisson",
			     responseVar = "Species_richness", 
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
			     keepVars = keepVars,
                       randomStruct =best.random.p$best.random,
			     otherRandoms=character(0),
                       verbose=TRUE)
# "Species_richness~poly(ag_suit,1)+poly(log_elevation,2)+(1+Within_PA|SS)+(1|SSBS)"

model.i <- model_select(all.data  = inverts, 
			     responseVar = "Species_richness", 
				fitFamily = "poisson",
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
			     keepVars = keepVars,
                       randomStruct =best.random.i$best.random,
			     otherRandoms=character(0),
                       verbose=TRUE)
#"Species_richness~poly(log_elevation,1)+(1+Within_PA|SS)+(1|SSBS)+(1|SSB)"


model.v <- model_select(all.data  = verts, 
			     responseVar = "Species_richness", 
				fitFamily = "poisson",
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
			     keepVars = keepVars,
                       randomStruct =best.random.v$best.random,
			     otherRandoms=character(0),
                       verbose=TRUE)
#"Species_richness~poly(ag_suit,3)+poly(log_elevation,3)+(1+Within_PA|SS)+(1|SSBS)+(1|SSB)"


# run models
data.p <- plants[,c("Within_PA", "ag_suit", "log_elevation", "log_slope", "SS", "SSB", "SSBS", "Species_richness")]
data.p <- na.omit(data.p)
m1txp <- glmer(Species_richness ~ Within_PA + poly(ag_suit,1)+poly(log_elevation,2)+poly(log_slope,1)
	+ (Within_PA|SS)+  (1|SSBS), family = "poisson", 
	 data = data.p)
m2txp <- glmer(Species_richness ~ 1+poly(ag_suit,1)+poly(log_elevation,2)+poly(log_slope,1)
	+ (Within_PA|SS)+ (1|SSBS), family = "poisson", 
	 data = data.p)
anova(m1txp , m2txp)
#0.3559      1     0.5508

data.i <- inverts[,c("Within_PA", "ag_suit", "log_elevation", "log_slope", "SS", "SSB", "SSBS","Species_richness")]
data.i <- na.omit(data.i)
m1txi <- glmer(Species_richness ~ Within_PA +poly(log_elevation,1)+poly(log_slope,1) +ag_suit
	+ (Within_PA|SS)+ (1|SSB)+ (1|SSBS), family = "poisson", 
	 data = data.i)
m2txi<- glmer(Species_richness ~ 1 +poly(log_elevation,1)+poly(log_slope,1) + ag_suit
	+ (Within_PA|SS)+ (1|SSB)+ (1|SSBS), family = "poisson", 
	 data = data.i)
anova(m1txi, m2txi)
#1.4926      1     0.2218

data.v <- verts[,c("Within_PA", "ag_suit", "log_elevation", "log_slope", "SS", "SSB","SSBS", "Species_richness")]
data.v <- na.omit(data.v)
m1txv <- glmer(Species_richness ~ Within_PA + poly(ag_suit,3)+poly(log_elevation,3)+ poly(log_slope,1)
	+ (Within_PA|SS)+ (1|SSB)+ (1|SSBS), family = "poisson", 
	 data = data.v)
m2txv <- glmer(Species_richness ~ 1 + poly(ag_suit,3)+poly(log_elevation,3)+ poly(log_slope,1)
	+ (Within_PA|SS)+ (1|SSB)+ (1|SSBS), family = "poisson", 
	 data = data.v)
anova(m1txv, m2txv)
#7e-04      1     0.9787

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
				n.site = c(nrow(data.p), nrow(data.i), nrow(data.v)))
sp.plot <- rbind(sp.plot3, tax)




# master plot



tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/02_15/simple models matchedlanduse sp rich.tif",
	width = 23, height = 16, units = "cm", pointsize = 12, res = 300)

trop.col <- rgb(0.9,0,0)
temp.col <- rgb(0,0.1,0.7)
p.col <- rgb(0.2,0.7,0.2)
i.col <- rgb(0,0.7,0.9)
v.col <- rgb(0.9,0.5,0)

par(mar = c(9,6,4,1))
plot(1,1, 
	ylim = c(65,190), xlim = c(0.5,nrow(sp.plot)+1),
	bty = "l", 
	axes = F,
	ylab = "Species richness difference (%)",
	cex.lab = 1.5,
	xlab = "")
arrows(1:nrow(sp.plot),sp.plot$upper,
	col = c(1,1,rep(rgb(0.5, 0.5, 0.5), 3), c(trop.col, temp.col, p.col, i.col, v.col)),
	lwd = 2,
	1:nrow(sp.plot),sp.plot$lower, code = 3, length = 0, angle = 90)
points(sp.plot$est ~ c(1:nrow(sp.plot)),
	pch = c(21, rep(16,4), rep(15,2),rep(17,3)), 
	lwd = 2,
	col = c(1,1,rep(rgb(0.5, 0.5, 0.5), 3), c(trop.col, temp.col, p.col, i.col, v.col)),
	bg = "white", 
	cex = 1.5)
abline(v = c(2.5,5.5,7.5), col = 8)
abline(h= 100, lty = 2)
text(1:nrow(sp.plot),65, sp.plot$n.site)
axis(1, c(1:nrow(sp.plot)), sp.plot$label, cex.axis = 1.5, las = 2, tick = 0)
axis(2, c(80,100,120,140,160,180), c(80,100,120,140,160,180))


dev.off()



