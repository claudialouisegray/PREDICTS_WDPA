
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

#load("\\\\smbhome.uscs.susx.ac.uk\\clg32\\Documents\\PREDICTS\\WDPA analysis\\RData files\\8 landuses\\simple models - sp rich.RData")



### model species richness

# check polynomials for confounding variables
fF <- c("Within_PA") 
fT <- list("ag_suit" = "3", "log_slope" = "3", "log_elevation" = "3")
keepVars <- list()
fI <- character(0)
RS <-  c("Within_PA")

Species_richness.best.random <- compare_randoms(multiple.taxa.PA_11_14, "Species_richness",
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


sp.model <- model_select(all.data  = multiple.taxa.PA_11_14, 
			     responseVar = "Species_richness", 
			     fitFamily = "poisson", 
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
			     keepVars = keepVars,
                       randomStruct = Species_richness.best.random$best.random,
			     otherRandoms=character(0),
                       verbose=TRUE)
sp.model$final.call
# Species_richness~poly(ag_suit,1)+poly(log_elevation,1)+Within_PA+(Within_PA|SS)+(1|SSB)+(1|SSBS)

data <- multiple.taxa.PA_11_14[,c("ag_suit", "log_elevation", "log_slope", "Within_PA", "SS", "SSB", "SSBS", "Species_richness")]
data <- na.omit(data)
m1 <- glmer(Species_richness ~ Within_PA + poly(log_elevation,1) + ag_suit
	+ (Within_PA|SS) + (1|SSB) + (1|SSBS), 
	family = "poisson", data = data)
m2 <- glmer(Species_richness ~ 1 + poly(log_elevation,1) + ag_suit
	+ (Within_PA|SS) + (1|SSB) + (1|SSBS), 
	family = "poisson", data = data)
anova(m1, m2)


# redo with non-orthogonal polynomials for predicting from rasters

m1.est <- glmer(Species_richness ~ Within_PA + log_slope + log_elevation + I(log_elevation^2) + ag_suit
	+ (Within_PA|SS) + (1|SSB) + (1|SSBS), 
	family = "poisson", data = data)
fixef(m1.est)
fixef(m1)

### plotting
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





# plot

tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/01_15/simple model sp rich.tif",
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


text(2,80, paste("n =", length(data$SS[which(data$Within_PA == "yes")]), sep = " "))
text(1,80, paste("n =", length(data$SS[which(data$Within_PA == "no")]), sep = " "))


axis(1, c(1,2), labels)
axis(2, c(80,100,120,140), c(80,100,120,140))
arrows(2,CI[1],2,CI[2], code = 3, length = 0.03, angle = 90)
abline(h = 100, lty = 2)
points(points ~ c(1,2), pch = 16, col = c(1,3), cex = 1.5)

dev.off()

#keep points for master plot
sp.plot1 <- data.frame(label = c("unprotected", "all protected"), est = points, 
		upper = c(100, CI[1]), lower = c(100,CI[2]),
		n.site = c(length(data$SS[which(data$Within_PA == "no")]), length(data$SS[which(data$Within_PA == "yes")])))



### simple species richness with IUCN cat 

multiple.taxa.PA_11_14$IUCN_CAT <- relevel(multiple.taxa.PA_11_14$IUCN_CAT, "0")

# check polynomials for confounding variables
fF <- c("IUCN_CAT") 
fT <- list("ag_suit" = "3", "log_slope" = "3", "log_elevation" = "3")
keepVars <- list()
fI <- character(0)
RS <-  c("IUCN_CAT")


Species_richness.best.random.IUCN <- compare_randoms(multiple.taxa.PA_11_14, "Species_richness",
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


sp.model.IUCN <- model_select(all.data  = multiple.taxa.PA_11_14, 
			     responseVar = "Species_richness", 
			     fitFamily = "poisson", 
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
			     keepVars = keepVars,
                       randomStruct = Species_richness.best.random.IUCN$best.random,
			     otherRandoms=character(0),
                       verbose=TRUE)
sp.model.IUCN$final.call
# "Species_richness~IUCN_CAT+poly(ag_suit,1)+poly(log_elevation,3)+(IUCN_CAT|SS)+(1|SSB)+(1|SSBS)"
# but, doesnt converge for elevation down to quadratic
# compare cubic and linear

data <- multiple.taxa.PA_11_14[,c("ag_suit", "log_elevation", "log_slope", "IUCN_CAT", "SS", "SSB", "SSBS", "Species_richness")]
data <- na.omit(data)

m2i <- glmer(Species_richness ~ 1 + poly(log_elevation,1) + ag_suit
	+ (IUCN_CAT|SS) + (1|SSB) + (1|SSBS), 
	family = "poisson", data = data)
#doesnt converge with cubic - use linear elevation relationship

m3ii <- glmer(Species_richness ~ IUCN_CAT + poly(log_elevation,3) + ag_suit
	+ (IUCN_CAT|SS) + (1|SSB) + (1|SSBS), 
	family = "poisson", data = data)
m3i <- glmer(Species_richness ~ IUCN_CAT  + poly(log_elevation,1) + ag_suit
	+ (IUCN_CAT|SS) + (1|SSB) + (1|SSBS), 
	family = "poisson", data = data)
fixef(m3i)
fixef(m3i)
#coefficients are similar for both models

anova(m2i, m3i) # 10.785      3    0.01294
summary(m3i)

# look at differences between strictest categories and others
data$IUCN_CAT <- relevel(data$IUCN_CAT, "1.5")
m3i <- glmer(Species_richness ~ IUCN_CAT  + poly(log_elevation,1) + ag_suit
	+ (IUCN_CAT|SS) + (1|SSB) + (1|SSBS), 
	family = "poisson", data = data)
summary(m3i)



# plot 

tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/01_15/simple model sp rich IUCN.tif",
	width = 15, height = 12, units = "cm", pointsize = 12, res = 300)


labels <- c("Unprotected", "III  - VI", "unknown",  "I & II" )

data <- multiple.taxa.PA_11_14[,c("ag_suit", "log_elevation", "log_slope", "IUCN_CAT", "SS", "SSB", "SSBS", "Species_richness")]
data <- na.omit(data)
data$IUCN_CAT <- relevel(data$IUCN_CAT, "0")

m1i <- glmer(Species_richness ~ IUCN_CAT + poly(log_elevation,1) + ag_suit
	+ (IUCN_CAT|SS) + (1|SSB) + (1|SSBS), 
	family = "poisson", data = data)

summary(m1i)

pos <- c(grep("4.5", names(fixef(m1i))),grep("7", names(fixef(m1i))),grep("1.5", names(fixef(m1i))))
y <- as.numeric(fixef(m1i)[pos])
se <- as.numeric(se.fixef(m1i)[pos])
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



text(1, 80, paste("n =", length(data$SS[which(data$IUCN_CAT == "0")]), sep = " "))
text(2, 80, paste("n =", length(data$SS[which(data$IUCN_CAT == "4.5")]), sep = " "))
text(3, 80, paste("n =", length(data$SS[which(data$IUCN_CAT == "7")]), sep = " "))
text(4, 80, paste("n =", length(data$SS[which(data$IUCN_CAT == "1.5")]), sep = " "))

arrows(seq(2,length(points),1),CI[,1],
	seq(2,length(points),1),CI[,2], code = 3, length = 0.03, angle = 90)
abline(h = 100, lty = 2)
points(points ~ c(1,2,3,4), pch = 16, col = c(1,3,3,3), cex = 1.5)

dev.off()


#add points for master plot

IUCN.plot <- data.frame(label = labels[2:4], est = points[2:4], 
		upper = CI[,1], lower = CI[,2],
		n.site = c(length(data$SS[which(data$IUCN_CAT == "4.5")]), 
			length(data$SS[which(data$IUCN_CAT == "7")]),
			length(data$SS[which(data$IUCN_CAT == "1.5")])))
sp.plot2 <- rbind(sp.plot1, IUCN.plot)





# simple species richness for Zone data


sp.tropical <- subset(multiple.taxa.PA_11_14, Zone == "Tropical")
sp.temperate <- subset(multiple.taxa.PA_11_14, Zone == "Temperate")

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
m1ztr <- glmer(Species_richness ~ Within_PA + poly(ag_suit,1)+poly(log_elevation,3)
	+ (1+Within_PA|SS)+ (1|SSBS)+ (1|SSB), family = "poisson",
	 data = data.trop)
m2ztr <- glmer(Species_richness ~ 1  +poly(ag_suit,1)+poly(log_elevation,3)
	+(1+Within_PA|SS)+ (1|SSBS)+ (1|SSB), family = "poisson", 
	 data = data.trop)
anova(m1ztr, m2ztr)
# 8.7687      1   0.003064

data.temp <- sp.temperate[,c("Within_PA", "ag_suit", "log_elevation", "log_slope", "SS", "SSB","SSBS", "Species_richness")] 
data.temp <- na.omit(data.temp)
m1zte <- glmer(Species_richness ~ Within_PA +poly(log_slope,2) + poly(log_elevation,2) 
	+ (1+Within_PA|SS)+ (1|SSBS)+ (1|SSB), family = "poisson", 
	 data = data.temp)
m2zte <- glmer(Species_richness ~ 1 + poly(log_slope,2) + poly(log_elevation,2)
	+ (1+Within_PA|SS)+ (1|SSBS)+ (1|SSB), family = "poisson", 
	 data = data.temp)
anova(m1zte, m2zte)
# 2.6042      1     0.1066


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
				n.site = c(nrow(data.trop[which(data.trop$Within_PA == "yes"),]), 
					nrow(data.temp[which(data.temp$Within_PA == "yes"),])))
sp.plot3 <- rbind(sp.plot2, a.zone)







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
best.random.p$best.random # "(1+Within_PA|SS)+ (1|SSB)"

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
best.random.i$best.random # "(1+Within_PA|SS)+ (1|SSB)"

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
#"Species_richness~poly(ag_suit,1)+poly(log_elevation,2)+poly(log_slope,3)+(1+Within_PA|SS)+(1|SSBS)+(1|SSB)"

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
#"Species_richness~poly(log_elevation,1)+poly(log_slope,1)+Within_PA+(1+Within_PA|SS)+(1|SSBS)+(1|SSB)"


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
#"Species_richness~poly(ag_suit,3)+poly(log_elevation,2)+(1+Within_PA|SS)+(1|SSBS)+(1|SSB)"


# run models
data.p <- plants[,c("Within_PA", "ag_suit", "log_elevation", "log_slope", "SS", "SSB", "SSBS", "Species_richness")]
data.p <- na.omit(data.p)
m1txp <- glmer(Species_richness ~ Within_PA + poly(ag_suit,1)+poly(log_elevation,2)+poly(log_slope,3)
	+ (Within_PA|SS)+ (1|SSB) + (1|SSBS), family = "poisson", 
	 data = data.p)
m2txp <- glmer(Species_richness ~ 1+poly(ag_suit,1)+poly(log_elevation,2)+poly(log_slope,3)
	+ (Within_PA|SS)+ (1|SSB)+ (1|SSBS), family = "poisson", 
	 data = data.p)
anova(m1txp , m2txp)
#1.761      1,13     0.1845

data.i <- inverts[,c("Within_PA", "ag_suit", "log_elevation", "log_slope", "SS", "SSB", "SSBS","Species_richness")]
data.i <- na.omit(data.i)
m1txi <- glmer(Species_richness ~ Within_PA +poly(log_elevation,1)+poly(log_slope,1) 
	+ (Within_PA|SS)+ (1|SSB)+ (1|SSBS), family = "poisson", 
	 data = data.i)
m2txi<- glmer(Species_richness ~ 1 +poly(log_elevation,1)+poly(log_slope,1) 
	+ (Within_PA|SS)+ (1|SSB)+ (1|SSBS), family = "poisson", 
	 data = data.i)
anova(m1txi, m2txi)
#4.0646      1    0.04379

data.v <- verts[,c("Within_PA", "ag_suit", "log_elevation", "log_slope", "SS", "SSB","SSBS", "Species_richness")]
data.v <- na.omit(data.v)
m1txv <- glmer(Species_richness ~ Within_PA + poly(ag_suit,3)+poly(log_elevation,2)
	+ (Within_PA|SS)+ (1|SSB)+ (1|SSBS), family = "poisson", 
	 data = data.v)
m2txv <- glmer(Species_richness ~ 1 + poly(ag_suit,3)+poly(log_elevation,2)
	+ (Within_PA|SS)+ (1|SSB)+ (1|SSBS), family = "poisson", 
	 data = data.v)
anova(m1txv, m2txv)
#3.3799      1      0.066

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
sp.plot <- rbind(sp.plot3, tax)




# master plot



tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/02_15/simple models sp rich.tif",
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




### effectiveness estimates
# benefit for species richness

b.sp <- exp(fixef(m1)[2]) -1
b.sp.max <- exp(fixef(m1)[2] + 1.96*se.fixef(m1)[2])-1
b.sp.min <- exp(fixef(m1)[2] - 1.96*se.fixef(m1)[2])-1


# if IUCN cat 1 or 2
pos <- which(names(fixef(m1i))== "IUCN_CAT1.5")
b.spIUCN1 <- exp(fixef(m1i)[pos])-1
b.spIUCN1.max <- exp(fixef(m1i)[pos] + 1.96*se.fixef(m1)[2])-1
b.spIUCN1.min <- exp(fixef(m1i)[pos] - 1.96*se.fixef(m1)[2])-1

# if IUCN cat  3 to 6
pos <- which(names(fixef(m1i))== "IUCN_CAT4.5")
b.spIUCN2 <- exp(fixef(m1i)[pos])-1
b.spIUCN2.max <- exp(fixef(m1i)[pos] + 1.96*se.fixef(m1)[2])-1
b.spIUCN2.min <- exp(fixef(m1i)[pos] - 1.96*se.fixef(m1)[2])-1


b.vals <- c(b.sp, b.sp.max, b.sp.min,
		b.spIUCN1, b.spIUCN1.max, b.spIUCN1.min,
		b.spIUCN2, b.spIUCN2.max, b.spIUCN2.min)

vals <- data.frame(name =  c("b.sp", "b.sp.max", "b.sp.min",
				"b.spIUCN1", "b.spIUCN1.max", "b.spIUCN1.min",
				"b.spIUCN2", "b.spIUCN2.max", "b.spIUCN2.min"),
			b.vals = b.vals, 
			metric = rep("sp.rich", 9), 
			NPA.abs = rep(NA,length(b.vals)),
			PA.abs = rep(NA,length(b.vals)),
			est = rep(NA,length(b.vals)))


#By 2005, we estimate that human impacts had reduced local richness by an average of 13.6% (95% CI: 9.1 – 17.8%) and 
#total abundance by 10.7% (95% CI: 3.8% gain – 23.7% reduction) compared with pre-impact times. 
#Approximately 60% of the decline in richness was independent of effects on abundance: 
#average rarefied richness has fallen by 8.1% (95% CI: 3.5 – 12.9%).
			
for(b in vals$b.vals){

benefit <- as.numeric(b)		# percentage increase in metric in PAs
PA.pct <- 14.7 				# percentage of total land area in PAs
# global loss of biodiversity (from Newbold et al)
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



### how much area needed to get effectiveness of 90 if benefit doesnt increase
### PA.abs and NPA.abs must change
### then if total loss is the same, what is the area needed

#given
#est = 1-(1-PA.abs)/(1-NPA.abs)
#NPA.abs = PA.abs*NPA.rel
#then
#est = 1-(1-PA.abs)/(1-PA.abs*NPA.rel)
#est*(1-PA.abs*NPA.rel) = (1-PA.abs*NPA.rel)-(1-PA.abs)
#est - est*PA.abs*NPA.rel = 1 - PA.abs*NPA.rel - 1 + PA.abs
#est - est*PA.abs*NPA.rel =  - PA.abs*NPA.rel  + PA.abs
#est =  - PA.abs*NPA.rel  + PA.abs + est*PA.abs*NPA.rel
#est = PA.abs(- NPA.rel  + 1 + est*NPA.rel)
#est/(- NPA.rel  + 1 + est*NPA.rel) = PA.abs

est <- vals$est[which(vals$name == "b.spIUCN1")]
benefit <- exp(fixef(m1)[2]) -1
NPA.rel <- 1-benefit
PA.abs <- est/(- NPA.rel  + 1 + est*NPA.rel)
PA.abs
PA.abs*NPA.rel

# so if higher effectiveness for same benefit, both PA and NPA abs must increase.
# how do we get current global loss given these scenarios?

#then given 
#global.int = (1 - PA.pct/100)*PA.abs*NPA.rel + (PA.pct/100)*PA.abs
#global.int == PA.abs*NPA.rel - (PA.pct/100)*PA.abs*NPA.rel  + (PA.pct/100)*PA.abs
#global.int = PA.abs*NPA.rel - (PA.pct/100)*(PA.abs*NPA.rel - PA.abs)
#(PA.pct/100)*(PA.abs*NPA.rel - PA.abs) = PA.abs*NPA.rel - global.int
#(PA.pct/100) = (PA.abs*NPA.rel - global.int)/(PA.abs*NPA.rel - PA.abs)

## nb checked that these work by using 
#est <- vals$est[which(vals$name == "b.sp")]

100*(PA.abs*NPA.rel - global.int)/(PA.abs*NPA.rel - PA.abs)
# we must unprotect everything, + another 10% of the planet. 
# this is not a useful figure. 

# with current equations, if area is increased, effectiveness goes down
# because PAs must have been doing less well if more of the planet is in them, given observed global loss 



# need to work out, if all PAs were IUCN cat 1 and still 14.7%, instead of what we have
# what would global loss be
# then how much area do we need to get this, given current benefit of PAs

# NPA abs remains the same as in normal benefit scenario
NPA.abs <- vals$NPA.abs[which(vals$name == "b.sp")]
# total area protected remains the same
PA.pct <- 14.7
# NPA.rel becomes 1-benefit of protection in IUCN category
pos <- which(names(fixef(m1i))== "IUCN_CAT1.5")
NPA.rel <- 1 - (exp(fixef(m1i)[pos])-1)

#given
#global.int = (1 - PA.pct/100)*NPA.abs + (PA.pct/100)*PA.abs
global.int = (1 - PA.pct/100)*NPA.abs + (PA.pct/100)*(NPA.abs/NPA.rel)

#now, given the normal benefit, how much additional PA area needed to get this
# need area in terms of global.int and NPA abs. 
#global.int = (1 - PA.pct/100)*NPA.abs + (PA.pct/100)*NPA.abs/NPA.rel
#global.int == NPA.abs - (PA.pct/100)*NPA.abs  + (PA.pct/100)*NPA.abs/NPA.rel
#global.int = NPA.abs - (PA.pct/100)*(NPA.abs - NPA.abs/NPA.rel)
#(PA.pct/100) = (NPA.abs - global.int)/(NPA.abs - NPA.abs/NPA.rel)

NPA.rel <- 1 - (exp(fixef(m1)[2]) -1)
100*(NPA.abs - global.int)/(NPA.abs - NPA.abs/NPA.rel)

# what if all PAs are lower management restriction (i.e. new PAs only 3 - 6 and old PAs degraded)
# ie NPA.rel is only that for III to VI

pos <- which(names(fixef(m1i))== "IUCN_CAT4.5")
NPA.rel <- 1 - (exp(fixef(m1i)[pos])-1)
100*(NPA.abs - global.int)/(NPA.abs - NPA.abs/NPA.rel)







save.image("N:\\Documents\\PREDICTS\\WDPA analysis\\RData files\\8 landuses\\simple models - sp rich.RData")

















