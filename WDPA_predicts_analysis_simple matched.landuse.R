#setwd("C:/Users/Claudia/Documents/PREDICTS/WDPA analysis")

rm(list=ls())




library(lme4)
library(yarg)
library(roquefort)

# load functions
setwd("R:/ecocon_d/clg32/GitHub/PREDICTS_WDPA")
#source("compare_randoms.R")
source("compare_randoms_lmer - with poly.R")
source("model_select.R")


# Load dataset 
setwd("N:/Documents/PREDICTS/WDPA analysis")
matched.landuse <- read.csv("matched.landuse_11_2014.csv")
nrow(matched.landuse )#5491





### LOAD MATCHING DATA

setwd("N:/Documents/PREDICTS/WDPA analysis/matching data")

access <- read.table("bfer_1km_acc50k_11_14.txt", header = T, sep = ",")
hpd <- read.table("bfer_1km_HPD_11_14.txt", header = T, sep = ",")
elevation <- read.table("bfer_1km_mn30_elevation_11_14.txt", header = T, sep = ",")
slope <- read.table("bfer_1km_mn30_slope_11_14.txt", header = T, sep = ",")

ag_suit <- read.csv("ag_suitability_11_2014_moll.csv") 

# this is taken from extract values to points function, saved as text file, then opened in excel
# FID column deleted and -9999 values in RASTERVALU switched to NA

ag.1 <- ag_suit[,c("SSS", "RASTERVALU")] 
colnames(ag.1) <- c("SSS", "ag_suit")

#these are water, must have fallen on the edge of lakes/rivers
ag.1$ag_suit[which(ag.1$ag_suit==9)] <- NA

# switch scale to be more intuitive - higher values more agriculturally suitable
ag.1$ag_suit <- 9 - ag.1$ag_suit

access.1 <- access[,c("SSS", "MEAN")]
hpd.1 <- hpd[,c("SSS", "MEAN")]
elevation.1 <- elevation[,c("SSS", "MEAN")]
slope.1 <- slope[,c("SSS", "MEAN")]

m <- merge(matched.landuse , access.1, "SSS", all.x = T)
colnames(m)[which(colnames(m) == "MEAN")] <- "access"
m <- merge(m, hpd.1, "SSS", all.x = T)
colnames(m)[which(colnames(m) == "MEAN")] <- "hpd"
m <- merge(m, elevation.1, "SSS", all.x = T)
colnames(m)[which(colnames(m) == "MEAN")] <- "elevation"
m <- merge(m, slope.1, "SSS", all.x = T)
colnames(m)[which(colnames(m) == "MEAN")] <- "slope"
m <- merge(m, ag.1, "SSS", all.x = T)
nrow(m)
matched.landuse  <- m

# create explanatory variables

matched.landuse$IUCN_CAT_number <- factor(matched.landuse$IUCN_CAT_number) # they arent really in an order
matched.landuse$log_slope <- log(matched.landuse$slope +1)
matched.landuse$log_elevation <- log(matched.landuse$elevation +1)
matched.landuse$log_hpd<- log(matched.landuse$hpd +1)
matched.landuse$log_access <- log(matched.landuse$access +1)
matched.landuse$log_GIS_AREA <- log(matched.landuse$GIS_AREA+1)

# make IUCN cat variable of I or II vs III to VI vs unknown vs unprotected
matched.landuse$IUCN_CAT <- matched.landuse$IUCN_CAT_number 
levels(matched.landuse$IUCN_CAT) <- c(levels(matched.landuse$IUCN_CAT), "0")
matched.landuse$IUCN_CAT[which(matched.landuse$Within_PA == "no")] <- 0

#make response variables

matched.landuse$range <- matched.landuse$CWM_Geographic_range_log10_square_km
matched.landuse$mass <- matched.landuse$CWM_Mass_log10_g  ### UPDATE
matched.landuse$veg <- matched.landuse$CWM_Vegetative_height_log10_m
matched.landuse$vol <- matched.landuse$CWM_Length_derived_volume_3log10_mm
matched.landuse$log_abundance <- log(matched.landuse$Total_abundance +1)


### CREATE MULTIPLE TAXA PER STUDY DATASET 
# create dataset for species richness analysis that without all studies that are only on one taxon

setwd("N:/Documents/PREDICTS/WDPA analysis")

studies.taxa <- read.csv("Number of taxa per study split_taxa_coarse 11_2014.csv")

which(studies.taxa$number.taxa == 1)
length(which(studies.taxa$number.taxa == 1)) #25

more.than.one.taxa <- studies.taxa$SS[which(studies.taxa$number.taxa != 1)]

multiple.taxa.matched.landuse  <- subset(matched.landuse , SS %in% more.than.one.taxa)
length(multiple.taxa.matched.landuse [,1]) #5311

multiple.taxa.matched.landuse  <- droplevels(multiple.taxa.matched.landuse)




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


s.model <- model_select(all.data  = multiple.taxa.matched.landuse, 
			     responseVar = "Species_richness", 
			     fitFamily = "poisson", 
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
			     keepVars = keepVars,
                       randomStruct = "(Within_PA|SS) + (1|SSB) + (1|SSBS)",
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

# Tims plots for landuse - the y axis is % difference
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





#simple species richness with  Zone

m0z <- glmer(Species_richness ~ Within_PA + Zone + Within_PA:Zone
	+ log_slope + log_elevation + ag_suit
	+ (Within_PA|SS) + (1|SSB) + (1|SSBS), 
	family = "poisson", data = multiple.taxa.matched.landuse)

m1z <- glmer(Species_richness ~ Within_PA + Zone + log_slope + log_elevation + ag_suit
	+ (Within_PA|SS) + (1|SSB) + (1|SSBS), 
	family = "poisson", data = multiple.taxa.matched.landuse)

anova(m0z,m1z)



m2z <- update(m1z, .~. -Within_PA)
m3z <- update(m1z, .~. -Zone)

anova(m1z, m2z)
anova(m1z, m3z)

summary(m1z)






#simple species richness with IUCN cat 


fF <- c("IUCN_CAT") 
fT <- list("ag_suit" = "3", "log_slope" = "3", "log_elevation" = "3")
keepVars <- list()
fI <- character(0)
RS <-  c("IUCN_CAT")
#cant converge with non-linear confounding variables


s.model <- model_select(all.data  = multiple.taxa.matched.landuse, 
			     responseVar = "Species_richness", 
			     fitFamily = "poisson", 
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
			     keepVars = keepVars,
                       randomStruct = "(IUCN_CAT|SS) + (1|SSB) + (1|SSBS)",
			     otherRandoms=character(0),
                       verbose=TRUE)


data <- multiple.taxa.matched.landuse[,c("ag_suit", "log_elevation", "log_slope", "Within_PA", "SS", "SSB", "SSBS", "Species_richness")]
data <- na.omit(data)


m0i <- glmer(Species_richness ~ Within_PA + log_slope + log_elevation + ag_suit
	+ (Within_PA|SS) + (1|SSB) + (1|SSBS), 
	family = "poisson", data = multiple.taxa.matched.landuse)

m1i <- glmer(Species_richness ~ IUCN_CAT + log_slope + log_elevation + ag_suit
	+ (Within_PA|SS) + (1|SSB) + (1|SSBS), 
	family = "poisson", data = multiple.taxa.matched.landuse)



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
m2i <- glmer(Species_richness ~ 1  + log_slope + log_elevation + ag_suit
	+ (IUCN_CAT|SS) + (1|SSB) + (1|SSBS), 
	family = "poisson", data = multiple.taxa.matched.landuse_ord,
	control= glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
m3i <- glmer(Species_richness ~ IUCN_CAT + log_slope + log_elevation + ag_suit
	+ (IUCN_CAT|SS) + (1|SSB) + (1|SSBS), 
	family = "poisson", data = multiple.taxa.matched.landuse_ord,
	control= glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))


anova(m1i, m0i)
anova(m2i, m3i)
summary(m1i)




# plot 

tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/01_15/simple model matched.landuse sp rich IUCN.tif",
	width = 15, height = 12, units = "cm", pointsize = 12, res = 300)


labels <- c("Unprotected", "IUCN I & II", "IUCN III  - VI", "unknown")

levels.IUCN <- levels(multiple.taxa.matched.landuse$IUCN_CAT)

multiple.taxa.matched.landuse$IUCN_CAT <- relevel(multiple.taxa.matched.landuse$IUCN_CAT, "0")

m4i <- glmer(Species_richness ~ IUCN_CAT + log_slope + log_elevation + ag_suit
	+ (IUCN_CAT|SS) + (1|SSB) + (1|SSBS), 
	family = "poisson", data = multiple.taxa.matched.landuse,
	control= glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))


y <- as.numeric(fixef(m4i)[2:4])
se <- as.numeric(se.fixef(m4i)[2:4])
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
text(2, 80, paste("n =", length(data$SSS[which(data$IUCN_CAT == "1.5")]), sep = " "))
text(3, 80, paste("n =", length(data$SSS[which(data$IUCN_CAT == "4.5")]), sep = " "))
text(4, 80, paste("n =", length(data$SSS[which(data$IUCN_CAT == "7")]), sep = " "))

arrows(seq(2,length(points),1),CI[,1],
	seq(2,length(points),1),CI[,2], code = 3, length = 0.03, angle = 90)
abline(h = 100, lty = 2)
points(points ~ c(1,2,3,4), pch = 16, col = c(1,3,3,3), cex = 1.5)

dev.off()






# rarefied richness

fF <- c("Within_PA") 
fT <- list("ag_suit" = "3", "log_slope" = "3", "log_elevation" = "3")
keepVars <- list()
fI <- character(0)
RS <-  c("Within_PA")
#Richness_rarefied~poly(ag_suit,3)+poly(log_elevation,3)+(Within_PA|SS)+(1|SSB)+(1|SSBS)

r.sp.model <- model_select(all.data  = multiple.taxa.matched.landuse, 
			     responseVar = "Richness_rarefied", 
			     fitFamily = "poisson", 
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
			     keepVars = keepVars,
                       randomStruct = "(Within_PA|SS) + (1|SSB) + (1|SSBS)",
			     otherRandoms=character(0),
                       verbose=TRUE)
r.sp.model$stats # p <0.005

data <- multiple.taxa.matched.landuse[,c("ag_suit", "log_elevation", "log_slope", "Within_PA", "SS", "SSB", "SSBS", "Richness_rarefied")]
data <- na.omit(data)


m3 <- glmer(Richness_rarefied ~ Within_PA + log_slope + poly(ag_suit,3)+poly(log_elevation,3)
	+ (Within_PA|SS)+ (1|SSB) + (1|SSBS), 
	family = "poisson", data = data)
m4 <- glmer(Richness_rarefied ~ 1 + log_slope + poly(ag_suit,3)+poly(log_elevation,3)
	+ (Within_PA|SS)+ (1|SSB) + (1|SSBS), 
	family = "poisson", data = data)
anova(m3, m4)


#no sig difference


tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/01_15/simple model matched.landuse rar rich.tif",
	width = 10, height = 10, units = "cm", pointsize = 12, res = 300)


labels <- c("Unprotected", "Protected")
y <- as.numeric(fixef(m3)[2])
se <- as.numeric(se.fixef(m3)[2])
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
data <- multiple.taxa.matched.landuse[,c("Within_PA", "SSS", "Richness_rarefied")]
data <- na.omit(data)
text(2,80, paste("n =", length(data$SSS[which(data$Within_PA == "yes")]), sep = " "))
text(1,80, paste("n =", length(data$SSS[which(data$Within_PA == "no")]), sep = " "))
axis(1, c(1,2), labels)
axis(2, c(80,100,120,140), c(80,100,120,140))
arrows(2,CI[1],2,CI[2], code = 3, length = 0.03, angle = 90)
abline(h = 100, lty = 2)
points(points ~ c(1,2), pch = 16, col = c(1,3), cex = 1.5)

dev.off()





# rarefied richness with IUCN cat


fF <- c("IUCN_CAT") 
fT <- list("ag_suit" = "3", "log_slope" = "3", "log_elevation" = "3")
keepVars <- list()
fI <- character(0)
RS <-  c("IUCN_CAT")
#Richness_rarefied~poly(ag_suit,3)+poly(log_elevation,3)+(IUCN_CAT|SS)+(1|SSB)+(1|SSBS)

r.sp.model.IUCN <- model_select(all.data  = multiple.taxa.matched.landuse, 
			     responseVar = "Richness_rarefied", 
			     fitFamily = "poisson", 
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
			     keepVars = keepVars,
                       randomStruct = "(IUCN_CAT|SS) + (1|SSB) + (1|SSBS)",
			     otherRandoms=character(0),
                       verbose=TRUE)
r.sp.model.IUCN$stats # p <0.015

data <- multiple.taxa.matched.landuse[,c("ag_suit", "log_elevation", "log_slope", "IUCN_CAT", "SS", "SSB", "SSBS", "Richness_rarefied")]
data <- na.omit(data)


m3i <- glmer(Richness_rarefied ~ Within_PA + log_slope + log_elevation + ag_suit
	+ (Within_PA|SS) + (1|SSB)+ (1|SSBS), 
	family = "poisson", data = multiple.taxa.matched.landuse)
m4i <- glmer(Richness_rarefied ~ IUCN_CAT + log_slope + log_elevation + ag_suit
	+ (Within_PA|SS) + (1|SSB)+ (1|SSBS), 
	family = "poisson", data = multiple.taxa.matched.landuse)

# this is the same as
#m4i <- glmer(Richness_rarefied ~ IUCN_CAT + Within_PA + log_slope + log_elevation + ag_suit
#	+ (Within_PA|SS)+ (1|SSB) + (1|SSBS), 
#	family = "poisson", data = multiple.taxa.matched.landuse)


m5i <- glmer(Richness_rarefied ~1 + log_slope + poly(ag_suit,3)+poly(log_elevation,3)
	+ (IUCN_CAT|SS) + (1|SSB)+ (1|SSBS), 
	family = "poisson", data = data,
	control= glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
m6i <- glmer(Richness_rarefied ~ IUCN_CAT + log_slope + poly(ag_suit,3)+poly(log_elevation,3)
	+ (IUCN_CAT|SS) + (1|SSB)+ (1|SSBS), 
	family = "poisson", data = data,
	control= glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
anova(m3i,m4i)
anova(m5i,m6i)
summary(m4i)




# plot 

tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/01_15/simple model matched.landuse rar rich IUCN.tif",
	width = 15, height = 12, units = "cm", pointsize = 12, res = 300)


labels <- c("Unprotected", "IUCN I & II", "IUCN III  - VI", "unknown")

levels.IUCN <- levels(multiple.taxa.matched.landuse$IUCN_CAT)

data <- multiple.taxa.matched.landuse[,c("ag_suit", "log_elevation", "log_slope", "IUCN_CAT", "SS", "SSB", "SSBS", "Richness_rarefied")]
data <- na.omit(data)
data$IUCN_CAT <- relevel(data$IUCN_CAT, "0")

m4i <- glmer(Richness_rarefied ~ IUCN_CAT + log_slope + poly(ag_suit,3)+poly(log_elevation,3)
	+ (IUCN_CAT|SS) + (1|SSB) + (1|SSBS), 
	family = "poisson", data = data,
	control= glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))


y <- as.numeric(fixef(m4i)[2:4])
se <- as.numeric(se.fixef(m4i)[2:4])
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

data <- multiple.taxa.matched.landuse[,c("IUCN_CAT", "SSS", "Richness_rarefied")]
data <- na.omit(data)
text(1, 80, paste("n =", length(data$SSS[which(data$IUCN_CAT == "0")]), sep = " "))
text(2, 80, paste("n =", length(data$SSS[which(data$IUCN_CAT == "1.5")]), sep = " "))
text(3, 80, paste("n =", length(data$SSS[which(data$IUCN_CAT == "4.5")]), sep = " "))
text(4, 80, paste("n =", length(data$SSS[which(data$IUCN_CAT == "7")]), sep = " "))

arrows(seq(2,length(points),1),CI[,1],
	seq(2,length(points),1),CI[,2], code = 3, length = 0.03, angle = 90)
abline(h = 100, lty = 2)
points(points ~ c(1,2,3,4), pch = 16, col = c(1,3,3,3), cex = 1.5)

dev.off()


# use points to see percentages
points








### model abundance

fF <- c("Within_PA") 
fT <- list("ag_suit" = "3", "log_slope" = "3", "log_elevation" = "3")
keepVars <- list()
fI <- character(0)
RS <-  c("Within_PA")
#"log_abundance~poly(log_elevation,1)+(Within_PA|SS)+(1|SSB)"

ab.model <- model_select(all.data  = matched.landuse, 
			     responseVar = "log_abundance", 
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
			     keepVars = keepVars,
                       randomStruct = "(Within_PA|SS) + (1|SSB)",
			     otherRandoms=character(0),
                       verbose=TRUE)
ab.model$stats # p <0.015

data <- matched.landuse[,c("ag_suit", "log_elevation", "log_slope", "IUCN_CAT", "SS", "SSB", "SSBS", "log_abundance")]
data <- na.omit(data)


m1a <- lmer(log_abundance ~ Within_PA + log_slope + log_elevation + ag_suit
	+ (Within_PA|SS)+ (1|SSB), 
	 data = matched.landuse)
m2a <- lmer(log_abundance ~ 1 + log_slope + log_elevation + ag_suit
	+ (Within_PA|SS)+ (1|SSB), 
	 data = matched.landuse)
anova(m1a, m2a)

summary(m1a)

exp(fixef(m1a)[2]) 





# plot


tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/01_15/simple model matched.landuse abundance.tif",
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

data <- matched.landuse[,c("Within_PA", "SSS", "log_abundance")]
data <- na.omit(data)
text(2,80, paste("n =", length(data$SSS[which(data$Within_PA == "yes")]), sep = " "))
text(1,80, paste("n =", length(data$SSS[which(data$Within_PA == "no")]), sep = " "))

dev.off()








# model abundance and zone

m0az <- lmer(log_abundance ~ Within_PA + Zone + Within_PA:Zone
	+log_slope + log_elevation + ag_suit
	+ (Within_PA|SS)+ (1|SSB) , 
	 data = matched.landuse)

m1az <- lmer(log_abundance ~ Within_PA + Zone +log_slope + log_elevation + ag_suit
	+ (Within_PA|SS)+ (1|SSB), 
	 data = matched.landuse)

anova(m0az, m1az)


m2az <- update(m1az, .~. - Within_PA)
m3az <- update(m1az, .~. - Zone)


anova(m1az, m2az)
anova(m1az, m3az)

summary(m1az)




# model abundance and IUCN_category



fF <- c("IUCN_CAT") 
fT <- list("ag_suit" = "3", "log_slope" = "3", "log_elevation" = "3")
keepVars <- list()
fI <- character(0)
RS <-  c("IUCN_CAT")
#doesnt converge with nonlinear terms

ab.model.IUCN <- model_select(all.data  = matched.landuse, 
			     responseVar = "log_abundance", 
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
			     keepVars = keepVars,
                       randomStruct = "(IUCN_CAT|SS) + (1|SSB)",
			     otherRandoms=character(0),
                       verbose=TRUE)
ab.model.IUCN$stats # p <0.015

data <- matched.landuse[,c("ag_suit", "log_elevation", "log_slope", "IUCN_CAT", "SS", "SSB", "SSBS", "log_abundance")]
data <- na.omit(data)


m1ai <- lmer(log_abundance ~ Within_PA +log_slope + log_elevation + ag_suit
	+ (Within_PA|SS)+ (1|SSB), 
	 data = matched.landuse)
m2ai <- lmer(log_abundance ~ IUCN_CAT +log_slope + log_elevation + ag_suit
	+ (Within_PA|SS)+ (1|SSB), 
	 data = matched.landuse)

anova(m1ai, m2ai)


#doesnt converge
m3ai <- lmer(log_abundance ~ 1 +log_slope + log_elevation + ag_suit
	+ (IUCN_CAT|SS)+ (1|SSB), 
	 data = matched.landuse,
	control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
m4ai <- lmer(log_abundance ~ IUCN_CAT +log_slope + log_elevation + ag_suit
	+ (IUCN_CAT|SS)+ (1|SSB), 
	 data = matched.landuse,
	control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

#try ordinal
m3ai <- lmer(log_abundance ~ 1 +log_slope + log_elevation + ag_suit
	+ (IUCN_CAT|SS)+ (1|SSB), 
	 data = matched.landuse_ord,
	control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
m4ai <- lmer(log_abundance ~ IUCN_CAT +log_slope + log_elevation + ag_suit
	+ (IUCN_CAT|SS)+ (1|SSB), 
	 data = matched.landuse_ord,
	control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
anova(m3ai, m4ai)


summary(m2ai)
validate(m2ai)


# plot 

tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/01_15/simple model matched.landuse abundance IUCN.tif",
	width = 15, height = 12, units = "cm", pointsize = 12, res = 300)


labels <- c("Unprotected", "IUCN I & II", "IUCN III  - VI", "unknown")

levels.IUCN <- levels(matched.landuse$IUCN_CAT)

#matched.landuse$IUCN_CAT <- relevel(matched.landuse$IUCN_CAT, "0")

#m2ai <- lmer(log_abundance ~ IUCN_CAT + log_slope + log_elevation + ag_suit
#	+ (Within_PA|SS) + (1|SSB), 
#	 data = matched.landuse)


y <- as.numeric(fixef(m4ai)[c(4,2,3)])
se <- as.numeric(se.fixef(m4ai)[c(4,2,3)])
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
	ylab = "Abundance difference (% ± 95%CI)",
	xlab = "")
axis(1,seq(1,length(points),1), labels)
axis(2, c(80,100,120,140), c(80,100,120,140))
arrows(seq(2,length(points),1),CI[,1],
	seq(2,length(points),1),CI[,2], code = 3, length = 0.03, angle = 90)
abline(h = 100, lty = 2)
points(points ~ c(1,2,3,4), pch = 16, col = c(1,3,3,3), cex = 1.5)

data <- matched.landuse[,c("IUCN_CAT", "SSS", "log_abundance")]
data <- na.omit(data)
text(1, 80, paste("n =", length(data$SSS[which(data$IUCN_CAT == "0")]), sep = " "))
text(2, 80, paste("n =", length(data$SSS[which(data$IUCN_CAT == "1.5")]), sep = " "))
text(3, 80, paste("n =", length(data$SSS[which(data$IUCN_CAT == "4.5")]), sep = " "))
text(4, 80, paste("n =", length(data$SSS[which(data$IUCN_CAT == "7")]), sep = " "))

dev.off()








### model range


fF <- c("Within_PA") 
fT <- list("ag_suit" = "3", "log_slope" = "3", "log_elevation" = "3")
keepVars <- list()
fI <- character(0)
RS <-  c("Within_PA")
#range~poly(ag_suit,3)+poly(log_elevation,2)+poly(log_slope,3)+(Within_PA|SS)+(1|SSB) 



range.model <- model_select(all.data  = matched.landuse, 
			     responseVar = "range", 
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
			     keepVars = keepVars,
                       randomStruct = "(Within_PA|SS) + (1|SSB)",
			     otherRandoms=character(0),
                       verbose=TRUE)
range.model$stats # p <0.015

data <- matched.landuse[,c("ag_suit", "log_elevation", "log_slope", "Within_PA", "SS", "SSB", "SSBS", "range")]
data <- na.omit(data)

m1r <- lmer(range ~ Within_PA + poly(ag_suit,3)+poly(log_elevation,2)+poly(log_slope,3)
	+ (Within_PA|SS)+ (1|SSB), 
	 data = data)
m2r <- lmer(range ~ 1 +poly(ag_suit,3)+poly(log_elevation,2)+poly(log_slope,3)
	+ (Within_PA|SS)+ (1|SSB), 
	 data = data)
anova(m1r, m2r)

summary(m1r)

exp(fixef(m1r)[2]) # 0.978

#convert to endemicity
# not this 1 - exp(fixef(m1r)[2])
#instead need difference between 1/range in protected and unprotected



# plot



#RANGE

tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/01_15/simple model matched.landuse range.tif",
	width = 10, height = 10, units = "cm", pointsize = 12, res = 300)

# no link function so plot as relative, not percentage

labels <- c("Unprotected", "Protected")
y <- as.numeric(fixef(m1r)[2])
se <- as.numeric(se.fixef(m1r)[2])
yplus <- y + se*1.96
yminus <- y - se*1.96
y <-(y + 1) # plot as relative to 1
yplus<-(yplus + 1)
yminus<-(yminus + 1)

points <- c(1, y)
CI <- c(yplus, yminus)

plot(points ~ c(1,2), ylim = c(0.8,1.5), xlim = c(0.5,2.5), 
	bty = "l", pch = 16, col = c(1,3), cex = 1.5,
	yaxt = "n", xaxt = "n",
	ylab = "Relative CWM range difference (± 95%CI)",
	xlab = "")
axis(1, c(1,2), labels)
axis(2, c(0.8,1,1.2,1.4), c(80,100,120,140))
arrows(2,CI[1],2,CI[2], code = 3, length = 0.03, angle = 90)
abline(h = 1, lty = 2)
points(points ~ c(1,2), pch = 16, col = c(1,3), cex = 1.5)




dev.off()





#ENDEMICITY

tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/01_15/simple model matched.landuse endemicity.tif",
	width = 10, height = 10, units = "cm", pointsize = 12, res = 300)


labels <- c("Unprotected", "Protected")
y0 <-as.numeric(fixef(m1r)[1])
e.y1 <- 1/y0
y <- as.numeric(fixef(m1r)[2])
y2 <- y0+y
e.y2 <- 1/y2

#as a percentage of outside 
e.relative <- e.y2/e.y1*100

se <- as.numeric(se.fixef(m1r)[2])
y2plus <- y0 + y + se*1.96
e.y2plus <- 1/y2plus
e.relative.plus <- e.y2plus/e.y1*100

y2minus <- y0 + y - se*1.96
e.y2minus <- 1/y2minus
e.relative.minus <- e.y2minus/e.y1*100




points <- c(100, e.relative)
CI <- c(e.relative.plus, e.relative.minus)

plot(points ~ c(1,2), ylim = c(98,102), xlim = c(0.5,2.5), 
	bty = "l", pch = 16, col = c(1,3), cex = 1.5,
	yaxt = "n", xaxt = "n",
	ylab = "Endemicity difference (% ± 95%CI)",
	xlab = "")
axis(1, c(1,2), labels)
axis(2, c(99,100,101,102), c(99,100,101,102))
arrows(2,CI[1],2,CI[2], code = 3, length = 0.03, angle = 90)
abline(h = 100, lty = 2)
points(points ~ c(1,2), pch = 16, col = c(1,3), cex = 1.5)

data <- matched.landuse[,c("Within_PA", "SSS", "range")]
data <- na.omit(data)
text(2,98, paste("n =", length(data$SSS[which(data$Within_PA == "yes")]), sep = " "))
text(1,98, paste("n =", length(data$SSS[which(data$Within_PA == "no")]), sep = " "))

dev.off()





### range and IUCN category


fF <- c("IUCN_CAT") 
fT <- list("ag_suit" = "3", "log_slope" = "3", "log_elevation" = "3")
keepVars <- list()
fI <- character(0)
RS <-  c("IUCN_CAT")
#doesnt converge with nonlinear terms

range.model.IUCN <- model_select(all.data  = matched.landuse, 
			     responseVar = "range", 
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
			     keepVars = keepVars,
                       randomStruct = "(IUCN_CAT|SS) + (1|SSB)",
			     otherRandoms=character(0),
                       verbose=TRUE)


m1ri <- lmer(range ~ Within_PA +log_slope + log_elevation + ag_suit
	+ (Within_PA|SS)+ (1|SSB), 
	 data = matched.landuse)
m2ri <- lmer(range ~ IUCN_CAT +log_slope + log_elevation + ag_suit
	+ (Within_PA|SS)+ (1|SSB), 
	 data = matched.landuse)

#doesnt converge
m3ri <- lmer(range ~ 1 +log_slope + log_elevation + ag_suit
	+ (IUCN_CAT|SS)+ (1|SSB), 
	 data = matched.landuse)
m4ri <- lmer(range ~ IUCN_CAT +log_slope + log_elevation + ag_suit
	+ (IUCN_CAT|SS)+ (1|SSB), 
	 data = matched.landuse,
	control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

#try ordinal
#also doesnt converge
m3ri <- lmer(range ~ 1 +log_slope + log_elevation + ag_suit
	+ (IUCN_CAT|SS)+ (1|SSB), 
	 data = matched.landuse_ord)
m4ri <- lmer(range ~ IUCN_CAT +log_slope + log_elevation + ag_suit
	+ (IUCN_CAT|SS)+ (1|SSB), 
	 data = matched.landuse_ord,
	control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))


anova(m1ri, m2ri)
anova(m3ri, m4ri)
summary(m4ri)


#PLOT

tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/01_15/simple model matched.landuse endemicity IUCN.tif",
	width = 15, height = 12, units = "cm", pointsize = 12, res = 300)

labels <- c("Unprotected", "IUCN I & II", "IUCN III  - VI", "unknown")

#matched.landuse$IUCN_CAT <- relevel( matched.landuse$IUCN_CAT, "0")

#m4ri <- lmer(range ~ IUCN_CAT +log_slope + log_elevation + ag_suit
#	+ (Within_PA|SS)+ (1|SSB), 
#	 data = matched.landuse)

y0 <-as.numeric(fixef(m4ri)[1])
e.y1 <- 1/y0
y <- as.numeric(fixef(m4ri)[c(4,2,3)])
y2 <- y0+y
e.y2 <- 1/y2

#as a percentage of outside 
e.relative <- e.y2/e.y1*100

se <- as.numeric(se.fixef(m4ri)[c(4,2,3)])
y2plus <- y0 + y + se*1.96
e.y2plus <- 1/y2plus
e.relative.plus <- e.y2plus/e.y1*100

y2minus <- y0 + y - se*1.96
e.y2minus <- 1/y2minus
e.relative.minus <- e.y2minus/e.y1*100




points <- c(100, e.relative)
CI <- cbind(e.relative.plus, e.relative.minus)

plot(points ~ c(1,2,3,4), ylim = c(98,102), xlim = c(0.5,4.5), 
	bty = "l", pch = 16, col = c(1,3,3,3), cex = 1.5,
	yaxt = "n", xaxt = "n",
	ylab = "Endemicity difference (% ± 95%CI)",
	xlab = "")
axis(1, c(1,2,3,4), labels)
axis(2, c(99,100,101,102), c(99,100,101,102))
arrows(c(2,3,4),CI[,1],c(2,3,4),CI[,2], code = 3, length = 0.03, angle = 90)
abline(h = 100, lty = 2)
points(points ~ c(1,2,3,4), pch = 16, col = c(1,3,3,3), cex = 1.5)


data <- matched.landuse[,c("IUCN_CAT", "SSS", "range")]
data <- na.omit(data)
text(1, 98, paste("n =", length(data$SSS[which(data$IUCN_CAT == "0")]), sep = " "))
text(2, 98, paste("n =", length(data$SSS[which(data$IUCN_CAT == "1.5")]), sep = " "))
text(3, 98, paste("n =", length(data$SSS[which(data$IUCN_CAT == "4.5")]), sep = " "))
text(4, 98, paste("n =", length(data$SSS[which(data$IUCN_CAT == "7")]), sep = " "))


dev.off()



### model proportion threatened ###

matched.landuse_a_m_b <- read.csv("matched.landuse_amph_mamm_bird.csv")
nrow(matched.landuse_a_m_b) #2695
names(matched.landuse_a_m_b)



m <- merge(matched.landuse_a_m_b , access.1, "SSS", all.x = T)
colnames(m)[which(colnames(m) == "MEAN")] <- "access"
m <- merge(m, hpd.1, "SSS", all.x = T)
colnames(m)[which(colnames(m) == "MEAN")] <- "hpd"
m <- merge(m, elevation.1, "SSS", all.x = T)
colnames(m)[which(colnames(m) == "MEAN")] <- "elevation"
m <- merge(m, slope.1, "SSS", all.x = T)
colnames(m)[which(colnames(m) == "MEAN")] <- "slope"
m <- merge(m, ag.1, "SSS", all.x = T)


nrow(m)
matched.landuse_a_m_b  <- m





# create explanatory variables

matched.landuse_a_m_b$IUCN_CAT_number <- factor(matched.landuse_a_m_b$IUCN_CAT_number) # they arent really in an order
matched.landuse_a_m_b$log_slope <- log(matched.landuse_a_m_b$slope +1)
matched.landuse_a_m_b$log_elevation <- log(matched.landuse_a_m_b$elevation +1)
matched.landuse_a_m_b$log_hpd<- log(matched.landuse_a_m_b$hpd +1)
matched.landuse_a_m_b$log_access <- log(matched.landuse_a_m_b$access +1)
matched.landuse_a_m_b$log_GIS_AREA <- log(matched.landuse_a_m_b$GIS_AREA+1)



matched.landuse_a_m_b$y <- cbind(matched.landuse_a_m_b$abundance_VU_EN_CR, matched.landuse_a_m_b$abundance_LC_NT)


m1t <- glmer(y ~ Within_PA + log_slope + log_elevation + ag_suit
	+ (Within_PA|SS)+ (1|SSB) + (1|SSBS), 
	family = "binomial", data = matched.landuse_a_m_b)
m2t <- glmer(y ~ 1 + log_slope + log_elevation + ag_suit
	+ (Within_PA|SS)+ (1|SSB) + (1|SSBS), 
	family = "binomial", data = matched.landuse_a_m_b)
anova(m1t, m2t)

summary(m1t)


#PLOT
tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/01_15/simple model prop threat.tif",
	width = 10, height = 10, units = "cm", pointsize = 12, res = 300)

labels <- c("Unprotected", "Protected")
y <- as.numeric(fixef(m1t)[2])
se <- as.numeric(se.fixef(m1t)[2])
yplus <- y + se*1.96
yminus <- y - se*1.96

intercept<-fixef(m1t)[1]
true.intercept <- 1/(1+exp(-(intercept)))
true.y <- 1/(1+exp(-(intercept+y)))
true.yplus <- 1/(1+exp(-(intercept+yplus)))
true.yminus <- 1/(1+exp(-(intercept+yminus)))

# this is code from plot_lu_effects -
# why is it backtranformed absolute value of non-ref level/back transformed value of ref variable?
# because!!!
# dont want the percentage of the difference, want it as 100+x percent. (Tims code then -100 at the end to get difference). 

y<-((true.y/true.intercept)*100)
yplus<-((true.yplus/true.intercept)*100)
yminus<-((true.yminus/true.intercept)*100)

points <- c(100, y)
CI <- c(yplus, yminus)

plot(points ~ c(1,2), ylim = c(0,500), xlim = c(0.5,2.5), 
	bty = "l", pch = 16, col = c(1,3), cex = 1.5,
	xaxt = "n", #yaxt = "n", 
	ylab = "Proportion threatened species difference (% ± 95%CI)",
	xlab = "")
axis(1, c(1,2), labels)
#axis(2, c(80,100,120,140), c(80,100,120,140))
arrows(2,CI[1],2,CI[2], code = 3, length = 0.03, angle = 90)
abline(h = 100, lty = 2)
points(points ~ c(1,2), pch = 16, col = c(1,3), cex = 1.5)

data <- matched.landuse_a_m_b[,c("Within_PA", "SSS", "y")]
data <- na.omit(data)
text(2,0, paste("n =", length(data$SSS[which(data$Within_PA == "yes")]), sep = " "))
text(1,0, paste("n =", length(data$SSS[which(data$Within_PA == "no")]), sep = " "))


dev.off()


# IUCN CAT
m1ti <- glmer(y ~ Within_PA + log_slope + log_elevation + ag_suit
	+ (Within_PA|SS)+ (1|SSB) + (1|SSBS), 
	family = "binomial", data = matched.landuse_a_m_b)
m2ti <- glmer(y ~ IUCN_CAT + log_slope + log_elevation + ag_suit
	+ (Within_PA|SS)+ (1|SSB) + (1|SSBS), 
	family = "binomial", data = matched.landuse_a_m_b)
anova(m1ti, m2ti)














#PA effectiveness estimates

# benefit for species richness
#
b.sp <- exp(fixef(m1)[2]) -1
b.sp.max <- (exp(fixef(m1)[2]) + 1.96*se.fixef(m1)[2])-1
b.sp.min <- (exp(fixef(m1)[2]) - 1.96*se.fixef(m1)[2])-1

b.a <- exp(fixef(m1a)[2]) -1
b.a.max <- (exp(fixef(m1a)[2]) + 1.96*se.fixef(m1)[2])-1
b.a.min <- (exp(fixef(m1a)[2]) - 1.96*se.fixef(m1)[2])-1

b.vals <- c(b.sp, b.sp.max, b.sp.min, b.a, b.a.max, b.a.min)
all.results <- vector()

b <- b.a.max
b <- b.a


for(b in b.vals){

benefit <- as.numeric(b)		# percentage increase in metric in PAs
PA.pct <- 13 				# percentage of total land area in PAs
global.loss <- 0.116			# global loss of biodiversity (from Newbold et al)

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

# if pristine is 1, then where between 1 and NPA.abs does PA.abs fall?
# get difference between PA.abs and NPA.abs as a percentage of NPA.abs
#((1-NPA.abs)-(1-PA.abs))/(1-NPA.abs)

est <- 1-(1-PA.abs)/(1-NPA.abs)

# ie
#est <- ((1-NPA.abs)-(1-PA.abs))/(1-NPA.abs)


all.results <- c(all.results, est)
}

d <- data.frame(estimates= all.results, 
	values = c("b.sp", "b.sp.max", "b.sp.min", "b.a", "b.a.max", "b.a.min"))

write.csv(t(d), "simple.effectiveness.estimates.csv")






# Getting andys numbers - my approach written before I got Andys code. 


# so, if PAs are 7.5% higher in species richness
# Sites outside PAs have 92.5% as many species as sites inside
# we also know 13% of land surface is in Protected Areas
# Globally PREDICTS (nature MS) estimates mean net loss of species from terrestrial sites 
	# across all landuses at 12.9% (relative to pristine)

# what are Protected and Unprotected relevant to pristine
# n = % sp  in unprotected relative to pristine
# y = % sp  in protected relative to pristine

# total loss  = 87% loss in unprotected and 13% loss in protected
# simultaneous equations are
#1 - 0.129 = 0.87*n + 0.13*y
#n * 1.075 = y


# so
#1 - 0.129 = 0.87*n + 0.13*1.075*n


#then
n <- 1 - 0.129/(0.87 + 0.13*1.075) 
y <- n * 1.075

n #0.872
y #0.927

loss of species in unprotected relative to pristine = 0.1277544
loss of species in protected relative to pristine = 

Simultaneous equations suggest:
#PA sites have 93.2% as many species as pristine sites; non-PA sites have 86.2%

# 0% effective is 86.2%, 100% effective is the same as pristine
100 - 86.2 # 13.8
93.2-86.2  #7
7/13.8 = 50.7%
# PAs are therefore 51% effective in retaining site-level species richness

# to vary parameters

total.loss <- 1 - 0.129
p.protect <- 0.13
p.unprotect <- 1- p.protect
x <- as.numeric(exp(fixef(m1)[2]))

n <- total.loss/(p.unprotect + p.protect*x) 
y <- n * x

a <- 1 - n
b <- y - n
b/a


