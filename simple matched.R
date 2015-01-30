#setwd("C:/Users/Claudia/Documents/PREDICTS/WDPA analysis jan 2015/WDPA analysis")

rm(list=ls())

library(lme4)
library(yarg)
library(roquefort)
library(influence.ME)

# load functions
setwd("R:/ecocon_d/clg32/GitHub/PREDICTS_WDPA")
#source("compare_randoms.R")
source("compare_randoms_lmer - with poly.R")
source("model_select.R")



source("WDPA_predicts_prep_PA_11_14_for_analysis.R")



### model species richness

# check is block useful
# yes definitely

# check polynomials for confounding variables
fF <- c("Within_PA") 
fT <- list("ag_suit" = "3", "log_slope" = "3", "log_elevation" = "3")
keepVars <- list()
fI <- character(0)
RS <-  c("Within_PA")


sp.model <- model_select(all.data  = multiple.taxa.PA_11_14, 
			     responseVar = "Species_richness", 
			     fitFamily = "poisson", 
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
			     keepVars = keepVars,
                       randomStruct = "(Within_PA|SS) + (1|SSB) + (1|SSBS)",
			     otherRandoms=character(0),
                       verbose=TRUE)

# Species_richness~poly(ag_suit,1)+poly(log_elevation,2)+Within_PA+(Within_PA|SS)+(1|SSB)+(1|SSBS)

data <- multiple.taxa.PA_11_14[,c("ag_suit", "log_elevation", "log_slope", "Within_PA", "SS", "SSB", "SSBS", "Species_richness")]
data <- na.omit(data)
m1 <- glmer(Species_richness ~ Within_PA + log_slope + poly(log_elevation,2) + ag_suit
	+ (Within_PA|SS) + (1|SSB) + (1|SSBS), 
	family = "poisson", data = data)
m2 <- glmer(Species_richness ~ 1 + log_slope + poly(log_elevation,2) + ag_suit
	+ (Within_PA|SS) + (1|SSB) + (1|SSBS), 
	family = "poisson", data = data)
anova(m1, m2)

# is ag suit significant
m1_ag <- glmer(Species_richness ~ Within_PA + log_slope + log_elevation
	+ (Within_PA|SS) + (1|SSB) + (1|SSBS), 
	family = "poisson", data = multiple.taxa.PA_11_14)

anova(m1, m1_ag)


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

data <- multiple.taxa.PA_11_14[,c("Within_PA", "SSS", "Species_richness")]
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
	family = "poisson", data = multiple.taxa.PA_11_14)

m1z <- glmer(Species_richness ~ Within_PA + Zone + log_slope + log_elevation + ag_suit
	+ (Within_PA|SS) + (1|SSB) + (1|SSBS), 
	family = "poisson", data = multiple.taxa.PA_11_14)

anova(m0z,m1z)



m2z <- update(m1z, .~. -Within_PA)
m3z <- update(m1z, .~. -Zone)

anova(m1z, m2z)
anova(m1z, m3z)

summary(m1z)






#simple species richness with IUCN cat 

multiple.taxa.PA_11_14$IUCN_CAT <- relevel(multiple.taxa.PA_11_14$IUCN_CAT, "0")

# check polynomials for confounding variables
fF <- c("IUCN_CAT") 
fT <- list("ag_suit" = "3", "log_slope" = "3", "log_elevation" = "3")
keepVars <- list()
fI <- character(0)
RS <-  c("IUCN_CAT")

# "Species_richness~IUCN_CAT+poly(ag_suit,1)+poly(log_elevation,3)+(IUCN_CAT|SS)+(1|SSB)+(1|SSBS)"
# but, doesnt converge for elevation down to quadratic

sp.model <- model_select(all.data  = multiple.taxa.PA_11_14, 
			     responseVar = "Species_richness", 
			     fitFamily = "poisson", 
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
			     keepVars = keepVars,
                       randomStruct = "(IUCN_CAT|SS) + (1|SSB) + (1|SSBS)",
			     otherRandoms=character(0),
                       verbose=TRUE)


data <- multiple.taxa.PA_11_14[,c("ag_suit", "log_elevation", "log_slope", "IUCN_CAT", "SS", "SSB", "SSBS", "Species_richness")]
data <- na.omit(data)
data$IUCN_CAT <- relevel(data$IUCN_CAT, "4.5")

m0i <- glmer(Species_richness ~ Within_PA + log_slope + log_elevation + ag_suit
	+ (Within_PA|SS) + (1|SSB) + (1|SSBS), 
	family = "poisson", data = multiple.taxa.PA_11_14)

m1i <- glmer(Species_richness ~ IUCN_CAT + log_slope + log_elevation + ag_suit
	+ (Within_PA|SS) + (1|SSB) + (1|SSBS), 
	family = "poisson", data = multiple.taxa.PA_11_14)


m2i <- glmer(Species_richness ~ 1 + log_slope + poly(log_elevation,3) + ag_suit
	+ (IUCN_CAT|SS) + (1|SSB) + (1|SSBS), 
	family = "poisson", data = data)

m3i <- glmer(Species_richness ~ IUCN_CAT + log_slope + poly(log_elevation,3) + ag_suit
	+ (IUCN_CAT|SS) + (1|SSB) + (1|SSBS), 
	family = "poisson", data = data)

anova(m1i, m0i) # ordinal doesnt make any difference
anova(m2i, m3i)
summary(m3i)

SS <- data.frame(ranef(m2i)$SS)
SS$SS <- rownames(SS)
View(SS)
pos <- order(SS[,1])
plot(SS[pos,1], ylim = c(-4,3), pch = 16)
points(SS[pos,1] +SS[pos,2], col = rgb(0.5,0.5,0.5,0.3), pch = 16)
points(SS[pos,1] +SS[pos,3], col = rgb(0.8,0.8,0.5,0.5), pch = 16)
points(SS[pos,1] +SS[pos,4], col = rgb(0.8,0.5,0.5,0.5), pch = 16)

mean(SS[,1]) ; sd(SS[,1])
mean(SS[,1] + SS[,2]) ; sd(SS[,1] + SS[,2])
mean(SS[,1] + SS[,3]) ; sd(SS[,1] + SS[,3])
mean(SS[,1] + SS[,4]) ; sd(SS[,1] + SS[,4])


# plot 

tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/01_15/simple model sp rich IUCN.tif",
	width = 15, height = 12, units = "cm", pointsize = 12, res = 300)


labels <- c("Unprotected", "IUCN I & II", "IUCN III  - VI", "unknown")

levels.IUCN <- levels(multiple.taxa.PA_11_14$IUCN_CAT)

data <- multiple.taxa.PA_11_14[,c("ag_suit", "log_elevation", "log_slope", "IUCN_CAT", "SS", "SSB", "SSBS", "Species_richness")]
data <- na.omit(data)
data$IUCN_CAT <- relevel(data$IUCN_CAT, "0")

m1i <- glmer(Species_richness ~ IUCN_CAT + log_slope + poly(log_elevation,3) + ag_suit
	+ (IUCN_CAT|SS) + (1|SSB) + (1|SSBS), 
	family = "poisson", data = data)

summary(m1i)

y <- as.numeric(fixef(m1i)[2:4])
se <- as.numeric(se.fixef(m1i)[2:4])
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


data <- multiple.taxa.PA_11_14[,c("IUCN_CAT", "SSS", "Species_richness")]
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

# check polynomials

# check polynomials for confounding variables
fF <- c("Within_PA") 
fT <- list("ag_suit" = "3", "log_slope" = "3", "log_elevation" = "3")
keepVars <- list()
fI <- character(0)
RS <-  c("Within_PA")
#Richness_rarefied~poly(ag_suit,3)+poly(log_elevation,3)+(Within_PA|SS)+(1|SSB)+(1|SSBS)

r.sp.model <- model_select(all.data  = multiple.taxa.PA_11_14, 
			     responseVar = "Richness_rarefied", 
			     fitFamily = "poisson", 
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
			     keepVars = keepVars,
                       randomStruct = "(Within_PA|SS) + (1|SSB) + (1|SSBS)",
			     otherRandoms=character(0),
                       verbose=TRUE)
r.sp.model$stats  # pvalues <0.008


data <- multiple.taxa.PA_11_14[,c("ag_suit", "log_elevation", "log_slope", "Within_PA", "SS", "SSB", "SSBS", "Richness_rarefied")]
data <- na.omit(data)
m3 <- glmer(Richness_rarefied ~ Within_PA + log_slope + poly(ag_suit,3)+poly(log_elevation,3)
	+ (Within_PA|SS)+ (1|SSB) + (1|SSBS), 
	family = "poisson", data = data)
m4 <- glmer(Richness_rarefied ~ 1 + log_slope + poly(ag_suit,3)+poly(log_elevation,3)
	+ (Within_PA|SS)+ (1|SSB) + (1|SSBS), 
	family = "poisson", data = data)
anova(m3, m4)


#no sig difference


tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/01_15/simple model rar rich.tif",
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
data <- multiple.taxa.PA_11_14[,c("Within_PA", "SSS", "Richness_rarefied")]
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



r.sp.model.IUCN <- model_select(all.data  = multiple.taxa.PA_11_14, 
			     responseVar = "Richness_rarefied", 
			     fitFamily = "poisson", 
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
			     keepVars = keepVars,
                       randomStruct = "(IUCN_CAT|SS) + (1|SSB) + (1|SSBS)",
			     otherRandoms=character(0),
                       verbose=TRUE)
r.sp.model.IUCN$stats 
# cant converge with slope less than quadratic.
# NB overall result the same even if all confounding variables non linear 

# the comparison against within_PA doesnt really make sense as theres no effect of protection
m3i <- glmer(Richness_rarefied ~ Within_PA + log_slope + log_elevation + ag_suit
	+ (Within_PA|SS) + (1|SSB)+ (1|SSBS), 
	family = "poisson", data = multiple.taxa.PA_11_14)
m4i <- glmer(Richness_rarefied ~ IUCN_CAT + log_slope + log_elevation + ag_suit
	+ (Within_PA|SS) + (1|SSB)+ (1|SSBS), 
	family = "poisson", data = multiple.taxa.PA_11_14)

# so compare against null
# null doesnt converge with the non linear confounding variables.  Test as non-linear.

m5i <- glmer(Richness_rarefied ~ 1 + log_slope + ag_suit + log_elevation
	+ (IUCN_CAT|SS) + (1|SSB)+ (1|SSBS), 
	family = "poisson", data = multiple.taxa.PA_11_14,
	control= glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
m6i <- glmer(Richness_rarefied ~ IUCN_CAT + log_slope + ag_suit + log_elevation 
	+ (IUCN_CAT|SS) + (1|SSB)+ (1|SSBS), 
	control= glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),
	family = "poisson", data = multiple.taxa.PA_11_14)

summary(m5i)

# try with ordinal data
m5i <- glmer(Richness_rarefied ~ 1 + log_slope + log_elevation + ag_suit
	+ (IUCN_CAT|SS) + (1|SSB)+ (1|SSBS), 
	family = "poisson", data = multiple.taxa.PA_11_14_ord,
	control= glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
m6i <- glmer(Richness_rarefied ~ IUCN_CAT + log_slope + log_elevation + ag_suit
	+ (IUCN_CAT|SS) + (1|SSB)+ (1|SSBS), 
	control= glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),
	family = "poisson", data = multiple.taxa.PA_11_14_ord)


# m4i above is the same as
#m4i <- glmer(Richness_rarefied ~ IUCN_CAT + Within_PA + log_slope + log_elevation + ag_suit
#	+ (Within_PA|SS)+ (1|SSB) + (1|SSBS), 
#	family = "poisson", data = multiple.taxa.PA_11_14)

anova(m3i,m4i)
anova(m5i, m6i)


summary(m6i)




# plot 

tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/01_15/simple model rar rich IUCN.tif",
	width = 15, height = 12, units = "cm", pointsize = 12, res = 300)


labels <- c("Unprotected", "IUCN I & II", "IUCN III  - VI", "unknown")

levels.IUCN <- levels(multiple.taxa.PA_11_14$IUCN_CAT)

multiple.taxa.PA_11_14$IUCN_CAT <- relevel(multiple.taxa.PA_11_14$IUCN_CAT, "0")

m6i <- glmer(Richness_rarefied ~ IUCN_CAT + log_slope + log_elevation + ag_suit
	+ (IUCN_CAT|SS) + (1|SSB) + (1|SSBS), 
	family = "poisson", data = multiple.taxa.PA_11_14,
	control= glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))


y <- as.numeric(fixef(m6i)[2:4])
se <- as.numeric(se.fixef(m6i)[2:4])
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

data <- multiple.taxa.PA_11_14[,c("IUCN_CAT", "SSS", "Richness_rarefied")]
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

# check polynomials

fF <- c("Within_PA") 
fT <- list("ag_suit" = "3", "log_slope" = "3", "log_elevation" = "3")
keepVars <- list()
fI <- character(0)
RS <-  c("Within_PA")
#log_abundance~Within_PA+(Within_PA|SS)+(1|SSB)

ab.model <- model_select(all.data  = PA_11_14, 
			     responseVar = "log_abundance", 
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
			     keepVars = keepVars,
                       randomStruct = "(Within_PA|SS) + (1|SSB)",
			     otherRandoms=character(0),
                       verbose=TRUE)

m1a <- lmer(log_abundance ~ Within_PA + log_slope + log_elevation + ag_suit
	+ (Within_PA|SS)+ (1|SSB), 
	 data =  PA_11_14)
m2a <- lmer(log_abundance ~ 1 + log_slope + log_elevation + ag_suit
	+ (Within_PA|SS)+ (1|SSB), 
	 data =  PA_11_14)
anova(m1a, m2a)

summary(m1a)

exp(fixef(m1a)[2]) # 1.127


m1a_ag <- lmer(log_abundance ~ Within_PA + log_slope + log_elevation
	+ (Within_PA|SS)+ (1|SSB), 
	 data = PA_11_14)
anova(m1a, m1a_ag)

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

dev.off()








# model abundance and zone

m0az <- lmer(log_abundance ~ Within_PA + Zone + Within_PA:Zone
	+log_slope + log_elevation + ag_suit
	+ (Within_PA|SS)+ (1|SSB) , 
	 data = PA_11_14)

m1az <- lmer(log_abundance ~ Within_PA + Zone +log_slope + log_elevation + ag_suit
	+ (Within_PA|SS)+ (1|SSB), 
	 data = PA_11_14)

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
# cant converge with nonlinear terms

ab.model.IUCN <- model_select(all.data  = PA_11_14, 
			     responseVar = "log_abundance", 
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
			     keepVars = keepVars,
                       randomStruct = "(IUCN_CAT|SS) + (1|SSB)",
			     otherRandoms=character(0),
                       verbose=TRUE)


m1ai <- lmer(log_abundance ~ Within_PA +log_slope + log_elevation + ag_suit
	+ (Within_PA|SS)+ (1|SSB), 
	 data = PA_11_14)
m2ai <- lmer(log_abundance ~ IUCN_CAT +log_slope + log_elevation + ag_suit
	+ (Within_PA|SS)+ (1|SSB), 
	 data = PA_11_14)

# m3ai doesnt converge - use ordinal stats but plot/coefficients from m4ai
m3ai <- lmer(log_abundance ~ 1 +log_slope + log_elevation + ag_suit
	+ (IUCN_CAT|SS)+ (1|SSB), 
	 data = PA_11_14,
	control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
m4ai <- lmer(log_abundance ~ IUCN_CAT +log_slope + log_elevation + ag_suit
	+ (IUCN_CAT|SS)+ (1|SSB), 
	 data = PA_11_14,
	control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

#(Intercept)    4.698465   0.221816  21.182
#IUCN_CAT1.5    0.198650   0.134044   1.482
#IUCN_CAT4.5    0.048272   0.063379   0.762
#IUCN_CAT7      0.122347   0.070013   1.747


#with ordinal - gives exactly the same results when compared ord vs not ordinal with other responses 
#but coefficients v different
m3ai <- lmer(log_abundance ~ 1 +log_slope + log_elevation + ag_suit
	+ (IUCN_CAT|SS)+ (1|SSB), 
	 data = PA_11_14_ord)
m4ai <- lmer(log_abundance ~ IUCN_CAT +log_slope + log_elevation + ag_suit
	+ (IUCN_CAT|SS)+ (1|SSB), 
	 data = PA_11_14_ord)

#(Intercept)    4.790783   0.225904  21.207
#IUCN_CAT.L     0.149823   0.090057   1.664
#IUCN_CAT.Q     0.014016   0.074037   0.189
#IUCN_CAT.C    -0.005272   0.066009  -0.080


anova(m1ai, m2ai)
anova(m3ai, m4ai)
summary(m4ai)
validate(m2ai)


# plot 

tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/01_15/simple model abundance IUCN.tif",
	width = 15, height = 12, units = "cm", pointsize = 12, res = 300)


labels <- c("Unprotected", "IUCN I & II", "IUCN III  - VI", "unknown")

levels.IUCN <- levels(PA_11_14$IUCN_CAT)

PA_11_14$IUCN_CAT <- relevel(PA_11_14$IUCN_CAT, "0")

m4ai <- lmer(log_abundance ~ IUCN_CAT + log_slope + log_elevation + ag_suit
	+ (IUCN_CAT|SS) + (1|SSB), 
	 data = PA_11_14)

#PA_11_14_ord$IUCN_CAT
#levels 0 < 4.5 < 7 < 1.5

y <- as.numeric(fixef(m4ai)[2:4])
se <- as.numeric(se.fixef(m4ai)[2:4])
yplus <- y + se*1.96
yminus <- y - se*1.96
y <-(exp(y)*100)
yplus<-(exp(yplus)*100)
yminus<-(exp(yminus)*100)

points <- c(100, y)
CI <- cbind(yplus, yminus)

plot(points ~ c(1,2,3,4), ylim = c(80,160), xlim = c(0.5,4.5),
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

data <- PA_11_14[,c("IUCN_CAT", "SSS", "log_abundance")]
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
# range~poly(ag_suit,3)+poly(log_elevation,2)+poly(log_slope,3)+(Within_PA|SS)+(1|SSB)"


r.model.IUCN <- model_select(all.data  = PA_11_14, 
			     responseVar = "range", 
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
			     keepVars = keepVars,
                       randomStruct = "(Within_PA|SS) + (1|SSB)",
			     otherRandoms=character(0),
                       verbose=TRUE)

r.model.IUCN$stats

data <- PA_11_14[,c("ag_suit", "log_elevation", "log_slope", "Within_PA", "SS", "SSB", "SSBS", "range")]
data <- na.omit(data)

m1r <- lmer(range ~ Within_PA +poly(ag_suit,3)+poly(log_elevation,2)+poly(log_slope,3)
	+ (Within_PA|SS)+ (1|SSB), 
	 data = data)
m2r <- lmer(range ~ 1 + poly(ag_suit,3)+poly(log_elevation,2)+poly(log_slope,3)
	+ (Within_PA|SS)+ (1|SSB), 
	 data = data)
anova(m1r, m2r)

summary(m1r)

exp(fixef(m1r)[2]) # 0.978


m1r_ag <- lmer(range ~ Within_PA + log_slope + log_elevation
	+ (Within_PA|SS)+ (1|SSB), 
	 data = PA_11_14)
anova(m1r, m1r_ag)


#convert to endemicity
# not this 1 - exp(fixef(m1r)[2])
#instead need difference between 1/range in protected and unprotected



# plot



#RANGE

tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/12_14/simple model range.tif",
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
axis(2, c(0.8,1,1.2,1.4), c(0.8,1,1.2,1.4))
arrows(2,CI[1],2,CI[2], code = 3, length = 0.03, angle = 90)
abline(h = 1, lty = 2)
points(points ~ c(1,2), pch = 16, col = c(1,3), cex = 1.5)




dev.off()





#ENDEMICITY

tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/01_15/simple model endemicity.tif",
	width = 10, height = 10, units = "cm", pointsize = 12, res = 300)


labels <- c("Unprotected", "Protected")
y0 <-as.numeric(fixef(m1r)[1])
e.y1 <- 1/y0
y <- as.numeric(fixef(m1r)[2])
y2 <- y0+y
e.y2 <- 1/y2

#as a percentage of outside 
e.relative <- e.y2/e.y1*100

# what does this mean?
# difference of 0.3% = e.y2 - e.y1 #0.00057
# or y2 - y0 #- 0.0022 log sq km

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

data <- PA_11_14[,c("Within_PA", "SSS", "range")]
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
# range~poly(ag_suit,3)+poly(log_elevation,2)+poly(log_slope,3)+(IUCN_CAT|SS)+(1|SSB)"


ab.model.IUCN <- model_select(all.data  = PA_11_14, 
			     responseVar = "range", 
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
			     keepVars = keepVars,
                       randomStruct = "(IUCN_CAT|SS) + (1|SSB)",
			     otherRandoms=character(0),
                       verbose=TRUE)

data <- PA_11_14[,c("ag_suit", "log_elevation", "log_slope", "IUCN_CAT", "SS", "SSB", "SSBS", "range")]
data <- na.omit(data)
data$IUCN_CAT <- relevel(data$IUCN_CAT, "0")

m1ri <- lmer(range ~ Within_PA +log_slope + log_elevation + ag_suit
	+ (Within_PA|SS)+ (1|SSB), 
	 data = PA_11_14)
m2ri <- lmer(range ~ IUCN_CAT +log_slope + log_elevation + ag_suit
	+ (Within_PA|SS)+ (1|SSB), 
	 data = PA_11_14)

m3ri <- lmer(range ~ 1 + poly(ag_suit,3)+poly(log_elevation,2)+poly(log_slope,3)
	+ (IUCN_CAT|SS)+ (1|SSB), 
	 data = data)
m4ri <- lmer(range ~ IUCN_CAT + poly(ag_suit,3)+poly(log_elevation,2)+poly(log_slope,3)
	+ (IUCN_CAT|SS)+ (1|SSB), 
	 data = data)

#ordinal - gives same significance
m3ri <- lmer(range ~ 1 +log_slope + log_elevation + ag_suit
	+ (IUCN_CAT|SS)+ (1|SSB), 
	 data = PA_11_14_ord)
m4ri <- lmer(range ~ IUCN_CAT +log_slope + log_elevation + ag_suit
	+ (IUCN_CAT|SS)+ (1|SSB), 
	 data = PA_11_14_ord)


anova(m1ri, m2ri)
anova(m3ri, m4ri)
summary(m4ri)


#PLOT

tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/01_15/simple model endemicity IUCN.tif",
	width = 15, height = 12, units = "cm", pointsize = 12, res = 300)

labels <- c("Unprotected", "IUCN I & II", "IUCN III  - VI", "unknown")

data <- PA_11_14[,c("ag_suit", "log_elevation", "log_slope", "IUCN_CAT", "SS", "SSB", "SSBS", "range")]
data <- na.omit(data)
data$IUCN_CAT <- relevel(data$IUCN_CAT, "0")

m4ri <- lmer(range ~ IUCN_CAT + poly(ag_suit,3)+poly(log_elevation,2)+poly(log_slope,3)
	+ (IUCN_CAT|SS)+ (1|SSB), 
	 data = data)

y0 <-as.numeric(fixef(m4ri)[1])
e.y1 <- 1/y0
y <- as.numeric(fixef(m4ri)[2:4])
y2 <- y0+y
e.y2 <- 1/y2

#as a percentage of outside 
e.relative <- e.y2/e.y1*100

se <- as.numeric(se.fixef(m4ri)[2:4])
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


data <- PA_11_14[,c("IUCN_CAT", "SSS", "range")]
data <- na.omit(data)
text(1, 98, paste("n =", length(data$SSS[which(data$IUCN_CAT == "0")]), sep = " "))
text(2, 98, paste("n =", length(data$SSS[which(data$IUCN_CAT == "1.5")]), sep = " "))
text(3, 98, paste("n =", length(data$SSS[which(data$IUCN_CAT == "4.5")]), sep = " "))
text(4, 98, paste("n =", length(data$SSS[which(data$IUCN_CAT == "7")]), sep = " "))


dev.off()





### model proportion threatened ###

PA_11_14_a_m_b <- read.csv("PA_11_14_amph_mamm_bird.csv")
nrow(PA_11_14_a_m_b) #2695
names(PA_11_14_a_m_b)



m <- merge(PA_11_14_a_m_b , access.1, "SSS", all.x = T)
colnames(m)[which(colnames(m) == "MEAN")] <- "access"
m <- merge(m, hpd.1, "SSS", all.x = T)
colnames(m)[which(colnames(m) == "MEAN")] <- "hpd"
m <- merge(m, elevation.1, "SSS", all.x = T)
colnames(m)[which(colnames(m) == "MEAN")] <- "elevation"
m <- merge(m, slope.1, "SSS", all.x = T)
colnames(m)[which(colnames(m) == "MEAN")] <- "slope"
m <- merge(m, ag.1, "SSS", all.x = T)


nrow(m)
PA_11_14_a_m_b  <- m





# create explanatory variables

PA_11_14_a_m_b$IUCN_CAT_number <- factor(PA_11_14_a_m_b$IUCN_CAT_number) # they arent really in an order
PA_11_14_a_m_b$log_slope <- log(PA_11_14_a_m_b$slope +1)
PA_11_14_a_m_b$log_elevation <- log(PA_11_14_a_m_b$elevation +1)
PA_11_14_a_m_b$log_hpd<- log(PA_11_14_a_m_b$hpd +1)
PA_11_14_a_m_b$log_access <- log(PA_11_14_a_m_b$access +1)
PA_11_14_a_m_b$log_GIS_AREA <- log(PA_11_14_a_m_b$GIS_AREA+1)



PA_11_14_a_m_b$y <- cbind(PA_11_14_a_m_b$abundance_VU_EN_CR, PA_11_14_a_m_b$abundance_LC_NT)


# test nonlinear confoundinf variables

fF <- c("Within_PA") 
fT <- list("ag_suit" = "3", "log_slope" = "3", "log_elevation" = "3")
keepVars <- list()
fI <- character(0)
RS <-  c("Within_PA")
# y~poly(log_elevation,3)+poly(log_slope,3)+(Within_PA|SS)+(1|SSB)+(1|SSBS)


RLS.model<- model_select(all.data  = PA_11_14_a_m_b, 
			     responseVar = "y", 
			     fitFamily = "binomial",
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
			     keepVars = keepVars,
                       randomStruct = "(Within_PA|SS) + (1|SSB)+ (1|SSBS)",
			     otherRandoms=character(0),
                       verbose=TRUE)

data <- PA_11_14_a_m_b[,c("ag_suit", "log_elevation", "log_slope", "Within_PA", "SS", "SSB", "SSBS", "y")]
data <- na.omit(data)

m1t <- glmer(y ~ Within_PA + poly(log_elevation,3)+poly(log_slope,3) + ag_suit
	+ (Within_PA|SS)+ (1|SSB) + (1|SSBS), 
	family = "binomial", data = data)
m2t <- glmer(y ~ 1 + poly(log_elevation,3)+poly(log_slope,3) + ag_suit
	+ (Within_PA|SS)+ (1|SSB) + (1|SSBS), 
	family = "binomial", data = data)
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

data <- PA_11_14_a_m_b[,c("Within_PA", "SSS", "y")]
data <- na.omit(data)
text(2,0, paste("n =", length(data$SSS[which(data$Within_PA == "yes")]), sep = " "))
text(1,0, paste("n =", length(data$SSS[which(data$Within_PA == "no")]), sep = " "))


dev.off()


# IUCN CAT
m1ti <- glmer(y ~ Within_PA + log_slope + log_elevation + ag_suit
	+ (Within_PA|SS)+ (1|SSB) + (1|SSBS), 
	family = "binomial", data = PA_11_14_a_m_b)
m2ti <- glmer(y ~ IUCN_CAT + log_slope + log_elevation + ag_suit
	+ (Within_PA|SS)+ (1|SSB) + (1|SSBS), 
	family = "binomial", data = PA_11_14_a_m_b)
anova(m1ti, m2ti)














#PA effectiveness estimates

# benefit for species richness
#
b.sp <- exp(fixef(m1)[2]) -1
b.sp.max <- (exp(fixef(m1)[2]) + 1.96*se.fixef(m1)[2])-1
b.sp.min <- (exp(fixef(m1)[2]) - 1.96*se.fixef(m1)[2])-1

b.a <- exp(fixef(m1a)[2]) -1
b.a.max <- (exp(fixef(m1a)[2]) + 1.96*se.fixef(m1a)[2])-1
b.a.min <- (exp(fixef(m1a)[2]) - 1.96*se.fixef(m1a)[2])-1

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


