

rm(list=ls()) 

library(yarg)
library(roquefort)
library(rgdal)
library(sp)
library(rgeos)
library(maptools)
library(gplots)
library(ggplot2)
library(scales)
library(gridExtra)
library(optimx)
library(SDMTools)
library(data.table)
library(gamm4)
library(scatterplot3d)
library(rgl)


#setwd("C:/Users/Claudia/Documents/PREDICTS/WDPA analysis")

setwd("R:/ecocon_d/clg32/GitHub/PREDICTS_WDPA")

# load functions

source("compare_randoms.R")
source("multiplot_pairwise_all.R")
source("bray_curtis_dissimilarity.R")
source("model_select.R")


#load data
source("prep_matched.landuse_for_analysis.R")

validate <- function(x) {
  par(mfrow = c(1,2))
  plot(resid(x)~ fitted(x))
  hist(resid(x))
  par(mfrow = c(1,1))
}



construct_call<-function(responseVar,fixedStruct,randomStruct){
  return(paste(responseVar,"~",fixedStruct,"+",randomStruct,sep=""))
}




names(matched.landuse)

plants <- subset(matched.landuse, taxon_of_interest == "Plants and Fungi")
inverts.data <- subset(matched.landuse, taxon_of_interest == "Invertebrates")
verts.data <- subset(matched.landuse, taxon_of_interest == "Vertebrates")






#########
VEG
#########

plot(veg ~ Within_PA, matched.landuse)
plot(veg ~ taxon_of_interest, matched.landuse)
plot(veg ~ jitter(ag_suit), matched.landuse)

xyplot(veg ~ Within_PA|Zone, matched.landuse)
xyplot(veg ~ Within_PA|Realm, matched.landuse)
xyplot(veg ~ Within_PA|taxon_of_interest, matched.landuse)


#check for nonlinear effects
fF <- c("Within_PA") 
fT <- list("ag_suit" = "3", "log_elevation" = "3", "log_slope" = "3")
fI <- character(0)
RS <-  c("Within_PA")

veg.best.random <- compare_randoms(matched.landuse, "veg",
				fixedFactors=fF,
                         fixedTerms=fT,
                       fixedInteractions=fI,
                         otherRandoms=c("Predominant_habitat"),
				fixed_RandomSlopes = RS,
                          fitInteractions=FALSE,
				verbose=TRUE)
veg.best.random$best.random


# model select
veg.model <- model_select(all.data  = matched.landuse, 
			     responseVar = "veg", 
                       fixedFactors= fF,
                       fixedTerms= fT,
                       fixedInteractions=fI,
                       randomStruct = veg.best.random$best.random,
			     otherRandoms=c("Predominant_habitat"),
                       verbose=TRUE)
veg.model$final.call
#"veg~poly(log_elevation,3)+poly(log_slope,3)+(1+Within_PA|SS)+(1|Predominant_habitat)"




# get only data that model runs on 

veg.data <- veg.model$data

veg1 <- lmer(veg~ Within_PA + poly(log_elevation,3) + poly(log_slope,3)+(1+Within_PA|SS)+(1|SSB)+(1|Predominant_habitat), data = veg.data )
veg2 <- lmer(veg~ 1 + poly(log_elevation,3) + poly(log_slope,3)+(1+Within_PA|SS)+(1|SSB)+(1|Predominant_habitat), data =veg.data )

anova(veg1, veg2)
#0.0205      1     0.8862


tiff("N:/Documents/PREDICTS/WDPA analysis/plots/02_15/simple models matchedlanduse veg height.tif",
	width = 10, height = 10, units = "cm", pointsize = 12, res = 300)

labels <- c("Unprotected", "Protected")
y <- as.numeric(fixef(veg1)[2])
se <- as.numeric(se.fixef(veg1)[2])
yplus <- y + se*1.96
yminus <- y - se*1.96
y <-100*(y + 1) # plot as relative to 100
yplus<-100*(yplus + 1)
yminus<-100*(yminus + 1)

points <- c(100, y)
CI <- c(yplus, yminus)

plot(points ~ c(1,2), ylim = c(80,140), xlim = c(0.5,2.5), 
	bty = "l", pch = 16, col = c(1,3), cex = 1.5,
	yaxt = "n", xaxt = "n",
	ylab = "Plant height difference (± 95%CI)",
	xlab = "")
axis(1, c(1,2), labels)
axis(2, c(80,100,120,140), c(80,100,120,140))
arrows(2,CI[1],2,CI[2], code = 3, length = 0, angle = 90)
abline(h = 100, lty = 2)
points(points ~ c(1,2), pch = c(21,16), bg = "white", cex = 1.5)


text(2,80, paste("n =", length(veg.data $SS[which(veg.data$Within_PA == "yes")]), sep = " "))
text(1,80, paste("n =", length(veg.data $SS[which(veg.data$Within_PA == "no")]), sep = " "))

dev.off()


# IUCN category
# Zone
# taxon










#####
INVERTS - length derived volume
#####



plot(vol ~ Within_PA, inverts.data)
plot(vol ~ taxon_of_interest, inverts.data)
plot(vol ~ jitter(ag_suit), inverts.data)

xyplot(vol ~ Within_PA|Zone, inverts.data)
xyplot(vol ~ Within_PA|Realm, inverts.data)
xyplot(vol ~ Within_PA|taxon_of_interest, inverts.data)


#check for nonlinear effects
fF <- c("Within_PA") 
fT <- list("ag_suit" = "3", "log_elevation" = "3", "log_slope" = "3")
fI <- character(0)
RS <-  c("Within_PA")

inverts.vol.best.random <- compare_randoms(inverts.data, "vol",
				fixedFactors=fF,
                         fixedTerms=fT,
                       fixedInteractions=fI,
                         otherRandoms=c("Predominant_habitat"),
				fixed_RandomSlopes = RS,
                          fitInteractions=FALSE,
				verbose=TRUE)
inverts.vol.best.random$best.random


# model select
inverts.vol.model <- model_select(all.data  = inverts.data, 
			     responseVar = "vol", 
                       fixedFactors= fF,
                       fixedTerms= fT,
                       fixedInteractions=fI,
                       randomStruct = inverts.vol.best.random$best.random,
			     otherRandoms=c("Predominant_habitat"),
                       verbose=TRUE)
validate(inverts.vol.model$model)
inverts.vol.model$final.call
# "vol~poly(ag_suit,1)+poly(log_elevation,1)+poly(log_slope,2)+(1+Within_PA|SS)"


# get only data that model runs on 

inverts.vol.data <- inverts.vol.model$data

inverts.vol1 <- lmer(vol~ Within_PA + poly(ag_suit,1)+poly(log_elevation,1)+poly(log_slope,2)+(1+Within_PA|SS), data = inverts.vol.data)
inverts.vol2 <- lmer(vol~ 1 + poly(ag_suit,1)+poly(log_elevation,1)+poly(log_slope,2)+(1+Within_PA|SS), data = inverts.vol.data)

anova(inverts.vol1, inverts.vol2)
#0.246      1     0.6199




tiff("N:/Documents/PREDICTS/WDPA analysis/plots/02_15/simple models matchedlanduse inverts vol.tif",
	width = 10, height = 10, units = "cm", pointsize = 12, res = 300)

labels <- c("Unprotected", "Protected")
y <- as.numeric(fixef(inverts.vol1)[2])
se <- as.numeric(se.fixef(inverts.vol1)[2])
yplus <- y + se*1.96
yminus <- y - se*1.96
y <-100*(y + 1) # plot as relative to 100
yplus<-100*(yplus + 1)
yminus<-100*(yminus + 1)

points <- c(100, y)
CI <- c(yplus, yminus)

plot(points ~ c(1,2), ylim = c(80,140), xlim = c(0.5,2.5), 
	bty = "l", pch = 16, col = c(1,3), cex = 1.5,
	yaxt = "n", xaxt = "n",
	ylab = "Invertebrate volume difference (± 95%CI)",
	xlab = "")
axis(1, c(1,2), labels)
axis(2, c(80,100,120,140), c(80,100,120,140))
arrows(2,CI[1],2,CI[2], code = 3, length = 0, angle = 90)
abline(h = 100, lty = 2)
points(points ~ c(1,2), pch = c(21,16), bg = "white", cex = 1.5)


text(2,80, paste("n =", length(inverts.vol.data $SS[which(inverts.vol.data$Within_PA == "yes")]), sep = " "))
text(1,80, paste("n =", length(inverts.vol.data $SS[which(inverts.vol.data$Within_PA == "no")]), sep = " "))

dev.off()






#########
MASS
#########

(just mammals and birds and reptiles)

hist(matched.landuse$mass)
hist(log(matched.landuse$mass+1)) 
#logging it doesnt make a massive difference

names(matched.landuse)
plot(mass ~ Within_PA, matched.landuse)
plot(mass ~ taxon_of_interest, matched.landuse)

xyplot(mass ~ Within_PA|Zone, matched.landuse)
xyplot(mass ~ Within_PA|Realm, matched.landuse)
xyplot(mass ~ Within_PA|taxon_of_interest, matched.landuse)


#check for nonlinear effects
fF <- c("Within_PA") 
fT <- list("ag_suit" = "3", "log_elevation" = "3", "log_slope" = "3")
fI <- character(0)
RS <-  c("Within_PA")

mass.best.random <- compare_randoms(matched.landuse, "mass",
				fixedFactors=fF,
                         fixedTerms=fT,
                       fixedInteractions=fI,
                         otherRandoms=c("Predominant_habitat"),
				fixed_RandomSlopes = RS,
                          fitInteractions=FALSE,
				verbose=TRUE)
mass.best.random$best.random


# model select
mass.model <- model_select(all.data  = matched.landuse, 
			     responseVar = "mass", 
                       fixedFactors= fF,
                       fixedTerms= fT,
                       fixedInteractions=fI,
                       randomStruct = mass.best.random$best.random,
			     otherRandoms=c("Predominant_habitat"),
                       verbose=TRUE)
mass.model$final.call
#"mass~poly(log_slope,3)+(1+Within_PA|SS)+(1|SSB)+(1|Predominant_habitat)"


# get only data that model runs on 

mass.data <- mass.model$data

mass1 <- lmer(mass~ Within_PA + poly(log_slope,3)+(1+Within_PA|SS)+(1|SSB)+(1|Predominant_habitat), data = mass.data)
mass2 <- lmer(mass~ 1 + poly(log_slope,3)+(1+Within_PA|SS)+(1|SSB)+(1|Predominant_habitat), data = mass.data)

anova(mass1, mass2)
# 0.4367      1     0.5087


tiff("N:/Documents/PREDICTS/WDPA analysis/plots/02_15/simple models matchedlanduse mass.tif",
	width = 10, height = 10, units = "cm", pointsize = 12, res = 300)

labels <- c("Unprotected", "Protected")
y <- as.numeric(fixef(mass1)[2])
se <- as.numeric(se.fixef(mass1)[2])
yplus <- y + se*1.96
yminus <- y - se*1.96
y <-100*(y + 1) # plot as relative to 100
yplus<-100*(yplus + 1)
yminus<-100*(yminus + 1)

points <- c(100, y)
CI <- c(yplus, yminus)

plot(points ~ c(1,2), ylim = c(80,140), xlim = c(0.5,2.5), 
	bty = "l", pch = 16, col = c(1,3), cex = 1.5,
	yaxt = "n", xaxt = "n",
	ylab = "Vert mass difference (± 95%CI)",
	xlab = "")
axis(1, c(1,2), labels)
axis(2, c(80,100,120,140), c(80,100,120,140))
arrows(2,CI[1],2,CI[2], code = 3, length = 0, angle = 90)
abline(h = 100, lty = 2)
points(points ~ c(1,2), pch = c(21,16), bg = "white", cex = 1.5)


text(2,80, paste("n =", length(mass.data$SS[which(mass.data$Within_PA == "yes")]), sep = " "))
text(1,80, paste("n =", length(mass.data$SS[which(mass.data$Within_PA == "no")]), sep = " "))

dev.off()



# IUCN category
# Zone
# taxon







#####
VERTS - length derived volume (amphibs)
#####


plot(vol ~ Within_PA, verts.data)
plot(vol ~ taxon_of_interest, verts.data)
plot(vol ~ jitter(ag_suit), verts.data)

xyplot(vol ~ Within_PA|Zone, verts.data) #no temperate data
xyplot(vol ~ Within_PA|Realm, verts.data)
xyplot(vol ~ Within_PA|taxon_of_interest, verts.data)


#check for nonlinear effects
fF <- c("Within_PA") 
fT <- list("ag_suit" = "3", "log_elevation" = "3", "log_slope" = "3")
fI <- character(0)
RS <-  c("Within_PA")

verts.vol.best.random <- compare_randoms(verts.data, "vol",
				fixedFactors=fF,
                         fixedTerms=fT,
                       fixedInteractions=fI,
                         otherRandoms=c("Predominant_habitat"),
				fixed_RandomSlopes = RS,
                          fitInteractions=FALSE,
				verbose=TRUE)
verts.vol.best.random$best.random


# model select
verts.vol.model <- model_select(all.data  = verts.data, 
			     responseVar = "vol", 
                       fixedFactors= fF,
                       fixedTerms= fT,
                       fixedInteractions=fI,
                       randomStruct = verts.vol.best.random$best.random,
			     otherRandoms=c("Predominant_habitat"),
                       verbose=TRUE)
validate(verts.vol.model$model)
verts.vol.model$final.call
#  "vol~poly(log_slope,1)+(1|SS)"


# get only data that model runs on 

verts.vol.data <- verts.vol.model$data

verts.vol1 <- lmer(vol~ Within_PA + poly(log_slope,1)+(1|SS), data = verts.vol.data)
verts.vol2 <- lmer(vol~ 1 +poly(log_slope,1)+(1|SS), data = verts.vol.data)

anova(verts.vol1, verts.vol2)
#0.459      1     0.4981



tiff("N:/Documents/PREDICTS/WDPA analysis/plots/02_15/simple models matchedlanduse verts vol.tif",
	width = 10, height = 10, units = "cm", pointsize = 12, res = 300)

labels <- c("Unprotected", "Protected")
y <- as.numeric(fixef(verts.vol1)[2])
se <- as.numeric(se.fixef(verts.vol1)[2])
yplus <- y + se*1.96
yminus <- y - se*1.96
y <-100*(y + 1) # plot as relative to 100
yplus<-100*(yplus + 1)
yminus<-100*(yminus + 1)

points <- c(100, y)
CI <- c(yplus, yminus)

plot(points ~ c(1,2), ylim = c(80,140), xlim = c(0.5,2.5), 
	bty = "l", pch = 16, col = c(1,3), cex = 1.5,
	yaxt = "n", xaxt = "n",
	ylab = "Vertebrate volume difference (± 95%CI)",
	xlab = "")
axis(1, c(1,2), labels)
axis(2, c(80,100,120,140), c(80,100,120,140))
arrows(2,CI[1],2,CI[2], code = 3, length = 0, angle = 90)
abline(h = 100, lty = 2)
points(points ~ c(1,2), pch = c(21,16), bg = "white", cex = 1.5)


text(2,80, paste("n =", length(verts.vol.data $SS[which(verts.vol.data$Within_PA == "yes")]), sep = " "))
text(1,80, paste("n =", length(verts.vol.data $SS[which(verts.vol.data$Within_PA == "no")]), sep = " "))

dev.off()






