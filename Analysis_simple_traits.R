

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
source("model_select.R")


#load data
source("prep_PA_11_14_for_analysis.R")

validate <- function(x) {
  par(mfrow = c(1,2))
  plot(resid(x)~ fitted(x))
  hist(resid(x))
  par(mfrow = c(1,1))
}



construct_call<-function(responseVar,fixedStruct,randomStruct){
  return(paste(responseVar,"~",fixedStruct,"+",randomStruct,sep=""))
}


load("N:\\Documents\\PREDICTS\\WDPA analysis\\RData files\\8 landuses\\simple models - traits.RData")
> 


names(PA_11_14)

plants <- subset(PA_11_14, taxon_of_interest == "Plants and Fungi")
inverts.data <- subset(PA_11_14, taxon_of_interest == "Invertebrates")
verts.data <- subset(PA_11_14, taxon_of_interest == "Vertebrates")






#########
VEG
#########

plot(veg ~ Within_PA, PA_11_14)
plot(veg ~ taxon_of_interest, PA_11_14)
plot(veg ~ jitter(ag_suit), PA_11_14)

xyplot(veg ~ Within_PA|Zone, PA_11_14)
xyplot(veg ~ Within_PA|Realm, PA_11_14)
xyplot(veg ~ Within_PA|taxon_of_interest, PA_11_14)


#check for nonlinear effects
fF <- c("Within_PA") 
fT <- list("ag_suit" = "3", "log_elevation" = "3", "log_slope" = "3")
fI <- character(0)
RS <-  c("Within_PA")

veg.best.random <- compare_randoms(PA_11_14, "veg",
				fixedFactors=fF,
                         fixedTerms=fT,
                       fixedInteractions=fI,
                         otherRandoms=c("Predominant_habitat"),
				fixed_RandomSlopes = RS,
                          fitInteractions=FALSE,
				verbose=TRUE)
veg.best.random$best.random


# model select
veg.model <- model_select(all.data  = PA_11_14, 
			     responseVar = "veg", 
                       fixedFactors= fF,
                       fixedTerms= fT,
                       fixedInteractions=fI,
                       randomStruct = veg.best.random$best.random,
			     otherRandoms=c("Predominant_habitat"),
                       verbose=TRUE)
veg.model$final.call
"veg~poly(ag_suit,2)+poly(log_elevation,3)+poly(log_slope,3)+(1+Within_PA|SS)+(1|SSB)+(1|Predominant_habitat)"
 




# get only data that model runs on 

veg.data <- veg.model$data

veg1 <- lmer(veg~ Within_PA + poly(ag_suit,2)+poly(log_elevation,3)+poly(log_slope,3)+(1+Within_PA|SS)+(1|SSB)+(1|Predominant_habitat), data = veg.data )
veg2 <- lmer(veg~ 1 + poly(ag_suit,2)+poly(log_elevation,3)+poly(log_slope,3)+(1+Within_PA|SS)+(1|SSB)+(1|Predominant_habitat), data =veg.data )

anova(veg1, veg2)
#4e-04      1     0.9845


tiff("N:/Documents/PREDICTS/WDPA analysis/plots/02_15/simple models veg height.tif",
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


### try with landuse
veg.data2 <- PA_11_14[,c("LUPA", "veg", "log_elevation", "log_slope", "SS", "SSB", 
		"Predominant_habitat", "Within_PA")]
veg.data2 <- na.omit(veg.data2)

table(veg.data2$Within_PA, veg.data2$Predominant_habitat)
#no urban outside
unique(veg.data2$LUPA)

#full interaction drops coeffcient
#veg_1 <- lmer(veg ~ Within_PA + Predominant_habitat +  Within_PA:Predominant_habitat + 
#	poly(log_elevation,3)+poly(log_slope,3)+
#	(1+Within_PA|SS)+(1|SSB), data = veg.data2)

veg_1 <- lmer(veg ~ LUPA + 
	poly(log_elevation,3)+poly(log_slope,3)+
	(1+Within_PA|SS)+(1|SSB), data = veg.data2)
veg_2 <- lmer(veg~  Within_PA + Predominant_habitat + 
	poly(log_elevation,3)+poly(log_slope,3)+
	(1+Within_PA|SS)+(1|SSB), data = veg.data2)

anova(veg_1, veg_2)

#for plotting
veg_LUPA <- lmer(veg ~ LUPA +
	poly(log_elevation,3)+poly(log_slope,3)+
	(1+Within_PA|SS)+(1|SSB), data = veg.data2)

summary(veg_LUPA)

tiff("N:/Documents/PREDICTS/WDPA analysis/plots/02_15/LUPA veg.tif",
	width = 15, height = 10, units = "cm", pointsize = 12, res = 300)

plotLUPA(model = veg_LUPA,
		responseVar = "plant height", 
		data = veg.data2,
		xvar = "LUPA",
		xvar.order = c("IN Primary Vegetation",
			"OUT Mature secondary vegetation",
			"IN Mature secondary vegetation",
			"OUT Intermediate secondary vegetation",
			"IN Intermediate secondary vegetation",
			"OUT Young secondary vegetation",
			"IN Young secondary vegetation",
			"OUT Plantation forest",
			"IN Plantation forest",
			"IN Cropland",
			"OUT Cropland",
			"IN Pasture",
			"OUT Pasture")
			#"OUT Urban",
			#"IN Urban"
			,
		col.key = NULL,
		logLink = "n",
		seMultiplier=1.96,
		forPaper = T,
		ylims = c(-1.2,1.2))

dev.off()










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
#vol~poly(ag_suit,1)+poly(log_elevation,3)+poly(log_slope,1)+(1+Within_PA|SS)+(1|SSB)+(1|Predominant_habitat)"




# get only data that model runs on 

inverts.vol.data <- inverts.vol.model$data

inverts.vol1 <- lmer(vol~ Within_PA +poly(ag_suit,1)+poly(log_elevation,3)+poly(log_slope,1)
	+(1+Within_PA|SS)+(1|SSB)+(1|Predominant_habitat), data = inverts.vol.data)
inverts.vol2 <- lmer(vol~ 1 + poly(ag_suit,1)+poly(log_elevation,3)+poly(log_slope,1)
	+(1+Within_PA|SS)+(1|SSB)+(1|Predominant_habitat), data = inverts.vol.data)

anova(inverts.vol1, inverts.vol2)
#0.0206      1     0.8859




tiff("N:/Documents/PREDICTS/WDPA analysis/plots/02_15/simple models inverts vol.tif",
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




### try with landuse
inverts.vol.data2 <- inverts.data[,c("LUPA", "vol", "log_elevation", "log_slope", "SS", "SSB", 
		"Predominant_habitat", "Within_PA")]
inverts.vol.data2 <- na.omit(inverts.vol.data2)

table(inverts.vol.data2$Within_PA, inverts.vol.data2$Predominant_habitat)
#no plantation inside
unique(inverts.vol.data2$LUPA)

#full interaction drops coeffcient
#vol_1 <- lmer(vol ~ Within_PA + Predominant_habitat +  Within_PA:Predominant_habitat + 
#	poly(log_elevation,3)+poly(log_slope,3)+
#	(1+Within_PA|SS)+(1|SSB), data = inverts.vol.data2)

inverts.vol_1 <- lmer(vol ~ LUPA + 
	poly(log_elevation,3)+poly(log_slope,3)+
	(1+Within_PA|SS)+(1|SSB), data = inverts.vol.data2)
inverts.vol_2 <- lmer(vol~  Within_PA + Predominant_habitat + 
	poly(log_elevation,3)+poly(log_slope,3)+
	(1+Within_PA|SS)+(1|SSB), data = inverts.vol.data2)

anova(inverts.vol_1, inverts.vol_2)

#for plotting
inverts.vol_LUPA <- lmer(vol ~ LUPA +
	poly(log_elevation,3)+poly(log_slope,3)+
	(1+Within_PA|SS)+(1|SSB), data = inverts.vol.data2)

summary(inverts.vol_LUPA)

tiff("N:/Documents/PREDICTS/WDPA analysis/plots/02_15/LUPA inverts vol.tif",
	width = 15, height = 10, units = "cm", pointsize = 12, res = 300)

plotLUPA(model = inverts.vol_LUPA,
		responseVar = "invertebrate volume", 
		data = inverts.vol.data2,
		xvar = "LUPA",
		xvar.order = c("IN Primary Vegetation",
			"OUT Mature secondary vegetation",
			"IN Mature secondary vegetation",
			"OUT Intermediate secondary vegetation",
			"IN Intermediate secondary vegetation",
			"OUT Young secondary vegetation",
			"IN Young secondary vegetation",
			#"OUT Plantation forest",
			#"IN Plantation forest",
			"IN Cropland",
			"OUT Cropland",
			"IN Pasture",
			"OUT Pasture",
			"OUT Urban",
			"IN Urban")
			,
		labels = rep(c("Primary Vegetation",
			"Mature secondary vegetation",
			"Intermediate secondary vegetation",
			"Young secondary vegetation",
			"Cropland",
			"Pasture",
			"Urban"), each = 2),
		col.key = NULL,
		logLink = "n",
		seMultiplier=1.96,
		forPaper = T,
		ylims = c(-1.2,1.2))

dev.off()


#########
MASS
#########

(just mammals and birds and reptiles)

hist(PA_11_14$mass)
hist(log(PA_11_14$mass+1)) 
#logging it doesnt make a massive difference


plot(mass ~ Within_PA, PA_11_14)
plot(mass ~ taxon_of_interest, PA_11_14)

xyplot(mass ~ Within_PA|Zone, PA_11_14)
xyplot(mass ~ Within_PA|Realm, PA_11_14)
xyplot(mass ~ Within_PA|taxon_of_interest, PA_11_14)


#check for nonlinear effects
fF <- c("Within_PA") 
fT <- list("ag_suit" = "3", "log_elevation" = "3", "log_slope" = "3")
fI <- character(0)
RS <-  c("Within_PA")

mass.best.random <- compare_randoms(PA_11_14, "mass",
				fixedFactors=fF,
                         fixedTerms=fT,
                       fixedInteractions=fI,
                         otherRandoms=c("Predominant_habitat"),
				fixed_RandomSlopes = RS,
                          fitInteractions=FALSE,
				verbose=TRUE)
mass.best.random$best.random


# model select
mass.model <- model_select(all.data  = PA_11_14, 
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
#0.5985      1     0.4391


tiff("N:/Documents/PREDICTS/WDPA analysis/plots/02_15/simple models mass.tif",
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


### try with landuse
mass.data2 <- PA_11_14[,c("LUPA", "mass", "log_elevation", "log_slope", "SS", "SSB", 
		"Predominant_habitat", "Within_PA")]
mass.data2 <- na.omit(mass.data2)

mass_1 <- lmer(mass ~ Within_PA + Predominant_habitat +  Within_PA:Predominant_habitat + 
	poly(log_elevation,3)+poly(log_slope,3)+
	(1+Within_PA|SS)+(1|SSB), data = mass.data2)
mass_1 <- lmer(mass ~ LUPA + 
	poly(log_elevation,3)+poly(log_slope,3)+
	(1+Within_PA|SS)+(1|SSB), data = mass.data2)
mass_2 <- lmer(mass~  Within_PA + Predominant_habitat + 
	poly(log_elevation,3)+poly(log_slope,3)+
	(1+Within_PA|SS)+(1|SSB), data = mass.data2)

anova(mass_ab1, mass_ab2)

#for plotting
mass_LUPA <- lmer(mass ~ LUPA +
	poly(log_elevation,3)+poly(log_slope,3)+
	(1+Within_PA|SS)+(1|SSB), data = mass.data2)
summary(mass_LUPA)

tiff("N:/Documents/PREDICTS/WDPA analysis/plots/02_15/LUPA mass.tif",
	width = 15, height = 10, units = "cm", pointsize = 12, res = 300)

plotLUPA(model = mass_LUPA,
		responseVar = " vertebrate mass", 
		data = mass.data2,
		xvar = "LUPA",
		xvar.order = c("IN Primary Vegetation",
			"OUT Mature secondary vegetation",
			"IN Mature secondary vegetation",
			"OUT Intermediate secondary vegetation",
			"IN Intermediate secondary vegetation",
			"OUT Young secondary vegetation",
			"IN Young secondary vegetation",
			"OUT Plantation forest",
			"IN Plantation forest",
			"IN Cropland",
			"OUT Cropland",
			"IN Pasture",
			"OUT Pasture",
			"OUT Urban",
			"IN Urban"),
		col.key = NULL,
		logLink = "n",
		seMultiplier=1.96,
		forPaper = T,
		ylims = c(-1.2,1.2))

dev.off()






### try just mammals and birds data

names(PA_11_14)

#check for nonlinear effects
fF <- c("Within_PA") 
fT <- list("ag_suit" = "3", "log_elevation" = "3", "log_slope" = "3")
fI <- character(0)
RS <-  c("Within_PA")

mass_ab.best.random <- compare_randoms(PA_11_14, "CWM_Adult_wet_mass_log10_g",
				fixedFactors=fF,
                         fixedTerms=fT,
                       fixedInteractions=fI,
                         otherRandoms=c("Predominant_habitat"),
				fixed_RandomSlopes = RS,
                          fitInteractions=FALSE,
				verbose=TRUE)
mass_ab.best.random$best.random


# model select
mass_ab.model <- model_select(all.data  = PA_11_14, 
			     responseVar = "CWM_Adult_wet_mass_log10_g", 
                       fixedFactors= fF,
                       fixedTerms= fT,
                       fixedInteractions=fI,
                       randomStruct = mass_ab.best.random$best.random,
			     otherRandoms=c("Predominant_habitat"),
                       verbose=TRUE)
mass_ab.model$final.call
validate(mass_ab.model$model)
#"CWM_Adult_wet_mass_log10_g~poly(log_elevation,3)+poly(log_slope,3)+(1+Within_PA|SS)+(1|SSB)+(1|Predominant_habitat)"

# get only data that model runs on 

mass1.data <- mass_ab.model$data

mass_ab1 <- lmer(CWM_Adult_wet_mass_log10_g ~ Within_PA + poly(log_elevation,3)+poly(log_slope,3)+
	(1+Within_PA|SS)+(1|SSB)+(1|Predominant_habitat), data = mass1.data)
mass_ab2 <- lmer(CWM_Adult_wet_mass_log10_g~ 1 + poly(log_elevation,3)+poly(log_slope,3)+
	(1+Within_PA|SS)+(1|SSB)+(1|Predominant_habitat), data = mass1.data)

anova(mass_ab1, mass_ab2)
#0.4518      1     0.5015


### try with landuse

mass_ab1 <- lmer(CWM_Adult_wet_mass_log10_g ~ Within_PA + Predominant_habitat +  Within_PA:Predominant_habitat + 
	poly(log_elevation,3)+poly(log_slope,3)+
	(1+Within_PA|SS)+(1|SSB), data = mass1.data)
mass_ab1 <- lmer(CWM_Adult_wet_mass_log10_g ~ LUPA + 
	poly(log_elevation,3)+poly(log_slope,3)+
	(1+Within_PA|SS)+(1|SSB), data = mass1.data)
mass_ab2 <- lmer(CWM_Adult_wet_mass_log10_g~  Within_PA + Predominant_habitat + 
	poly(log_elevation,3)+poly(log_slope,3)+
	(1+Within_PA|SS)+(1|SSB), data = mass1.data)

anova(mass_ab1, mass_ab2)

#for plotting
mass.data2 <- PA_11_14[,c("LUPA", "CWM_Adult_wet_mass_log10_g", "log_elevation", "log_slope", "SS", "SSB", 
		"Predominant_habitat", "Within_PA")]
mass.data2 <- na.omit(mass.data2)
mass_ab_LUPA <- lmer(CWM_Adult_wet_mass_log10_g ~ LUPA +
	poly(log_elevation,3)+poly(log_slope,3)+
	(1+Within_PA|SS)+(1|SSB), data = mass.data2)

summary(mass_ab_LUPA)

tiff("N:/Documents/PREDICTS/WDPA analysis/plots/02_15/LUPA mass birds and mammals.tif",
	width = 10, height = 10, units = "cm", pointsize = 12, res = 300)

plotLUPA(model = mass_ab_LUPA,
		responseVar = "CWM_Adult_wet_mass_log10_g", 
		data = mass.data2,
		xvar = "LUPA",
		xvar.order = c("IN Primary Vegetation",
			"OUT Mature secondary vegetation",
			"IN Mature secondary vegetation",
			"OUT Intermediate secondary vegetation",
			"IN Intermediate secondary vegetation",
			"OUT Young secondary vegetation",
			"IN Young secondary vegetation",
			"OUT Plantation forest",
			"IN Plantation forest",
			"IN Cropland",
			"OUT Cropland",
			"IN Pasture",
			"OUT Pasture",
			"OUT Urban",
			"IN Urban"),
		col.key = NULL,
		logLink = "n",
		seMultiplier=1.96,
		forPaper = T,
		ylims = c(-1.2,1.2))

dev.off()



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
#0.0131      1     0.9087



tiff("N:/Documents/PREDICTS/WDPA analysis/plots/02_15/simple models verts vol.tif",
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


### try with landuse
verts.vol.data2 <- verts.data[,c("LUPA", "vol", "log_elevation", "log_slope", "SS", "SSB", 
		"Predominant_habitat", "Within_PA")]
verts.vol.data2 <- na.omit(verts.vol.data2)

table(verts.vol.data2$Within_PA, verts.vol.data2$Predominant_habitat)
#no cropland, mature secondary, urban, young secondary
unique(verts.vol.data2$LUPA)

#full interaction drops coeffcient
#vol_1 <- lmer(vol ~ Within_PA + Predominant_habitat +  Within_PA:Predominant_habitat + 
#	poly(log_elevation,3)+poly(log_slope,3)+
#	(1+Within_PA|SS)+(1|SSB), data = verts.vol.data2)

verts.vol_1 <- lmer(vol ~ LUPA + 
	poly(log_elevation,3)+poly(log_slope,3)+
	(1+Within_PA|SS)+(1|SSB), data = verts.vol.data2)
verts.vol_2 <- lmer(vol~  Within_PA + Predominant_habitat + 
	poly(log_elevation,3)+poly(log_slope,3)+
	(1+Within_PA|SS)+(1|SSB), data = verts.vol.data2)

anova(verts.vol_1, verts.vol_2)

#for plotting
verts.vol_LUPA <- lmer(vol ~ LUPA +
	poly(log_elevation,3)+poly(log_slope,3)+
	(1+Within_PA|SS)+(1|SSB), data = verts.vol.data2)

summary(verts.vol_LUPA)

tiff("N:/Documents/PREDICTS/WDPA analysis/plots/02_15/LUPA verts vol.tif",
	width = 10, height = 10, units = "cm", pointsize = 12, res = 300)

plotLUPA(model = verts.vol_LUPA,
		responseVar = "vertebrate volume", 
		data = verts.vol.data2,
		xvar = "LUPA",
		xvar.order = c("IN Primary Vegetation",
			#"OUT Mature secondary vegetation",
			#"IN Mature secondary vegetation",
			"OUT Intermediate secondary vegetation",
			"IN Intermediate secondary vegetation",
			#"OUT Young secondary vegetation",
			#"IN Young secondary vegetation",
			"OUT Plantation forest",
			"IN Plantation forest",
			#"IN Cropland",
			#"OUT Cropland",
			"IN Pasture",
			"OUT Pasture"),
			#"OUT Urban",
			#"IN Urban")
			,
		labels = rep(c("Primary Vegetation",
			"Intermediate secondary vegetation",
			"Plantation forest",
			"Pasture"), each = 2),
		col.key = NULL,
		logLink = "n",
		seMultiplier=1.96,
		forPaper = T,
		ylims = c(-1.2,1.2))

dev.off()




