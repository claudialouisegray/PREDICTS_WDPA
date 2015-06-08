
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

setwd("R:/ecocon_d/clg32/GitHub/PREDICTS_WDPA")

# load functions
source("compare_randoms.R")
source("model_select.R")
source("plotFactorInteraction.R")

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


### Size and Age analysis

xyplot(Species_richness ~ taxon_of_interest|AREA_DoP, multiple.taxa.matched.landuse)

fF <- c("Zone", "taxon_of_interest", "AREA_DoP") 
fT <- list("ag_suit" = "3", "log_slope" = "3", "log_elevation" = "3")
keepVars <- character(0)
fI <- c("AREA_DoP:taxon_of_interest", "AREA_DoP:Zone")
RS <-  character(0)
Species_richness.best.random <- compare_randoms(multiple.taxa.matched.landuse, "Species_richness",
				fitFamily = "poisson",
				siteRandom = TRUE,
				fixedFactors=fF,
                        fixedTerms=fT,
			     	keepVars = keepVars,
                       	fixedInteractions=fI,
                        otherRandoms=c("Predominant_habitat"),
				fixed_RandomSlopes = RS,
                        fitInteractions=FALSE,
				verbose=TRUE)

Species_richness.best.random$best.random #"(1|SS)+ (1|SSBS)+ (1|SSB)+(1|Predominant_habitat)"

Species_richness.model <- model_select(all.data  = multiple.taxa.matched.landuse, 
			     responseVar = "Species_richness",
			     fitFamily = "poisson",
			     siteRandom = TRUE, 
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
			     keepVars = keepVars,
                       fixedInteractions=fI,
                       randomStruct =Species_richness.best.random$best.random ,
			     otherRandoms=c("Predominant_habitat"),
                       verbose=TRUE)
Species_richness.model$final.call
#Species_richness~AREA_DoP+taxon_of_interest+Zone+poly(ag_suit,3)+poly(log_slope,1)+poly(log_elevation,2)+(1|SS)+(1|SSBS)+(1|SSB)+(1|Predominant_habitat)"

summary(Species_richness.model$model)
Species_richness.model$stats


### plot

tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/02_15/sp rich vs size and age.tif",
	width = 9, height = 15, units = "cm", pointsize = 12, res = 300)

plotFactorInteraction(model = Species_richness.model$model,
responseVar = "Species_richness",
data = Species_richness.model$data,
xvar = "AREA_DoP",
xvar.order = c("small_young", "small_old", "large_young", "large_old"), #this must be the levels of the factor in the order to be plotted, not including reference level
xvar.labels = c("Small, Young", "Small, Old", "Large, Young", "Large, Old"),
logLink = "e",
xlab = "PA size and age class",
ylab = "Relative species richness %")

#text(1,26, "Young = 0 - 20 yrs \nOld = 20 - 85 yrs \nSmall = 0 - 400 km2 \nLarge = 400 - 12000km2", 
#	adj = 0, cex = 0.8)

dev.off()


# test interaction

m1 <- glmer(Species_richness  ~ taxon_of_interest + Zone + AREA_DoP + 
		(1|SS)+ (1|SSBS)+ (1|SSB)+(1|Predominant_habitat), multiple.taxa.matched.landuse,
		family = "poisson")
m2 <- glmer(Species_richness  ~ taxon_of_interest + Zone + AREA.PA.f + 
		(1|SS)+ (1|SSBS)+ (1|SSB)+(1|Predominant_habitat), multiple.taxa.matched.landuse,
		family = "poisson")
m3 <- glmer(Species_richness  ~ taxon_of_interest + Zone + DoP.PA.f +
		(1|SS)+ (1|SSBS)+ (1|SSB)+(1|Predominant_habitat), multiple.taxa.matched.landuse,
		family = "poisson")
m4 <- glmer(Species_richness  ~ taxon_of_interest + Zone + 
		(1|SS)+ (1|SSBS)+ (1|SSB)+(1|Predominant_habitat), multiple.taxa.matched.landuse,
		family = "poisson")
anova(m1, m2)#interaction term better than model with just area
anova(m1, m3)#not better than model with just age
anova(m2, m4) #no effect of area alone
anova(m3, m4) #effect of age alone

