
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
#setwd("N:/Documents/PREDICTS/WDPA analysis")

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



### age/size analysis

xyplot(range ~ taxon_of_interest|AREA_DoP, matched.landuse)

fF <- c("Zone", "taxon_of_interest", "AREA_DoP") 
fT <- list("ag_suit" = "3", "log_slope" = "3", "log_elevation" = "3")
keepVars <- character(0)
fI <- c("AREA_DoP:taxon_of_interest", "AREA_DoP:Zone")
RS <-  character(0)

range.best.random <- compare_randoms(matched.landuse, "range",
				fixedFactors=fF,
                         fixedTerms=fT,
				keepVars = keepVars,
                       fixedInteractions=fI,
                         otherRandoms=c("Predominant_habitat"),
				fixed_RandomSlopes = RS,
                          fitInteractions=FALSE,
				verbose=TRUE)

range.best.random$best.random #"(1+Within_PA|SS)+ (1|SSB)+(1|Predominant_habitat)"

# model select
range.model <- model_select(all.data  = matched.landuse, 
			     responseVar = "range", 
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
			     keepVars = keepVars,
                       fixedInteractions=fI,
                       randomStruct = range.best.random$best.random,
			     otherRandoms=c("Predominant_habitat"),
                       verbose=TRUE)
#range~taxon_of_interest+Zone+AREA_DoP:taxon_of_interest+AREA_DoP+poly(ag_suit,3)+poly(log_slope,3)+poly(log_elevation,2)+(1|SS)+(1|SSB)+(1|Predominant_habitat)"
range.model$stats

summary(range.model2$model)

XV <- CrossValidate(range.model2$model,10, divFactor = "SS", fitFamily = "gaussian", data= range.model2$data, 
	fixedFactors = fF, fixedTerms = fT, fixedInteractions = fI, randomStruct = "(1+Within_PA|SS)+ (1|SSB) + (1|Predominant_habitat)")


tiff( "R:/ecocon_d/clg32/PREDICTS/WDPA analysis/plots/02_15/range vs size and age vs taxon.tif",
	width = 20, height = 15, units = "cm", pointsize = 12, res = 300)

plotFactorInteraction(model = range.model$model,
responseVar = "Endemicity",
data = range.model$data,
xvar = "AREA_DoP",
xvar.order = c("small_young", "small_old", "large_young", "large_old"),  #this must be the levels of the factor in the order to be plotted, not including reference level
xvar.labels = c("Small, Young", "Small, Old", "Large, Young", "Large, Old"),
intvar = "taxon_of_interest",
intvar.order = c("Invertebrates", "Vertebrates"),
logLink = "n",
xlab = "PA size and age class",
ylab = "Relative endemicity %")

#text(0.7,8, "Young = 0 - 20 yrs \nOld = 20 - 85 yrs \nSmall = 0 - 400 km2 \nLarge = 400 - 12000km2", 
#	adj = 0, cex = 0.8)


dev.off()

 save.image("R:/ecocon_d/clg32/PREDICTS/WDPA analysis/RData files\\8 landuses\\age_size range.RData")


