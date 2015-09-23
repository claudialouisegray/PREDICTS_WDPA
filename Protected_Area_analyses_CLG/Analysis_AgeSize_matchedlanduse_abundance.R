

rm(list=ls()) 

library(yarg)
library(roquefort)
library(gamm4)


setwd("R:/ecocon_d/clg32/GitHub/PREDICTS_WDPA")
source("compare_randoms.R")
source("model_select.R")
setwd("R:/ecocon_d/clg32/PREDICTS/WDPA analysis")
PREDICTS_WDPA <- read.csv("PREDICTS_WDPA.csv")

validate <- function(x) {
  par(mfrow = c(1,2))
  plot(resid(x)~ fitted(x))
  hist(resid(x))
  par(mfrow = c(1,1))
}


matched.landuse <- subset(PREDICTS_WDPA, matched.landuse == "yes")
nrow(matched.landuse) #5015


### age/size analysis

xyplot(log_abundance ~ taxon_of_interest|AREA_DoP, matched.landuse)

fF <- c("Zone", "taxon_of_interest", "AREA_DoP") 
fT <- list("ag_suit" = "3", "log_slope" = "3", "log_elevation" = "3")
keepVars <- character(0)
fI <- c("AREA_DoP:taxon_of_interest", "AREA_DoP:Zone")
RS <-  character(0)

log_abundance.best.random <- compare_randoms(matched.landuse, 
				responseVar = "log_abundance",
				keepVars = keepVars,
				fixedFactors=fF,
        fixedTerms=fT,
        fixedInteractions=fI,
        otherRandoms=c("Predominant_habitat"),
				fixed_RandomSlopes = RS,
        fitInteractions=FALSE,
				verbose=TRUE)

log_abundance.best.random$best.random

# model select
log_abundance.model<- model_select(all.data  = matched.landuse, 
        responseVar = "log_abundance", 
        keepVars = keepVars,
        alpha = 0.05,
        fixedFactors= fF,
        fixedTerms= fT,
        fixedInteractions=fI,
        randomStruct = log_abundance.best.random$best.random,
        otherRandoms=c("Predominant_habitat"),
        verbose=TRUE)

log_abundance.model$stats
log_abundance.model$final.call
#"log_abundance~poly(ag_suit,1)+taxon_of_interest+AREA_DoP:taxon_of_interest+AREA_DoP+(1|SS)+(1|SSB)+(1|Predominant_habitat)"


validate(log_abundance.model$model) #ok
summary(log_abundance.model$model)
log_abundance.model$warnings

summary(log_abundance.model2$model)
fT <- list("DoP.PA" = "1", "log_AREA.PA" = "3", "ag_suit" = "1", log_elevation = "1")

XV <- CrossValidate(log_abundance.model2$model,10, divFactor = "SS", fitFamily = "gaussian", data= log_abundance.model2$data, 
	fixedFactors = fF, fixedTerms = fT, fixedInteractions = fI, randomStruct = "(1+Within_PA|SS)+ (1|SSB) + (1|Predominant_habitat)")




tiff( "abundance vs size and age vs taxon.tif",
	width = 20, height = 15, units = "cm", pointsize = 12, res = 300)

plotFactorInteraction(model = log_abundance.model$model,
responseVar = "log_abundance",
data = log_abundance.model$data,
xvar = "AREA_DoP",
xvar.order = c("small_young", "small_old", "large_young", "large_old"), 
#this must be the levels of the factor in the order to be plotted, not including reference level
intvar = "taxon_of_interest",
intvar.order = c("Invertebrates", "Vertebrates"),
logLink = "e",
xlab = "PA size and age class",
ylab = "Relative abundance %")

#text(0.7,90, "Young = 0 - 20 yrs \nOld = 20 - 85 yrs \nSmall = 0 - 400 km2 \nLarge = 400 - 12000km2", 
#	adj = 0, cex = 0.8)

dev.off()

