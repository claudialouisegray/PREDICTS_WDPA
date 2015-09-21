
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
multiple.taxa.matched.landuse <- subset(matched.landuse, multiple_taxa == "yes")
nrow(matched.landuse) #5015


### Size and Age analysis

xyplot(Richness_rarefied ~ taxon_of_interest|AREA_DoP, multiple.taxa.matched.landuse)

fF <- c("Zone", "taxon_of_interest", "AREA_DoP") 
fT <- list("ag_suit" = "3", "log_slope" = "3", "log_elevation" = "3")
keepVars <- character(0)
fI <- c("AREA_DoP:taxon_of_interest", "AREA_DoP:Zone")
RS <-  character(0)

Richness_rarefied.best.random <- compare_randoms(multiple.taxa.matched.landuse, "Richness_rarefied",
				fitFamily = "poisson",
				siteRandom = TRUE,
				fixedFactors=fF,
        fixedTerms=fT,
        fixedInteractions=fI,
        otherRandoms=c("Predominant_habitat"),
				fixed_RandomSlopes = RS,
        fitInteractions=FALSE,
				verbose=TRUE)

Richness_rarefied.best.random$best.random #"(1|SS)+ (1|SSBS)+ (1|SSB)+(1|Predominant_habitat)"

# model select
Richness_rarefied.model <- model_select(all.data  = multiple.taxa.matched.landuse, 
			 responseVar = "Richness_rarefied",
			 fitFamily = "poisson",
			 siteRandom = TRUE, 
			 alpha = 0.05,
       fixedFactors= fF,
       fixedTerms= fT,
       keepVars = keepVars,
       fixedInteractions=fI,
       randomStruct = Richness_rarefied.best.random$best.random,
       otherRandoms=c("Predominant_habitat"),
       
       verbose=TRUE)
Richness_rarefied.model$final.call
Richness_rarefied.model$stats
#Richness_rarefied~Zone+poly(ag_suit,3)+poly(log_elevation,3)+(1|SS)+(1|SSBS)+(1|SSB)+(1|Predominant_habitat)"




