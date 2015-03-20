

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

xyplot(Richness_rarefied ~ taxon_of_interest|AREA_DoP, multiple.taxa.matched.landuse)

#run with no interactions first to see which variables have nonlinear relationships
fF <- c("Zone", "taxon_of_interest", "AREA_DoP") 
fT <- list("ag_suit" = "3", "log_slope" = "3", "log_elevation" = "3")
keepVars <- character(0)
fI <- character(0)
RS <-  character(0)
Richness_rarefied.best.random2 <- compare_randoms(multiple.taxa.matched.landuse, "Richness_rarefied",
				fitFamily = "poisson",
				siteRandom = TRUE,
				fixedFactors=fF,
                        fixedTerms=fT,
                       	fixedInteractions=fI,
                        otherRandoms=c("Predominant_habitat"),
				fixed_RandomSlopes = RS,
                        fitInteractions=FALSE,
				verbose=TRUE)


Richness_rarefied.best.random2$best.random #"(1|SS)+ (1|SSBS)+ (1|SSB)+(1|Predominant_habitat)"

Richness_rarefied.poly <- model_select(all.data  = multiple.taxa.matched.landuse, 
			     responseVar = "Richness_rarefied",
			     fitFamily = "poisson",
			     siteRandom = TRUE, 
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
			     keepVars = keepVars,
                       fixedInteractions=fI,
                       randomStruct = Richness_rarefied.best.random2$best.random,
			     otherRandoms=c("Predominant_habitat"),
                       verbose=TRUE)
Richness_rarefied.poly$final.call
#Richness_rarefied~poly(ag_suit,3)+poly(log_elevation,3)+Zone+(1|SS)+(1|SSBS)+(1|SSB)+(1|Predominant_habitat)"



# add interactions
fF <- c("Zone", "taxon_of_interest", "AREA_DoP") 
fT <- list()
keepVars <- list("ag_suit" = "3",  "log_elevation" = "3")
fI <- c("AREA_DoP:taxon_of_interest", "AREA_DoP:Zone")
RS <-  character(0)

# model select
Richness_rarefied.model2 <- model_select(all.data  = multiple.taxa.matched.landuse, 
			     responseVar = "Richness_rarefied",
			     fitFamily = "poisson",
			     siteRandom = TRUE, 
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
			     keepVars = keepVars,
                       fixedInteractions=fI,
                       randomStruct = Richness_rarefied.best.random2$best.random,
			     otherRandoms=c("Predominant_habitat"),
                       verbose=TRUE)
Richness_rarefied.model2$final.call
#Richness_rarefied~Zone+poly(ag_suit,3)+poly(log_elevation,3)+(1|SS)+(1|SSBS)+(1|SSB)+(1|Predominant_habitat)"

write.csv(Richness_rarefied.model2$stats, "N:/Documents/PREDICTS/WDPA analysis/stats tables all/age size bins/rarrich.model.age.size.stats.20.03.2014.csv")

save.image("N:\\Documents\\PREDICTS\\WDPA analysis\\RData files\\8 landuses\\age_size rar rich.RData")
