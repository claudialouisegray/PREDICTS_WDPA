


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

source("compare_randoms_lmer - with poly.R")
source("model_select.R")


#load data
source("WDPA_predicts_prep_PA_11_14_for_analysis.R")

validate <- function(x) {
  par(mfrow = c(1,2))
  plot(resid(x)~ fitted(x))
  hist(resid(x))
  par(mfrow = c(1,1))
}



construct_call<-function(responseVar,fixedStruct,randomStruct){
  return(paste(responseVar,"~",fixedStruct,"+",randomStruct,sep=""))
}





names(PA_11_14)

m <- lmer(log_abundance ~ Within_PA + Predominant_habitat +
	Within_PA:Predominant_habitat + log_slope + log_elevation + ag_suit
	+ Zone + taxon_of_interest + 
	(1|SS) + (1|SSB), PA_11_14)

m1 <- lmer(log_abundance ~ Within_PA + Predominant_habitat +
	 log_slope + log_elevation + ag_suit
	+ Zone + taxon_of_interest + 
	(1|SS) + (1|SSB), PA_11_14)

anova(m1, m)


m <- lmer(log_abundance ~ Within_PA + Predominant_habitat +
	Within_PA:Predominant_habitat + log_slope + log_elevation + ag_suit
	+ Zone + taxon_of_interest + 
	(Within_PA|SS) + (1|SSB), matched.landuse)

m1 <- lmer(log_abundance ~ Within_PA + Predominant_habitat +
	 log_slope + log_elevation + ag_suit
	+ Zone + taxon_of_interest + 
	(Within_PA|SS) + (1|SSB), matched.landuse)

anova(m1, m)



m <- lmer(Species_richness ~ Within_PA + Predominant_habitat +
	Within_PA:Predominant_habitat + log_slope + log_elevation + ag_suit
	+ Zone + taxon_of_interest + 
	(Within_PA|SS) + (1|SSB) + (1|SSBS), family = "poisson", 
	matched.landuse)
m <- lmer(Species_richness ~ Within_PA + Predominant_habitat +
	Within_PA:Predominant_habitat + log_slope + log_elevation + ag_suit
	+ Zone + taxon_of_interest + 
	(Within_PA|SS) + (1|SSB) + (1|SSBS), family = "poisson", 
	matched.landuse)


m <- glmer(Species_richness ~ Within_PA + Predominant_habitat +
	Within_PA:Predominant_habitat + log_slope + log_elevation + ag_suit
	+ Zone + taxon_of_interest + 
	(Within_PA|SS) + (1|SSB) + (1|SSBS), family = "poisson", 
	data = PA_11_14)