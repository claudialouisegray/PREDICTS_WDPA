

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
source("WDPA_predicts_prep_matched.landuse_for_analysis.R")

validate <- function(x) {
  par(mfrow = c(1,2))
  plot(resid(x)~ fitted(x))
  hist(resid(x))
  par(mfrow = c(1,1))
}



construct_call<-function(responseVar,fixedStruct,randomStruct){
  return(paste(responseVar,"~",fixedStruct,"+",randomStruct,sep=""))
}



nrow(matched.landuse_VU_EN_CR)






######
Species_richness
######

title(main = "all matched data")
hist(multiple.taxa.matched.landuse_VU_EN_CR.s$CWM_Geographic_Species_richness_log10_square_km)

plot(log(hpd+1) ~ ag_suit, multiple.taxa.matched.landuse_VU_EN_CR)


plot(Species_richness ~ Within_PA, multiple.taxa.matched.landuse_VU_EN_CR)
plot(Species_richness ~ Zone, multiple.taxa.matched.landuse_VU_EN_CR)
plot(Species_richness ~ taxon_of_interest, multiple.taxa.matched.landuse_VU_EN_CR)
plot(Species_richness ~ elevation, multiple.taxa.matched.landuse_VU_EN_CR)
plot(Species_richness ~ log_bound_dist_km_PA_neg, multiple.taxa.matched.landuse_VU_EN_CR)
plot(Species_richness ~ bound_dist_km_PA_neg, multiple.taxa.matched.landuse_VU_EN_CR)


plot(Species_richness ~ DoP.PA, multiple.taxa.matched.landuse_VU_EN_CR)


xyplot(Species_richness ~ Within_PA|Zone, multiple.taxa.matched.landuse_VU_EN_CR)
xyplot(Species_richness ~ Within_PA|Realm, multiple.taxa.matched.landuse_VU_EN_CR)
xyplot(Species_richness ~ Within_PA|taxon_of_interest, multiple.taxa.matched.landuse_VU_EN_CR)

xyplot(Species_richness ~ log_bound_dist_km_PA_neg|Zone, multiple.taxa.matched.landuse_VU_EN_CR)
xyplot(Species_richness ~ log_bound_dist_km_PA_neg|taxon_of_interest, multiple.taxa.matched.landuse_VU_EN_CR)
xyplot(Species_richness ~ log_bound_dist_km_PA_neg|ag_suit, multiple.taxa.matched.landuse_VU_EN_CR) # could be tricky to fit, few points in ag_suit = 1
plot3d( multiple.taxa.matched.landuse_VU_EN_CR$log_bound_dist_km_PA_neg, multiple.taxa.matched.landuse_VU_EN_CR$log_access, multiple.taxa.matched.landuse_VU_EN_CR$Species_richness) #reasonable
plot3d( multiple.taxa.matched.landuse_VU_EN_CR$log_bound_dist_km_PA_neg, multiple.taxa.matched.landuse_VU_EN_CR$log_hpd, multiple.taxa.matched.landuse_VU_EN_CR$Species_richness)	#reasonable

xyplot(Species_richness ~ hpd|Within_PA, multiple.taxa.matched.landuse_VU_EN_CR)


xyplot(Species_richness ~ elevation|Zone, multiple.taxa.matched.landuse_VU_EN_CR)
xyplot(Species_richness ~ access|Zone, multiple.taxa.matched.landuse_VU_EN_CR)
xyplot(Species_richness ~ ag_suit|Zone, multiple.taxa.matched.landuse_VU_EN_CR)
xyplot(Species_richness ~ hpd|Zone, multiple.taxa.matched.landuse_VU_EN_CR)
xyplot(Species_richness ~ slope|Zone, multiple.taxa.matched.landuse_VU_EN_CR)



# get only data that model runs on 

data <- multiple.taxa.matched.landuse_VU_EN_CR[,c("Species_richness", "Zone", "taxon_of_interest", "ag_suit", "log_access",
	 "log_hpd", "log_bound_dist_km", "log_bound_dist_km_PA_neg", "Within_PA", "DoP.PA", "IUCN.PA", "log_AREA.PA", "Predominant_habitat", "SS", "SSBS")]

data <- na.omit(data) #have to get rid of NA to try poly on ag_suit







### test with cubic to see what works


fF <- c("Zone", "Within_PA") 
fT <- list("ag_suit" = "3", "log_hpd" = "3", "log_access" = "3", "log_bound_dist_km_PA_neg" = "3",  "log_AREA.PA" = "3")
fI <- character(0)
RS <-  c("Within_PA")

# 07_14 Species_richness~poly(ag_suit,2)+poly(log_hpd,3)+(1+Within_PA|SS)+(1|Predominant_habitat)+(1|SSBS)
# 11_14 Species_richness~poly(ag_suit,2)+poly(log_access,1)+poly(log_hpd,3)+(1+Within_PA|SS)+(1|Predominant_habitat)+(1|SSBS)"




# add interactions
#doesnt work
fF <- c("Zone", "Within_PA") 
fT <- list("ag_suit" = "2", "log_hpd" = "3", "log_access" = "1", "log_bound_dist_km_PA_neg" = "1", "DoP.PA" = "1", "log_AREA.PA" = "1")
fI <- c("Within_PA:poly(ag_suit,2)", "Within_PA:poly(log_hpd,3)", "Within_PA:poly(log_access,1)",
	"Within_PA:Zone",
	"Zone:poly(log_bound_dist_km_PA_neg,1)","Zone:poly(DoP.PA,1)", "Zone:poly(log_AREA.PA,1)")
RS <-  c("Within_PA")



# works ok no interactions significant
fF <- c("Zone", "Within_PA") 
fT <- list("ag_suit" = "2", "log_hpd" = "3", "log_access" = "1", "log_bound_dist_km_PA_neg" = "1", "log_AREA.PA" = "1")
fI <- c("Within_PA:poly(ag_suit,2)", "Within_PA:poly(log_hpd,3)", "Within_PA:poly(log_access,1)",
	"Within_PA:Zone",
	"Zone:poly(log_bound_dist_km_PA_neg,1)", "Zone:poly(log_AREA.PA,1)")
RS <-  c("Within_PA")
#11_14 Species_richness~poly(ag_suit,2)+poly(log_access,1)+poly(log_hpd,3)+(1+Within_PA|SS)+(1|Predominant_habitat)+(1|SSBS)"
 


Species_richness.best.random <- compare_randoms(multiple.taxa.matched.landuse_VU_EN_CR, "Species_richness",
				fitFamily = "poisson",
				siteRandom = TRUE,
				fixedFactors=fF,
                         fixedTerms=fT,
                       fixedInteractions=fI,
                         otherRandoms=c("Predominant_habitat"),
				fixed_RandomSlopes = RS,
                          fitInteractions=FALSE,
				verbose=TRUE)


Species_richness.best.random$best.random
 


# model select
#Species_richness.model.int now has the model without interactions
Species_richness.model <- model_select(all.data  = multiple.taxa.matched.landuse_VU_EN_CR, 
			     responseVar = "Species_richness",
			     fitFamily = "poisson",
			     siteRandom = TRUE, 
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
                       fixedInteractions=fI,
                       randomStruct = "(1+Within_PA|SS)+(1|Predominant_habitat) + (1|SSBS)",
			     otherRandoms=c("Predominant_habitat"),
                       verbose=TRUE)


validate(Species_richness.model$model) #ok
Species_richness.model$warnings
Species_richness.model$final.call
Species_richness.model$stats


summary(Species_richness.model$model)
XV <- CrossValidate(Species_richness.model$model,50, divFactor = "SS", fitFamily = "gaussian", data= Species_richness.model$data, 
	fixedFactors = fF, fixedTerms = fT, fixedInteractions = fI, randomStruct = Species_richness.best.random$best.random)
# all significant estimates seem acceptably similar


XVr <- CrossValidate(Species_richness.model$model,-1, divFactor = "Predominant_habitat", fitFamily = "gaussian", data= Species_richness.model$data, 
	fixedFactors = fF, fixedTerms = fT, fixedInteractions = fI, randomStruct = Species_richness.best.random$best.random)
# but if you dropped a whole land use it wouldnt be that surprising if the results changed? 

#glmmR.wikidot/faq

#




write.csv(Species_richness.model$stats, "Species_richness.model.stats.15.10.2014.csv")








