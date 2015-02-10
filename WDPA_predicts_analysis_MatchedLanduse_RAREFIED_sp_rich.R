

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









######
Richness_rarefied
######

title(main = "all matched data")
hist(multiple.taxa.matched.landuse.s$CWM_Geographic_Richness_rarefied_log10_square_km)

plot(log(hpd+1) ~ ag_suit, multiple.taxa.matched.landuse)


plot(Richness_rarefied ~ Within_PA, multiple.taxa.matched.landuse)
plot(Richness_rarefied ~ Zone, multiple.taxa.matched.landuse)
plot(Richness_rarefied ~ taxon_of_interest, multiple.taxa.matched.landuse)
plot(Richness_rarefied ~ elevation, multiple.taxa.matched.landuse)
plot(Richness_rarefied ~ log_dist_in, multiple.taxa.matched.landuse)
plot(Richness_rarefied ~ log_dist_out, multiple.taxa.matched.landuse)
plot(Richness_rarefied ~ log_bound_dist_km_PA_neg, multiple.taxa.matched.landuse)
plot(Richness_rarefied ~ bound_dist_km_PA_neg, multiple.taxa.matched.landuse)


plot(Richness_rarefied ~ DoP.PA, multiple.taxa.matched.landuse)


xyplot(Richness_rarefied ~ Within_PA|Zone, multiple.taxa.matched.landuse)
xyplot(Richness_rarefied ~ Within_PA|Realm, multiple.taxa.matched.landuse)
xyplot(Richness_rarefied ~ Within_PA|taxon_of_interest, multiple.taxa.matched.landuse)

xyplot(Richness_rarefied ~ log_bound_dist_km_PA_neg|Zone, multiple.taxa.matched.landuse)
xyplot(Richness_rarefied ~ log_bound_dist_km_PA_neg|taxon_of_interest, multiple.taxa.matched.landuse)
xyplot(Richness_rarefied ~ log_bound_dist_km_PA_neg|ag_suit, multiple.taxa.matched.landuse) # could be tricky to fit, few points in ag_suit = 1
plot3d( multiple.taxa.matched.landuse$log_bound_dist_km_PA_neg, multiple.taxa.matched.landuse$log_access, multiple.taxa.matched.landuse$Richness_rarefied) #reasonable
plot3d( multiple.taxa.matched.landuse$log_bound_dist_km_PA_neg, multiple.taxa.matched.landuse$log_hpd, multiple.taxa.matched.landuse$Richness_rarefied)	#reasonable

xyplot(Richness_rarefied ~ hpd|Within_PA, multiple.taxa.matched.landuse)


xyplot(Richness_rarefied ~ elevation|Zone, multiple.taxa.matched.landuse)
xyplot(Richness_rarefied ~ access|Zone, multiple.taxa.matched.landuse)
xyplot(Richness_rarefied ~ ag_suit|Zone, multiple.taxa.matched.landuse)
xyplot(Richness_rarefied ~ hpd|Zone, multiple.taxa.matched.landuse)
xyplot(Richness_rarefied ~ slope|Zone, multiple.taxa.matched.landuse)




hist(multiple.taxa.matched.landuse$elevation)
hist(multiple.taxa.matched.landuse$ag_suit)
hist(multiple.taxa.matched.landuse$hpd)
hist(multiple.taxa.matched.landuse$DoP)
hist(log(multiple.taxa.matched.landuse$elevation+1))
hist(log(multiple.taxa.matched.landuse$hpd+1))
hist(log(PA.data$DoP+1))

plot(Richness_rarefied ~ elevation, multiple.taxa.matched.landuse.s)
plot(Richness_rarefied ~ log(elevation+1), multiple.taxa.matched.landuse.s)







### distance to boundary analysis

fF <- c("Zone", "taxon_of_interest") 
fT <- list("ag_suit" = "3", "log_slope" = "3", "log_elevation" = "3", "log_bound_dist_km_PA_neg" = "3", "DoP.PA" = "3", "log_AREA.PA" = "3")
keepVars <- character(0)
fI <- character(0)
RS <-  c("log_bound_dist_km_PA_neg")
# without block "Richness_rarefied~poly(ag_suit,3)+poly(log_AREA.PA,2)+poly(log_elevation,3)+Zone+(1+log_bound_dist_km_PA_neg|SS)+(1|SSBS)+(1|Predominant_habitat)"
# with block: "Richness_rarefied~poly(ag_suit,3)+poly(log_elevation,3)+Zone+(1+log_bound_dist_km_PA_neg|SS)+(1|SSBS)+(1|SSB)+(1|Predominant_habitat)
# 8 landuses: "Richness_rarefied~poly(ag_suit,3)+poly(log_AREA.PA,2)+poly(log_elevation,3)+Zone+(1|SS)+(1|SSBS)+(1|SSB)+(1|Predominant_habitat)"


# add interactions

fF <- c("Zone", "taxon_of_interest") 
fT <- list("log_bound_dist_km_PA_neg" = "1", "DoP.PA" = "1", "log_AREA.PA" = "2")
keepVars <- list("ag_suit" = "3", "log_slope" = "1", "log_elevation" = "3")
fI <- c("Zone:poly(log_bound_dist_km_PA_neg,1)", "taxon_of_interest:poly(log_bound_dist_km_PA_neg,1)")
RS <-  c("log_bound_dist_km_PA_neg")


#without block - but elevation is wrong here, should have been cubic. Richness_rarefied~poly(ag_suit,3)+poly(log_AREA.PA,2)+poly(log_elevation,1)+Zone
#+(1+log_bound_dist_km_PA_neg|SS)+(1|SSBS)+(1|Predominant_habitat)"
# with block: "Richness_rarefied~Zone+poly(ag_suit,3)+poly(log_slope,1)+poly(log_elevation,3)+(1+log_bound_dist_km_PA_neg|SS)+(1|SSBS)+(1|SSB)+(1|Predominant_habitat)"
# 8 landuses: "Richness_rarefied~Zone+poly(ag_suit,3)+poly(log_slope,1)+poly(log_elevation,3)+(1|SS)+(1|SSBS)+(1|SSB)+(1|Predominant_habitat)

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
# use this - keep to PREDICTS approach

best.random <- "(1|SS)+ (1|SSBS) + (1|SSB) + (1|Predominant_habitat)"
 


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
                       randomStruct = best.random,
			     otherRandoms=c("Predominant_habitat"),
                       verbose=TRUE)


validate(Richness_rarefied.model$model) #ok
write.csv(Richness_rarefied.model$stats, 
	"N:/Documents/PREDICTS/WDPA analysis/stats tables all/8 landuses/Richness_rarefied.model.stats.23.01.2015.csv")






### within PA analysis


fF <- c("Zone", "taxon_of_interest", "Within_PA") 
fT <- list("ag_suit" = "3", "log_slope" = "3", "log_elevation" = "3", "DoP.PA" = "3", "log_AREA.PA" = "3")
keepVars <- character(0)
fI <- character(0)
RS <-  c("Within_PA")

#without block: Richness_rarefied~poly(ag_suit,3)+poly(log_AREA.PA,2)+poly(log_elevation,3)+Zone+(1+Within_PA|SS)+(1|SSBS)+(1|Predominant_habitat)"
# with block: Richness_rarefied~poly(ag_suit,3)+poly(log_AREA.PA,2)+poly(log_elevation,3)+Zone+(1+Within_PA|SS)+(1|SSBS)+(1|SSB)+(1|Predominant_habitat)"
# 8 landuses: Richness_rarefied~poly(ag_suit,3)+poly(log_AREA.PA,2)+poly(log_elevation,3)+Zone+(1+Within_PA|SS)+(1|SSBS)+(1|SSB)+(1|Predominant_habitat)


# add interactions
fF <- c("Zone", "taxon_of_interest", "Within_PA") 
fT <- list("DoP.PA" = "1", "log_AREA.PA" = "2")
keepVars <- list("ag_suit" = "3", "log_slope" = "1", "log_elevation" = "3")
fI <- c("Within_PA:poly(ag_suit,3)", "Within_PA:poly(log_slope,1)", "Within_PA:poly(log_elevation,3)",
	"Within_PA:taxon_of_interest", 
	"taxon_of_interest:poly(DoP.PA,1)", "taxon_of_interest:poly(log_AREA.PA,2)",
	"Within_PA:Zone",
	"Zone:poly(DoP.PA,1)", "Zone:poly(log_AREA.PA,2)")
RS <-  c("Within_PA")

#without block: Richness_rarefied~poly(ag_suit,3)+poly(log_AREA.PA,2)+poly(log_elevation,3)+Zone
#+Within_PA:poly(ag_suit,3)
#+Within_PA:taxon_of_interest
#+Within_PA+taxon_of_interest
#+(1+Within_PA|SS)+(1|SSBS)+(1|Predominant_habitat)

# with block: Richness_rarefied~poly(ag_suit,3)+poly(log_AREA.PA,2)+poly(log_elevation,3)+Zone
#+(1+Within_PA|SS)+(1|SSBS)+(1|SSB)+(1|Predominant_habitat)

# 8 landuses: Richness_rarefied~poly(log_AREA.PA,2)+Zone+poly(ag_suit,3)+poly(log_slope,1)+poly(log_elevation,3)+(1+Within_PA|SS)+(1|SSBS)+(1|SSB)+(1|Predominant_habitat)

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


Richness_rarefied.best.random2$best.random #"(1+Within_PA|SS)+ (1|SSBS)+ (1|SSB)+(1|Predominant_habitat)"


best.random <- "(1+Within_PA|SS)+ (1|SSBS)+ (1|SSB) + (1|Predominant_habitat)"
 


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
                       randomStruct = best.random,
			     otherRandoms=c("Predominant_habitat"),
                       verbose=TRUE)


validate(Richness_rarefied.model2$model) #ok
write.csv(Richness_rarefied.model2$stats, "Richness_rarefied.model2.stats.23.01.2015.csv")

Richness_rarefied.model$warnings
Richness_rarefied.model$final.call
Richness_rarefied.model$stats




summary(Richness_rarefied.model.int$model)
XV <- CrossValidate(Richness_rarefied.model$model,50, divFactor = "SS", fitFamily = "gaussian", data= Richness_rarefied.model$data, 
	fixedFactors = fF, fixedTerms = fT, fixedInteractions = fI, randomStruct = Richness_rarefied.best.random$best.random)
# all significant estimates seem acceptably similar


XVr <- CrossValidate(Richness_rarefied.model$model,-1, divFactor = "Predominant_habitat", fitFamily = "gaussian", data= Richness_rarefied.model$data, 
	fixedFactors = fF, fixedTerms = fT, fixedInteractions = fI, randomStruct = Richness_rarefied.best.random$best.random)
# but if you dropped a whole land use it wouldnt be that surprising if the results changed? 

#glmmR.wikidot/faq










