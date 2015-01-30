

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
Species_richness
######

title(main = "all matched data")
hist(multiple.taxa.matched.landuse.s$CWM_Geographic_Species_richness_log10_square_km)

plot(log(hpd+1) ~ ag_suit, multiple.taxa.matched.landuse)


plot(Species_richness ~ Within_PA, multiple.taxa.matched.landuse)
plot(Species_richness ~ Zone, multiple.taxa.matched.landuse)
plot(Species_richness ~ taxon_of_interest, multiple.taxa.matched.landuse)
plot(Species_richness ~ elevation, multiple.taxa.matched.landuse)
plot(Species_richness ~ log_bound_dist_km_PA_neg, multiple.taxa.matched.landuse)
plot(Species_richness ~ bound_dist_km_PA_neg, multiple.taxa.matched.landuse)


plot(Species_richness ~ DoP.PA, multiple.taxa.matched.landuse)


xyplot(Species_richness ~ Within_PA|Zone, multiple.taxa.matched.landuse)
xyplot(Species_richness ~ Within_PA|Realm, multiple.taxa.matched.landuse)
xyplot(Species_richness ~ Within_PA|taxon_of_interest, multiple.taxa.matched.landuse)

xyplot(Species_richness ~ log_bound_dist_km_PA_neg|Zone, multiple.taxa.matched.landuse)
xyplot(Species_richness ~ log_bound_dist_km_PA_neg|taxon_of_interest, multiple.taxa.matched.landuse)
xyplot(Species_richness ~ log_bound_dist_km_PA_neg|ag_suit, multiple.taxa.matched.landuse) # could be tricky to fit, few points in ag_suit = 1
plot3d( multiple.taxa.matched.landuse$log_bound_dist_km_PA_neg, multiple.taxa.matched.landuse$log_access, multiple.taxa.matched.landuse$Species_richness) #reasonable
plot3d( multiple.taxa.matched.landuse$log_bound_dist_km_PA_neg, multiple.taxa.matched.landuse$log_hpd, multiple.taxa.matched.landuse$Species_richness)	#reasonable

xyplot(Species_richness ~ hpd|Within_PA, multiple.taxa.matched.landuse)


xyplot(Species_richness ~ elevation|Zone, multiple.taxa.matched.landuse)
xyplot(Species_richness ~ access|Zone, multiple.taxa.matched.landuse)
xyplot(Species_richness ~ ag_suit|Zone, multiple.taxa.matched.landuse)
xyplot(Species_richness ~ hpd|Zone, multiple.taxa.matched.landuse)
xyplot(Species_richness ~ slope|Zone, multiple.taxa.matched.landuse)




hist(multiple.taxa.matched.landuse$elevation)
hist(multiple.taxa.matched.landuse$ag_suit)
hist(multiple.taxa.matched.landuse$hpd)
hist(multiple.taxa.matched.landuse$DoP)
hist(log(multiple.taxa.matched.landuse$elevation+1))
hist(log(multiple.taxa.matched.landuse$hpd+1))
hist(log(PA.data$DoP+1))

plot(Species_richness ~ elevation, multiple.taxa.matched.landuse.s)
plot(Species_richness ~ log(elevation+1), multiple.taxa.matched.landuse.s)






### log distance to boundary analysis

fF <- c("Zone", "taxon_of_interest") 
fT <- list("ag_suit" = "3", "log_elevation" = "3", "log_slope" = "3", "log_bound_dist_km_PA_neg" = "3", "DoP.PA" = "3", "log_AREA.PA" = "3")
keepVars <- character(0)
fI <- character(0)
RS <-  c("log_bound_dist_km_PA_neg")

# without block: Species_richness~poly(ag_suit,3)+poly(log_AREA.PA,3)+poly(log_bound_dist_km_PA_neg,3)+poly(log_elevation,2)+poly(log_slope,2)+taxon_of_interest+Zone+(1+log_bound_dist_km_PA_neg|SS)+(1|SSBS)+(1|Predominant_habitat)"
# with block: Species_richness~poly(log_AREA.PA,3)+poly(log_elevation,1)+taxon_of_interest+Zone+(1+log_bound_dist_km_PA_neg|SS)+(1|SSBS)+(1|SSB)+(1|Predominant_habitat)
# 8 landuses: Species_richness~poly(ag_suit,3)+poly(log_bound_dist_km_PA_neg,2)+poly(log_elevation,2)+taxon_of_interest+Zone+(1+log_bound_dist_km_PA_neg|SS)+(1|SSBS)+(1|SSB)+(1|Predominant_habitat)




# add interactions
# and now keep confounding variables
fF <- c("Zone", "taxon_of_interest") 
fT <- list("log_bound_dist_km_PA_neg" = "2", "DoP.PA" = "1", "log_AREA.PA" = "1")
keepVars <- c("ag_suit" = "3", "log_elevation" = "2", "log_slope" = "1")
fI <- c("Zone:poly(log_bound_dist_km_PA_neg,2)", "taxon_of_interest:poly(log_bound_dist_km_PA_neg,2)")
RS <-  c("log_bound_dist_km_PA_neg")

# without block: Species_richness~poly(ag_suit,3)+poly(log_AREA.PA,3)+poly(log_bound_dist_km_PA_neg,3)+poly(log_elevation,2)+poly(log_slope,2)+taxon_of_interest+Zone+taxon_of_interest:poly(log_bound_dist_km_PA_neg,3)+(1+log_bound_dist_km_PA_neg|SS)+(1|SSBS)+(1|Predominant_habitat)


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


Species_richness.best.random$best.random #includes block and land use
 
best.random <- "(1+log_bound_dist_km_PA_neg|SS)+ (1|SSBS)+ (1|SSB)+(1|Predominant_habitat)"


# model select
Species_richness.model <- model_select(all.data  = multiple.taxa.matched.landuse, 
			     responseVar = "Species_richness",
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


validate(Species_richness.model$model) #ok

write.csv(Species_richness.model$stats, "Species_richness.model.stats.23.01.2015.csv")











### Within PA analysis


#run with no interactions first to see which variables have nonlinear relationships


fF <- c("Zone", "taxon_of_interest", "Within_PA") 
fT <- list("ag_suit" = "3", "log_slope" = "3", "log_elevation" = "3", "DoP.PA" = "3", "log_AREA.PA" = "3")
keepVars <- character(0)
fI <- character(0)
RS <-  c("Within_PA")

# without block: Species_richness~poly(ag_suit,3)+poly(log_AREA.PA,3)+poly(log_elevation,2)+poly(log_slope,2)+taxon_of_interest+Zone+(1+Within_PA|SS)+(1|SSBS)+(1|Predominant_habitat)
# with block: Species_richness~poly(log_AREA.PA,3)+poly(log_elevation,2)+taxon_of_interest+Zone+(1+Within_PA|SS)+(1|SSBS)+(1|SSB)+(1|Predominant_habitat)
# 8 landuse: Species_richness~poly(ag_suit,3)+poly(log_elevation,2)+taxon_of_interest+Zone+(1+Within_PA|SS)+(1|SSBS)+(1|SSB)+(1|Predominant_habitat)



# try all - doesnt work
fF <- c("Zone", "taxon_of_interest", "Within_PA") 
fT <- list("ag_suit" = "3", "log_elevation" = "2", "log_slope" = "2", "DoP.PA" = "1", "log_AREA.PA" = "3")
fI <- c("Within_PA:poly(ag_suit,3)", "Within_PA:poly(log_elevation,2)", "Within_PA:poly(log_slope,2)",
	"Within_PA:taxon_of_interest", "Within_PA:Zone",
	"taxon_of_interest:poly(DoP.PA,1)","taxon_of_interest:poly(log_AREA.PA,3)",
	"Zone:poly(DoP.PA,1)","Zone:poly(log_AREA.PA,3)")
RS <-  c("Within_PA")

#remove 2 PA area interactions - problems before with these. now works. 
fF <- c("Zone", "taxon_of_interest", "Within_PA") 
fT <- list("DoP.PA" = "1", "log_AREA.PA" = "1")
keepVars <-list("ag_suit" = "3", "log_elevation" = "2", "log_slope" = "1")
fI <- c("Within_PA:poly(ag_suit,3)", "Within_PA:poly(log_elevation,2)", "Within_PA:poly(log_slope,1)",
	"Within_PA:taxon_of_interest", "Within_PA:Zone",
	"taxon_of_interest:poly(DoP.PA,1)",
	"Zone:poly(DoP.PA,1)")
RS <-  c("Within_PA")


# without block: Species_richness~poly(ag_suit,3)+poly(log_AREA.PA,3)+poly(log_elevation,2)+poly(log_slope,2)+taxon_of_interest+Zone
#+Within_PA:poly(ag_suit,3)
#+Within_PA:poly(log_slope,2)
#+Within_PA+(1+Within_PA|SS)+(1|SSBS)+(1|Predominant_habitat)

#with block: 
#Species_richness~poly(log_AREA.PA,3)+taxon_of_interest+Zone+Within_PA:poly(ag_suit,1)+Within_PA:poly(log_slope,1)+Within_PA+poly(ag_suit,1)+poly(log_elevation,2)+poly(log_slope,1)
#+(1+Within_PA|SS)+(1|SSBS)+(1|SSB)+(1|Predominant_habitat)

#8 landuses: Species_richness~taxon_of_interest+Zone+Within_PA:poly(ag_suit,3)+Within_PA:poly(log_slope,1)+Within_PA+poly(ag_suit,3)+poly(log_elevation,2)+poly(log_slope,1)+(1+Within_PA|SS)+(1|SSBS)+(1|SSB)+(1|Predominant_habitat)

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


Species_richness.best.random$best.random #
 
best.random <- "(1+Within_PA|SS)+ (1|SSBS)+ (1|SSB) + (1|Predominant_habitat)"


# model select

Species_richness.model2 <- model_select(all.data  = multiple.taxa.matched.landuse, 
			     responseVar = "Species_richness",
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


validate(Species_richness.model2$model) #ok

write.csv(Species_richness.model2$stats, "Species_richness.model2.stats.23.01.2015.csv")



Species_richness.model2$warnings
Species_richness.model2$final.call
Species_richness.model2$stats


summary(Species_richness.model$model)
XV.10 <- CrossValidate(Species_richness.model$model,10, divFactor = "SS", fitFamily = "poisson", data= Species_richness.model$data, 
	fixedFactors = fF, fixedTerms = fT, fixedInteractions = fI, randomStruct = "(1+Within_PA|SS)+(1|SSBS)+(1|Predominant_habitat)")
# some convergence issues

XV.20 <- CrossValidate(Species_richness.model$model,20, divFactor = "SS", fitFamily = "poisson", data= Species_richness.model$data, 
	fixedFactors = fF, fixedTerms = fT, fixedInteractions = fI, randomStruct = "(1+Within_PA|SS)+(1|SSBS)+(1|Predominant_habitat)")
# DoP and log_bound_dist,2 are worrying as upper and lower estimates either side of 0
# many convergence issues

XV.50 <- CrossValidate(Species_richness.model$model,50, divFactor = "SS", fitFamily = "poisson", data= Species_richness.model$data, 
	fixedFactors = fF, fixedTerms = fT, fixedInteractions = fI, randomStruct = "(1+Within_PA|SS)+(1|SSBS)+(1|Predominant_habitat)")




XVr <- CrossValidate(Species_richness.model$model,-1, divFactor = "Predominant_habitat", fitFamily = "poisson", data= Species_richness.model$data, 
	fixedFactors = fF, fixedTerms = fT, fixedInteractions = fI, randomStruct = "(1+Within_PA|SS)+(1|SSBS)+(1|Predominant_habitat)")
# but if you dropped a whole land use it wouldnt be that surprising if the results changed? 

#glmmR.wikidot/faq

#













