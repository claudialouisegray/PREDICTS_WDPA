

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
source("plotFactorInteraction.R")

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

nrow(matched.landuse) #5015

######
log_abundance
######

title(main = "all matched data")
hist(matched.landuse.s$CWM_Geographic_log_abundance_log10_square_km)

plot(log(slope+1) ~ ag_suit, matched.landuse)

plot(log_abundance ~ Within_PA, matched.landuse)
plot(log_abundance ~ IUCN.PA, matched.landuse)
plot(log_abundance ~ Zone, matched.landuse)
plot(log_abundance ~ taxon_of_interest, matched.landuse)
plot(log_abundance ~ elevation, matched.landuse)
plot(log_abundance ~ log_dist_in, matched.landuse)
plot(log_abundance ~ log_dist_out, matched.landuse)
plot(log_abundance ~ log_bound_dist_km_PA_neg, matched.landuse)

plot(log_abundance ~ DoP.PA, matched.landuse)

xyplot(log_abundance ~ Within_PA|Zone, matched.landuse)
xyplot(log_abundance ~ Within_PA|Realm, matched.landuse)
xyplot(log_abundance ~ Within_PA|taxon_of_interest, matched.landuse)

xyplot(log_abundance ~ log_bound_dist_km_PA_neg|Zone, matched.landuse)
xyplot(log_abundance ~ log_bound_dist_km_PA_neg|taxon_of_interest, matched.landuse)
xyplot(log_abundance ~ log_bound_dist_km_PA_neg|ag_suit, matched.landuse) # could be tricky to fit, few points in ag_suit = 1
plot3d( matched.landuse$log_bound_dist_km_PA_neg, matched.landuse$log_elevation, matched.landuse$log_abundance) #reasonable
plot3d( matched.landuse$log_bound_dist_km_PA_neg, matched.landuse$log_slope, matched.landuse$log_abundance)	#reasonable

xyplot(log_abundance ~ slope|Within_PA, matched.landuse)

xyplot(log_abundance ~ elevation|Zone, matched.landuse)
xyplot(log_abundance ~ elevation|Zone, matched.landuse)
xyplot(log_abundance ~ ag_suit|Zone, matched.landuse)
xyplot(log_abundance ~ slope|Zone, matched.landuse)
xyplot(log_abundance ~ slope|Zone, matched.landuse)

hist(matched.landuse$elevation)
hist(matched.landuse$ag_suit)
hist(matched.landuse$slope)
hist(matched.landuse$DoP)
hist(log(matched.landuse$elevation+1))
hist(log(matched.landuse$slope+1))
hist(log(PA.data$DoP+1))

plot(log_abundance ~ elevation, matched.landuse.s)
plot(log_abundance ~ log(elevation+1), matched.landuse.s)


tropical <- subset(matched.landuse, Zone == "Tropical")
temperate <- subset(matched.landuse, Zone == "Temperate")

table(matched.landuse$IUCN.PA, matched.landuse$taxon_of_interest)
table(tropical$IUCN.PA, tropical$taxon_of_interest)
table(temperate$IUCN.PA, temperate$taxon_of_interest)







### distance to boundary analysis
# check polynomials
fF <- c("Zone", "taxon_of_interest") 
fT <- list("ag_suit" = "3", "log_slope" = "3", "log_elevation" = "3", "log_bound_dist_km_PA_neg" = "3", "DoP.PA" = "3", "log_AREA.PA" = "3")
keepVars <- character(0)
fI <- character(0)
RS <-  c("log_bound_dist_km_PA_neg")

# without block:log_abundance~poly(log_AREA.PA,3)+poly(log_elevation,2)+poly(log_slope,2)+taxon_of_interest+(1+log_bound_dist_km_PA_neg|SS)+(1|Predominant_habitat)
# with block: log_abundance~poly(log_bound_dist_km_PA_neg,3)+poly(log_elevation,2)+poly(log_slope,2)+taxon_of_interest+(1+log_bound_dist_km_PA_neg|SS)+(1|SSB)+(1|Predominant_habitat)"
# all 8 landuses in random:  log_abundance~poly(ag_suit,3)+poly(log_bound_dist_km_PA_neg,3)+poly(log_elevation,2)+poly(log_slope,2)+taxon_of_interest+(1+log_bound_dist_km_PA_neg|SS)+(1|SSB)+(1|Predominant_habitat)

# look at gamm

#doesnt run well
gamm.model <- gamm4(log_abundance ~ Zone + taxon_of_interest + ag_suit + s(log_elevation) + s(log_slope)+ s(log_bound_dist_km_PA_neg),
# + s(DoP.PA) + s(log_AREA.PA), 
	random = ~(1+log_bound_dist_km_PA_neg|SS)+ (1|SSB) + (1|Predominant_habitat), data =  matched.landuse)

gamm.model.dist <- gamm4(log_abundance ~ s(log_bound_dist_km_PA_neg), 
	random = ~(1+log_bound_dist_km_PA_neg|SS)+ (1|SSB) + (1|Predominant_habitat), data =  matched.landuse)

# the output of this is the same as the linear model below. 
gamm.model.dist.conf <- gamm4(log_abundance ~ s(log_bound_dist_km_PA_neg)+ ag_suit + s(log_elevation) + s(log_slope), 
	random = ~(1+log_bound_dist_km_PA_neg|SS)+ (1|SSB) + (1|Predominant_habitat), data =  matched.landuse)

#doesnt test for this properly
gamm.model.dist.conf.int <- gamm4(log_abundance ~ s(log_bound_dist_km_PA_neg)+ ag_suit + s(log_elevation) + s(log_slope)
	+ s(log_bound_dist_km_PA_neg, by = taxon_of_interest), 
	random = ~(1+log_bound_dist_km_PA_neg|SS)+ (1|SSB) + (1|Predominant_habitat), data =  matched.landuse)


# look at gamm output
anova(gamm.model$gam)
plot(gamm.model$gam)



# add interactions
# and now keep confounding variables
fF <- c("Zone", "taxon_of_interest") 
fT <- list("log_bound_dist_km_PA_neg" = "3", "DoP.PA" = "1", "log_AREA.PA" = "1")
keepVars <- c("ag_suit" = "3", "log_slope" = "2", "log_elevation" = "2")
fI <- c("Zone:poly(log_bound_dist_km_PA_neg,3)", "taxon_of_interest:poly(log_bound_dist_km_PA_neg,3)")
RS <-  c("log_bound_dist_km_PA_neg")

#without block: log_abundance~poly(log_AREA.PA,3)+poly(log_elevation,2)+poly(log_slope,2)+taxon_of_interest+(1+log_bound_dist_km_PA_neg|SS)+(1|Predominant_habitat)"

# with block : "log_abundance~poly(log_bound_dist_km_PA_neg,3)+taxon_of_interest
# taxon_of_interest:poly(log_bound_dist_km_PA_neg,3)+poly(ag_suit,1)+poly(log_slope,2)+poly(log_elevation,2)+(1+log_bound_dist_km_PA_neg|SS)+(1|SSB)+(1|Predominant_habitat)

# 8 landuses:log_abundance~poly(log_bound_dist_km_PA_neg,3)+taxon_of_interest
#+taxon_of_interest:poly(log_bound_dist_km_PA_neg,3)+poly(ag_suit,3)+poly(log_slope,2)+poly(log_elevation,2)+(1+log_bound_dist_km_PA_neg|SS)+(1|SSB)+(1|Predominant_habitat)"




log_abundance.best.random <- compare_randoms(matched.landuse, "log_abundance",
				fixedFactors=fF,
                         fixedTerms=fT,
                       fixedInteractions=fI,
                         otherRandoms=c("Predominant_habitat"),
				fixed_RandomSlopes = RS,
                          fitInteractions=FALSE,
				verbose=TRUE)


log_abundance.best.random$best.random 
 


# model select
log_abundance.model <- model_select(all.data  = matched.landuse, 
			     responseVar = "log_abundance", 
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
			     keepVars = keepVars,
                       fixedInteractions=fI,
                       randomStruct = log_abundance.best.random$best.random,
			     otherRandoms=c("Predominant_habitat"),
                       verbose=TRUE)


write.csv(log_abundance.model$stats, "ab.model.stats.23.01.2014.csv")








### age/size analysis



xyplot(log_abundance ~ taxon_of_interest|AREA_DoP, matched.landuse)


fF <- c("Zone", "taxon_of_interest", "AREA_DoP") 
fT <- list("ag_suit" = "3", "log_slope" = "3", "log_elevation" = "3")
keepVars <- character(0)
fI <- character(0)
RS <-  character(0)
# "log_abundance~poly(ag_suit,1)+taxon_of_interest+(1|SS)+(1|SSB)+(1|Predominant_habitat)"


# add interactions 
fF <- c("Zone", "taxon_of_interest", "AREA_DoP") 
fT <- list()
keepVars <- list("ag_suit" = "1", "log_slope" = "1", "log_elevation" = "1")
fI <- c("AREA_DoP:taxon_of_interest", "AREA_DoP:Zone")
RS <-  character(0)
#log_abundance~taxon_of_interest+AREA_DoP:taxon_of_interest+AREA_DoP+poly(ag_suit,1)+poly(log_slope,1)+poly(log_elevation,1)+(1|SS)+(1|SSB)+(1|Predominant_habitat)"




# check gamm
gamm.model <- gamm4(log_abundance ~ Zone + taxon_of_interest + ag_suit + s(log_elevation) + s(log_slope)+ Within_PA,
	random = ~(1+Within_PA|SS)+ (1|SSB) + (1|Predominant_habitat), data =  matched.landuse)
anova(gamm.model$gam)
plot(gamm.model$gam)

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
log_abundance.model2 <- model_select(all.data  = matched.landuse, 
			     responseVar = "log_abundance", 
			     keepVars = keepVars,
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
                       fixedInteractions=fI,
                       randomStruct = log_abundance.best.random$best.random,
			     otherRandoms=c("Predominant_habitat"),
                       verbose=TRUE)


write.csv(log_abundance.model2$stats, "N:/Documents/PREDICTS/WDPA analysis/stats tables all/age size bins/ab.model.age.size.stats.10.02.2014.csv")

log_abundance.model2$stats
validate(log_abundance.model2$model) #ok
summary(log_abundance.model2$model)
log_abundance.model2$warnings



summary(log_abundance.model2$model)
fT <- list("DoP.PA" = "1", "log_AREA.PA" = "3", "ag_suit" = "1", log_elevation = "1")

XV <- CrossValidate(log_abundance.model2$model,10, divFactor = "SS", fitFamily = "gaussian", data= log_abundance.model2$data, 
	fixedFactors = fF, fixedTerms = fT, fixedInteractions = fI, randomStruct = "(1+Within_PA|SS)+ (1|SSB) + (1|Predominant_habitat)")
# all significant estimates seem acceptably similar
XV
#isnt what I really want this process but the whole model select, and the upper and lower p values? 

XVr <- CrossValidate(log_abundance.model$model,-1, divFactor = "Predominant_habitat", fitFamily = "gaussian", data= log_abundance.model$data, 
	fixedFactors = fF, fixedTerms = fT, fixedInteractions = fI, randomStruct = "(1+Within_PA|SS)+ (1|SSB) + (1|Predominant_habitat)")
# but if you dropped a whole land use it wouldnt be that surprising if the results changed? 

#glmmR.wikidot/faq



tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/02_15/abundance vs size and age vs taxon.tif",
	width = 25, height =20, units = "cm", pointsize = 12, res = 300)

plotFactorInteraction(model = log_abundance.model2$model,
responseVar = "log_abundance",
data = log_abundance.model2$data,
xvar = "AREA_DoP",
xvar.order = c("small_young", "small_old", "large_young", "large_old"), 
#this must be the levels of the factor in the order to be plotted, not including reference level
intvar = "taxon_of_interest",
intvar.order = c("Invertebrates", "Vertebrates"),
logLink = "e",
xlab = "PA size and age class",
ylab = "Relative abundance %")

text(0.7,90, "Young = 0 - 20 yrs \nOld = 20 - 85 yrs \nSmall = 0 - 400 km2 \nLarge = 400 - 12000km2", 
	adj = 0, cex = 0.8)

dev.off()


table(data$AREA_DoP, data$taxon_of_interest)

### Within_PA analysis


fF <- c("Zone", "taxon_of_interest", "Within_PA") 
fT <- list("ag_suit" = "3", "log_slope" = "3", "log_elevation" = "3", "DoP.PA" = "3", "log_AREA.PA" = "3")
keepVars <- character(0)
fI <- character(0)
RS <-  c("Within_PA")


#without block: log_abundance~poly(ag_suit,1)+poly(log_AREA.PA,3)+poly(log_elevation,3)+taxon_of_interest+(1+Within_PA|SS)+(1|Predominant_habitat)"
#with block: log_abundance~poly(ag_suit,1)+poly(log_AREA.PA,3)+poly(log_elevation,1)+taxon_of_interest+(1+Within_PA|SS)+(1|SSB)+(1|Predominant_habitat)"
# 8 landuses:  log_abundance~poly(ag_suit,3)+poly(log_elevation,2)+taxon_of_interest+(1+Within_PA|SS)+(1|SSB)+(1|Predominant_habitat)

# check gamm
gamm.model <- gamm4(log_abundance ~ Zone + taxon_of_interest + ag_suit + s(log_elevation) + s(log_slope)+ Within_PA,
# + s(DoP.PA) + s(log_AREA.PA), 
	random = ~(1+Within_PA|SS)+ (1|SSB) + (1|Predominant_habitat), data =  matched.landuse)
anova(gamm.model$gam)
plot(gamm.model$gam)



# add interactions
fF <- c("Zone", "taxon_of_interest", "Within_PA") 
fT <- list("DoP.PA" = "1", "log_AREA.PA" = "1")
keepVars <- c("ag_suit" = "3", "log_slope" = "1", "log_elevation" = "2")
fI <- c("Within_PA:poly(ag_suit,3)", "Within_PA:poly(log_slope,1)", "Within_PA:poly(log_elevation,2)",
	"Within_PA:taxon_of_interest", 
	"taxon_of_interest:poly(DoP.PA,1)", "taxon_of_interest:poly(log_AREA.PA,1)",
	"Within_PA:Zone",
	"Zone:poly(DoP.PA,1)", "Zone:poly(log_AREA.PA,1)")
RS <-  c("Within_PA")

#without block: log_abundance~poly(ag_suit,1)+poly(log_AREA.PA,3)+poly(log_elevation,3)+taxon_of_interest
#+Within_PA:poly(ag_suit,1)
#+Within_PA:poly(log_slope,1)
#+Within_PA+poly(log_slope,1)+(1+Within_PA|SS)+(1|Predominant_habitat)"

#with block and keeping in confounding variables
#log_abundance~poly(log_AREA.PA,3)+taxon_of_interest
#+Within_PA:poly(ag_suit,1)
#+Within_PA:poly(log_slope,1)
#+Zone:poly(DoP.PA,1)
#+Zone:poly(log_AREA.PA,3)+Within_PA+Zone+poly(DoP.PA,1)+poly(ag_suit,1)+poly(log_slope,1)+poly(log_elevation,1)
#+(1+Within_PA|SS)+(1|SSB)+(1|Predominant_habitat)"

# 8 landuses
# log_abundance~taxon_of_interest+poly(ag_suit,3)+poly(log_slope,1)+poly(log_elevation,2)+(1+Within_PA|SS)+(1|SSB)+(1|Predominant_habitat)

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
log_abundance.model2 <- model_select(all.data  = matched.landuse, 
			     responseVar = "log_abundance", 
			     keepVars = keepVars,
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
                       fixedInteractions=fI,
                       randomStruct = log_abundance.best.random$best.random,
			     otherRandoms=c("Predominant_habitat"),
                       verbose=TRUE)

log_abundance.model2_NM <- model_select(all.data  = matched.landuse, 
			     responseVar = "log_abundance", 
			     keepVars = keepVars,
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
                       fixedInteractions=fI,
                       randomStruct = "(1+Within_PA|SS)+ (1|SSB) + (1|Predominant_habitat)",
			     otherRandoms=c("Predominant_habitat"),
			     optimizer="Nelder_Mead",
                       verbose=TRUE)

# same result both optimisers

write.csv(log_abundance.model2$stats, "ab.model2.stats.23.01.2014.csv")


validate(log_abundance.model$model) #ok
log_abundance.model2$warnings
log_abundance.model2$final.call
log_abundance.model2$stats


summary(log_abundance.model2$model)
fT <- list("DoP.PA" = "1", "log_AREA.PA" = "3", "ag_suit" = "1", log_elevation = "1")

XV <- CrossValidate(log_abundance.model2$model,10, divFactor = "SS", fitFamily = "gaussian", data= log_abundance.model2$data, 
	fixedFactors = fF, fixedTerms = fT, fixedInteractions = fI, randomStruct = "(1+Within_PA|SS)+ (1|SSB) + (1|Predominant_habitat)")
# all significant estimates seem acceptably similar
XV
#isnt what I really want this process but the whole model select, and the upper and lower p values? 

XVr <- CrossValidate(log_abundance.model$model,-1, divFactor = "Predominant_habitat", fitFamily = "gaussian", data= log_abundance.model$data, 
	fixedFactors = fF, fixedTerms = fT, fixedInteractions = fI, randomStruct = "(1+Within_PA|SS)+ (1|SSB) + (1|Predominant_habitat)")
# but if you dropped a whole land use it wouldnt be that surprising if the results changed? 

#glmmR.wikidot/faq




