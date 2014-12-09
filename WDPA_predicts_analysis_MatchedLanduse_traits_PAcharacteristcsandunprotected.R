

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


setwd("N:/Documents/PREDICTS/WDPA analysis")



# load functions

#source("compare_randoms.R")
source("compare_randoms_lmer - with poly.R")
source("multiplot_pairwise_all.R")
source("bray_curtis_dissimilarity.R")
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








tapply(matched.landuse$mass, matched.landuse$taxon_of_interest, length)

plants <- subset(matched.landuse, taxon_of_interest == "Plants and Fungi")
inverts.data <- subset(matched.landuse, taxon_of_interest == "Invertebrates")






#########
MASS
#########

(just mammals and birds at the moment)

hist(matched.landuse$mass)
hist(log(matched.landuse$mass+1)) 
#logging it doesnt make a massive difference

names(matched.landuse)
plot(mass ~ Within_PA, matched.landuse)
plot(mass ~ taxon_of_interest, matched.landuse)
plot(mass ~ jitter(ag_suit), matched.landuse)
plot(mass ~ log_bound_dist_km, matched.landuse)
plot(mass ~ log_bound_dist_km_PA_neg, matched.landuse)


xyplot(mass ~ Within_PA|Zone, matched.landuse)
xyplot(mass ~ Within_PA|Realm, matched.landuse)
xyplot(mass ~ Within_PA|taxon_of_interest, matched.landuse)
xyplot(mass ~ log_bound_dist_km_PA_neg|Zone, matched.landuse)
xyplot(mass ~ log_bound_dist_km_PA_neg|taxon_of_interest, matched.landuse)





plot(elevation ~ access, matched.landuse.s) 	# not massively correlated
plot(elevation ~ slope, matched.landuse.s)	# not massively correlated
plot(access ~ slope, matched.landuse.s)		# not massively correlated

# nb 2/3 of this data is just one landuse - so only comparisons within these boxes are valid
# not the relative positions between them 
xyplot(mass ~ Within_PA|Predominant_habitat, matched.landuse)




# get only data that model runs on 

data <- matched.landuse[,c("mass", "Zone", "taxon_of_interest", "ag_suit", "log_access",
	 "log_hpd", "log_bound_dist_km", "log_dist_in", "log_dist_out", "Within_PA", "DoP.PA", "IUCN.PA", "log_AREA.PA", "Predominant_habitat", "SS")]

data <- na.omit(data) #have to get rid of NA to try poly on ag_suit





# just inc PA characteristics



unique(data$Zone)
unique(data$taxon_of_interest) # only verts have this data, cant include taxon of interest

#check for nonlinear effects
fF <- c("Zone", "Within_PA") 
fT <- list("ag_suit" = "3", "log_hpd" = "3", "log_access" = "3", 
		"log_bound_dist_km_PA_neg" = "3", "DoP.PA" = "3", "log_AREA.PA" = "3")
fI <- character(0)
RS <-  c("Within_PA")


# add interactions for non-linear terms
#mass~poly(ag_suit,2)+poly(log_AREA.PA,1)+poly(log_bound_dist_km_PA_neg,2)+poly(log_hpd,2)+(1+Within_PA|SS)+(1|Predominant_habitat)"


fF <- c("Zone", "Within_PA") 
fT <- list("ag_suit" = "2", "log_hpd" = "2", "log_access" = "1", "log_bound_dist_km_PA_neg" = "2", "DoP.PA" = "1", "log_AREA.PA" = "1")
fI <- c("Within_PA:poly(ag_suit,2)", "Within_PA:poly(log_hpd,2)", "Within_PA:poly(log_access,1)",
	"Within_PA:Zone",
	"Zone:poly(log_bound_dist_km_PA_neg,2)","Zone:poly(DoP.PA,1)", "Zone:poly(log_AREA.PA,1)")
RS <-  c("Within_PA")

#mass~poly(ag_suit,2)+poly(log_AREA.PA,1)+poly(log_bound_dist_km_PA_neg,2)+poly(log_hpd,2)+(1+Within_PA|SS)+(1|Predominant_habitat)"





mass.best.random <- compare_randoms(matched.landuse, "mass",
				fixedFactors=fF,
                         fixedTerms=fT,
                       fixedInteractions=fI,
                         otherRandoms=c("Predominant_habitat"),
				fixed_RandomSlopes = RS,
                          fitInteractions=FALSE,
				verbose=TRUE)


mass.best.random$best.random
#"(1|SS)+(1|Predominant_habitat)"
 


# model select
mass.model <- model_select(all.data  = matched.landuse, 
			     responseVar = "mass", 
                       fixedFactors= fF,
                       fixedTerms= fT,
                       fixedInteractions=fI,
                       randomStruct = "(1+Within_PA|SS)+(1|Predominant_habitat)",
			     otherRandoms=c("Predominant_habitat"),
                       verbose=TRUE)


validate(mass.model$model) #doesnt look ideal, slightly better logged...
mass.model$warnings


write.csv(mass.model$stats, "mass.model.stats.15.09.2012.csv")
















#########
VEG
#########

(just mammals and birds at the moment)

hist(matched.landuse$veg)
# cant log it, already logged, negative values

names(matched.landuse)
plot(veg ~ Within_PA, matched.landuse)
plot(veg ~ taxon_of_interest, matched.landuse)
plot(veg ~ jitter(ag_suit), matched.landuse)
plot(veg ~ log_bound_dist_km, matched.landuse)
plot(veg ~ log_bound_dist_km_PA_neg, matched.landuse)


xyplot(veg ~ Within_PA|Zone, matched.landuse)
xyplot(veg ~ Within_PA|Realm, matched.landuse)
xyplot(veg ~ Within_PA|taxon_of_interest, matched.landuse)
xyplot(veg ~ log_bound_dist_km_PA_neg|Zone, matched.landuse)
xyplot(veg ~ log_bound_dist_km_PA_neg|taxon_of_interest, matched.landuse)





# nb 2/3 of this data is just one landuse - so only comparisons within these boxes are valid
# not the relative positions between them 
xyplot(veg ~ Within_PA|Predominant_habitat, matched.landuse)


# just inc PA characteristics


#check for nonlinear effects
fF <- c("Zone", "Within_PA") 
fT <- list("ag_suit" = "3", "log_hpd" = "3", "log_access" = "3", 
		"log_bound_dist_km_PA_neg" = "3", "DoP.PA" = "3", "log_AREA.PA" = "3")
fI <- character(0)
RS <-  c("Within_PA")

#veg~poly(ag_suit,3)+poly(log_AREA.PA,3)+Zone+(1+Within_PA|SS)+(1|Predominant_habitat)



# add interactions for non-linear terms
fF <- c("Zone", "Within_PA") 
fT <- list("ag_suit" = "3", "log_hpd" = "1", "log_access" = "1", "log_bound_dist_km_PA_neg" = "1", "DoP.PA" = "1", "log_AREA.PA" = "3")
fI <- c("Within_PA:poly(ag_suit,3)", "Within_PA:poly(log_hpd,1)", "Within_PA:poly(log_access,1)",
	"Within_PA:Zone",
	"Zone:poly(log_bound_dist_km_PA_neg,1)","Zone:poly(DoP.PA,1)", "Zone:poly(log_AREA.PA,3)")
RS <-  c("Within_PA")

#veg~poly(ag_suit,3)+poly(log_AREA.PA,3)+Zone+(1+Within_PA|SS)+(1|Predominant_habitat)






veg.best.random <- compare_randoms(matched.landuse, "veg",
				fixedFactors=fF,
                         fixedTerms=fT,
                       fixedInteractions=fI,
                         otherRandoms=c("Predominant_habitat"),
				fixed_RandomSlopes = RS,
                          fitInteractions=FALSE,
				verbose=TRUE)


veg.best.random$best.random
#"(1+Within_PA|SS)+(1|Predominant_habitat)"


# model select
veg.model <- model_select(all.data  = matched.landuse, 
			     responseVar = "veg", 
                       fixedFactors= fF,
                       fixedTerms= fT,
                       fixedInteractions=fI,
                       randomStruct = "(1+Within_PA|SS)+(1|Predominant_habitat)",
			     otherRandoms=c("Predominant_habitat"),
                       verbose=TRUE)


validate(veg.model$model) #doesnt look ideal, slightly better logged...
veg.model$stats
veg.model$warnings


length(which(matched.landuse$veg > 0))
write.csv(veg.model$stats, "veg.model.stats.17.09.2014.csv")













#####
INVERTS - length derived volume
#####









data <- inverts.data[,c("vol", "Zone", "taxon_of_interest", "ag_suit", "log_access",
	 "log_hpd", "log_bound_dist_km", "Within_PA", "Predominant_habitat", "SS")]



nrow(data)
length(unique(data$SS))


#check for nonlinear effects
fF <- c("Zone", "Within_PA") 
fT <- list("ag_suit" = "3", "log_hpd" = "3", "log_access" = "3", 
		"log_bound_dist_km_PA_neg" = "3", "DoP.PA" = "3", "log_AREA.PA" = "3")
fI <- character(0)
RS <-  c("Within_PA")

#"vol~poly(ag_suit,1)+poly(log_access,3)+poly(log_AREA.PA,3)+poly(log_hpd,2)+Within_PA+(1+Within_PA|SS)+(1|Predominant_habitat)"





# add interactions for non-linear terms
fF <- c("Zone", "Within_PA") 
fT <- list("ag_suit" = "1", "log_hpd" = "2", "log_access" = "3", "log_bound_dist_km_PA_neg" = "1", "DoP.PA" = "1", "log_AREA.PA" = "3")
fI <- c("Within_PA:poly(ag_suit,1)", "Within_PA:poly(log_hpd,2)", "Within_PA:poly(log_access,3)",
	"Within_PA:Zone",
	"Zone:poly(log_bound_dist_km_PA_neg,1)","Zone:poly(DoP.PA,1)", "Zone:poly(log_AREA.PA,3)")
RS <-  c("Within_PA")

# vol~poly(ag_suit,1)+poly(log_access,3)+poly(log_AREA.PA,3)+poly(log_hpd,2)+Within_PA
#+Zone:poly(log_AREA.PA,3)
#+Zone +(1+Within_PA|SS)+(1|Predominant_habitat)"
 





vol.best.random <- compare_randoms(inverts.data, "vol",
				fixedFactors=fF,
                         fixedTerms=fT,
                       fixedInteractions=fI,
                         otherRandoms=c("Predominant_habitat"),
				fixed_RandomSlopes = RS,
                          fitInteractions=FALSE,
				verbose=TRUE)


vol.best.random$best.random
#(1|SS)
 


# model select
vol.model <- model_select(all.data  = inverts.data, 
			     responseVar = "vol", 
                       fixedFactors= fF,
                       fixedTerms= fT,
                       fixedInteractions=fI,
                       randomStruct = "(1+Within_PA|SS)+(1|Predominant_habitat)",
			     otherRandoms=c("Predominant_habitat"),
                       verbose=TRUE)


validate(vol.model$model) #doesnt look ideal, slightly better logged...
summary(vol.model$model)
vol.model$stats
vol.model$warnings


write.csv(vol.model$stats, "vol.model.stats.17.09.2014.csv")

nrow(veg.model$data)












#####
VERTS - length derived volume (amphibs)
#####






verts.data <- subset(matched.landuse, taxon_of_interest == "Vertebrates")



data <- verts.data[,c("vol", "Zone", "taxon_of_interest", "ag_suit", "log_access",
	 "log_hpd", "log_bound_dist_km", "Within_PA", "Predominant_habitat", "SS")]

data <- na.omit(data) #have to get rid of NA to try poly on ag_suit

nrow(data)




gamm.model <- gamm4(vol ~ ag_suit + s(log_access) + s(log_hpd)+ s(log_bound_dist_km), 
	random = ~(1+Within_PA+log_access+log_hpd|SS)+(1|Predominant_habitat), data =  data)
anova(gamm.model$gam)
plot(gamm.model$gam) = access linear, hpd linear, lbd quad
 
# but the model_select tests and reduces if not significant
# still need to know what the fit should be to specify for interaction terms


#all points are tropical
testmodel <- lmer(vol~ Within_PA  + ag_suit+ log_access + log_hpd + log_bound_dist_km
	+ (1+Within_PA|SS)+(1|Predominant_habitat), data)
AIC(testmodel) #275


testmodel1 <- lmer(vol~ Within_PA + poly(ag_suit,1)+ poly(log_access,1) + poly(log_hpd,1) + poly(log_bound_dist_km,1)
	+ (1+Within_PA|SS)+(1|Predominant_habitat) ,data)
AIC(testmodel1) #252

testmodel2 <- lmer(vol~ Within_PA +poly(ag_suit,2)+ poly(log_access,2) + poly(log_hpd,2) + poly(log_bound_dist_km,2)
	+ (1+Within_PA|SS)+(1|Predominant_habitat) , data)
AIC(testmodel2) #248

testmodel3 <- lmer(vol~ Within_PA + poly(ag_suit,3)+ poly(log_access,3) + poly(log_hpd,3) + poly(log_bound_dist_km,3)
	+ (1+Within_PA|SS)+(1|Predominant_habitat) , data)
AIC(testmodel3) #243

#based on gamm probably log hpd just linear
testmodel4 <- lmer(vol~ Within_PA + poly(ag_suit,3)+ poly(log_access,1) + poly(log_hpd,1) + poly(log_bound_dist_km,2)
	+ (1+Within_PA|SS)+(1|Predominant_habitat) ,data)
AIC(testmodel4) #249




# start with all cubic to test main terms
fF <- c("Within_PA") 
fT <- list("log_hpd" = "3", "log_access" = "3", "log_bound_dist_km" = "3", "ag_suit" = "3")
fI <- character(0)
RS <-  c("Within_PA", "log_access", "log_hpd", "log_bound_dist_km")
#2 warnings but ok as significance tested without error in subsequent iterations


# only keep what main terms only test suggests are needed
fF <- c("Within_PA") 
fT <- list("log_hpd" = "1", "log_access" = "1", "log_bound_dist_km" = "2", "ag_suit" = "2")
fI <- c("Within_PA:poly(log_bound_dist_km,2)", "Within_PA:poly(ag_suit,2)", "Within_PA:poly(log_hpd,1)", "Within_PA:poly(log_access,1)")
RS <-  c("Within_PA", "log_access", "log_hpd", "log_bound_dist_km")


# final call
# "vol~poly(ag_suit,2)+poly(log_bound_dist_km,2)+(1+log_access|SS)"


verts.data$log_vol <- log(verts.data$vol+1)


vol.best.random <- compare_randoms(verts.data, "vol",
				fixedFactors=fF,
                         fixedTerms=fT,
                       fixedInteractions=fI,
                         otherRandoms=c("Predominant_habitat"),
				fixed_RandomSlopes = RS,
                          fitInteractions=FALSE,
				verbose=TRUE)


vol.best.random$best.random
#(1+log_access|SS)
 


# model select
vol.model <- model_select(all.data  = verts.data, 
			     responseVar = "vol", 
                       fixedFactors= fF,
                       fixedTerms= fT,
                       fixedInteractions=fI,
                       randomStruct = vol.best.random$best.random,
			     otherRandoms=c("Predominant_habitat"),
                       verbose=TRUE)


validate(vol.model$model) #doesnt look ideal, slightly better logged...
vol.model$stats
vol.model$warnings


write.csv(vol.model$stats, "vol.amphib.model.stats.17.09.2014.csv")

nrow(veg.model$data)











