

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



# get only data that model runs on 

data <- multiple.taxa.matched.landuse[,c("Richness_rarefied", "Zone", "taxon_of_interest", "ag_suit", "log_access",
	 "log_hpd", "log_bound_dist_km", "log_bound_dist_km_PA_neg", "Within_PA", "DoP.PA", "IUCN.PA", "log_AREA.PA", "Predominant_habitat", "SS", "SSBS")]

data <- na.omit(data) #have to get rid of NA to try poly on ag_suit




##########
# inc PA characteristics
############



# simplest model

m1 <- glmer(Richness_rarefied~Within_PA 
	+ (1+Within_PA|SS)+(1|Predominant_habitat)+ (1|SSBS), family = "poisson", data)
m2 <- glmer(Richness_rarefied~ 1
	+ (1+Within_PA|SS)+(1|Predominant_habitat)+ (1|SSBS), family = "poisson", data)




m1 <- glmer(Richness_rarefied~Within_PA 
	+ (1|SS)+(1|Predominant_habitat)+ (1|SSBS), family = "poisson", data)
m2 <- glmer(Richness_rarefied~ 1
	+ (1|SS)+(1|Predominant_habitat)+ (1|SSBS), family = "poisson", data)


m1 <- glmer(Richness_rarefied~Within_PA 
	+ (1+Within_PA|SS)+(1|SSBS), family = "poisson", data)
m2 <- glmer(Richness_rarefied~ 1
	+ (1+Within_PA|SS)+ (1|SSBS), family = "poisson", data)


m1 <- glmer(Richness_rarefied~Within_PA 
	+ (1|SS)+ (1|SSBS), family = "poisson", data)
m2 <- glmer(Richness_rarefied~ 1
	+ (1|SS)+ (1|SSBS), family = "poisson", data)







# simplest model

m1 <- glmer(Richness_rarefied~Within_PA + IUCN.PA + log_AREA.PA + poly(DoP.PA,1)
	+ (1+Within_PA|SS)+(1|Predominant_habitat)+ (1|SSBS), family = "poisson", data)


m2 <- glmer(Richness_rarefied~Within_PA + log_AREA.PA + poly(DoP.PA,1)
	+ (1+Within_PA|SS)+(1|Predominant_habitat)+ (1|SSBS), family = "poisson", data)

m3 <- glmer(Richness_rarefied~Within_PA + IUCN.PA  + poly(DoP.PA,1)
	+ (1+Within_PA|SS)+(1|Predominant_habitat) + (1|SSBS), family = "poisson", data)

m4 <- glmer(Richness_rarefied~Within_PA + IUCN.PA + log_AREA.PA
	+ (1+Within_PA|SS)+(1|Predominant_habitat)+ (1|SSBS), family = "poisson", data)

anova(m1, m2)
anova(m1, m3)
anova(m1, m4)
#only loss of PA size is significant #chisq = 7.8585, df = 1, p =0.005058







testmodel <- glmer(Richness_rarefied~Within_PA+Zone+taxon_of_interest +IUCN_CAT_number
	+ ag_suit + log_access + log_hpd 
	+ log_dist_in + log_dist_out + log_GIS_AREA + DoP
	+ (1+Within_PA|SS)+(1|Predominant_habitat) , multiple.taxa.matched.landuse)


testmodel1 <- glmer(Richness_rarefied~Zone+taxon_of_interest + IUCN.PA
	+ poly(ag_suit,1)+ poly(log_access,1) + poly(log_hpd,1) 
	+ poly(log_bound_dist_km,1)+ poly(log_AREA.PA,1) + poly(DoP.PA,1)
	+ (1+Within_PA+DoP.PA+log_AREA.PA|SS)+(1|Predominant_habitat) , data)
AIC(testmodel1) #-1951

testmodel2 <- glmer(Richness_rarefied~Zone+taxon_of_interest + IUCN.PA
	+ poly(ag_suit,2)+ poly(log_access,2) + poly(log_hpd,2) 
	+ poly(log_bound_dist_km,2)+ poly(log_AREA.PA,2) + poly(DoP.PA,2)
	+ (1+Within_PA+DoP.PA+log_AREA.PA|SS)+(1|Predominant_habitat) , data)
AIC(testmodel2) #-1973

testmodel3 <- glmer(Richness_rarefied~Zone+taxon_of_interest + IUCN.PA
	+ poly(ag_suit,3)+ poly(log_access,3) + poly(log_hpd,3) 
	+ poly(log_bound_dist_km,3)+ poly(log_AREA.PA,3) + poly(DoP.PA,3)
	+ (1+Within_PA+DoP.PA+log_AREA.PA|SS)+(1|Predominant_habitat) , data)
AIC(testmodel3) #-1981





### try with IUCN.PA

# start with all cubic to test main terms
fF <- c("Zone", "taxon_of_interest", "IUCN.PA") 
fT <- list("ag_suit" = "3", "log_hpd" = "3", "log_access" = "3","log_bound_dist_km_PA_neg" = "3", "DoP.PA" = "3", "log_AREA.PA" = "3")
fI <- character(0)
RS <-  character(0)
# Richness_rarefied~poly(ag_suit,3)+poly(log_AREA.PA,2)+Zone+(1|SS)+(1|SSBS)+(1|Predominant_habitat)"



# keep only significant nonlinear terms
# cant do them all 
fF <- c("Zone", "taxon_of_interest", "IUCN.PA") 
fT <- list("ag_suit" = "3", "log_hpd" = "1", "log_access" = "1", "log_bound_dist_km_PA_neg" = "1", "DoP.PA" = "1", "log_AREA.PA" = "2")
fI <- c("IUCN.PA:poly(ag_suit,3)", "IUCN.PA:poly(log_hpd,1)", "IUCN.PA:poly(log_access,1)",
	"IUCN.PA:taxon_of_interest", 
	"taxon_of_interest:poly(log_bound_dist_km_PA_neg,1)","taxon_of_interest:poly(DoP.PA,1)", "taxon_of_interest:poly(log_AREA.PA,2)",
	"IUCN.PA:Zone",
	"Zone:poly(log_bound_dist_km_PA_neg,1)","Zone:poly(DoP.PA,1)", "Zone:poly(log_AREA.PA,2)")
RS <-  character(0)




#doesnt work
fF <- c("Zone", "taxon_of_interest", "IUCN.PA") 
fT <- list("ag_suit" = "3", "log_hpd" = "1", "log_access" = "1", "log_bound_dist_km_PA_neg" = "1", "DoP.PA" = "1", "log_AREA.PA" = "2")
fI <- c("IUCN.PA:poly(ag_suit,3)", "IUCN.PA:poly(log_hpd,1)", "IUCN.PA:poly(log_access,1)",
	"taxon_of_interest:poly(log_bound_dist_km_PA_neg,1)")
RS <-  character(0)



# doesnt run properly either
fF <- c("Zone", "taxon_of_interest", "IUCN.PA") 
fT <- list("ag_suit" = "3", "log_hpd" = "1", "log_access" = "1", "log_bound_dist_km_PA_neg" = "1", "DoP.PA" = "1", "log_AREA.PA" = "2")
fI <- c("IUCN.PA:poly(ag_suit,3)", "IUCN.PA:poly(log_hpd,1)", "IUCN.PA:poly(log_access,1)")
RS <-  character(0)


#works
fF <- c("Zone", "taxon_of_interest", "IUCN.PA") 
fT <- list("ag_suit" = "3", "log_hpd" = "1", "log_access" = "1", "log_bound_dist_km_PA_neg" = "1", "DoP.PA" = "1", "log_AREA.PA" = "2")
fI <- c("taxon_of_interest:poly(log_bound_dist_km_PA_neg,1)","taxon_of_interest:poly(DoP.PA,1)", "taxon_of_interest:poly(log_AREA.PA,2)")
RS <-  character(0)


#works
fF <- c("Zone", "taxon_of_interest", "IUCN.PA") 
fT <- list("ag_suit" = "3", "log_hpd" = "1", "log_access" = "1", "log_bound_dist_km_PA_neg" = "1", "DoP.PA" = "1", "log_AREA.PA" = "2")
fI <- c("taxon_of_interest:poly(log_bound_dist_km_PA_neg,1)","taxon_of_interest:poly(DoP.PA,1)", "taxon_of_interest:poly(log_AREA.PA,2)",
	"Zone:poly(log_bound_dist_km_PA_neg,1)","Zone:poly(DoP.PA,1)", "Zone:poly(log_AREA.PA,2)")
RS <-  character(0)



#Richness_rarefied~poly(ag_suit,3)+poly(log_AREA.PA,2)+Zone+
#taxon_of_interest:poly(log_AREA.PA,2)
#+Zone:poly(log_AREA.PA,2)
#+taxon_of_interest+(1|SS)+(1|SSBS)+(1|Predominant_habitat)"







######
### put in everything apart from IUCN_CAT
######

fF <- c("Zone", "taxon_of_interest", "Within_PA") 
fT <- list("ag_suit" = "3", "log_hpd" = "3", "log_access" = "3", "log_bound_dist_km_PA_neg" = "3", "DoP.PA" = "3", "log_AREA.PA" = "3")
fI <- character(0)
RS <-  c("Within_PA")

# 07_14 "Richness_rarefied~poly(ag_suit,3)+Zone+(1+Within_PA|SS)+(1|SSBS)+(1|Predominant_habitat)"
# 11_14"Richness_rarefied~poly(ag_suit,3)+poly(log_AREA.PA,2)+(1+Within_PA|SS)+(1|SSBS)+(1|Predominant_habitat)"





# add interactions
fF <- c("Zone", "taxon_of_interest", "Within_PA") 
fT <- list("ag_suit" = "3", "log_hpd" = "1", "log_access" = "1", "log_bound_dist_km_PA_neg" = "1", "DoP.PA" = "1", "log_AREA.PA" = "2")
fI <- c("Within_PA:poly(ag_suit,3)", "Within_PA:poly(log_hpd,1)", "Within_PA:poly(log_access,1)",
	"Within_PA:taxon_of_interest", 
	"taxon_of_interest:poly(log_bound_dist_km_PA_neg,1)","taxon_of_interest:poly(DoP.PA,1)", "taxon_of_interest:poly(log_AREA.PA,2)",
	"Within_PA:Zone",
	"Zone:poly(log_bound_dist_km_PA_neg,1)","Zone:poly(DoP.PA,1)", "Zone:poly(log_AREA.PA,2)")
RS <-  c("Within_PA")



#07_14 Richness_rarefied~poly(ag_suit,3)+poly(DoP.PA,1)+poly(log_AREA.PA,1)+Zone
#+Within_PA:taxon_of_interest
#+Within_PA+taxon_of_interest+(1+Within_PA|SS)+(1|SSBS)+(1|Predominant_habitat)"

# 11_14 "Richness_rarefied~poly(ag_suit,3)+poly(log_AREA.PA,2)
#+Within_PA:poly(ag_suit,3)
#+taxon_of_interest:poly(DoP.PA,1)
#+Within_PA+taxon_of_interest+poly(DoP.PA,1)+(1+Within_PA|SS)+(1|SSBS)+(1|Predominant_habitat)"


# the main effects that are now retained were dropped at a stage earlier in the polynomial checking round
# as other polynomial terms were still in the model at that point.



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


Richness_rarefied.best.random$best.random
#"(1|SS)+ (1|SSBS)+(1|Predominant_habitat)"
#doesnt actually choose Within_PA RS - force this for comparison
 


# model select
#Richness_rarefied.model has the simple model
Richness_rarefied.model.int <- model_select(all.data  = multiple.taxa.matched.landuse, 
			     responseVar = "Richness_rarefied",
			     fitFamily = "poisson",
			     siteRandom = TRUE, 
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
                       fixedInteractions=fI,
                       randomStruct = "(1 + Within_PA|SS)+ (1|SSBS)+(1|Predominant_habitat)",
			     otherRandoms=c("Predominant_habitat"),
                       verbose=TRUE)


validate(Richness_rarefied.model$model) #ok
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

#write.csv(Richness_rarefied.model.int$stats, "Richness_rarefied.model.stats.15.10.2014.csv")








