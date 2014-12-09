

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




names(matched.landuse)




######
log_abundance
######

title(main = "all matched data")
hist(matched.landuse.s$CWM_Geographic_log_abundance_log10_square_km)

plot(log(hpd+1) ~ ag_suit, matched.landuse)


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
plot3d( matched.landuse$log_bound_dist_km_PA_neg, matched.landuse$log_access, matched.landuse$log_abundance) #reasonable
plot3d( matched.landuse$log_bound_dist_km_PA_neg, matched.landuse$log_hpd, matched.landuse$log_abundance)	#reasonable

xyplot(log_abundance ~ hpd|Within_PA, matched.landuse)


xyplot(log_abundance ~ elevation|Zone, matched.landuse)
xyplot(log_abundance ~ access|Zone, matched.landuse)
xyplot(log_abundance ~ ag_suit|Zone, matched.landuse)
xyplot(log_abundance ~ hpd|Zone, matched.landuse)
xyplot(log_abundance ~ slope|Zone, matched.landuse)



xyplot(log_abundance ~ ag_suit|IUCN.PA, matched.landuse)
xyplot(log_abundance ~ log_hpd|IUCN.PA, matched.landuse)

hist(matched.landuse$elevation)
hist(matched.landuse$ag_suit)
hist(matched.landuse$hpd)
hist(matched.landuse$DoP)
hist(log(matched.landuse$elevation+1))
hist(log(matched.landuse$hpd+1))
hist(log(PA.data$DoP+1))

plot(log_abundance ~ elevation, matched.landuse.s)
plot(log_abundance ~ log(elevation+1), matched.landuse.s)


tropical <- subset(matched.landuse, Zone == "Tropical")
temperate <- subset(matched.landuse, Zone == "Temperate")

table(matched.landuse$IUCN.PA, matched.landuse$taxon_of_interest)
table(tropical$IUCN.PA, tropical$taxon_of_interest)
table(temperate$IUCN.PA, temperate$taxon_of_interest)



low.protection <- subset(data, IUCN.PA == "4.5")
hist(low.protection$hpd)

# get only data that model runs on 

data <- matched.landuse[,c("log_abundance", "Zone", "taxon_of_interest", "ag_suit", "log_access",
	 "log_hpd", "log_bound_dist_km", "log_dist_in", "log_dist_out", "Within_PA", "DoP.PA", "IUCN.PA", "log_AREA.PA", "Predominant_habitat", "SS")]

data <- na.omit(data) #have to get rid of NA to try poly on ag_suit




##########
# inc PA characteristics
############


#simplest model



# simplest model

m1 <- lmer(log_abundance ~ Within_PA 
	+ (1+Within_PA|SS)+(1|Predominant_habitat), data)
m2 <- lmer(log_abundance~ 1
	+ (1+Within_PA|SS)+(1|Predominant_habitat), data)


m1 <- lmer(log_abundance~Within_PA 
	+ (1|SS)+(1|Predominant_habitat), data)
m2 <- lmer(log_abundance~ 1
	+ (1|SS)+(1|Predominant_habitat), data)



m1 <- lmer(log_abundance ~ Within_PA 
	+ (Within_PA|SS), data)
m2 <- lmer(log_abundance ~ 1
	+ (Within_PA|SS), data)

m1 <- lmer(log_abundance ~ Within_PA 
	+ (1|SS), data)
m2 <- lmer(log_abundance ~ 1
	+ (1|SS), data)

anova(m1, m2)

data$Within_PA <- relevel(data$Within_PA, "no")
m1 <- lmer(log_abundance~Within_PA 
	+ (1|SS), data)
z.outside <- fixef(m1)[1]
se.outside <- se.fixef(m1)[2]
upper.out <- z.outside + 1.96*se.outside
lower.out <- z.outside - 1.96*se.outside


data$Within_PA <- relevel(data$Within_PA, "yes")
m1 <- lmer(log_abundance~Within_PA 
	+ (1|SS), data)
z.inside <- fixef(m1)[1]
se.inside <- se.fixef(m1)[2]
upper.in <- z.inside + 1.96*se.inside 
lower.in = z.inside - 1.96*se.inside 


plot(c(z.outside, z.inside) ~ c(1,2), xlim = c(0,3), ylim = c(4.3,4.8), ylab = "Log abundance",
	pch = 16, bty = "l", xlab = "", xaxt = "n")
arrows(1,lower.out,1,upper.out, length = 0.1, angle = 90, code = 3)
arrows(2,lower.in,2,upper.in, length = 0.1, angle = 90, code = 3)
axis(1, c(1,2), c("Unprotected", "Protected"))






#  fit gamm to model without interactions to see what shape of model might be

gamm.model <- gamm4(log_abundance ~ Within_PA + Zone + taxon_of_interest + IUCN_CAT_number + ag_suit + s(log_access) + s(log_hpd)+ s(log_dist_in)+ s(log_dist_out) + s(DoP.PA) + s(log_AREA.PA), 
	random = ~(1+Within_PA|SS)+(1|Predominant_habitat), data =  matched.landuse)


# all look pretty linear
anova(gamm.model$gam)
plot(gamm.model$gam)



matched.landuse$IUCN.PA <- relevel(matched.landuse$IUCN.PA, "7")
# the model drops IUCN CAT = 0 anyway if within_PA is in the model 
# if we have IUCN with PA out = NA, model select will drop all the unprotected sites. 

testmodel <- lmer(log_abundance~ Zone + taxon_of_interest + IUCN.PA
	+ ag_suit + log_access + log_hpd 
	+ log_bound_dist_km_PA_neg + log_AREA.PA + DoP.PA
	+ (1|SS)+(1|Predominant_habitat) , matched.landuse)
summary(testmodel)

testmodel1 <- lmer(log_abundance~Zone+taxon_of_interest + IUCN.PA
	+ poly(ag_suit,1)+ poly(log_access,1) + poly(log_hpd,1) 
	+ poly(log_bound_dist_km,1)+ poly(log_AREA.PA,1) + poly(DoP.PA,1)
	+ (1+Within_PA+DoP.PA+log_AREA.PA|SS)+(1|Predominant_habitat) , data)
AIC(testmodel1) #-1951

testmodel2 <- lmer(log_abundance~Zone+taxon_of_interest + IUCN.PA
	+ poly(ag_suit,2)+ poly(log_access,2) + poly(log_hpd,2) 
	+ poly(log_bound_dist_km,2)+ poly(log_AREA.PA,2) + poly(DoP.PA,2)
	+ (1+Within_PA+DoP.PA+log_AREA.PA|SS)+(1|Predominant_habitat) , data)
AIC(testmodel2) #-1973

testmodel3 <- lmer(log_abundance~Zone+taxon_of_interest + IUCN.PA
	+ poly(ag_suit,3)+ poly(log_access,3) + poly(log_hpd,3) 
	+ poly(log_bound_dist_km,3)+ poly(log_AREA.PA,3) + poly(DoP.PA,3)
	+ (1+Within_PA+DoP.PA+log_AREA.PA|SS)+(1|Predominant_habitat) , data)
AIC(testmodel3) #-1981


#based on gamm 




# start with all cubic to test main terms

fF <- c("Zone", "taxon_of_interest", "IUCN.PA") 
fT <- list("ag_suit" = "3", "log_hpd" = "3", "log_access" = "3", "log_bound_dist_km_PA_neg" = "3", "DoP.PA" = "3", "log_AREA.PA" = "3")
fI <- character(0)
RS <-  character(0) #cant have IUCN_PA as not enough levels per site

#log_abundance~IUCN.PA+poly(log_AREA.PA,3)+poly(log_hpd,2)+taxon_of_interest+Zone+(1|SS)+(1|Predominant_habitat)




# keep only significant nonlinear terms
fF <- c("Zone", "taxon_of_interest", "IUCN.PA") 
fT <- list("ag_suit" = "1", "log_hpd" = "2", "log_access" = "1", "log_bound_dist_km_PA_neg" = "1", "DoP.PA" = "1", "log_AREA.PA" = "3")
fI <- c("IUCN.PA:poly(ag_suit,1)","IUCN.PA:poly(log_hpd,2)", "IUCN.PA:poly(log_access,1)",
	"IUCN.PA:taxon_of_interest",  "IUCN.PA:Zone",
	"taxon_of_interest:poly(log_bound_dist_km_PA_neg,1)","taxon_of_interest:poly(DoP.PA,1)", "taxon_of_interest:poly(log_AREA.PA,3)", 
	"Zone:poly(log_bound_dist_km_PA_neg,1)","Zone:poly(DoP.PA,1)", "Zone:poly(log_AREA.PA,3)")
RS <-  character(0)




#log_abundance~IUCN.PA+poly(log_AREA.PA,3)+poly(log_hpd,2)+taxon_of_interest+Zone
##+IUCN.PA:poly(log_hpd,2)
##+IUCN.PA:Zone
##+taxon_of_interest:poly(log_bound_dist_km_PA_neg,1)
#+taxon_of_interest:poly(log_AREA.PA,3)
##+Zone:poly(DoP.PA,1)
#+poly(log_bound_dist_km_PA_neg,1)+poly(DoP.PA,1)+(1|SS)+(1|Predominant_habitat)"
 




######
### put in everything apart from IUCN_CAT
######

fF <- c("Zone", "taxon_of_interest", "Within_PA") 
fT <- list("ag_suit" = "3", "log_hpd" = "3", "log_access" = "3", "log_bound_dist_km_PA_neg" = "3", "DoP.PA" = "3", "log_AREA.PA" = "3")
fI <- character(0)
RS <-  c("Within_PA")

# 07_14 log_abundance~poly(DoP.PA,1)+poly(log_AREA.PA,2)+poly(log_hpd,2)+taxon_of_interest+Zone
#+(1+Within_PA|SS)+(1|Predominant_habitat)"

# 11_14 log_abundance~poly(log_access,1)+poly(log_bound_dist_km_PA_neg,3)+poly(log_hpd,2)+taxon_of_interest+(1+Within_PA|SS)+(1|Predominant_habitat)"




# add interactions
fF <- c("Zone", "taxon_of_interest", "Within_PA") 
fT <- list("ag_suit" = "1", "log_hpd" = "2", "log_access" = "1", "log_bound_dist_km_PA_neg" = "3", "DoP.PA" = "1", "log_AREA.PA" = "1")
fI <- c("Within_PA:poly(ag_suit,1)", "Within_PA:poly(log_hpd,2)", "Within_PA:poly(log_access,1)",
	"Within_PA:taxon_of_interest", 
	"taxon_of_interest:poly(log_bound_dist_km_PA_neg,3)","taxon_of_interest:poly(DoP.PA,1)", "taxon_of_interest:poly(log_AREA.PA,1)",
	"Within_PA:Zone",
	"Zone:poly(log_bound_dist_km_PA_neg,3)","Zone:poly(DoP.PA,1)", "Zone:poly(log_AREA.PA,1)")
RS <-  c("Within_PA")


# 07_14 log_abundance~poly(DoP.PA,1)+poly(log_AREA.PA,2)+poly(log_hpd,2)+taxon_of_interest+Zone
#+Within_PA:poly(log_hpd,2)
#+taxon_of_interest:poly(log_bound_dist_km_PA_neg,1)
#+Within_PA:Zone
#+Zone:poly(DoP.PA,1)
#+Within_PA+poly(log_bound_dist_km_PA_neg,1)+(1+Within_PA|SS)+(1|Predominant_habitat)"
 

#11 _14 
#log_abundance~poly(log_access,1)+poly(log_bound_dist_km_PA_neg,3)+poly(log_hpd,2)+taxon_of_interest
#+Within_PA:poly(log_hpd,2)
#+taxon_of_interest:poly(log_bound_dist_km_PA_neg,3)
#+Zone:poly(DoP.PA,1)
#+Within_PA+Zone+poly(DoP.PA,1)+(1+Within_PA|SS)+(1|Predominant_habitat)"
 






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
                       fixedInteractions=fI,
                       randomStruct = log_abundance.best.random$best.random,
			     otherRandoms=c("Predominant_habitat"),
                       verbose=TRUE)




validate(log_abundance.model$model) #ok
log_abundance.model$warnings
log_abundance.model$final.call
log_abundance.model$stats


summary(log_abundance.model$model)
XV <- CrossValidate(log_abundance.model$model,50, divFactor = "SS", fitFamily = "gaussian", data= log_abundance.model$data, 
	fixedFactors = fF, fixedTerms = fT, fixedInteractions = fI, randomStruct = log_abundance.best.random$best.random)
# all significant estimates seem acceptably similar


XVr <- CrossValidate(log_abundance.model$model,-1, divFactor = "Predominant_habitat", fitFamily = "gaussian", data= log_abundance.model$data, 
	fixedFactors = fF, fixedTerms = fT, fixedInteractions = fI, randomStruct = log_abundance.best.random$best.random)
# but if you dropped a whole land use it wouldnt be that surprising if the results changed? 

#glmmR.wikidot/faq

write.csv(log_abundance.model$stats, "ab.model.stats.18.11.2014.csv")








