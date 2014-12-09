

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



#min(data$bound_dist_km_PA_neg)

#data$Source_ID[which(data$bound_dist_km_PA_neg < -50)]

#matched.landuse <- subset(matched.landuse, Source_ID != "SE2_2010__McCarthy")


nrow(matched.landuse) #5471




######
range
######

title(main = "all matched data")
hist(matched.landuse.s$CWM_Geographic_range_log10_square_km)

plot(log(hpd+1) ~ ag_suit, matched.landuse)


plot(range ~ Within_PA, matched.landuse)
plot(range ~ IUCN.PA, matched.landuse)
plot(range ~ Zone, matched.landuse)
plot(range ~ taxon_of_interest, matched.landuse)
plot(range ~ elevation, matched.landuse)
plot(range ~ log_bound_dist_km_PA_neg, matched.landuse)
plot(range ~ bound_dist_km_PA_neg, matched.landuse)


plot(range ~ DoP.PA, matched.landuse)


xyplot(range ~ Within_PA|Zone, matched.landuse)
xyplot(range ~ Within_PA|Realm, matched.landuse)
xyplot(range ~ Within_PA|taxon_of_interest, matched.landuse)

xyplot(range ~ log_bound_dist_km_PA_neg|Zone, matched.landuse)
xyplot(range ~ log_bound_dist_km_PA_neg|taxon_of_interest, matched.landuse)
xyplot(range ~ log_bound_dist_km_PA_neg|ag_suit, matched.landuse) # could be tricky to fit, few points in ag_suit = 1
plot3d( matched.landuse$log_bound_dist_km_PA_neg, matched.landuse$log_access, matched.landuse$range) #reasonable
plot3d( matched.landuse$log_bound_dist_km_PA_neg, matched.landuse$log_hpd, matched.landuse$range)	#reasonable

xyplot(range ~ hpd|Within_PA, matched.landuse)


xyplot(range ~ elevation|Zone, matched.landuse)
xyplot(range ~ access|Zone, matched.landuse)
xyplot(range ~ ag_suit|Zone, matched.landuse)
xyplot(range ~ hpd|Zone, matched.landuse)
xyplot(range ~ slope|Zone, matched.landuse)




hist(matched.landuse$elevation)
hist(matched.landuse$ag_suit)
hist(matched.landuse$hpd)
hist(matched.landuse$DoP)
hist(log(matched.landuse$elevation+1))
hist(log(matched.landuse$hpd+1))
hist(log(PA.data$DoP+1))

plot(range ~ elevation, matched.landuse.s)
plot(range ~ log(elevation+1), matched.landuse.s)



# get only data that model runs on 

data <- matched.landuse[,c("range", "Zone", "taxon_of_interest", "ag_suit", "log_access",
	 "log_hpd", "log_bound_dist_km", "Within_PA", "DoP.PA", "IUCN.PA", "log_AREA.PA", "Predominant_habitat", "SS")]

data <- na.omit(data) #have to get rid of NA to try poly on ag_suit




#### IUCN only model

m1 <- lmer(range~ Zone +taxon_of_interest + IUCN.PA
	+ ag_suit + log_slope + log_elevation 
	+ DoP.PA + log_AREA.PA 
	+ (1|SS)+(1|Predominant_habitat), matched.landuse)

m2 <- update(m1, .~. -IUCN.PA)
m3 <- update(m1, .~. -log_AREA.PA)
m4 <- update(m1, .~. -DoP.PA)

anova(m1, m2)
anova(m1, m3)
anova(m1, m4)


m1 <- lmer(range~ Zone +taxon_of_interest + IUCN.PA
	+ ag_suit + log_slope + log_elevation 
	+ (1|SS)+(1|Predominant_habitat), matched.landuse)


#adding Within_PA is weird - gives coefficents so that site can be either IUCN cats III to VI or "yes"
summary(m1)





### Log dist to boundary model 

#  fit gamm to model without interactions to see what shape of model might be

gamm.model <- gamm4(range ~ Within_PA + Zone + taxon_of_interest + IUCN_CAT_number + ag_suit + s(log_access) + s(log_hpd)+ s(log_dist_in)+ s(log_dist_out) + s(DoP.PA) + s(log_AREA.PA), 
	random = ~(1+Within_PA|SS)+(1|Predominant_habitat), data =  matched.landuse)


# all look pretty linear
anova(gamm.model$gam)
plot(gamm.model$gam)





testmodel <- lmer(range~Within_PA+Zone+taxon_of_interest +IUCN_CAT_number
	+ ag_suit + log_access + log_hpd 
	+ log_dist_in + log_dist_out + log_GIS_AREA + DoP
	+ (1+Within_PA|SS)+(1|Predominant_habitat) , matched.landuse)


testmodel1 <- lmer(range~Zone+taxon_of_interest + IUCN.PA
	+ poly(ag_suit,1)+ poly(log_access,1) + poly(log_hpd,1) 
	+ poly(log_bound_dist_km,1)+ poly(log_AREA.PA,1) + poly(DoP.PA,1)
	+ (1+Within_PA+DoP.PA+log_AREA.PA|SS)+(1|Predominant_habitat) , data)
AIC(testmodel1) #-1951

testmodel2 <- lmer(range~Zone+taxon_of_interest + IUCN.PA
	+ poly(ag_suit,2)+ poly(log_access,2) + poly(log_hpd,2) 
	+ poly(log_bound_dist_km,2)+ poly(log_AREA.PA,2) + poly(DoP.PA,2)
	+ (1+Within_PA+DoP.PA+log_AREA.PA|SS)+(1|Predominant_habitat) , data)
AIC(testmodel2) #-1973

testmodel3 <- lmer(range~Zone+taxon_of_interest + IUCN.PA
	+ poly(ag_suit,3)+ poly(log_access,3) + poly(log_hpd,3) 
	+ poly(log_bound_dist_km,3)+ poly(log_AREA.PA,3) + poly(DoP.PA,3)
	+ (1+Within_PA+DoP.PA+log_AREA.PA|SS)+(1|Predominant_habitat) , data)
AIC(testmodel3) #-1981


#based on gamm 




# USING IUCN instead of within_PA

# start with all cubic to test main terms

# log neg dist inside PA
fF <- c("Zone", "taxon_of_interest", "IUCN.PA") 
fT <- list("ag_suit" = "3", "log_hpd" = "3", "log_access" = "3", "log_bound_dist_km_PA_neg" = "3", "DoP.PA" = "3", "log_AREA.PA" = "3")
fI <- character(0)
RS <-  character(0)


# with no Within_PA, just IUCN cat in the model
# weird dist to boundary results.... but if IUCN cat has an effect, wouldnt it make more sense to do dist to boundary for different IUCN cats...
# how to allocate?

#range~IUCN.PA+poly(ag_suit,3)+poly(log_access,2)+poly(log_AREA.PA,3)+poly(log_hpd,3)+taxon_of_interest+Zone
#+(1|SS)+(1|Predominant_habitat)"

# with no Within_PA, just IUCN cat in the model, and IUCN RS - but equally doesnt make sense to do that, as few IUCN vals per study
# convergence issues too
#range~poly(ag_suit,3)+poly(log_access,2)+poly(log_hpd,3)+taxon_of_interest+Zone+(1+IUCN.PA|SS)+(1|Predominant_habitat)"

#log dist to boundary doesnt work as RS


# within_PA and IUCN.PA in the model - same if RS the same
# BUT THIS DOESNT MAKE SENSE TO SET AS A MODEL - rank deficiencies
#range~IUCN.PA+poly(ag_suit,3)+poly(log_access,2)+poly(log_AREA.PA,3)+poly(log_hpd,3)+taxon_of_interest+Zone+(1|SS)+(1|Predominant_habitat)

# within_PA and IUCN.PA in the model - if within_PA in RS
#range~poly(ag_suit,3)+poly(log_access,2)+poly(log_AREA.PA,3)+poly(log_bound_dist_km_PA_neg,3)+poly(log_hpd,3)+taxon_of_interest+Within_PA+Zone+(1+Within_PA|SS)+(1|Predominant_habitat)"




# keep only significant nonlinear terms
fF <- c("Zone", "taxon_of_interest", "IUCN.PA") 
fT <- list("ag_suit" = "3", "log_hpd" = "3", "log_access" = "2", "log_bound_dist_km_PA_neg" = "1", "DoP.PA" = "1", "log_AREA.PA" = "3")
fI <- c("IUCN.PA:poly(ag_suit,3)", "IUCN.PA:poly(log_hpd,3)", "IUCN.PA:poly(log_access,2)",
	"IUCN.PA:taxon_of_interest", 	"IUCN.PA:Zone",
	"taxon_of_interest:poly(log_bound_dist_km_PA_neg,1)","taxon_of_interest:poly(DoP.PA,1)", "taxon_of_interest:poly(log_AREA.PA,3)", 
	"Zone:poly(log_bound_dist_km_PA_neg,1)","Zone:poly(DoP.PA,1)", "Zone:poly(log_AREA.PA,3)")
RS <-  character(0)


#range~IUCN.PA+poly(ag_suit,3)+poly(log_access,2)+poly(log_AREA.PA,3)+poly(log_hpd,3)+taxon_of_interest+Zone
##+IUCN.PA:poly(ag_suit,3)
##+IUCN.PA:poly(log_hpd,3)
##+taxon_of_interest:poly(log_bound_dist_km_PA_neg,1)
##+Zone:poly(DoP.PA,1)
##+Zone:poly(log_AREA.PA,3)
#+poly(log_bound_dist_km_PA_neg,1)+poly(DoP.PA,1)+(1|SS)+(1|Predominant_habitat)"






#just look at interactions with IUCN
#range~IUCN.PA+poly(ag_suit,3)+poly(log_access,2)+poly(log_AREA.PA,3)+poly(log_hpd,3)+taxon_of_interest+Zone
#+IUCN.PA:poly(ag_suit,3)+IUCN.PA:poly(log_hpd,3)+(1|SS)+(1|Predominant_habitat)"
# seems to be very different to what I had before, and dont really trust it as now so dependent on fewer IUCN cat I & II points 
# theres very little going on thats different in different IUCN cats in the plot - says its significant but doesnt look it
# same with PA AREA - says significant but plot unconvincing
# interactions with ag_suit and hpd also look susceptible to outliers or lack of data spread




#its still more interesting to look at PA vs Not overall, then PA traits later without interactions? 


fF <- c("Zone", "taxon_of_interest", "Within_PA") 
fT <- list("ag_suit" = "3", "log_hpd" = "3", "log_access" = "3", "log_bound_dist_km_PA_neg" = "3")
fI <- character(0)
RS <-  c("Within_PA", "log_access", "log_hpd")

#only within_PA RS
#range~poly(ag_suit,2)+poly(log_access,2)+poly(log_bound_dist_km_PA_neg,2)+poly(log_hpd,3)+taxon_of_interest+Zone+(1+Within_PA|SS)+(1|Predominant_habitat)"
# Within PA and log hpd and log access
#range~poly(ag_suit,3)+poly(log_access,2)+poly(log_bound_dist_km_PA_neg,2)+taxon_of_interest+Zone+(1+Within_PA+log_access+log_hpd|SS)+(1|Predominant_habitat)"


# retain sig non-linear terms
fF <- c("Zone", "taxon_of_interest", "Within_PA") 
fT <- list("ag_suit" = "3", "log_hpd" = "1", "log_access" = "2", "log_bound_dist_km_PA_neg" = "2")
fI <- c("Within_PA:poly(ag_suit,3)", "Within_PA:poly(log_hpd,1)", "Within_PA:poly(log_access,2)",
	"Within_PA:taxon_of_interest", 	"Within_PA:Zone", 
	"taxon_of_interest:poly(log_bound_dist_km_PA_neg,2)", "Zone:poly(log_bound_dist_km_PA_neg,2)")


#range~poly(ag_suit,3)+poly(log_access,2)+poly(log_bound_dist_km_PA_neg,2)+taxon_of_interest+Zone
#+Within_PA:poly(ag_suit,3)
#+Within_PA:poly(log_hpd,1)
#+taxon_of_interest:poly(log_bound_dist_km_PA_neg,2)
#+Zone:poly(log_bound_dist_km_PA_neg,2)
#+Within_PA+poly(log_hpd,1)+(1+Within_PA+log_access+log_hpd|SS)+(1|Predominant_habitat)"
 





# then do PA characteristics seperately???


fF <- c("Zone", "taxon_of_interest", "IUCN.PA") 
fT <- list("ag_suit" = "3", "log_hpd" = "3", "log_access" = "3", "log_bound_dist_km_PA_neg" = "3", "DoP.PA" = "3", "log_AREA.PA" = "3")
fI <- character(0)
RS <-  character(0)

# range~IUCN.PA+poly(ag_suit,3)+poly(log_access,2)+poly(log_AREA.PA,3)+poly(log_hpd,3)+taxon_of_interest+Zone+(1|SS)+(1|Predominant_habitat)"
# but when plotted, IUCN CAT and AREA show no real difference??? this doesnt seem right







######
### put in everything apart from IUCN_CAT
######



fF <- c("Zone", "taxon_of_interest", "Within_PA") 
fT <- list("ag_suit" = "3", "log_hpd" = "3", "log_access" = "3", "log_bound_dist_km_PA_neg" = "3", "DoP.PA" = "3", "log_AREA.PA" = "3")
fI <- character(0)
RS <-  c("Within_PA")

# 07_14 range~poly(ag_suit,3)+poly(log_access,2)+poly(log_AREA.PA,3)+poly(log_bound_dist_km_PA_neg,3)+poly(log_hpd,3)+taxon_of_interest+Within_PA+Zone+(1+Within_PA|SS)+(1|Predominant_habitat)"
# 11_14 range~poly(ag_suit,3)+poly(DoP.PA,3)+poly(log_access,2)+poly(log_AREA.PA,3)+poly(log_bound_dist_km_PA_neg,1)+poly(log_hpd,3)+taxon_of_interest+Within_PA+Zone+(1+Within_PA|SS)+(1|Predominant_habitat)"




# add interactions
fF <- c("Zone", "taxon_of_interest", "Within_PA") 
fT <- list("ag_suit" = "3", "log_hpd" = "3", "log_access" = "2", "log_bound_dist_km_PA_neg" = "1", "DoP.PA" = "3", "log_AREA.PA" = "3")
fI <- c("Within_PA:poly(ag_suit,3)", "Within_PA:poly(log_hpd,3)", "Within_PA:poly(log_access,2)",
	"Within_PA:taxon_of_interest", 
	"taxon_of_interest:poly(log_bound_dist_km_PA_neg,1)","taxon_of_interest:poly(DoP.PA,3)", "taxon_of_interest:poly(log_AREA.PA,3)",
	"Within_PA:Zone",
	"Zone:poly(log_bound_dist_km_PA_neg,1)","Zone:poly(DoP.PA,3)", "Zone:poly(log_AREA.PA,3)")
RS <-  c("Within_PA")


#07_11
#"range~poly(ag_suit,3)+poly(log_access,2)+poly(log_AREA.PA,3)+poly(log_bound_dist_km_PA_neg,3)+poly(log_hpd,3)+taxon_of_interest+Within_PA+Zone
#+Within_PA:poly(ag_suit,3)
#+Within_PA:poly(log_access,2)
#+taxon_of_interest:poly(log_bound_dist_km_PA_neg,3)
#+Zone:poly(log_bound_dist_km_PA_neg,3)
#+Zone:poly(DoP.PA,1)
#+Zone:poly(log_AREA.PA,3)
#+poly(DoP.PA,1)+(1+Within_PA|SS)+(1|Predominant_habitat)"

#11_14
#range~poly(ag_suit,3)+poly(DoP.PA,3)+poly(log_access,2)+poly(log_AREA.PA,3)+poly(log_bound_dist_km_PA_neg,1)+poly(log_hpd,3)+taxon_of_interest+Within_PA+Zone
#+Within_PA:poly(ag_suit,3)
#taxon_of_interest:poly(log_bound_dist_km_PA_neg,1)
#Zone:poly(DoP.PA,3)
#Zone:poly(log_AREA.PA,3)
#+(1+Within_PA|SS)+(1|Predominant_habitat)"
 





range.best.random <- compare_randoms(matched.landuse, "range",
				fixedFactors=fF,
                         fixedTerms=fT,
                       fixedInteractions=fI,
                         otherRandoms=c("Predominant_habitat"),
				fixed_RandomSlopes = RS,
                          fitInteractions=FALSE,
				verbose=TRUE)


range.best.random$best.random
 


# model select
range.model <- model_select(all.data  = matched.landuse, 
			     responseVar = "range", 
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
                       fixedInteractions=fI,
                       randomStruct = range.best.random$best.random,
			     otherRandoms=c("Predominant_habitat"),
                       verbose=TRUE)






validate(range.model$model) #ok
range.model$warnings
range.model$final.call
range.model$stats


summary(range.model$model)


XV.10 <- CrossValidate(range.model$model,10, divFactor = "SS", fitFamily = "gaussian", data= range.model$data, 
	fixedFactors = fF, fixedTerms = fT, fixedInteractions = fI, randomStruct = range.best.random$best.random)
# all significant estimates seem acceptably similar

XV.20 <- CrossValidate(range.model$model,20, divFactor = "SS", fitFamily = "gaussian", data= range.model$data, 
	fixedFactors = fF, fixedTerms = fT, fixedInteractions = fI, randomStruct = range.best.random$best.random)
# all significant estimates seem acceptably similar


XV.50 <- CrossValidate(range.model$model,50, divFactor = "SS", fitFamily = "gaussian", data= range.model$data, 
	fixedFactors = fF, fixedTerms = fT, fixedInteractions = fI, randomStruct = range.best.random$best.random)
# all significant estimates seem acceptably similar




XVr <- CrossValidate(range.model$model,-1, divFactor = "Predominant_habitat", fitFamily = "gaussian", data= range.model$data, 
	fixedFactors = fF, fixedTerms = fT, fixedInteractions = fI, randomStruct = range.best.random$best.random)
# but if you dropped a whole land use it wouldnt be that surprising if the results changed? 

#glmmR.wikidot/faq

write.csv(range.model$stats, "range.model.stats.18.11.2014.csv")








