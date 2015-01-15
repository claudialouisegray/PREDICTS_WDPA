

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




names(matched.landuse_amp_mam_bir)

#previous analysis attempts that do not work
#matched.landuse_amp_mam_bir$RLS <- log10(matched.landuse_amp_mam_bir$CWM_IUCN_Red_List_Score+1)
#matched.landuse_amp_mam_bir$RLS <- matched.landuse_amp_mam_bir$CWM_IUCN_Red_List_Score
#matched.landuse_amp_mam_bir$RLS <- matched.landuse_amp_mam_bir$CWM_endangered


matched.landuse_amp_mam_bir$RLS.y <- cbind(matched.landuse_amp_mam_bir$abundance_VU_EN_CR,
									matched.landuse_amp_mam_bir$abundance_LC_NT)  



# most of the sites the CWM value of the red list score is least concern 
# if LC is set to 0, still doesnt make sense to use zero-inflated poisson model as they are not real zeros



######
RLS
######

title(main = "all matched data")
hist(matched.landuse_amp_mam_bir$RLS)

hist(log(matched.landuse_amp_mam_bir$RLS+1))


plot(RLS ~ Within_PA, matched.landuse_amp_mam_bir)
plot(RLS ~ IUCN.PA, matched.landuse_amp_mam_bir)
plot(RLS ~ Zone, matched.landuse_amp_mam_bir)
plot(RLS ~ elevation, matched.landuse_amp_mam_bir)
plot(RLS ~ log_dist_in, matched.landuse_amp_mam_bir)
plot(RLS ~ log_dist_out, matched.landuse_amp_mam_bir)
plot(RLS ~ log_bound_dist_km_PA_neg, matched.landuse_amp_mam_bir)
plot(RLS ~ bound_dist_km_PA_neg, matched.landuse_amp_mam_bir)


plot(RLS ~ DoP.PA, matched.landuse_amp_mam_bir)


xyplot(RLS ~ Within_PA|Zone, matched.landuse_amp_mam_bir)
xyplot(RLS ~ Within_PA|Realm, matched.landuse_amp_mam_bir)
xyplot(RLS ~ Within_PA|taxon_of_interest, matched.landuse_amp_mam_bir)

xyplot(RLS ~ log_bound_dist_km_PA_neg|Zone, matched.landuse_amp_mam_bir)
xyplot(RLS ~ log_bound_dist_km_PA_neg|taxon_of_interest, matched.landuse_amp_mam_bir)
xyplot(RLS ~ log_bound_dist_km_PA_neg|ag_suit, matched.landuse_amp_mam_bir) # could be tricky to fit, few points in ag_suit = 1
plot3d( matched.landuse_amp_mam_bir$log_bound_dist_km_PA_neg, matched.landuse_amp_mam_bir$log_access, matched.landuse_amp_mam_bir$RLS) #reasonable
plot3d( matched.landuse_amp_mam_bir$log_bound_dist_km_PA_neg, matched.landuse_amp_mam_bir$log_hpd, matched.landuse_amp_mam_bir$RLS)	#reasonable

xyplot(RLS ~ hpd|Within_PA, matched.landuse_amp_mam_bir)

xyplot(RLS ~ elevation|Zone, matched.landuse_amp_mam_bir)
xyplot(RLS ~ access|Zone, matched.landuse_amp_mam_bir)
xyplot(RLS ~ ag_suit|Zone, matched.landuse_amp_mam_bir)
xyplot(RLS ~ hpd|Zone, matched.landuse_amp_mam_bir)
xyplot(RLS ~ slope|Zone, matched.landuse_amp_mam_bir)

hist(matched.landuse_amp_mam_bir$elevation)
hist(matched.landuse_amp_mam_bir$ag_suit)
hist(matched.landuse_amp_mam_bir$hpd)
hist(matched.landuse_amp_mam_bir$DoP)
hist(log(matched.landuse_amp_mam_bir$elevation+1))
hist(log(matched.landuse_amp_mam_bir$hpd+1))
hist(log(PA.data$DoP+1))

plot(RLS ~ elevation, matched.landuse_amp_mam_bir.s)
plot(RLS ~ log(elevation+1), matched.landuse_amp_mam_bir.s)




### distance to boundary model

fF <- c("Zone") 
fT <- list("ag_suit" = "3", "log_slope" = "3", "log_elevation" = "3", "log_bound_dist_km_PA_neg" = "3", "DoP.PA" = "3", "log_AREA.PA" = "3")
keepVars <- character(0)
fI <- character(0)
RS <-  c("log_bound_dist_km_PA_neg")

# without block: RLS.y~poly(ag_suit,2)+poly(DoP.PA,1)+poly(log_AREA.PA,1)+poly(log_bound_dist_km_PA_neg,2)+poly(log_elevation,3)+(1+log_bound_dist_km_PA_neg|SS)+(1|SSBS)+(1|Predominant_habitat)
# with block: "RLS.y~poly(ag_suit,1)+poly(DoP.PA,1)+poly(log_AREA.PA,1)+poly(log_elevation,3)+(1+log_bound_dist_km_PA_neg|SS)+(1|SSBS)+(1|SSB)+(1|Predominant_habitat)

#add interactions
fF <- c("Zone") 
fT <- list( "log_bound_dist_km_PA_neg" = "1", "DoP.PA" = "1", "log_AREA.PA" = "1")
keepVars <- list("ag_suit" = "1", "log_slope" = "1", "log_elevation" = "3")
fI <- c("Zone:poly(log_bound_dist_km_PA_neg,1)")
RS <-  c("log_bound_dist_km_PA_neg")

#without block: RLS.y~poly(ag_suit,2)+poly(DoP.PA,1)+poly(log_AREA.PA,1)+poly(log_bound_dist_km_PA_neg,2)+poly(log_elevation,3)+(1+log_bound_dist_km_PA_neg|SS)+(1|SSBS)+(1|Predominant_habitat)
# with block: "RLS.y~poly(DoP.PA,1)+poly(log_AREA.PA,1)+poly(ag_suit,1)+poly(log_slope,1)+poly(log_elevation,3)+(1+log_bound_dist_km_PA_neg|SS)+(1|SSBS)+(1|SSB)+(1|Predominant_habitat)


RLS.best.random <- compare_randoms(matched.landuse_amp_mam_bir, "RLS.y",
				fitFamily = "binomial",
				siteRandom = T,
				fixedFactors=fF,
                         	fixedTerms=fT,
                       	fixedInteractions=fI,
                        	otherRandoms=c("Predominant_habitat"),
				fixed_RandomSlopes = RS,
                          	fitInteractions=FALSE,
				verbose=TRUE)


RLS.best.random$best.random #"(1+log_bound_dist_km_PA_neg|SS)+ (1|SSBS)+ (1|SSB)"
 
best.random <- "(1+log_bound_dist_km_PA_neg|SS)+ (1|SSBS)+ (1|SSB)+(1|Predominant_habitat)"



# model select
RLS.model <- model_select(all.data  = matched.landuse_amp_mam_bir, 
			     	fitFamily = "binomial",
			     	siteRandom = T,
			     	responseVar = "RLS.y", 
			     	alpha = 0.05,
                       	fixedFactors= fF,
                       	fixedTerms= fT,
			     	keepVars = keepVars,
                       	fixedInteractions=fI,
                       	randomStruct = best.random,
			     	otherRandoms=c("Predominant_habitat"),
                       	verbose=TRUE)


validate(RLS.model$model)
write.csv(RLS.model$stats, "RLS.model.stats.05.01.2015.csv")




### within PA model


fF <- c("Zone", "Within_PA") 
fT <- list("ag_suit" = "1", "log_slope" = "1", "log_elevation" = "1", "DoP.PA" = "1", "log_AREA.PA" = "1")
keepVars <- character(0)
fI <- character(0)
RS <-  c("Within_PA")

### first round doesnt converge with all terms having cubic polynomials
# hence try just with PA variables cubic - most interested in whether non linear effects of these
# doesnt work - all linear
# (previously this couldn't estimate impact of dropping ag suit until DoP is reduced, but subsequently estimated successfully.

# without block: RLS.y~poly(ag_suit,1)+poly(log_elevation,1)+(1+Within_PA|SS)+(1|SSBS)+(1|Predominant_habitat)
# with block: RLS.y~poly(ag_suit,1)+poly(log_elevation,1)+(1+Within_PA|SS)+(1|SSBS)+(1|SSB)+(1|Predominant_habitat)


# add interactions
# doesnt work try removing PA.area interaction, as that didnt work before, neither did within_PA:Zone
fF <- c("Zone","Within_PA") 
fT <- list("ag_suit" = "1", "log_slope" = "1", "log_elevation" = "1", "DoP.PA" = "1", "log_AREA.PA" = "1")
fI <- c("Within_PA:poly(ag_suit,1)", "Within_PA:poly(log_slope,1)", "Within_PA:poly(log_elevation,1)",
	"Within_PA:Zone",
	"Zone:poly(DoP.PA,1)", "Zone:poly(log_AREA.PA,1)")
RS <-  c("Within_PA")

#doesnt work
fF <- c("Zone","Within_PA") 
fT <- list("ag_suit" = "1", "log_slope" = "1", "log_elevation" = "1", "DoP.PA" = "1", "log_AREA.PA" = "1")
fI <- c("Within_PA:poly(ag_suit,1)", "Within_PA:poly(log_slope,1)", "Within_PA:poly(log_elevation,1)",
	"Within_PA:Zone",
	"Zone:poly(DoP.PA,1)")
RS <-  c("Within_PA")

#works
fF <- c("Zone","Within_PA") 
fT <- list("DoP.PA" = "1", "log_AREA.PA" = "1")
keepVars <- list("ag_suit" = "1", "log_slope" = "1", "log_elevation" = "1")
fI <- c("Within_PA:poly(ag_suit,1)", "Within_PA:poly(log_slope,1)", "Within_PA:poly(log_elevation,1)",
	"Zone:poly(DoP.PA,1)")
RS <-  c("Within_PA")
# without block: RLS.y~poly(ag_suit,1)+poly(log_elevation,1)+Within_PA:poly(log_slope,1)+Within_PA+poly(log_slope,1)
#+(1+Within_PA|SS)+(1|SSBS)+(1|Predominant_habitat)

#with block: RLS.y~Within_PA:poly(log_slope,1)+Within_PA+poly(ag_suit,1)+poly(log_slope,1)+poly(log_elevation,1)
#+(1+Within_PA|SS)+(1|SSBS)+(1|SSB)+(1|Predominant_habitat)"

# problem interactions
#"Within_PA:Zone"
#"Zone:poly(log_AREA.PA,1)"




RLS.best.random <- compare_randoms(matched.landuse_amp_mam_bir, "RLS.y",
				fitFamily = "binomial",
				siteRandom = T,
				fixedFactors=fF,
                         fixedTerms=fT,
                       fixedInteractions=fI,
                         otherRandoms=c("Predominant_habitat"),
				fixed_RandomSlopes = RS,
                          fitInteractions=FALSE,
				verbose=TRUE)


RLS.best.random$best.random #"(1+log_bound_dist_km_PA_neg|SS)+ (1|SSBS)"
 
best.random <- "(1+Within_PA|SS)+ (1|SSBS)+ (1|SSB)+(1|Predominant_habitat)"



# model select
RLS.model2 <- model_select(all.data  = matched.landuse_amp_mam_bir, 
			     fitFamily = "binomial",
			     siteRandom = T,
			     responseVar = "RLS.y", 
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
			     keepVars = keepVars,
                       fixedInteractions=fI,
                       randomStruct = best.random,
			     otherRandoms=c("Predominant_habitat"),
                       verbose=TRUE)


validate(RLS.model2$model) 
write.csv(RLS.model2$stats, "RLS.model2.stats.05.01.2015.csv")





RLS.model$warnings
RLS.model$final.call
RLS.model$stats


summary(RLS.model$model)
XV <- CrossValidate(RLS.model$model,50, divFactor = "SS", fitFamily = "gaussian", data= RLS.model$data, 
	fixedFactors = fF, fixedTerms = fT, fixedInteractions = fI, randomStruct = RLS.best.random$best.random)
# all significant estimates seem acceptably similar


XVr <- CrossValidate(RLS.model$model,-1, divFactor = "Predominant_habitat", fitFamily = "gaussian", data= RLS.model$data, 
	fixedFactors = fF, fixedTerms = fT, fixedInteractions = fI, randomStruct = RLS.best.random$best.random)
# but if you dropped a whole land use it wouldnt be that surprising if the results changed? 

#glmmR.wikidot/faq










