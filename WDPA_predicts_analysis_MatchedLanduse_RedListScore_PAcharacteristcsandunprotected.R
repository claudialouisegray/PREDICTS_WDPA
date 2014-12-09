

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




names(matched.landuse_amp_mam_bir)

#previous analysis attempts
#matched.landuse_amp_mam_bir$RLS <- log10(matched.landuse_amp_mam_bir$CWM_IUCN_Red_List_Score+1)
#matched.landuse_amp_mam_bir$RLS <- matched.landuse_amp_mam_bir$CWM_IUCN_Red_List_Score

#matched.landuse_amp_mam_bir$RLS <- matched.landuse_amp_mam_bir$CWM_endangered


matched.landuse_amp_mam_bir$RLS.y <- cbind(matched.landuse_amp_mam_bir$abundance_VU_EN_CR,
									matched.landuse_amp_mam_bir$abundance_LC_NT)  



# most of the sites the CWM value of the red list score is least concern 
# if LC is set to 0, still doesnt make sense to use zero-inflated poisson model as they are not real zeros

str(matched.landuse_amp_mam_bir)

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



# get only data that model runs on 

data <- matched.landuse_amp_mam_bir[,c("RLS", "Zone", "taxon_of_interest", "ag_suit", "log_access",
	 "log_hpd", "log_bound_dist_km", "log_dist_in", "log_dist_out", "Within_PA", "DoP.PA", "IUCN.PA", "log_AREA.PA", "Predominant_habitat", "SS")]

data <- na.omit(data) #have to get rid of NA to try poly on ag_suit


nrow(data)



##########
# inc PA characteristics
############





### run model, without IUCN_CAT


names(matched.landuse_amp_mam_bir)

#
fF <- c("Zone", "Within_PA") 
fT <- list("ag_suit" = "3", "log_hpd" = "3", "log_access" = "3", "log_bound_dist_km_PA_neg" = "3", "DoP.PA" = "3", "log_AREA.PA" = "3")
fI <- character(0)
RS <-  c("Within_PA")


# with two-column response variable - doesnt converge with polynomials
# try linear for all except dist to edge, works #RLS.y~poly(log_access,1)+Zone+(1+Within_PA|SS)+(1|Predominant_habitat)+(1|SSBS)
# try with dist to edge, DoP.PA and log_AREA.PA cubic - too many convergence errors
# try with dist to edge, DoP.PA and log_AREA.PA quad - too many convergence errors
# with ag suit, log hpd and log_access cubic, log bound dist cubic, DoP and log_AREA linear - #RLS.y~poly(log_access,1)+Zone+(1+Within_PA|SS)+(1|Predominant_habitat)+(1|SSBS)
fF <- c("Zone", "Within_PA") 
fT <- list("ag_suit" = "1", "log_hpd" = "1", "log_access" = "1", "log_bound_dist_km_PA_neg" = "3", "DoP.PA" = "1", "log_AREA.PA" = "1")
fI <- character(0)
RS <-  c("Within_PA")

#RLS.y~poly(log_access,1)+(1+Within_PA|SS)+(1|Predominant_habitat)+(1|SSBS)



# add interactions
#doesnt work
fF <- c("Zone","Within_PA") 
fT <- list("ag_suit" = "1", "log_hpd" = "1", "log_access" = "1", "log_bound_dist_km_PA_neg" = "1", "DoP.PA" = "1", "log_AREA.PA" = "1")
fI <- c("Within_PA:poly(ag_suit,1)", "Within_PA:poly(log_hpd,1)", "Within_PA:poly(log_access,1)",
	"Within_PA:Zone",
	"Zone:poly(log_bound_dist_km_PA_neg,1)","Zone:poly(DoP.PA,1)", "Zone:poly(log_AREA.PA,1)")
RS <-  c("Within_PA")

#07_14 RLS.y~poly(log_access,1)+Zone
#+Within_PA:poly(log_hpd,1)
#+Zone:poly(log_bound_dist_km_PA_neg,1)
#+Within_PA+poly(log_hpd,1)+poly(log_bound_dist_km_PA_neg,1)
#+(1+Within_PA|SS)+(1|Predominant_habitat)+(1|SSBS)




# works
fF <- c("Zone","Within_PA") 
fT <- list("ag_suit" = "1", "log_hpd" = "1", "log_access" = "1", "log_bound_dist_km_PA_neg" = "1", "DoP.PA" = "1", "log_AREA.PA" = "1")
fI <- c("Within_PA:poly(ag_suit,1)", "Within_PA:poly(log_hpd,1)", "Within_PA:poly(log_access,1)",
	"Zone:poly(log_bound_dist_km_PA_neg,1)", "Zone:poly(DoP.PA,1)")
RS <-  c("Within_PA")
# RLS.y~poly(log_access,1)+Within_PA:poly(log_hpd,1)+Within_PA+poly(log_hpd,1)+(1+Within_PA|SS)+(1|Predominant_habitat)+(1|SSBS)

# problem interactions to add
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


RLS.best.random$best.random
 


# model select
RLS.model <- model_select(all.data  = matched.landuse_amp_mam_bir, 
				fitFamily = "binomial",
				siteRandom = T,
			     responseVar = "RLS.y", 
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
                       fixedInteractions=fI,
                       randomStruct = "(1+Within_PA|SS)+(1|Predominant_habitat) + (1|SSBS)",
			     otherRandoms=c("Predominant_habitat"),
                       verbose=TRUE)








validate(RLS.model$model) 
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

write.csv(RLS.model$stats, "RLS.model.stats.combined.18.11.2014.csv")








