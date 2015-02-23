

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


nrow(matched.landuse) #5491 down to 5015




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






### log dist to boundary analysis

fF <- c("Zone", "taxon_of_interest") 
fT <- list("ag_suit" = "3", "log_slope" = "3", "log_elevation" = "3", "log_bound_dist_km_PA_neg" = "3", "DoP.PA" = "3", "log_AREA.PA" = "3")
keepVars <- character(0)
fI <- character(0)
RS <-  c("log_bound_dist_km_PA_neg")

#without block: range~poly(ag_suit,3)+poly(log_elevation,2)+poly(log_slope,3)+taxon_of_interest+Zone+(1+log_bound_dist_km_PA_neg|SS)+(1|Predominant_habitat)
#with block:  range~poly(ag_suit,3)+poly(log_elevation,2)+poly(log_slope,3)+taxon_of_interest+Zone+(1+log_bound_dist_km_PA_neg|SS)+(1|SSB)+(1|Predominant_habitat)
# 8 landuses: "range~poly(ag_suit,3)+poly(log_bound_dist_km_PA_neg,2)+poly(log_elevation,2)+poly(log_slope,3)+taxon_of_interest+Zone+(1+log_bound_dist_km_PA_neg|SS)+(1|SSB)+(1|Predominant_habitat)


#check gamm
#doesnt converge
gamm.model <- gamm4(range ~ Zone + taxon_of_interest # + ag_suit + s(log_elevation) + s(log_slope),
	+ s(log_bound_dist_km_PA_neg),
 	#+ s(DoP.PA) + s(log_AREA.PA), 
	random = ~(1+log_bound_dist_km_PA_neg|SS)+ (1|SSB) + (1|Predominant_habitat), data =  matched.landuse)
anova(gamm.model$gam)
plot(gamm.model$gam)


#add interactions
# non linear dist to boundary relationship not strongly supported, keep linear

fF <- c("Zone", "taxon_of_interest") 
fT <- list("log_bound_dist_km_PA_neg" = "1", "DoP.PA" = "1", "log_AREA.PA" = "1")
keepVars <- c("ag_suit" = "3", "log_slope" = "3", "log_elevation" = "2")
fI <- c("Zone:poly(log_bound_dist_km_PA_neg,1)", "taxon_of_interest:poly(log_bound_dist_km_PA_neg,1)") 
RS <-  c("log_bound_dist_km_PA_neg")

#without block: range~poly(ag_suit,3)+poly(log_elevation,2)+poly(log_slope,3)+taxon_of_interest+Zone+(1+log_bound_dist_km_PA_neg|SS)+(1|Predominant_habitat)
# with block: "range~taxon_of_interest+Zone+poly(ag_suit,3)+poly(log_slope,3)+poly(log_elevation,2)+(1+log_bound_dist_km_PA_neg|SS)+(1|SSB)+(1|Predominant_habitat)
# 8 landuses: range~taxon_of_interest+Zone+poly(ag_suit,3)+poly(log_slope,3)+poly(log_elevation,2)+(1+log_bound_dist_km_PA_neg|SS)+(1|SSB)+(1|Predominant_habitat)

range.best.random <- compare_randoms(matched.landuse, "range",
				fixedFactors=fF,
                         fixedTerms=fT,
                       fixedInteractions=fI,
                         otherRandoms=c("Predominant_habitat"),
				fixed_RandomSlopes = RS,
                          fitInteractions=FALSE,
				verbose=TRUE)


range.best.random$best.random  #"(1+log_bound_dist_km_PA_neg|SS)+ (1|SSB)+(1|Predominant_habitat)"
 


# model select
range.model <- model_select(all.data  = matched.landuse, 
			     responseVar = "range", 
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
			     keepVars = keepVars,
                       fixedInteractions=fI,
                       randomStruct = range.best.random$best.random,
			     otherRandoms=c("Predominant_habitat"),
                       verbose=TRUE)


write.csv(range.model$stats, "range.model.stats.23.01.2015.csv")




### IUCN cat analysis
# this doesnt converge - return to including it in simple models



### age/size analysis

xyplot(range ~ taxon_of_interest|AREA_DoP, matched.landuse)

fF <- c("Zone", "taxon_of_interest", "AREA_DoP") 
fT <- list("ag_suit" = "3", "log_slope" = "3", "log_elevation" = "3")
keepVars <- character(0)
fI <- character(0)
RS <-  character(0)
#range~poly(ag_suit,3)+poly(log_elevation,2)+poly(log_slope,3)+taxon_of_interest+Zone+(1|SS)+(1|SSB)+(1|Predominant_habitat)"
> 


# add interactions
fF <- c("Zone", "taxon_of_interest", "AREA_DoP") 
fT <- list()
keepVars <- list("ag_suit" = "3", "log_slope" = "3", "log_elevation" = "2")
fI <- c("AREA_DoP:taxon_of_interest", "AREA_DoP:Zone")
RS <-  character(0)
#range~taxon_of_interest+Zone+AREA_DoP:taxon_of_interest+AREA_DoP+poly(ag_suit,3)+poly(log_slope,3)+poly(log_elevation,2)+(1|SS)+(1|SSB)+(1|Predominant_habitat)"




range.best.random2 <- compare_randoms(matched.landuse, "range",
				fixedFactors=fF,
                         fixedTerms=fT,
				keepVars = keepVars,
                       fixedInteractions=fI,
                         otherRandoms=c("Predominant_habitat"),
				fixed_RandomSlopes = RS,
                          fitInteractions=FALSE,
				verbose=TRUE)


range.best.random2$best.random #"(1+Within_PA|SS)+ (1|SSB)+(1|Predominant_habitat)"
 


# model select
range.model2 <- model_select(all.data  = matched.landuse, 
			     responseVar = "range", 
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
			     keepVars = keepVars,
                       fixedInteractions=fI,
                       randomStruct = range.best.random2$best.random,
			     otherRandoms=c("Predominant_habitat"),
                       verbose=TRUE)



write.csv(range.model2$stats, "N:/Documents/PREDICTS/WDPA analysis/stats tables all/age size bins/range.model.age.size.stats.10.02.2014.csv")




summary(range.model2$model)
fT <- list("DoP.PA" = "1", "log_AREA.PA" = "3", "ag_suit" = "1", log_elevation = "1")

XV <- CrossValidate(range.model2$model,10, divFactor = "SS", fitFamily = "gaussian", data= range.model2$data, 
	fixedFactors = fF, fixedTerms = fT, fixedInteractions = fI, randomStruct = "(1+Within_PA|SS)+ (1|SSB) + (1|Predominant_habitat)")
# all significant estimates seem acceptably similar
XV
#isnt what I really want this process but the whole model select, and the upper and lower p values? 

XVr <- CrossValidate(range.model$model,-1, divFactor = "Predominant_habitat", fitFamily = "gaussian", data= range.model$data, 
	fixedFactors = fF, fixedTerms = fT, fixedInteractions = fI, randomStruct = "(1+Within_PA|SS)+ (1|SSB) + (1|Predominant_habitat)")
# but if you dropped a whole land use it wouldnt be that surprising if the results changed? 

#glmmR.wikidot/faq



tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/02_15/range vs size and age vs taxon.tif",
	width = 25, height = 20, units = "cm", pointsize = 12, res = 300)

plotFactorInteraction(model = range.model2$model,
responseVar = "Endemicity",
data = range.model2$data,
xvar = "AREA_DoP",
xvar.order = c("small_young", "small_old", "large_young", "large_old"),  #this must be the levels of the factor in the order to be plotted, not including reference level
intvar = "taxon_of_interest",
intvar.order = c("Invertebrates", "Vertebrates"),
logLink = "n",
xlab = "PA size and age class",
ylab = "Relative endemicity %")

text(0.7,8, "Young = 0 - 20 yrs \nOld = 20 - 85 yrs \nSmall = 0 - 400 km2 \nLarge = 400 - 12000km2", 
	adj = 0, cex = 0.8)


dev.off()






### Within_PA analysis

fF <- c("Zone", "taxon_of_interest", "Within_PA") 
fT <- list("ag_suit" = "3", "log_slope" = "3", "log_elevation" = "3", "DoP.PA" = "1", "log_AREA.PA" = "1")
keepVars <- character(0)
fI <- character(0)
RS <-  c("Within_PA")


# without block: range~poly(ag_suit,3)+poly(DoP.PA,3)+poly(log_AREA.PA,3)+poly(log_elevation,3)+poly(log_slope,3)+taxon_of_interest+Within_PA+Zone+(1+Within_PA|SS)+(1|Predominant_habitat)"
# with block: range~poly(ag_suit,3)+poly(log_AREA.PA,3)+poly(log_elevation,2)+poly(log_slope,3)+taxon_of_interest+Within_PA+Zone+(1+Within_PA|SS)+(1|SSB)+(1|Predominant_habitat)
# 8 landuses: range~poly(ag_suit,3)+poly(log_AREA.PA,1)+poly(log_elevation,3)+poly(log_slope,3)+taxon_of_interest+Within_PA+Zone+(1+Within_PA|SS)+(1|SSB)+(1|Predominant_habitat)
# compare to gamm


gamm.model <- gamm4(range ~ Within_PA + Zone + taxon_of_interest + ag_suit + s(log_elevation) + s(log_slope)+ s(DoP.PA) + s(log_AREA.PA), 
	random = ~(1+ Within_PA|SS)+ (1|SSB) + (1|Predominant_habitat), data =  matched.landuse)
plot(gamm.model$gam)
summary(gamm.model$gam)


# add interactions
fF <- c("Zone", "taxon_of_interest", "Within_PA") 
fT <- list("DoP.PA" = "1", "log_AREA.PA" = "1")
keepVars <- list("ag_suit" = "3", "log_slope" = "3", "log_elevation" = "3")
fI <- c("Within_PA:poly(ag_suit,3)", "Within_PA:poly(log_slope,3)", "Within_PA:poly(log_elevation,3)",
	"Within_PA:taxon_of_interest", 
	"taxon_of_interest:poly(DoP.PA,1)", "taxon_of_interest:poly(log_AREA.PA,1)",
	"Within_PA:Zone",
	"Zone:poly(DoP.PA,1)", "Zone:poly(log_AREA.PA,1)")
RS <-  c("Within_PA")

#without block
#range~poly(ag_suit,3)+poly(DoP.PA,3)+poly(log_AREA.PA,3)+poly(log_elevation,3)+poly(log_slope,3)+taxon_of_interest+Within_PA+Zone
#+Within_PA:poly(ag_suit,3)
#+Within_PA:poly(log_elevation,3)
#+Zone:poly(log_AREA.PA,3)
#+(1+Within_PA|SS)+(1|Predominant_habitat)"
 
#with block:
#range~poly(log_AREA.PA,3)+taxon_of_interest+Within_PA+Zone
#+Within_PA:poly(ag_suit,3)
#+Within_PA:poly(log_slope,3)
#+Zone:poly(DoP.PA,1)
#+Zone:poly(log_AREA.PA,3)
#+poly(DoP.PA,1)+poly(ag_suit,3)+poly(log_slope,3)+poly(log_elevation,2)
#+(1+Within_PA|SS)+(1|SSB)+(1|Predominant_habitat)"


#8landuse:
#range~poly(log_AREA.PA,1)+taxon_of_interest+Within_PA+Zone
#+Within_PA:poly(ag_suit,3)
#+Within_PA:poly(log_slope,3)
#+Within_PA:taxon_of_interest
#+taxon_of_interest:poly(DoP.PA,1)
#+taxon_of_interest:poly(log_AREA.PA,1)
#+Zone:poly(DoP.PA,1)
#+Zone:poly(log_AREA.PA,1)
#+poly(DoP.PA,1)+poly(ag_suit,3)+poly(log_slope,3)+poly(log_elevation,3)+(1+Within_PA|SS)+(1|SSB)+(1|Predominant_habitat)"





# 8 landuse, limited to quadratic
# not that different, as p value is 0.0489 for taxon:PA_area.
"range~poly(log_AREA.PA,1)+taxon_of_interest+Within_PA+Zone
#+Within_PA:poly(ag_suit,2)
#+Within_PA:poly(log_slope,2)
#+Within_PA:taxon_of_interest
#+taxon_of_interest:poly(DoP.PA,1)
#+Zone:poly(DoP.PA,1)
#+Zone:poly(log_AREA.PA,1)+poly(DoP.PA,1)+poly(ag_suit,2)+poly(log_slope,2)+poly(log_elevation,2)
#+(1+Within_PA|SS)+(1|SSB)+(1|Predominant_habitat)

#same results with other optimizer too



range.best.random2 <- compare_randoms(matched.landuse, "range",
				fixedFactors=fF,
                         fixedTerms=fT,
				keepVars = keepVars,
                       fixedInteractions=fI,
                         otherRandoms=c("Predominant_habitat"),
				fixed_RandomSlopes = RS,
                          fitInteractions=FALSE,
				verbose=TRUE)


range.best.random2$best.random #"(1+Within_PA|SS)+ (1|SSB)+(1|Predominant_habitat)"
 


# model select
range.model2 <- model_select(all.data  = matched.landuse, 
			     responseVar = "range", 
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
			     keepVars = keepVars,
                       fixedInteractions=fI,
                       randomStruct = range.best.random$best.random,
			     otherRandoms=c("Predominant_habitat"),
                       verbose=TRUE)

range.model2_NM <- model_select(all.data  = matched.landuse, 
			     responseVar = "range", 
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
			     keepVars = keepVars,
                       fixedInteractions=fI,
                       randomStruct = range.best.random$best.random,
			     otherRandoms=c("Predominant_habitat"),
			     optimizer = "Nelder_Mead",
                       verbose=TRUE)

# check if result is the same with ag suit and elevation and slope quadratic
range.model2.quad <- model_select(all.data  = matched.landuse, 
			     responseVar = "range", 
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
			     keepVars = keepVars,
                       fixedInteractions=fI,
                       randomStruct = range.best.random$best.random,
			     otherRandoms=c("Predominant_habitat"),
                       verbose=TRUE)





validate(range.model2$model) #ok
write.csv(range.model2$stats, "range.model2.stats.23.01.2015.csv")



range.model$warnings
range.model$final.call
range.model$stats


summary(range.model$model)


XV.10 <- CrossValidate(range.model2$model,10, divFactor = "SS", fitFamily = "gaussian", data= range.model2$data, 
	fixedFactors = fF, fixedTerms = fT, fixedInteractions = fI, randomStruct = range.best.random$best.random)
# all significant estimates seem acceptably similar

XV.20 <- CrossValidate(range.model2$model,20, divFactor = "SS", fitFamily = "gaussian", data= range.model2$data, 
	fixedFactors = fF, fixedTerms = fT, fixedInteractions = fI, randomStruct = range.best.random$best.random)
# all significant estimates seem acceptably similar


XV.50 <- CrossValidate(range.model$model,50, divFactor = "SS", fitFamily = "gaussian", data= range.model$data, 
	fixedFactors = fF, fixedTerms = fT, fixedInteractions = fI, randomStruct = range.best.random$best.random)
# all significant estimates seem acceptably similar




XVr <- CrossValidate(range.model$model,-1, divFactor = "Predominant_habitat", fitFamily = "gaussian", data= range.model$data, 
	fixedFactors = fF, fixedTerms = fT, fixedInteractions = fI, randomStruct = range.best.random$best.random)
# but if you dropped a whole land use it wouldnt be that surprising if the results changed? 

#glmmR.wikidot/faq









