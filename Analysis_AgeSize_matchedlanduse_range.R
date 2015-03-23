
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

source("compare_randoms.R")
source("model_select.R")
source("plotFactorInteraction.R")

#load data
source("prep_matched.landuse_for_analysis.R")

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



### age/size analysis

xyplot(range ~ taxon_of_interest|AREA_DoP, matched.landuse)

fF <- c("Zone", "taxon_of_interest", "AREA_DoP") 
fT <- list("ag_suit" = "3", "log_slope" = "3", "log_elevation" = "3")
keepVars <- character(0)
fI <- character(0)
RS <-  character(0)
#range~poly(ag_suit,3)+poly(log_elevation,2)+poly(log_slope,3)+taxon_of_interest+Zone+(1|SS)+(1|SSB)+(1|Predominant_habitat)"

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

range.poly <- model_select(all.data  = matched.landuse, 
			     responseVar = "range", 
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
			     keepVars = keepVars,
                       fixedInteractions=fI,
                       randomStruct = range.best.random2$best.random,
			     otherRandoms=c("Predominant_habitat"),
                       verbose=TRUE)
#range~poly(ag_suit,3)+poly(log_elevation,2)+poly(log_slope,3)+taxon_of_interest+Zone+(1|SS)+(1|SSB)+(1|Predominant_habitat)"



# add interactions
fF <- c("Zone", "taxon_of_interest", "AREA_DoP") 
fT <- list()
keepVars <- list("ag_suit" = "3", "log_slope" = "3", "log_elevation" = "2")
fI <- c("AREA_DoP:taxon_of_interest", "AREA_DoP:Zone")
RS <-  character(0)

# model select
range.model <- model_select(all.data  = matched.landuse, 
			     responseVar = "range", 
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
			     keepVars = keepVars,
                       fixedInteractions=fI,
                       randomStruct = range.best.random2$best.random,
			     otherRandoms=c("Predominant_habitat"),
                       verbose=TRUE)
#range~taxon_of_interest+Zone+AREA_DoP:taxon_of_interest+AREA_DoP+poly(ag_suit,3)+poly(log_slope,3)+poly(log_elevation,2)+(1|SS)+(1|SSB)+(1|Predominant_habitat)"
range.model$stats


write.csv(range.model$stats, "N:/Documents/PREDICTS/WDPA analysis/stats tables all/age size bins/range.model.age.size.stats.20.03.2014.csv")




summary(range.model2$model)

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
	width = 20, height = 15, units = "cm", pointsize = 12, res = 300)

plotFactorInteraction(model = range.model$model,
responseVar = "Endemicity",
data = range.model$data,
xvar = "AREA_DoP",
xvar.order = c("small_young", "small_old", "large_young", "large_old"),  #this must be the levels of the factor in the order to be plotted, not including reference level
xvar.labels = c("Small, Young", "Small, Old", "Large, Young", "Large, Old"),
intvar = "taxon_of_interest",
intvar.order = c("Invertebrates", "Vertebrates"),
logLink = "n",
xlab = "PA size and age class",
ylab = "Relative endemicity %")

#text(0.7,8, "Young = 0 - 20 yrs \nOld = 20 - 85 yrs \nSmall = 0 - 400 km2 \nLarge = 400 - 12000km2", 
#	adj = 0, cex = 0.8)


dev.off()

 save.image("N:\\Documents\\PREDICTS\\WDPA analysis\\RData files\\8 landuses\\age_size range.RData")


