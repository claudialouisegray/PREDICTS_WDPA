

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
source("plotLU.R")
source("plotFactorInteraction.R")

#load data
source("prep_PA_11_14_for_analysis.R")

validate <- function(x) {
  par(mfrow = c(1,2))
  plot(resid(x)~ fitted(x))
  hist(resid(x))
  par(mfrow = c(1,1))
}



construct_call<-function(responseVar,fixedStruct,randomStruct){
  return(paste(responseVar,"~",fixedStruct,"+",randomStruct,sep=""))
}


length(unique(PA_11_14$SSS[which(PA_11_14$Within_PA == "yes")]))
length(unique(PA_11_14$SSS[which(PA_11_14$Within_PA == "no")]))

### age/size analysis

 PA_11_14$taxon_of_interest <- relevel(PA_11_14$taxon_of_interest, "Invertebrates")

xyplot(log_abundance ~ taxon_of_interest|AREA_DoP, PA_11_14)


fF <- c("Zone", "taxon_of_interest", "AREA_DoP") 
fT <- list("ag_suit" = "3", "log_slope" = "3", "log_elevation" = "3")
keepVars <- character(0)
fI <- character(0)
RS <-  character(0)
# "log_abundance~AREA_DoP+poly(log_elevation,1)+taxon_of_interest+Zone+(1|SS)+(1|SSB)+(1|Predominant_habitat)"


# add interactions 
fF <- c("Zone", "taxon_of_interest", "AREA_DoP") 
fT <- list()
keepVars <- list("ag_suit" = "1", "log_slope" = "1", "log_elevation" = "1")
fI <- c("AREA_DoP:taxon_of_interest", "AREA_DoP:Zone")
RS <-  character(0)
#log_abundance~taxon_of_interest+AREA_DoP:taxon_of_interest+AREA_DoP+poly(ag_suit,1)+poly(log_slope,1)+poly(log_elevation,1)+(1|SS)+(1|SSB)+(1|Predominant_habitat)"




# check gamm
gamm.model <- gamm4(log_abundance ~ Zone + taxon_of_interest + ag_suit + s(log_elevation) + s(log_slope)+ Within_PA,
	random = ~(1+Within_PA|SS)+ (1|SSB) + (1|Predominant_habitat), data =  PA_11_14)
anova(gamm.model$gam)
plot(gamm.model$gam)

log_abundance.best.random <- compare_randoms(PA_11_14, 
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
log_abundance.model2 <- model_select(all.data  = PA_11_14, 
			     responseVar = "log_abundance", 
			     keepVars = keepVars,
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
                       fixedInteractions=fI,
                       randomStruct = log_abundance.best.random$best.random,
			     otherRandoms=c("Predominant_habitat"),
                       verbose=TRUE)


write.csv(log_abundance.model2$stats, "N:/Documents/PREDICTS/WDPA analysis/stats tables all/age size bins/ab.model.LUPA.age.size.stats.10.02.2014.csv")

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



tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/02_15/PA_11_14 abundance vs size and age vs taxon.tif",
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


tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/02_15/PA_11_14 abundance vs size and age vs zone.tif",
	width = 25, height =20, units = "cm", pointsize = 12, res = 300)

plotFactorInteraction(model = log_abundance.model2$model,
responseVar = "log_abundance",
data = log_abundance.model2$data,
xvar = "AREA_DoP",
xvar.order = c("small_young", "small_old", "large_young", "large_old"), 
#this must be the levels of the factor in the order to be plotted, not including reference level
intvar = "Zone",
intvar.order = c("Tropical"),
logLink = "e",
xlab = "PA size and age class",
ylab = "Relative abundance %")

text(0.7,90, "Young = 0 - 20 yrs \nOld = 20 - 85 yrs \nSmall = 0 - 400 km2 \nLarge = 400 - 12000km2", 
	adj = 0, cex = 0.8)

dev.off()


table(data$AREA_DoP, data$taxon_of_interest)


