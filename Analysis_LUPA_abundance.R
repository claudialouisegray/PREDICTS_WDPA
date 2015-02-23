


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
source("plotLU.R")

#load data
source("WDPA_predicts_prep_PA_11_14_for_analysis.R")

validate <- function(x) {
  par(mfrow = c(1,2))
  plot(resid(x)~ fitted(x))
  hist(resid(x))
  par(mfrow = c(1,1))
}



construct_call<-function(responseVar,fixedStruct,randomStruct){
  return(paste(responseVar,"~",fixedStruct,"+",randomStruct,sep=""))
}




#check non linear relationships

fF <- c("Zone", "taxon_of_interest", "Within_PA", "Predominant_habitat") 
fT <- list("ag_suit" = "3", "log_slope" = "3", "log_elevation" = "3")
keepVars <- character(0)
fI <- character(0)
RS <-  c("Within_PA")

#log_abundance~Predominant_habitat+taxon_of_interest+Zone+(Within_PA|SS)+(1|SSB)

# add interactions
fF <- c("Zone", "taxon_of_interest", "Within_PA", "Predominant_habitat") 
fT <- character(0)
keepVars <- list("ag_suit" = "1", "log_elevation" = "1", "log_slope" = "1") 
fI <- c("Within_PA:Predominant_habitat")
RS <-  c("Within_PA")
#"log_abundance~taxon_of_interest+Within_PA+Within_PA:Predominant_habitat+Predominant_habitat+poly(ag_suit,1)+poly(log_elevation,1)+poly(log_slope,1)+(Within_PA|SS)+(1|SSB)"



#other interactions to add
#fI <- c("Within_PA:poly(ag_suit,3)", "Within_PA:poly(log_elevation,2)", "Within_PA:poly(log_slope,1)",
#	"Within_PA:taxon_of_interest", "Within_PA:Zone")



log_abundance.best.random <- compare_randoms(PA_11_14, "log_abundance",
				siteRandom = TRUE,
				fixedFactors=fF,
                        fixedTerms=fT,
			     	keepVars = keepVars,
                       	fixedInteractions=fI,
                        otherRandoms=c("Predominant_habitat"),
				fixed_RandomSlopes = RS,
                        fitInteractions=FALSE,
				verbose=TRUE)


log_abundance.best.random$best.random #
 
best.random <- "(Within_PA|SS)+(1|SSB)"


# model select

log_abundance.model <- model_select(all.data  = PA_11_14, 
			     responseVar = "log_abundance",
			     siteRandom = TRUE, 
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
			     keepVars = keepVars,
                       fixedInteractions=fI,
                       randomStruct = best.random,
			     otherRandoms=c("Predominant_habitat"),
                       verbose=TRUE)


validate(log_abundance.model$model) #ok

write.csv(log_abundance.model$stats, "N:/Documents/PREDICTS/WDPA analysis/stats tables all/LUPA/log_abundance.model.LUPA.stats.02.02.2015.csv")


tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/01_15/LUPA_abundance.tif",
	width = 20, height = 12, units = "cm", pointsize = 12, res = 300)

plotLU (responseVar = "Log abundance", 
				xvar = "Predominant_habitat",
				intVar = "Within_PA",
				level2 = "yes",
				model = log_abundance.model$model,
				col.key = NULL,
				logLink = "e",
				seMultiplier=1.96,
				forPaper = FALSE,
				cex.txt = 0.5
				)

dev.off()



multiple.taxa.PA_11_14$Predominant_habitat <- relevel(multiple.taxa.PA_11_14$Predominant_habitat, "Primary Vegetation")
multiple.taxa.PA_11_14$Within_PA <- relevel(multiple.taxa.PA_11_14$Within_PA, "yes")

m1 <- glmer(log_abundance ~ Within_PA + Predominant_habitat +
	Within_PA:Predominant_habitat + log_slope + log_elevation + ag_suit
	+ Zone + taxon_of_interest + 
	(Within_PA|SS) + (1|SSB) + (1|SSBS), family = "poisson", 
#	control= glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000),
	data = multiple.taxa.PA_11_14)

m1 <- glmer(log_abundance ~ Within_PA + Predominant_habitat +
	Within_PA:Predominant_habitat + log_slope + log_elevation + ag_suit
	+ Zone + taxon_of_interest + 
	(Within_PA|SS) + (1|SSB) + (1|SSBS), family = "poisson", 
#	control= glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000),
	data = multiple.taxa.PA_11_14)

#no convergence warnings when cropland is reference though
res <- data.frame(estimate = fixef(m1), se = se.fixef(m1))
write.csv(res, "log_abundance.LUPA.CLG.30.01.2015.csv")

m <- glmer(log_abundance ~ Within_PA + Predominant_habitat +
	Within_PA:Predominant_habitat + log_slope + log_elevation + ag_suit
	+ Zone + taxon_of_interest + 
	(Within_PA|SS) + (1|SSB) + (1|SSBS), family = "poisson", 
	matched.landuse)
m <- glmer(log_abundance ~ Within_PA + Predominant_habitat +
	Within_PA:Predominant_habitat + log_slope + log_elevation + ag_suit
	+ Zone + taxon_of_interest + 
	(Within_PA|SS) + (1|SSB) + (1|SSBS), family = "poisson", 
	matched.landuse)



