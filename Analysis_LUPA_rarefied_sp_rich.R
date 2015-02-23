


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
source("plotFactorInteraction.R")

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




# test polynomials

fF <- c("Zone", "taxon_of_interest", "Within_PA", "Predominant_habitat") 
fT <- list("ag_suit" = "3", "log_slope" = "3", "log_elevation" = "3")
keepVars <- character(0)
fI <- character(0)
RS <-  c("Within_PA")

# "Richness_rarefied~poly(ag_suit,3)+poly(log_elevation,3)+Predominant_habitat+Zone+(1+Within_PA|SS)+(1|SSBS)+(1|SSB)"


# add interactions
fF <- c("Zone", "taxon_of_interest", "Within_PA", "Predominant_habitat") 
fT <-character(0)
keepVars <- list("ag_suit" = "3", "log_elevation" = "3", "log_slope" = "1")
fI <- c("Within_PA:Predominant_habitat")
RS <-  c("Within_PA")
# Richness_rarefied~Predominant_habitat+Zone+poly(ag_suit,3)+poly(log_elevation,3)+poly(log_slope,1)+(Within_PA|SS)+(1|SSBS)+(1|SSB)



fI <- c("Within_PA:poly(ag_suit,2)", "Within_PA:poly(log_elevation,2)", "Within_PA:poly(log_slope,1)",
	"Within_PA:taxon_of_interest", "Within_PA:Zone")



Richness_rarefied.best.random <- compare_randoms(multiple.taxa.PA_11_14, "Richness_rarefied",
				fitFamily = "poisson",
				siteRandom = TRUE,
				fixedFactors=fF,
                        fixedTerms=fT,
			     	keepVars = keepVars,
                       	fixedInteractions=fI,
                        otherRandoms=character(0),
				fixed_RandomSlopes = RS,
                        fitInteractions=FALSE,
				verbose=TRUE)


Richness_rarefied.best.random$best.random #
 
best.random <- "(Within_PA|SS)+ (1|SSBS)+ (1|SSB)"


# model select

Richness_rarefied.model <- model_select(all.data  = multiple.taxa.PA_11_14, 
			     responseVar = "Richness_rarefied",
			     fitFamily = "poisson",
			     siteRandom = TRUE, 
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
			     keepVars = keepVars,
                       fixedInteractions=fI,
                       randomStruct = Richness_rarefied.best.random$best.random,
			     otherRandoms=character(0),
                       verbose=TRUE)


validate(Richness_rarefied.model$model) #ok
res <- data.frame(est = fixef(Richness_rarefied.model$model), se = se.fixef(Richness_rarefied.model$model))
write.csv(res, "Richness_rarefied.model.LUPA.estimates.30.01.2015.csv")

write.csv(Richness_rarefied.model$stats, "Richness_rarefied.model.LUPA.stats.30.01.2015.csv")

#Plots
#PlotErrBar
PlotErrBar (responseVar = "Rarefied richness", 
			model = Richness_rarefied.model$model,
			data = Richness_rarefied.model$data,
			logLink = "n",
			forPaper = F,
			catEffects = "Predominant_habitat")



#simpler plotting function derived from PlotErrBar
tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/01_15/LUPA_rar_rich.tif",
	width = 20, height = 12, units = "cm", pointsize = 12, res = 300)

plotLU (responseVar = "Rarefied richness", 
				xvar = "Predominant_habitat",
				model = Richness_rarefied.model$model,
				col.key = NULL,
				logLink = "e",
				seMultiplier=1.96,
				forPaper = FALSE,
				cex.txt = 0.5)
dev.off()






### comparison to a model where all confounding variables are linear

data <- multiple.taxa.PA_11_14[,c("Richness_rarefied", "Within_PA", "Predominant_habitat", "log_slope", "log_elevation", "ag_suit",
	"Zone", "taxon_of_interest", "SS", "SSB", "SSBS")]
data <- na.omit(data)

m1 <- glmer(Richness_rarefied ~ Within_PA + Predominant_habitat +
	Within_PA:Predominant_habitat + poly(log_slope,1) + poly(log_elevation,1) + poly(ag_suit,1)
	+ Zone + taxon_of_interest + 
	(Within_PA|SS) + (1|SSB) + (1|SSBS), family = "poisson", 
#	control= glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000),
	data = data)
#no convergence warnings when no poly specified at all
#Model failed to converge with max|grad| = 0.00148691 (tol = 0.001, component 15)


res <- data.frame(estimate = fixef(m1), se = se.fixef(m1))
write.csv(res, "Richness_rarefied.LUPA.CLG.30.01.2015.withpolylinear.csv")





### LUPA with 3 landuses


#test best random and polynomials
fF <- c("Zone", "taxon_of_interest", "Within_PA", "LU_3") 
fT <- list("ag_suit" = "3", "log_slope" = "3", "log_elevation" = "3")
keepVars <- character(0)
fI <- character(0)
RS <-  c("Within_PA")
# 

# add interactions
fF <- c("Zone", "taxon_of_interest", "Within_PA", "LU_3") 
fT <-character(0)
keepVars <- list("ag_suit" = "3", "log_elevation" = "3", "log_slope" = "1")
fI <- c("Within_PA:Predominant_habitat")
RS <-  c("Within_PA")
# Richness_rarefied~Predominant_habitat+Zone+poly(ag_suit,3)+poly(log_elevation,3)+poly(log_slope,1)+(Within_PA|SS)+(1|SSBS)+(1|SSB)



fI <- c("Within_PA:poly(ag_suit,2)", "Within_PA:poly(log_elevation,2)", "Within_PA:poly(log_slope,1)",
	"Within_PA:taxon_of_interest", "Within_PA:Zone")



Richness_rarefied.best.random <- compare_randoms(multiple.taxa.PA_11_14, "Richness_rarefied",
				fitFamily = "poisson",
				siteRandom = TRUE,
				fixedFactors=fF,
                        fixedTerms=fT,
			     	keepVars = keepVars,
                       	fixedInteractions=fI,
                        otherRandoms=character(0),
				fixed_RandomSlopes = RS,
                        fitInteractions=FALSE,
				verbose=TRUE)

Richness_rarefied.best.random$best.random #"(1+Within_PA|SS)+ (1|SSBS)+ (1|SSB)"

# model select
Richness_rarefied.model <- model_select(all.data  = multiple.taxa.PA_11_14, 
			     responseVar = "Richness_rarefied",
			     fitFamily = "poisson",
			     siteRandom = TRUE, 
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
			     keepVars = keepVars,
                       fixedInteractions=fI,
                       randomStruct = Richness_rarefied.best.random$best.random,
			     otherRandoms=character(0),
                       verbose=TRUE)











