


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




# Species richness

fF <- c("Zone", "taxon_of_interest", "Within_PA", "Predominant_habitat") 
fT <- list("ag_suit" = "3", "log_slope" = "3", "log_elevation" = "3")
keepVars <- character(0)
fI <- character(0)
RS <-  c("Within_PA")

# "Species_richness~poly(ag_suit,1)+poly(log_elevation,2)+Predominant_habitat+taxon_of_interest+Within_PA+Zone+(Within_PA|SS)+(1|SSBS)+(1|SSB)"


# add interactions
fF <- c("Zone", "taxon_of_interest", "Within_PA", "Predominant_habitat") 
fT <-character(0)
keepVars <- list("ag_suit" = "1", "log_elevation" = "2", "log_slope" = "1")
fI <- c("Within_PA:Predominant_habitat")
RS <-  c("Within_PA")
#Species_richness~Predominant_habitat+taxon_of_interest+Within_PA+Zone+Within_PA:Predominant_habitat+poly(ag_suit,1)+poly(log_elevation,2)+poly(log_slope,1)+(Within_PA|SS)+(1|SSBS)+(1|SSB)"



fI <- c("Within_PA:poly(ag_suit,2)", "Within_PA:poly(log_elevation,2)", "Within_PA:poly(log_slope,1)",
	"Within_PA:taxon_of_interest", "Within_PA:Zone")



Species_richness.best.random <- compare_randoms(multiple.taxa.PA_11_14, "Species_richness",
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


Species_richness.best.random$best.random #
 
best.random <- "(Within_PA|SS)+ (1|SSBS)+ (1|SSB)"


# model select

Species_richness.model <- model_select(all.data  = multiple.taxa.PA_11_14, 
			     responseVar = "Species_richness",
			     fitFamily = "poisson",
			     siteRandom = TRUE, 
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
			     keepVars = keepVars,
                       fixedInteractions=fI,
                       randomStruct = best.random,
			     otherRandoms=character(0),
                       verbose=TRUE)


validate(Species_richness.model$model) #ok
res <- data.frame(est = fixef(Species_richness.model$model), se = se.fixef(Species_richness.model$model))
write.csv(res, "Species_richness.model.LUPA.estimates.02.02.2015.csv")

write.csv(Species_richness.model$stats, "N:/Documents/PREDICTS/WDPA analysis/stats tables all/LUPA/Species_richness.model.LUPA.stats.30.01.2015.csv")



tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/01_15/LUPA_sp_rich.tif",
	width = 20, height = 12, units = "cm", pointsize = 12, res = 300)

plotLU (responseVar = "Species richness", 
				xvar = "Predominant_habitat",
				intVar = "Within_PA",
				level2 = "yes",
				model = Species_richness.model$model,
				col.key = NULL,
				logLink = "e",
				seMultiplier=1.96,
				forPaper = FALSE,
				cex.txt = 0.5
				)

dev.off()


PlotErrBar(model = Species_richness.model$model, data = multiple.taxa.PA_11_14, responseVar = "Species richness",
		logLink = "e", 
		secdAge= T,
		catEffects = c("Within_PA", "Predominant_habitat"))



# try with Within_PA as random slope

multiple.taxa.PA_11_14$Predominant_habitat <- relevel(multiple.taxa.PA_11_14$Predominant_habitat, "Primary Vegetation")
multiple.taxa.PA_11_14$Within_PA <- relevel(multiple.taxa.PA_11_14$Within_PA, "no")


data <- multiple.taxa.PA_11_14[,c("Species_richness", "Within_PA", "Predominant_habitat", "log_slope", "log_elevation", "ag_suit",
	"Zone", "taxon_of_interest", "SS", "SSB", "SSBS")]
data <- na.omit(data)

m1 <- glmer(Species_richness ~ Within_PA + Predominant_habitat +
	Within_PA:Predominant_habitat + poly(log_slope,1) + poly(log_elevation,2) + poly(ag_suit,2)
	+ Zone + taxon_of_interest + 
	(Within_PA|SS) + (1|SSB) + (1|SSBS), family = "poisson", 
#	control= glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000),
	data = data)
#no convergence warnings when no poly specified at all
#Model failed to converge with max|grad| = 0.00148691 (tol = 0.001, component 15)


m2 <- glmer(Species_richness ~ Within_PA + Predominant_habitat +
	poly(log_slope,1) + poly(log_elevation,2) + poly(ag_suit,2)
	+ Zone + taxon_of_interest + 
	(Within_PA|SS) + (1|SSB) + (1|SSBS), family = "poisson", 
#	control= glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000),
	data = data)

anova(m1,m2)

res <- data.frame(estimate = fixef(m1), se = se.fixef(m1))
write.csv(res, "Species_richness.LUPA.CLG.30.01.2015.withpolylinear.csv")





#combine secondary for feeding to land use

#make new dataset
PA_11_14_sec <- PA_11_14
PA_11_14_sec$Predominant_habitat <- gsub("Young secondary vegetation", "Secondary vegetation", PA_11_14_sec$Predominant_habitat)
PA_11_14_sec$Predominant_habitat <- gsub("Intermediate secondary vegetation", "Secondary vegetation", PA_11_14_sec$Predominant_habitat)
PA_11_14_sec$Predominant_habitat <- gsub("Mature secondary vegetation", "Secondary vegetation", PA_11_14_sec$Predominant_habitat)

multiple.taxa.PA_11_14_sec <- multiple.taxa.PA_11_14
multiple.taxa.PA_11_14_sec$Predominant_habitat <- gsub("Young secondary vegetation", "Secondary vegetation", multiple.taxa.PA_11_14_sec$Predominant_habitat)
multiple.taxa.PA_11_14_sec$Predominant_habitat <- gsub("Intermediate secondary vegetation", "Secondary vegetation", multiple.taxa.PA_11_14_sec$Predominant_habitat)
multiple.taxa.PA_11_14_sec$Predominant_habitat <- gsub("Mature secondary vegetation", "Secondary vegetation", multiple.taxa.PA_11_14_sec$Predominant_habitat)

multiple.taxa.PA_11_14_sec$Predominant_habitat <- factor(multiple.taxa.PA_11_14_sec$Predominant_habitat)
multiple.taxa.PA_11_14_sec$Predominant_habitat <- relevel(multiple.taxa.PA_11_14_sec$Predominant_habitat, "Primary Vegetation")



fF <- c("Zone", "taxon_of_interest", "Within_PA", "Predominant_habitat") 
fT <- list("ag_suit" = "3", "log_slope" = "3", "log_elevation" = "3")
keepVars <- character(0)
fI <- character(0)
RS <-  c("Within_PA")

# "Species_richness~poly(ag_suit,3)+poly(log_elevation,2)+Predominant_habitat+taxon_of_interest+Within_PA+Zone+(1+Within_PA|SS)+(1|SSBS)+(1|SSB)"


# add interactions
fF <- c("Zone", "taxon_of_interest", "Within_PA", "Predominant_habitat") 
fT <-character(0)
keepVars <- list("ag_suit" = "3", "log_elevation" = "2", "log_slope" = "1")
fI <- c("Within_PA:Predominant_habitat")
RS <-  c("Within_PA")
#"Species_richness~Predominant_habitat+taxon_of_interest+Within_PA+Zone+Within_PA:Predominant_habitat+poly(ag_suit,3)+poly(log_elevation,2)+poly(log_slope,1)+(1+Within_PA|SS)+(1|SSBS)+(1|SSB)"

Species_richness.best.random <- compare_randoms(multiple.taxa.PA_11_14_sec, "Species_richness",
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


Species_richness.best.random$best.random #(1+Within_PA|SS)+ (1|SSBS)+ (1|SSB)
 
best.random <- "(Within_PA|SS)+ (1|SSBS)+ (1|SSB)"


# model select

Species_richness.model <- model_select(all.data  = multiple.taxa.PA_11_14_sec, 
			     responseVar = "Species_richness",
			     fitFamily = "poisson",
			     siteRandom = TRUE, 
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
			     keepVars = keepVars,
                       fixedInteractions=fI,
                       randomStruct = Species_richness.best.random$best.random ,
			     otherRandoms=character(0),
                       verbose=TRUE)

Species_richness.model$stats
summary(Species_richness.model$model)






