


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




# test polynomials

fF <- c("Zone", "taxon_of_interest", "Within_PA", "Predominant_habitat") 
fT <- list("ag_suit" = "3", "log_slope" = "3", "log_elevation" = "3")
keepVars <- character(0)
fI <- character(0)
RS <-  c("Within_PA")

# "Richness_rarefied~poly(ag_suit,3)+poly(log_elevation,3)+Predominant_habitat+Zone+(1+Within_PA|SS)+(1|SSBS)+(1|SSB)"


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

Richness_rarefied.poly <- model_select(all.data  = multiple.taxa.PA_11_14, 
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

Richness_rarefied.poly$stats
#sams approach would suggest cubic

data <- Richness_rarefied.poly$data

m0 <- glmer(Richness_rarefied~Predominant_habitat+taxon_of_interest+Zone+Within_PA+ (1+Within_PA|SS)+(1|SSBS)+(1|SSB), data, family = "poisson")
m1 <- glmer(Richness_rarefied~Predominant_habitat+taxon_of_interest+Zone+Within_PA+ poly(log_slope,1)+(1+Within_PA|SS)+(1|SSBS)+(1|SSB), data, family = "poisson")
m2 <- glmer(Richness_rarefied~Predominant_habitat+taxon_of_interest+Zone+Within_PA + poly(log_slope,2)+(1+Within_PA|SS)+(1|SSBS)+(1|SSB), data, family = "poisson")
m3 <- glmer(Richness_rarefied~Predominant_habitat+taxon_of_interest+Zone+Within_PA + poly(log_slope,3)+(1+Within_PA|SS)+(1|SSBS)+(1|SSB), data, family = "poisson")

#looking at what to do if not significant
anova(m1, m0)
anova(m2, m0)
anova(m3, m0)
#this suggests quadratic



# add interactions
fF <- c("Zone", "taxon_of_interest", "Within_PA", "Predominant_habitat") 
fT <-character(0)
keepVars <- list("ag_suit" = "3", "log_elevation" = "3", "log_slope" = "1")
fI <- c("Within_PA:Predominant_habitat")
RS <-  c("Within_PA")
# Richness_rarefied~Predominant_habitat+Zone+poly(ag_suit,3)+poly(log_elevation,3)+poly(log_slope,1)+(Within_PA|SS)+(1|SSBS)+(1|SSB)


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

Richness_rarefied.model$final.call
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
Richness_rarefied.model.LU3 <- model_select(all.data  = multiple.taxa.PA_11_14, 
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







###combine secondary for land use estimate

#make new datasets
PA_11_14_sec <- PA_11_14
PA_11_14_sec$Predominant_habitat <- gsub("Young secondary vegetation", "Secondary vegetation", PA_11_14_sec$Predominant_habitat)
PA_11_14_sec$Predominant_habitat <- gsub("Intermediate secondary vegetation", "Secondary vegetation", PA_11_14_sec$Predominant_habitat)
PA_11_14_sec$Predominant_habitat <- gsub("Mature secondary vegetation", "Secondary vegetation", PA_11_14_sec$Predominant_habitat)

multiple.taxa.PA_11_14_sec <- multiple.taxa.PA_11_14
multiple.taxa.PA_11_14_sec$Predominant_habitat <- gsub("Young secondary vegetation", "Secondary vegetation", multiple.taxa.PA_11_14_sec$Predominant_habitat)
multiple.taxa.PA_11_14_sec$Predominant_habitat <- gsub("Intermediate secondary vegetation", "Secondary vegetation", multiple.taxa.PA_11_14_sec$Predominant_habitat)
multiple.taxa.PA_11_14_sec$Predominant_habitat <- gsub("Mature secondary vegetation", "Secondary vegetation", multiple.taxa.PA_11_14_sec$Predominant_habitat)

#make a factor and set primary as reference
multiple.taxa.PA_11_14_sec$Predominant_habitat <- factor(multiple.taxa.PA_11_14_sec$Predominant_habitat)
multiple.taxa.PA_11_14_sec$Predominant_habitat <- relevel(multiple.taxa.PA_11_14_sec$Predominant_habitat, "Primary Vegetation")

#also need to make new LUPA
multiple.taxa.PA_11_14_sec$LUPA <- factor(paste(multiple.taxa.PA_11_14_sec$PA, multiple.taxa.PA_11_14_sec$Predominant_habitat))



fF <- c("Zone", "taxon_of_interest", "Within_PA", "Predominant_habitat") 
fT <- list("ag_suit" = "3", "log_slope" = "3", "log_elevation" = "3")
keepVars <- character(0)
fI <- character(0)
RS <-  c("Within_PA")
Richness_rarefied.best.random <- compare_randoms(multiple.taxa.PA_11_14_sec, "Richness_rarefied",
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


Richness_rarefied.best.random$best.random #(1+Within_PA|SS)+ (1|SSBS)+ (1|SSB)
 
best.random <- "(Within_PA|SS)+ (1|SSBS)+ (1|SSB)"


Richness_rarefied.poly <- model_select(all.data  = multiple.taxa.PA_11_14_sec, 
			     responseVar = "Richness_rarefied",
			     fitFamily = "poisson",
			     siteRandom = TRUE, 
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
			     keepVars = keepVars,
                       fixedInteractions=fI,
                       randomStruct = Richness_rarefied.best.random$best.random ,
			     otherRandoms=character(0),
                       verbose=TRUE)
Richness_rarefied.poly$final.call
# "Richness_rarefied~poly(ag_suit,3)+poly(log_elevation,3)+Predominant_habitat+Zone+(1+Within_PA|SS)+(1|SSBS)+(1|SSB)"


# add interactions
fF <- c("Zone", "taxon_of_interest", "Within_PA", "Predominant_habitat") 
fT <-character(0)
keepVars <- list("ag_suit" = "3", "log_elevation" = "3")
fI <- c("Within_PA:Predominant_habitat")
RS <-  c("Within_PA")


Richness_rarefied.model <- model_select(all.data  = multiple.taxa.PA_11_14_sec, 
			     responseVar = "Richness_rarefied",
			     fitFamily = "poisson",
			     siteRandom = TRUE, 
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
			     keepVars = keepVars,
                       fixedInteractions=fI,
                       randomStruct = Richness_rarefied.best.random$best.random ,
			     otherRandoms=character(0),
                       verbose=TRUE)
Richness_rarefied.model$final.call
"Richness_rarefied~Predominant_habitat+Zone+poly(ag_suit,3)+poly(log_elevation,3)+(1+Within_PA|SS)+(1|SSBS)+(1|SSB)"


data <- multiple.taxa.PA_11_14_sec[,c("Richness_rarefied", "Within_PA", "Predominant_habitat", "LUPA",
	"log_slope", "log_elevation", "ag_suit",
	"Zone", "taxon_of_interest", "SS", "SSB", "SSBS")]
data <- na.omit(data)


#to plot need LUPA model as easy way to get standard errors for each landuse
# LUPA is the same as LU + PA + LU:PA
# LUPA vs just LU and PA

LUPA.sp3 <- glmer(Richness_rarefied~taxon_of_interest+Zone+LUPA
		+poly(ag_suit,3)+poly(log_elevation,2)+poly(log_slope,1)+(1+PA|SS)+(1|SSBS)+(1|SSB),
		data = data, family = "poisson")
LUPA.sp4 <- glmer(Richness_rarefied~Predominant_habitat+taxon_of_interest+PA+Zone
		+poly(ag_suit,3)+poly(log_elevation,2)+poly(log_slope,1)+(1+PA|SS)+(1|SSBS)+(1|SSB),
		data = data, family = "poisson")
anova(LUPA.sp3, LUPA.sp4)



save.image("N:\\Documents\\PREDICTS\\WDPA analysis\\RData files\\8 landuses\\LUPA - rar rich sec_combined.RData")








