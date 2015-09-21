


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
fT <- list("ag_suit" = "3", "log_elevation" = "3", "log_slope" = "1")
fI <- c("Within_PA:Predominant_habitat")
RS <-  c("Within_PA")

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

Richness_rarefied.poly$stats



Richness_rarefied.model$final.call
validate(Richness_rarefied.model$model) #ok
res <- data.frame(est = fixef(Richness_rarefied.model$model), se = se.fixef(Richness_rarefied.model$model))



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




### LUPA with 3 landuses
### to compare to 8 land uses above


# add interactions
fF <- c("Zone", "taxon_of_interest", "Within_PA", "LU_3") 
fT <- list("ag_suit" = "3", "log_elevation" = "3", "log_slope" = "1")
fI <- c("Within_PA:Predominant_habitat")
RS <-  c("Within_PA")
# Richness_rarefied~Predominant_habitat+Zone+poly(ag_suit,3)+poly(log_elevation,3)+poly(log_slope,1)+(Within_PA|SS)+(1|SSBS)+(1|SSB)

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


# LUPA TAX




fF <- c("taxon_of_interest", "Within_PA", "LU_3") 
fT <- list("ag_suit" = "3", "log_elevation" = "3", "log_slope" = "3")
fI <- c("Within_PA:LU_3", "taxon_of_interest:LU_3", "Within_PA:taxon_of_interest", "Within_PA:taxon_of_interest:LU_3")
RS <-  c("Within_PA")


Richness_rarefied.LUPAtax.best.random <- compare_randoms(multiple.taxa.PA_11_14, "Richness_rarefied",
				fitFamily = "poisson",
				siteRandom = TRUE,
				fixedFactors=fF,
                        fixedTerms=fT,			   
                       	fixedInteractions=fI,
                        otherRandoms=character(0),
				fixed_RandomSlopes = RS,
                        fitInteractions=FALSE,
				verbose=TRUE)

Richness_rarefied.LUPAtax.best.random$best.random #"(1+Within_PA|SS)+ (1|SSBS)+ (1|SSB)"


# model select
Richness_rarefied.LUPAtax.model <- model_select(all.data  = multiple.taxa.PA_11_14, 
			     responseVar = "Richness_rarefied",
			     fitFamily = "poisson",
			     siteRandom = TRUE, 
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
                       fixedInteractions=fI,
                       randomStruct = Richness_rarefied.LUPAtax.best.random$best.random,
			     otherRandoms=character(0),
                       verbose=TRUE)
Richness_rarefied.LUPAtax.model$warnings
Richness_rarefied.LUPAtax.model$stats
Richness_rarefied.LUPAtax.model$final.call
# yes, the df = 0 is there




fF <- c("taxon_of_interest", "LUPA") 
fT <- list("ag_suit" = "3", "log_elevation" = "3", "log_slope" = "3")
fI <- c("LUPA:taxon_of_interest")
RS <-  c("Within_PA")

Richness_rarefied.LUPAtax.model.LUPA <- model_select(all.data  = multiple.taxa.PA_11_14, 
			     responseVar = "Richness_rarefied",
			     fitFamily = "poisson",
			     siteRandom = TRUE, 
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
                       fixedInteractions=fI,
                       randomStruct = Richness_rarefied.LUPAtax.best.random$best.random,
			     otherRandoms=character(0),
                       verbose=TRUE)


#try by hand, taking confounding vars from above
multiple.taxa.PA_11_14$LUPA <- paste(multiple.taxa.PA_11_14$LU_3 , multiple.taxa.PA_11_14$Within_PA, sep = "")
data <- multiple.taxa.PA_11_14[,c("LU_3", "Within_PA", "log_slope", "log_elevation", "ag_suit", "Richness_rarefied", 
	"taxon_of_interest", "Zone", "SS", "SSB", "SSBS", "LUPA")]
data <- na.omit(data)
data$LU_3 <- factor(data$LU_3)
xyplot(Richness_rarefied ~ LU_3|taxon_of_interest, data)

M1 <- glmer(Richness_rarefied ~ LUPA + taxon_of_interest + LUPA:taxon_of_interest +Zone
	+ poly(log_elevation,3) + poly(ag_suit,3)
	+(1+Within_PA|SS)+ (1|SSBS)+ (1|SSB), data, family = "poisson")
M2 <- glmer(Richness_rarefied ~ LUPA + taxon_of_interest +Zone
	+ poly(log_elevation,3) + poly(ag_suit,3)
	+(1+Within_PA|SS)+ (1|SSBS)+ (1|SSB), data, family = "poisson")

M3 <- glmer(Richness_rarefied ~ taxon_of_interest +Zone
	+ poly(log_elevation,3) + poly(ag_suit,3)
	+(1+Within_PA|SS)+ (1|SSBS)+ (1|SSB), data, family = "poisson")
anova(M2, M3)

summary(M1)











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








