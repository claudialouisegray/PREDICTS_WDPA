
rm(list=ls()) 

library(yarg)
library(roquefort)
library(gamm4)

setwd("R:/ecocon_d/clg32/GitHub/PREDICTS_WDPA")
source("compare_randoms.R")
source("model_select.R")
setwd("R:/ecocon_d/clg32/PREDICTS/WDPA analysis")
PREDICTS_WDPA <- read.csv("PREDICTS_WDPA.csv")

validate <- function(x) {
  par(mfrow = c(1,2))
  plot(resid(x)~ fitted(x))
  hist(resid(x))
  par(mfrow = c(1,1))
}

matched.landuse <- subset(PREDICTS_WDPA, matched.landuse == "yes")
multiple.taxa.matched.landuse <- subset(matched.landuse, multiple_taxa == "yes")
nrow(matched.landuse) #5015


### Size and Age analysis

xyplot(Species_richness ~ taxon_of_interest|AREA_DoP, multiple.taxa.matched.landuse)

fF <- c("Zone", "taxon_of_interest", "AREA_DoP") 
fT <- list("ag_suit" = "3", "log_slope" = "3", "log_elevation" = "3")
keepVars <- character(0)
fI <- c("AREA_DoP:taxon_of_interest", "AREA_DoP:Zone")
RS <-  character(0)
Species_richness.best.random <- compare_randoms(multiple.taxa.matched.landuse, "Species_richness",
				fitFamily = "poisson",
				siteRandom = TRUE,
				fixedFactors=fF,
        fixedTerms=fT,
        keepVars = keepVars,
        fixedInteractions=fI,
        otherRandoms=c("Predominant_habitat"),
				fixed_RandomSlopes = RS,
        fitInteractions=FALSE,
				verbose=TRUE)

Species_richness.best.random$best.random #"(1|SS)+ (1|SSBS)+ (1|SSB)+(1|Predominant_habitat)"

Species_richness.model <- model_select(all.data  = multiple.taxa.matched.landuse, 
			     responseVar = "Species_richness",
			     fitFamily = "poisson",
			     siteRandom = TRUE, 
			     alpha = 0.05,
           fixedFactors= fF,
           fixedTerms= fT,
			     keepVars = keepVars,
           fixedInteractions=fI,
           randomStruct =Species_richness.best.random$best.random ,
			     otherRandoms=c("Predominant_habitat"),
           verbose=TRUE)
Species_richness.model$final.call
#Species_richness~AREA_DoP+taxon_of_interest+Zone+poly(ag_suit,3)+poly(log_slope,1)+poly(log_elevation,2)+(1|SS)+(1|SSBS)+(1|SSB)+(1|Predominant_habitat)"

summary(Species_richness.model$model)
Species_richness.model$stats


### plot

tiff( "sp rich vs size and age.tif",
	width = 17.4, height = 7, units = "cm", pointsize = 12, res = 300)

par(mar = c(4,4,0,0), mgp = c(2,0.5,0))

plotFactorInteraction(model = Species_richness.model$model,
responseVar = "Species_richness",
data = Species_richness.model$data,
xvar = "AREA_DoP",
xvar.order = c("small_young", "small_old", "large_young", "large_old"), #this must be the levels of the factor in the order to be plotted, not including reference level
xvar.labels = c("small, young", "small, old", "large, young", "large, old"),
logLink = "e",
xlab = "PA size and age class",
ylab = "Species richness \ndifference (%)", 
ylim = c(-15,21),
ylab.axis = c(-10,0,10,20),
cex = 1.5,
las = 1,
seperate.line = F)

#mtext(side = 1, text = "PA size and age class", line = 6, cex = 1)

#text(1,26, "Young = 0 - 20 yrs \nOld = 20 - 85 yrs \nSmall = 0 - 400 km2 \nLarge = 400 - 12000km2", 
#	adj = 0, cex = 0.8)

dev.off()


# test interaction

m1 <- glmer(Species_richness  ~ taxon_of_interest + Zone + AREA_DoP + 
		(1|SS)+ (1|SSBS)+ (1|SSB)+(1|Predominant_habitat), multiple.taxa.matched.landuse,
		family = "poisson")
m2 <- glmer(Species_richness  ~ taxon_of_interest + Zone + AREA.PA.f + 
		(1|SS)+ (1|SSBS)+ (1|SSB)+(1|Predominant_habitat), multiple.taxa.matched.landuse,
		family = "poisson")
m3 <- glmer(Species_richness  ~ taxon_of_interest + Zone + DoP.PA.f +
		(1|SS)+ (1|SSBS)+ (1|SSB)+(1|Predominant_habitat), multiple.taxa.matched.landuse,
		family = "poisson")
m4 <- glmer(Species_richness  ~ taxon_of_interest + Zone + 
		(1|SS)+ (1|SSBS)+ (1|SSB)+(1|Predominant_habitat), multiple.taxa.matched.landuse,
		family = "poisson")
anova(m1, m2)#interaction term better than model with just area
anova(m1, m3)#not better than model with just age
anova(m2, m4) #no effect of area alone
anova(m3, m4) #effect of age alone

