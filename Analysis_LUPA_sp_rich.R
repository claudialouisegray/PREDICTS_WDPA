


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




# Species richness

fF <- c("Zone", "taxon_of_interest", "Within_PA", "Predominant_habitat") 
fT <- list("ag_suit" = "3", "log_slope" = "3", "log_elevation" = "3")
keepVars <- character(0)
fI <- character(0)
RS <-  c("Within_PA")

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


Species_richness.best.random$best.random #"(Within_PA|SS)+ (1|SSBS)+ (1|SSB)"
 
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

# "Species_richness~poly(ag_suit,1)+poly(log_elevation,2)+Predominant_habitat+taxon_of_interest+Within_PA+Zone+(Within_PA|SS)+(1|SSBS)+(1|SSB)"



# add interactions
fF <- c("Zone", "taxon_of_interest", "Within_PA", "Predominant_habitat") 
fT <-character(0)
keepVars <- list("ag_suit" = "1", "log_elevation" = "2", "log_slope" = "1")
fI <- c("Within_PA:Predominant_habitat")
RS <-  c("Within_PA")
#Species_richness~Predominant_habitat+taxon_of_interest+Within_PA+Zone+Within_PA:Predominant_habitat+poly(ag_suit,1)+poly(log_elevation,2)+poly(log_slope,1)+(Within_PA|SS)+(1|SSBS)+(1|SSB)"

#other ints I could add
#fI <- c("Within_PA:poly(ag_suit,2)", "Within_PA:poly(log_elevation,2)", "Within_PA:poly(log_slope,1)",
#	"Within_PA:taxon_of_interest", "Within_PA:Zone")


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
multiple.taxa.PA_11_14_sec$LUPA <- factor(paste(multiple.taxa.PA_11_14_sec$Within_PA, multiple.taxa.PA_11_14_sec$Predominant_habitat))

### test polynomials
fF <- c("Zone", "taxon_of_interest", "Within_PA", "Predominant_habitat") 
fT <- list("ag_suit" = "3", "log_slope" = "3", "log_elevation" = "3")
keepVars <- character(0)
fI <- character(0)
RS <-  c("Within_PA")
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

Species_richness.poly <- model_select(all.data  = multiple.taxa.PA_11_14_sec, 
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
Species_richness.poly$final.call
# "Species_richness~poly(ag_suit,3)+poly(log_elevation,2)+Predominant_habitat+taxon_of_interest+Within_PA+Zone+(1+Within_PA|SS)+(1|SSBS)+(1|SSB)"

Species_richness.model$stats
summary(Species_richness.model$model)

# add interactions
fF <- c("Zone", "taxon_of_interest", "Within_PA", "Predominant_habitat") 
fT <-character(0)
keepVars <- list("ag_suit" = "3", "log_elevation" = "2")
fI <- c("Within_PA:Predominant_habitat")
RS <-  c("Within_PA")

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

Species_richness.model$final.call


### get model with reference = in for comparison to Sams coefficients
data <- Species_richness.model$data
nrow(data)
data$Within_PA <- relevel(data$Within_PA, "yes")

m <- glmer(Species_richness.model$final.call, data = data, family = "poisson",
	control= glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
summary(m)
vcov(m)


#to plot need LUPA model as easy way to get standard errors for each landuse for plotting
#Tim Newbold, 12/02/15
#In other words, your two models should contain either:
#LU + PA + LU:PA, or
#LUPA

# LUPA is the same as LU + PA + LU:PA
# LUPA vs just LU and PA
data <- multiple.taxa.PA_11_14_sec[,c("Species_richness", "Within_PA", "Predominant_habitat", "LUPA",
	"log_slope", "log_elevation", "ag_suit",
	"Zone", "taxon_of_interest", "SS", "SSB", "SSBS")]
data <- na.omit(data)
nrow(data)

data$LUPA <- relevel(data$LUPA, "no Primary Vegetation")
m2 <- glmer(Species_richness~+taxon_of_interest+Zone
	+LUPA+poly(ag_suit,3)+poly(log_elevation,2)+poly(log_slope,1)
	+(1+Within_PA|SS)+(1|SSBS)+(1|SSB), family = "poisson", data = data)
summary(m2)




# recreate model without orthogonal polynomials so that values can be used in prediction

Species_richness.model$final.call
data <- multiple.taxa.PA_11_14_sec[,c("Species_richness", "Within_PA", "Predominant_habitat", "LUPA",
	"log_slope", "log_elevation", "ag_suit",
	"Zone", "taxon_of_interest", "SS", "SSB", "SSBS")]
data <- na.omit(data)

sp.m <- glmer(Species_richness~Predominant_habitat+Within_PA+Zone+taxon_of_interest
	+ Within_PA:Predominant_habitat
	+ ag_suit + I(ag_suit^2) + I(ag_suit^3)
	+ log_elevation + I(log_elevation^2)
	+ log_slope 
	+ (1+Within_PA|SS)+(1|SSBS)+(1|SSB), data = data, family = poisson,
	control= glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
summary(sp.m)


#check results the same
fixef(sp.m)
fixef(Species_richness.model$model)







### get effectiveness estimate using proportion land use in each land use type

crop.PAs <- read.csv("N:/Documents/PREDICTS/WDPA analysis/effectiveness estimates/proportion land use/crop_WDPA.txt", header = T)
pasture.PAs <- read.csv("N:/Documents/PREDICTS/WDPA analysis/effectiveness estimates/proportion land use/pasture_WDPA.txt", header = T)
primary.PAs <- read.csv("N:/Documents/PREDICTS/WDPA analysis/effectiveness estimates/proportion land use/primary_WDPA.txt", header = T)
secondary.PAs <- read.csv("N:/Documents/PREDICTS/WDPA analysis/effectiveness estimates/proportion land use/secondary_WDPA.txt", header = T)
urban.PAs <- read.csv("N:/Documents/PREDICTS/WDPA analysis/effectiveness estimates/proportion land use/urban_WDPA.txt", header = T)

crop.out <- read.csv("N:/Documents/PREDICTS/WDPA analysis/effectiveness estimates/proportion land use/crop_outside.txt", header = T)
pasture.out <- read.csv("N:/Documents/PREDICTS/WDPA analysis/effectiveness estimates/proportion land use/pasture_outside.txt", header = T)
primary.out <- read.csv("N:/Documents/PREDICTS/WDPA analysis/effectiveness estimates/proportion land use/primary_outside.txt", header = T)
secondary.out <- read.csv("N:/Documents/PREDICTS/WDPA analysis/effectiveness estimates/proportion land use/secondary_outside.txt", header = T)
urban.out <- read.csv("N:/Documents/PREDICTS/WDPA analysis/effectiveness estimates/proportion land use/urban_outside.txt", header = T)

landuses <- c("Primary", "Secondary", "Pasture", "Cropland", "Urban")
LU.props <- expand.grid(landuse = landuses, outside = NA, PAs = NA)

# get proportion land use
# this equals sum of number of cells with each proportion divided by total number of cells
# total number of cells multiplied by 1000 as proportions out of 1000, not 1. 

LU.props[which(landuses == "Cropland"), which(names(LU.props) == "PAs")] <- 
	sum(crop.PAs$COUNT_*crop.PAs$VALUE_)/(sum(crop.PAs$COUNT_)*1000)
LU.props[which(landuses == "Primary"), which(names(LU.props) == "PAs")] <- 
	sum(as.numeric(primary.PAs$COUNT_*primary.PAs$VALUE_))/(sum(primary.PAs$COUNT_)*1000)
LU.props[which(landuses == "Secondary"), which(names(LU.props) == "PAs")] <- 
	sum(as.numeric(secondary.PAs$COUNT_*secondary.PAs$VALUE_))/(sum(secondary.PAs$COUNT_)*1000)
LU.props[which(landuses == "Pasture"), which(names(LU.props) == "PAs")] <- 
	sum(as.numeric(pasture.PAs$COUNT_*pasture.PAs$VALUE_))/(sum(pasture.PAs$COUNT_)*1000)
LU.props[which(landuses == "Urban"), which(names(LU.props) == "PAs")] <- 
	sum(urban.PAs$COUNT_*urban.PAs$VALUE_)/(sum(urban.PAs$COUNT_)*1000)

LU.props[which(landuses == "Cropland"), which(names(LU.props) == "outside")] <- 
	sum(as.numeric(crop.out$COUNT_*crop.out$VALUE_))/(sum(crop.out$COUNT_)*1000)
LU.props[which(landuses == "Primary"), which(names(LU.props) == "outside")] <- 
	sum(as.numeric(primary.out$COUNT_*primary.out$VALUE_))/(sum(primary.out$COUNT_)*1000)
LU.props[which(landuses == "Secondary"), which(names(LU.props) == "outside")] <- 
	sum(as.numeric(secondary.out$COUNT_*secondary.out$VALUE_))/(sum(secondary.out$COUNT_)*1000)
LU.props[which(landuses == "Pasture"), which(names(LU.props) == "outside")] <- 
	sum(as.numeric(pasture.out$COUNT_*pasture.out$VALUE_))/(sum(pasture.out$COUNT_)*1000)
LU.props[which(landuses == "Urban"), which(names(LU.props) == "outside")] <- 
	sum(urban.out$COUNT_*urban.out$VALUE_)/(sum(urban.out$COUNT_)*1000)

LU.props

sec_out <- fixef(Species_richness.model$model)[
	which(names(fixef(Species_richness.model$model)) == "Predominant_habitatSecondary vegetation")]
crop_out <- fixef(Species_richness.model$model)[
	which(names(fixef(Species_richness.model$model)) == "Predominant_habitatCropland")]
pas_out <- fixef(Species_richness.model$model)[
	which(names(fixef(Species_richness.model$model)) == "Predominant_habitatPasture")]
urb_out <- fixef(Species_richness.model$model)[
	which(names(fixef(Species_richness.model$model)) == "Predominant_habitatUrban")]

pri_in <- fixef(Species_richness.model$model)[
	which(names(fixef(Species_richness.model$model)) == "Within_PAyes")]
sec_in <- fixef(Species_richness.model$model)[
	which(names(fixef(Species_richness.model$model)) == "Within_PAyes")] +
	fixef(Species_richness.model$model)[
	which(names(fixef(Species_richness.model$model)) == "Predominant_habitatSecondary vegetation")] +
	fixef(Species_richness.model$model)[
	which(names(fixef(Species_richness.model$model)) == "Predominant_habitatSecondary vegetation:Within_PAyes")]
crop_in <- fixef(Species_richness.model$model)[
	which(names(fixef(Species_richness.model$model)) == "Within_PAyes")] +
	fixef(Species_richness.model$model)[
	which(names(fixef(Species_richness.model$model)) == "Predominant_habitatCropland")] +
	fixef(Species_richness.model$model)[
	which(names(fixef(Species_richness.model$model)) == "Predominant_habitatCropland:Within_PAyes")]
pas_in <- fixef(Species_richness.model$model)[
	which(names(fixef(Species_richness.model$model)) == "Within_PAyes")] +
	fixef(Species_richness.model$model)[
	which(names(fixef(Species_richness.model$model)) == "Predominant_habitatPasture")] +
	fixef(Species_richness.model$model)[
	which(names(fixef(Species_richness.model$model)) == "Predominant_habitatPasture:Within_PAyes")]
urb_in <- fixef(Species_richness.model$model)[
	which(names(fixef(Species_richness.model$model)) == "Within_PAyes")] +
	fixef(Species_richness.model$model)[
	which(names(fixef(Species_richness.model$model)) == "Predominant_habitatUrban")] +
	fixef(Species_richness.model$model)[
	which(names(fixef(Species_richness.model$model)) == "Predominant_habitatUrban:Within_PAyes")]


response_out <- exp(0 + 
	 	sec_out*LU.props[which(landuses == "Secondary"), which(names(LU.props) == "outside")] +
		crop_out*LU.props[which(landuses == "Cropland"), which(names(LU.props) == "outside")] +
		pas_out*LU.props[which(landuses == "Pasture"), which(names(LU.props) == "outside")]  + 
		urb_out*LU.props[which(landuses == "Urban"), which(names(LU.props) == "outside")])

response_in <- exp(pri_in*LU.props[which(landuses == "Primary"), which(names(LU.props) == "PAs")] + 
	 	sec_in*LU.props[which(landuses == "Secondary"), which(names(LU.props) == "PAs")] +
		crop_in*LU.props[which(landuses == "Cropland"), which(names(LU.props) == "PAs")] +
		pas_in*LU.props[which(landuses == "Pasture"), which(names(LU.props) == "PAs")]  + 
		urb_in*LU.props[which(landuses == "Urban"), which(names(LU.props) == "PAs")])


b.sp <-  response_in/response_out -1



benefit <- as.numeric(b.sp)		# percentage increase in metric in PAs
PA.pct <- 14.6 				# percentage of total land area in PAs
# global loss of biodiversity (from Newbold et al)
global.loss <- 0.136


NPA.rel <- 1-benefit			# relative biodiversity in unprotected sites
PA.rel <- 1					# biodiversity in protected sites
NPA.pct <- 100-PA.pct			# land area unprotected 
global.int <- 1-global.loss		# overall status of biodiversity relative to pristine


# we want NPA.abs and PA.abs - where these are the biodiv metrics in unprotected and protected relative to pristine
# simultaneous equations are

# global.int = NPA.pct/100*NPA.abs + PA.pct/100*PA.abs #overall loss is loss in PAs and nonPAs relative to pristine
# NPA.abs = PA.abs*NPA.rel			     


# so
#global.int = (NPA.pct/100)*PA.abs*NPA.rel + (PA.pct/100)*PA.abs
# which is the same as
#global.int = PA.abs*(NPA.pct/100*(NPA.rel) + PA.pct/100)

# then
PA.abs <- global.int/(PA.pct/100 + (NPA.pct/100 * (NPA.rel)))
NPA.abs <- PA.abs*NPA.rel

# if pristine is 1, then where between 1 and NPA.abs does PA.abs fall?
# get difference between PA.abs and NPA.abs as a percentage of NPA.abs
#((1-NPA.abs)-(1-PA.abs))/(1-NPA.abs)

est <- 1-(1-PA.abs)/(1-NPA.abs)

# ie
#est <- ((1-NPA.abs)-(1-PA.abs))/(1-NPA.abs)

est



