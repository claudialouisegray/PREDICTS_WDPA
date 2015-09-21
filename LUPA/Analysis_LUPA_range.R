


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



# first get best random
fF <- c("Zone", "taxon_of_interest", "Within_PA", "Predominant_habitat") 
fT <- list("ag_suit" = "3", "log_slope" = "3", "log_elevation" = "3")
keepVars <- character(0)
fI <- c("Within_PA:Predominant_habitat")
RS <- c("Within_PA")
#"range~Predominant_habitat+taxon_of_interest+Zone+Within_PA:Predominant_habitat+Within_PA+poly(ag_suit,3)+poly(log_elevation,2)+poly(log_slope,3)+(Within_PA|SS)+(1|SSB)"

range.best.random <- compare_randoms(PA_11_14, "range",
				siteRandom = FALSE,
				fixedFactors=fF,
                        fixedTerms=fT,
			     	keepVars = keepVars,
                       	fixedInteractions=fI,
                        otherRandoms=c("Predominant_habitat"),
				fixed_RandomSlopes = RS,
                        fitInteractions=FALSE,
				verbose=TRUE)


range.best.random$best.random 

# model select
range.model <- model_select(all.data  = PA_11_14, 
			     responseVar = "range",
			     siteRandom = TRUE, 
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
			     keepVars = keepVars,
                       fixedInteractions=fI,
                       randomStruct = range.best.random$best.random,
			     otherRandoms=c("Predominant_habitat"),
                       verbose=TRUE)
range.model$final.call

validate(range.model$model) #ok

tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/01_15/LUPA_range.tif",
	width = 20, height = 12, units = "cm", pointsize = 12, res = 300)

plotLU (responseVar = "Endemicity", 
				xvar = "Predominant_habitat",
				intVar = "Within_PA",
				level2 = "yes",
				model = range.model$model,
				col.key = NULL,
				logLink = "n",
				seMultiplier = 1.96,
				forPaper = FALSE,
				cex.txt = 0.5
				)

dev.off()


plotLU (responseVar = "Range", 
				xvar = "Predominant_habitat",
				intVar = "Within_PA",
				level2 = "yes",
				model = range.model$model,
				col.key = NULL,
				logLink = "n",
				seMultiplier = 1.96,
				forPaper = FALSE,
				cex.txt = 0.5
				)



### get R2 values for increasingly complex models
# LU + PA + LU:PA
M1 <- lmer(range ~ Predominant_habitat+Within_PA:Predominant_habitat+Within_PA
	#+poly(ag_suit,3)+poly(log_elevation,2)+poly(log_slope,3)
		+(1+Within_PA|SS)+(1|SSB),
	data = range.model$data, REML = F)

# LU + PA 
M2<- lmer(range ~ Predominant_habitat+Within_PA 
	#+poly(ag_suit,3)+poly(log_elevation,2)+poly(log_slope,3)
	 +(1+Within_PA|SS)+(1|SSB),
	data = range.model$data, REML = F)

# LU
M3 <- lmer(range ~ Predominant_habitat+ 
	#+poly(ag_suit,3)+poly(log_elevation,2)+poly(log_slope,3)
	+(1+Within_PA|SS)+(1|SSB),
	data = range.model$data, REML = F)

# PA 
M4 <- lmer(range ~ Within_PA+ 
	#+poly(ag_suit,3)+poly(log_elevation,2)+poly(log_slope,3)
	+(1+Within_PA|SS)+(1|SSB),
	data = range.model$data, REML = F)

# No land use or PA
M5 <- lmer(range ~ 1
	#+poly(ag_suit,3)+poly(log_elevation,2)+poly(log_slope,3)
	+(1+Within_PA|SS)+(1|SSB),
	data = range.model$data, REML = F)


res <- data.frame(models  =  c("LU + PA + LU:PA", "LU + PA", "LU", "PA", "1"),
		AIC = c(AIC(M1), AIC(M2), AIC(M3), AIC(M4), AIC(M5)),
		r2_conditional = c(R2GLMER(M1)$conditional, R2GLMER(M2)$conditional, R2GLMER(M3)$conditional, R2GLMER(M4)$conditional, R2GLMER(M5)$conditional),
		r2_marginal = c(R2GLMER(M1)$marginal, R2GLMER(M2)$marginal, R2GLMER(M3)$marginal, R2GLMER(M4)$marginal, R2GLMER(M5)$marginal))
res$dAIC = res$AIC - AIC(M5)
res$d_r2_conditional = res$r2_conditional - res$r2_conditional[which(res$models == 1)]
res$d_r2_marginal = res$r2_marginal - res$r2_marginal[which(res$models == 1)]
res <- res[,c(1,2,5,3,6,4,7)]
res

# what percentage of total explanatory power is due to having land use in model only
res$d_r2_marginal[3]/res$d_r2_marginal[1]



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
PA_11_14_sec$Predominant_habitat <- factor(PA_11_14_sec$Predominant_habitat)
PA_11_14_sec$Predominant_habitat <- relevel(PA_11_14_sec$Predominant_habitat, "Primary Vegetation")

#also need to make new LUPA
PA_11_14_sec$LUPA <- factor(paste(PA_11_14_sec$PA, PA_11_14_sec$Predominant_habitat))


#check non linear relationships

fF <- c("Zone", "taxon_of_interest", "Within_PA", "Predominant_habitat") 
fT <- list("ag_suit" = "3", "log_slope" = "3", "log_elevation" = "3")
keepVars <- character(0)
fI <- character(0)
RS <-  c("Within_PA")
#"range~poly(ag_suit,3)+poly(log_elevation,2)+poly(log_slope,3)+Predominant_habitat+taxon_of_interest+Zone+(1+Within_PA|SS)+(1|SSB)"


# add interactions
fF <- c("Zone", "taxon_of_interest", "Within_PA", "Predominant_habitat") 
fT <- character(0)
keepVars <- list("ag_suit" = "3", "log_elevation" = "2", "log_slope" = "3") 
fI <- c("Within_PA:Predominant_habitat")
RS <-  c("Within_PA")
#"range~Predominant_habitat+taxon_of_interest+Zone+Within_PA:Predominant_habitat+Within_PA+poly(ag_suit,3)+poly(log_elevation,2)+poly(log_slope,3)+(1+Within_PA|SS)+(1|SSB)"
 

#other interactions to add
#fI <- c("Within_PA:poly(ag_suit,3)", "Within_PA:poly(log_elevation,2)", "Within_PA:poly(log_slope,1)",
#	"Within_PA:taxon_of_interest", "Within_PA:Zone")



range.best.random <- compare_randoms(PA_11_14_sec, "range",
				fixedFactors=fF,
                        fixedTerms=fT,
			     	keepVars = keepVars,
                       	fixedInteractions=fI,
                        otherRandoms=c("Predominant_habitat"),
				fixed_RandomSlopes = RS,
                        fitInteractions=FALSE,
				verbose=TRUE)


range.best.random$best.random #
 

# model select

range.model <- model_select(all.data  = PA_11_14_sec, 
			     responseVar = "range",
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
			     keepVars = keepVars,
                       fixedInteractions=fI,
                       randomStruct = range.best.random$best.random,
			     otherRandoms=c("Predominant_habitat"),
                       verbose=TRUE)


validate(range.model$model) #ok


# recreate model without orthogonal polynomials so that values can be used in prediction

range.model$final.call
data <- PA_11_14_sec[,c("range", "Within_PA", "Predominant_habitat", "LUPA",
	"log_slope", "log_elevation", "ag_suit",
	"Zone", "taxon_of_interest", "SS", "SSB", "SSBS")]
data <- na.omit(data)

r.m <- lmer(range~Predominant_habitat+Within_PA+Zone+taxon_of_interest
	+ Within_PA:Predominant_habitat
	+ ag_suit + I(ag_suit^2) + I(ag_suit^3)
	+ log_elevation+ I(log_elevation^2)
	+ log_slope + I(log_slope^2) + I(log_slope^3)
	+ (1+Within_PA|SS)+(1|SSB), data = data, 
	control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
summary(r.m)
length(fixef(r.m))
length(fixef(range.model$model))


# get endemicity estimates relative to reference

y1 <- 10 ^(as.numeric(fixef(r.m)[1]))
e.y1 <- 1/y1
y2 <- 10^(as.numeric(fixef(r.m)[pos]) + as.numeric(fixef(r.m)[1]))
e.y2 <- 1/y2

#as a percentage of outside 
e.relative <- e.y2/e.y1*100

se <- as.numeric(se.fixef(r.m)[pos])
y2plus <- 10^(as.numeric(fixef(r.m)[pos]) + as.numeric(fixef(r.m)[1])+ se*1.96)
e.y2plus <- 1/y2plus
e.relative.plus <- e.y2plus/e.y1*100

y2minus <- 10^(as.numeric(fixef(r.m)[pos]) + as.numeric(fixef(r.m)[1])- se*1.96)
e.y2minus <- 1/y2minus
e.relative.minus <- e.y2minus/e.y1*100

points <- c(100, e.relative)
CI <- cbind(e.relative.plus, e.relative.minus)


























