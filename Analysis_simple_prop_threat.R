### model proportion threatened ###

rm(list=ls())

library(lme4)
library(yarg)
library(roquefort)
library(influence.ME)

# load functions
setwd("R:/ecocon_d/clg32/GitHub/PREDICTS_WDPA")
source("compare_randoms.R")
source("model_select.R")


source("prep_PA_11_14_for_analysis.R")

setwd("N:/Documents/PREDICTS/WDPA analysis")


PA_11_14_a_m_b <- read.csv("PA_11_2014_amph_mamm_bird.csv")
nrow(PA_11_14_a_m_b) #2695
names(PA_11_14_a_m_b)



m <- merge(PA_11_14_a_m_b , elevation.1, "SSS") 
colnames(m)[which(colnames(m) == "MEAN")] <- "elevation"
m <- merge(m, slope.1, "SSS") 
colnames(m)[which(colnames(m) == "MEAN")] <- "slope"
m <- merge(m, ag.1, "SSS") 

nrow(m)
PA_11_14_a_m_b  <- m





# create explanatory variables

PA_11_14_a_m_b$IUCN_CAT_number <- factor(PA_11_14_a_m_b$IUCN_CAT_number) # they arent really in an order
PA_11_14_a_m_b$log_slope <- log(PA_11_14_a_m_b$slope +1)
PA_11_14_a_m_b$log_elevation <- log(PA_11_14_a_m_b$elevation +1)
PA_11_14_a_m_b$log_GIS_AREA <- log(PA_11_14_a_m_b$GIS_AREA+1)



PA_11_14_a_m_b$y <- cbind(PA_11_14_a_m_b$abundance_VU_EN_CR, PA_11_14_a_m_b$abundance_LC_NT)


# test nonlinear confoundinf variables

fF <- c("Within_PA") 
fT <- list("ag_suit" = "3", "log_slope" = "3", "log_elevation" = "3")
keepVars <- list()
fI <- character(0)
RS <-  c("Within_PA")

RLS.model<- model_select(all.data  = PA_11_14_a_m_b, 
			     responseVar = "y", 
			     fitFamily = "binomial",
			     alpha = 0.05,
                       fixedFactors= fF,
                       fixedTerms= fT,
			     keepVars = keepVars,
                       randomStruct = "(Within_PA|SS) + (1|SSB)+ (1|SSBS)",
			     otherRandoms=character(0),
                       verbose=TRUE)
RLS.model$final.call
# y~poly(log_elevation,3)+poly(log_slope,3)+(Within_PA|SS)+(1|SSB)+(1|SSBS)

data <- PA_11_14_a_m_b[,c("ag_suit", "log_elevation", "log_slope", "Within_PA", "SS", "SSB", "SSBS", "y")]
data <- na.omit(data)

m1t <- glmer(y ~ Within_PA + poly(log_elevation,3)+poly(log_slope,3) + ag_suit
	+ (Within_PA|SS)+ (1|SSB) + (1|SSBS), 
	family = "binomial", data = data)
m2t <- glmer(y ~ 1 + poly(log_elevation,3)+poly(log_slope,3) + ag_suit
	+ (Within_PA|SS)+ (1|SSB) + (1|SSBS), 
	family = "binomial", data = data)
anova(m1t, m2t)

summary(m1t)


#PLOT
tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/02_15/simple model prop threat.tif",
	width = 10, height = 10, units = "cm", pointsize = 12, res = 300)

labels <- c("Unprotected", "Protected")
y <- as.numeric(fixef(m1t)[2])
se <- as.numeric(se.fixef(m1t)[2])
yplus <- y + se*1.96
yminus <- y - se*1.96

intercept <-fixef(m1t)[1]
true.intercept <- 1/(1+exp(-(intercept)))
true.y <- 1/(1+exp(-(intercept+y)))
true.yplus <- 1/(1+exp(-(intercept+yplus)))
true.yminus <- 1/(1+exp(-(intercept+yminus)))

# this is code from plot_lu_effects -
# why is it backtranformed absolute value of non-ref level/back transformed value of ref variable?
# because!!!
# dont want the percentage of the difference, want it as 100+x percent. (Tims code then -100 at the end to get difference). 

y<-((true.y/true.intercept)*100)
yplus<-((true.yplus/true.intercept)*100)
yminus<-((true.yminus/true.intercept)*100)

points <- c(100, y)
CI <- c(yplus, yminus)

par(mar = c(4,6,1,1))
plot(points ~ c(1,2), ylim = c(0,350), xlim = c(0.5,2.5), 
	bty = "l", pch = c(21,16), bg = "white", col = 1, cex = 1.5,
	xaxt = "n", yaxt = "n", 
	ylab = "Proportion threatened species\ndifference (% ± 95%CI)",
	xlab = "")
axis(1, c(1,2), labels)
axis(2, c(0,50, 100,200,300), c(0,50,100,200,300))
arrows(2,CI[1],2,CI[2], code = 3, length = 0, angle = 90)
abline(h = 100, lty = 2)
points(points ~ c(1,2), pch = c(21,16), bg = "white", cex = 1.5)

data <- PA_11_14_a_m_b[,c("Within_PA", "SSS", "y")]
data <- na.omit(data)
text(2,0, paste("n =", length(data$SSS[which(data$Within_PA == "yes")]), sep = " "))
text(1,0, paste("n =", length(data$SSS[which(data$Within_PA == "no")]), sep = " "))


dev.off()
