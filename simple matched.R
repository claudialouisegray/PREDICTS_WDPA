#setwd("C:/Users/Claudia/Documents/PREDICTS/WDPA analysis")

#rm(list=ls())

setwd("N:/Documents/PREDICTS/WDPA analysis")


library(lme4)
library(yarg)
library(roquefort)

# Load dataset on taxa_split matched all


PA_11_2014 <- read.csv("PA_11_2014.csv")
nrow(PA_11_2014 )#7077





### LOAD MATCHING DATA

setwd("N:/Documents/PREDICTS/WDPA analysis/matching data")

access <- read.table("bfer_1km_acc50k_11_14.txt", header = T, sep = ",")
hpd <- read.table("bfer_1km_HPD_11_14.txt", header = T, sep = ",")
elevation <- read.table("bfer_1km_mn30_elevation_11_14.txt", header = T, sep = ",")
slope <- read.table("bfer_1km_mn30_slope_11_14.txt", header = T, sep = ",")

ag_suit <- read.csv("ag_suitability_11_2014_moll.csv") 

# this is taken from extract values to points function, saved as text file, then opened in excel
# FID column deleted and -9999 values in RASTERVALU switched to NA

ag.1 <- ag_suit[,c("SSS", "RASTERVALU")] 
colnames(ag.1) <- c("SSS", "ag_suit")

#these are water, must have fallen on the edge of lakes/rivers
ag.1$ag_suit[which(ag.1$ag_suit==9)] <- NA

# switch scale to be more intuitive - higher values more agriculturally suitable
ag.1$ag_suit <- 9 - ag.1$ag_suit

access.1 <- access[,c("SSS", "MEAN")]
hpd.1 <- hpd[,c("SSS", "MEAN")]
elevation.1 <- elevation[,c("SSS", "MEAN")]
slope.1 <- slope[,c("SSS", "MEAN")]

m <- merge(PA_11_2014 , access.1, "SSS", all.x = T)
colnames(m)[which(colnames(m) == "MEAN")] <- "access"
m <- merge(m, hpd.1, "SSS", all.x = T)
colnames(m)[which(colnames(m) == "MEAN")] <- "hpd"
m <- merge(m, elevation.1, "SSS", all.x = T)
colnames(m)[which(colnames(m) == "MEAN")] <- "elevation"
m <- merge(m, slope.1, "SSS", all.x = T)
colnames(m)[which(colnames(m) == "MEAN")] <- "slope"
m <- merge(m, ag.1, "SSS", all.x = T)
nrow(m)
PA_11_2014  <- m

# create explanatory variables

PA_11_2014$IUCN_CAT_number <- factor(PA_11_2014$IUCN_CAT_number) # they arent really in an order
PA_11_2014$log_slope <- log(PA_11_2014$slope +1)
PA_11_2014$log_elevation <- log(PA_11_2014$elevation +1)
PA_11_2014$log_hpd<- log(PA_11_2014$hpd +1)
PA_11_2014$log_access <- log(PA_11_2014$access +1)
PA_11_2014$log_GIS_AREA <- log(PA_11_2014$GIS_AREA+1)

# make IUCN cat variable of I or II vs III to VI vs unknown vs unprotected
PA_11_2014$IUCN_CAT <- PA_11_2014$IUCN_CAT_number 
levels(PA_11_2014$IUCN_CAT) <- c(levels(PA_11_2014$IUCN_CAT), "0")
PA_11_2014$IUCN_CAT[which(PA_11_2014$Within_PA == "no")] <- 0

#make response variables

PA_11_2014$range <- PA_11_2014$CWM_Geographic_range_log10_square_km
PA_11_2014$mass <- PA_11_2014$CWM_Mass_log10_g  ### UPDATE
PA_11_2014$veg <- PA_11_2014$CWM_Vegetative_height_log10_m
PA_11_2014$vol <- PA_11_2014$CWM_Length_derived_volume_3log10_mm
PA_11_2014$log_abundance <- log(PA_11_2014$Total_abundance +1)


### CREATE MULTIPLE TAXA PER STUDY DATASET 
# create dataset for species richness analysis that without all studies that are only on one taxon

setwd("N:/Documents/PREDICTS/WDPA analysis")

studies.taxa <- read.csv("Number of taxa per study split_taxa_coarse 11_2014.csv")

which(studies.taxa$number.taxa == 1)
length(which(studies.taxa$number.taxa == 1)) #25

more.than.one.taxa <- studies.taxa$SS[which(studies.taxa$number.taxa != 1)]

multiple.taxa.PA_11_2014  <- subset(PA_11_2014 , SS %in% more.than.one.taxa)
length(multiple.taxa.PA_11_2014 [,1]) #6874

multiple.taxa.PA_11_2014  <- droplevels(multiple.taxa.PA_11_2014)







### model species richness

# check is block useful
# yes definitely

m1 <- glmer(Species_richness ~ Within_PA + log_slope + log_elevation + ag_suit
	+ (Within_PA|SS) + (1|SSBS), 
	family = "poisson", data = multiple.taxa.PA_11_2014)
m1b <- glmer(Species_richness ~ Within_PA + log_slope + log_elevation + ag_suit
	+ (Within_PA|SS) + (1|SSB)+ (1|SSBS), 
	family = "poisson", data = multiple.taxa.PA_11_2014)
anova(m1, m1b)



m1 <- glmer(Species_richness ~ Within_PA + log_slope + log_elevation + ag_suit
	+ (Within_PA|SS) + (1|SSB) + (1|SSBS), 
	family = "poisson", data = multiple.taxa.PA_11_2014)
m2 <- glmer(Species_richness ~ 1 + log_slope + log_elevation + ag_suit
	+ (Within_PA|SS) + (1|SSB) + (1|SSBS), 
	family = "poisson", data = multiple.taxa.PA_11_2014)
anova(m1, m2)

summary(m1)


fixef(m1)[2] # within_Pa yes = 0.0852
exp(fixef(m1)[2]) # 1.089


x <- exp(fixef(m1)[2]) - 1


# so if intercept is 100% (unprotected)
# protected is +8.9%



# plot


tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/12_14/simple model sp rich.tif",
	width = 10, height = 10, units = "cm", pointsize = 12, res = 300)


labels <- c("Unprotected", "Protected")
est.protect <- 1 + as.numeric(fixef(m1)[2])
se.protect <- as.numeric(se.fixef(m1)[2])
points <- c(1, est.protect)
CI <- c(est.protect - 1.96*se.protect, est.protect + 1.96*se.protect)


plot(points ~ c(1,2), ylim = c(0.8,1.3), xlim = c(0.5,2.5),
	bty = "l", pch = 16, col = c(1,3), cex = 1.5,
	yaxt = "n", xaxt = "n",
	ylab = "Species richness difference (% ± 95%CI)",
	xlab = "")
axis(1, c(1,2), labels)
axis(2, c(0.8,0.9, 1, 1.1, 1.2), c(80,90,100,110,120))
arrows(2,CI[1],2,CI[2], code = 3, length = 0.03, angle = 90)
abline(h = 1, lty = 2)
points(points ~ c(1,2), pch = 16, col = c(1,3), cex = 1.5)

dev.off()





#simple species richness with  Zone

m0z <- glmer(Species_richness ~ Within_PA + Zone + Within_PA:Zone
	+ log_slope + log_elevation + ag_suit
	+ (Within_PA|SS) + (1|SSB) + (1|SSBS), 
	family = "poisson", data = multiple.taxa.PA_11_2014)

m1z <- glmer(Species_richness ~ Within_PA + Zone + log_slope + log_elevation + ag_suit
	+ (Within_PA|SS) + (1|SSB) + (1|SSBS), 
	family = "poisson", data = multiple.taxa.PA_11_2014)

anova(m0z,m1z)



m2z <- update(m1z, .~. -Within_PA)
m3z <- update(m1z, .~. -Zone)

anova(m1z, m2z)
anova(m1z, m3z)

summary(m1z)






#simple species richness with IUCN cat 

m0i <- glmer(Species_richness ~ Within_PA + log_slope + log_elevation + ag_suit
	+ (Within_PA|SS) + (1|SSB) + (1|SSBS), 
	family = "poisson", data = multiple.taxa.PA_11_2014)

m1i <- glmer(Species_richness ~ IUCN_CAT + log_slope + log_elevation + ag_suit
	+ (Within_PA|SS) + (1|SSB) + (1|SSBS), 
	family = "poisson", data = multiple.taxa.PA_11_2014)


anova(m1i, m0i)

summary(m1i)






# rarefied richness

m3 <- glmer(Richness_rarefied ~ Within_PA + log_slope + log_elevation + ag_suit
	+ (Within_PA|SS)+ (1|SSB) + (1|SSBS), 
	family = "poisson", data = multiple.taxa.PA_11_2014)
m4 <- glmer(Richness_rarefied ~ 1 + log_slope + log_elevation + ag_suit
	+ (Within_PA|SS)+ (1|SSB) + (1|SSBS), 
	family = "poisson", data = multiple.taxa.PA_11_2014)
anova(m3, m4)


#no sig difference


tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/12_14/simple model rar rich.tif",
	width = 10, height = 10, units = "cm", pointsize = 12, res = 300)


labels <- c("Unprotected", "Protected")
est.protect <- 1 + as.numeric(fixef(m3)[2])
se.protect <- as.numeric(se.fixef(m3)[2])
points <- c(1, est.protect)
CI <- c(est.protect - 1.96*se.protect, est.protect + 1.96*se.protect)


plot(points ~ c(1,2), ylim = c(0.8,1.3), xlim = c(0.5,2.5),
	bty = "l", pch = 16, col = c(1,3), cex = 1.5,
	yaxt = "n", xaxt = "n",
	ylab = "Rarefied species richness difference (% ± 95%CI)",
	xlab = "")
axis(1, c(1,2), labels)
axis(2, c(0.8,0.9, 1, 1.1, 1.2), c(80,90,100,110,120))
arrows(2,CI[1],2,CI[2], code = 3, length = 0.03, angle = 90)
abline(h = 1, lty = 2)
points(points ~ c(1,2), pch = 16, col = c(1,3), cex = 1.5)

dev.off()





# rarefied richness with IUCN cat

m3i <- glmer(Richness_rarefied ~ Within_PA + log_slope + log_elevation + ag_suit
	+ (Within_PA|SS) + (1|SSB)+ (1|SSBS), 
	family = "poisson", data = multiple.taxa.PA_11_2014)
m4i <- glmer(Richness_rarefied ~ IUCN_CAT + log_slope + log_elevation + ag_suit
	+ (Within_PA|SS) + (1|SSB)+ (1|SSBS), 
	family = "poisson", data = multiple.taxa.PA_11_2014)

# this is the same as
#m4i <- glmer(Richness_rarefied ~ IUCN_CAT + Within_PA + log_slope + log_elevation + ag_suit
#	+ (Within_PA|SS)+ (1|SSB) + (1|SSBS), 
#	family = "poisson", data = multiple.taxa.PA_11_2014)

anova(m3i,m4i)

summary(m4i)




# plot 

tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/12_14/simple model rar rich IUCN.tif",
	width = 12, height = 10, units = "cm", pointsize = 12, res = 300)


labels <- c("Unprotected", "IUCN I & II", "IUCN III  - VI", "unknown")

levels.IUCN <- levels(multiple.taxa.PA_11_2014$IUCN_CAT)

multiple.taxa.PA_11_2014$IUCN_CAT <- relevel(multiple.taxa.PA_11_2014$IUCN_CAT, "0")

m <- glmer(Species_richness ~ IUCN_CAT + log_slope + log_elevation + ag_suit
	+ (Within_PA|SS) + (1|SSBS), 
	family = "poisson", data = multiple.taxa.PA_11_2014)

	
est.protect <- 1 + as.numeric(fixef(m)[2:4])
se.protect <- as.numeric(se.fixef(m)[2:4])
points <- c(1, est.protect)
CI <- cbind(est.protect - 1.96*se.protect, est.protect + 1.96*se.protect)

plot(points ~ seq(1,length(points),1), ylim = c(0.9,1.3), #xlim = c(0.5,2.5),
	bty = "l", pch = 16, col = c(1,3,3,3), cex = 1.5,
	yaxt = "n", xaxt = "n",
	ylab = "Rarefied species richness difference (% ± 95%CI)",
	xlab = "")
axis(1,seq(1,length(points),1), labels)
axis(2, c(0.8,0.9, 1, 1.1, 1.2, 1.3), c(80,90,100,110,120,130))
arrows(seq(2,length(points),1),CI[,1],
	seq(2,length(points),1),CI[,2], code = 3, length = 0.03, angle = 90)
abline(h = 1, lty = 2)
points(points ~ c(1,2), pch = 16, col = c(1,3), cex = 1.5)

dev.off()





















### model abundance
# block also very useful here
m1a <- lmer(log_abundance ~ Within_PA + log_slope + log_elevation + ag_suit
	+ (Within_PA|SS)+ (1|SSB), 
	 data = PA_11_2014)
m2a <- lmer(log_abundance ~ 1 + log_slope + log_elevation + ag_suit
	+ (Within_PA|SS)+ (1|SSB), 
	 data = PA_11_2014)
anova(m1a, m2a)

summary(m1a)

exp(fixef(m1a)[2]) # 1.127





# plot


tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/12_14/simple model abundance.tif",
	width = 10, height = 10, units = "cm", pointsize = 12, res = 300)

labels <- c("Unprotected", "Protected")
est.protect <- 1 + as.numeric(fixef(m1a)[2])
se.protect <- as.numeric(se.fixef(m1a)[2])
points <- c(1, est.protect)
CI <- c(est.protect - 1.96*se.protect, est.protect + 1.96*se.protect)


plot(points ~ c(1,2), ylim = c(0.8,1.3), xlim = c(0.5,2.5),
	bty = "l", pch = 16, col = c(1,3), cex = 1.5,
	yaxt = "n", xaxt = "n",
	ylab = "Abundance difference (% ± 95%CI)",
	xlab = "")
axis(1, c(1,2), labels)
axis(2, c(0.8,0.9, 1, 1.1, 1.2), c(80,90,100,110,120))
arrows(2,CI[1],2,CI[2], code = 3, length = 0.03, angle = 90)
abline(h = 1, lty = 2)
points(points ~ c(1,2), pch = 16, col = c(1,3), cex = 1.5)


dev.off()








# model abundance and zone

m0az <- lmer(log_abundance ~ Within_PA + Zone + Within_PA:Zone
	+log_slope + log_elevation + ag_suit
	+ (Within_PA|SS)+ (1|SSB) , 
	 data = PA_11_2014)

m1az <- lmer(log_abundance ~ Within_PA + Zone +log_slope + log_elevation + ag_suit
	+ (Within_PA|SS)+ (1|SSB), 
	 data = PA_11_2014)

anova(m0az, m1az)


m2az <- update(m1az, .~. - Within_PA)
m3az <- update(m1az, .~. - Zone)


anova(m1az, m2az)
anova(m1az, m3az)

summary(m1az)




# model abundance and IUCN_category
m1ai <- lmer(log_abundance ~ Within_PA +log_slope + log_elevation + ag_suit
	+ (Within_PA|SS)+ (1|SSB), 
	 data = PA_11_2014)
m2ai <- lmer(log_abundance ~ IUCN_CAT +log_slope + log_elevation + ag_suit
	+ (Within_PA|SS)+ (1|SSB), 
	 data = PA_11_2014)

anova(m1ai, m2ai)

summary(m2ai)
validate(m2ai)










### model range
m1r <- lmer(range ~ Within_PA + log_slope + log_elevation + ag_suit
	+ (Within_PA|SS)+ (1|SSB), 
	 data = PA_11_2014)
m2r <- lmer(range ~ 1 + log_slope + log_elevation + ag_suit
	+ (Within_PA|SS)+ (1|SSB), 
	 data = PA_11_2014)
anova(m1r, m2r)

summary(m1r)

exp(fixef(m1r)[2]) # 0.976

#convert to endemicity
1 - exp(fixef(m1r)[2])



# plot
#RANGE

tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/simple model PA_11_2014 range.tif",
	width = 10, height = 10, units = "cm", pointsize = 12, res = 300)

labels <- c("Unprotected", "Protected")
est.protect <- 1 + as.numeric(fixef(m1r)[2])
se.protect <- as.numeric(se.fixef(m1r)[2])
points <- c(1, est.protect)
CI <- c(est.protect - 1.96*se.protect, est.protect + 1.96*se.protect)


plot(points ~ c(1,2), ylim = c(0.8,1.3), xlim = c(0.5,2.5),
	bty = "l", pch = 16, col = c(1,3), cex = 1.5,
	yaxt = "n", xaxt = "n",
	ylab = "Range difference (% ± 95%CI)",
	xlab = "")
axis(1, c(1,2), labels)
axis(2, c(0.8,0.9, 1, 1.1, 1.2), c(80,90,100,110,120))
arrows(2,CI[1],2,CI[2], code = 3, length = 0.03, angle = 90)
abline(h = 1, lty = 2)
points(points ~ c(1,2), pch = 16, col = c(1,3), cex = 1.5)


dev.off()


#ENDEMICITY

tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/12_14/simple model endemicity.tif",
	width = 10, height = 10, units = "cm", pointsize = 12, res = 300)

labels <- c("Unprotected", "Protected")
est.protect <- 1 - as.numeric(fixef(m1r)[2]) #so that difference is actually above 1
se.protect <- as.numeric(se.fixef(m1r)[2])
points <- c(1, est.protect)
CI <- c(est.protect - 1.96*se.protect, est.protect + 1.96*se.protect)


plot(points ~ c(1,2), ylim = c(0.8,1.3), xlim = c(0.5,2.5),
	bty = "l", pch = 16, col = c(1,3), cex = 1.5,
	yaxt = "n", xaxt = "n",
	ylab = "Endemicity difference (% ± 95%CI)",
	xlab = "")
axis(1, c(1,2), labels)
axis(2, c(0.8,0.9, 1, 1.1, 1.2), c(80,90,100,110,120))
arrows(2,CI[1],2,CI[2], code = 3, length = 0.03, angle = 90)
abline(h = 1, lty = 2)
points(points ~ c(1,2), pch = 16, col = c(1,3), cex = 1.5)


dev.off()





### range and IUCN category

m1ri <- lmer(range ~ Within_PA +log_slope + log_elevation + ag_suit
	+ (Within_PA|SS)+ (1|SSB), 
	 data = PA_11_2014)
m2ri <- lmer(range ~ IUCN_CAT +log_slope + log_elevation + ag_suit
	+ (Within_PA|SS)+ (1|SSB), 
	 data = PA_11_2014)

anova(m1ri, m2ri)

summary(m2ri)









### model proportion threatened ###

PA_11_2014_a_m_b <- read.csv("PA_11_2014_amph_mamm_bird.csv")
nrow(PA_11_2014_a_m_b) #2695
names(PA_11_2014_a_m_b)



m <- merge(PA_11_2014_a_m_b , access.1, "SSS", all.x = T)
colnames(m)[which(colnames(m) == "MEAN")] <- "access"
m <- merge(m, hpd.1, "SSS", all.x = T)
colnames(m)[which(colnames(m) == "MEAN")] <- "hpd"
m <- merge(m, elevation.1, "SSS", all.x = T)
colnames(m)[which(colnames(m) == "MEAN")] <- "elevation"
m <- merge(m, slope.1, "SSS", all.x = T)
colnames(m)[which(colnames(m) == "MEAN")] <- "slope"
m <- merge(m, ag.1, "SSS", all.x = T)


nrow(m)
PA_11_2014_a_m_b  <- m





# create explanatory variables

PA_11_2014_a_m_b$IUCN_CAT_number <- factor(PA_11_2014_a_m_b$IUCN_CAT_number) # they arent really in an order
PA_11_2014_a_m_b$log_slope <- log(PA_11_2014_a_m_b$slope +1)
PA_11_2014_a_m_b$log_elevation <- log(PA_11_2014_a_m_b$elevation +1)
PA_11_2014_a_m_b$log_hpd<- log(PA_11_2014_a_m_b$hpd +1)
PA_11_2014_a_m_b$log_access <- log(PA_11_2014_a_m_b$access +1)
PA_11_2014_a_m_b$log_GIS_AREA <- log(PA_11_2014_a_m_b$GIS_AREA+1)



PA_11_2014_a_m_b$y <- cbind(PA_11_2014_a_m_b$abundance_VU_EN_CR, PA_11_2014_a_m_b$abundance_LC_NT)


m1t <- glmer(y ~ Within_PA + log_slope + log_elevation + ag_suit
	+ (Within_PA|SS)+ (1|SSB) + (1|SSBS), 
	family = "binomial", data = PA_11_2014_a_m_b)
m2t <- glmer(y ~ 1 + log_slope + log_elevation + ag_suit
	+ (Within_PA|SS)+ (1|SSB) + (1|SSBS), 
	family = "binomial", data = PA_11_2014_a_m_b)
anova(m1t, m2t)


#PLOT
tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/12_14/simple model prop threat.tif",
	width = 10, height = 10, units = "cm", pointsize = 12, res = 300)

labels <- c("Unprotected", "Protected")
est.protect <- 1 - as.numeric(fixef(m1t)[2]) #so that difference is actually above 1
se.protect <- as.numeric(se.fixef(m1t)[2])
points <- c(1, est.protect)
CI <- c(est.protect - 1.96*se.protect, est.protect + 1.96*se.protect)


plot(points ~ c(1,2), ylim = c(0.8,1.3), xlim = c(0.5,2.5),
	bty = "l", pch = 16, col = c(1,3), cex = 1.5,
	yaxt = "n", xaxt = "n",
	ylab = "Proportion threatened species difference (% ± 95%CI)",
	xlab = "")
axis(1, c(1,2), labels)
axis(2, c(0.9, 1, 1.1, 1.2, 1.3), c(90, 100, 110, 120, 130))
arrows(2,CI[1],2,CI[2], code = 3, length = 0.03, angle = 90)
abline(h = 1, lty = 2)
points(points ~ c(1,2), pch = 16, col = c(1,3), cex = 1.5)


dev.off()







# IUCN CAT
m1ti <- glmer(y ~ Within_PA + log_slope + log_elevation + ag_suit
	+ (Within_PA|SS)+ (1|SSB) + (1|SSBS), 
	family = "binomial", data = PA_11_2014_a_m_b)
m2ti <- glmer(y ~ IUCN_CAT + log_slope + log_elevation + ag_suit
	+ (Within_PA|SS)+ (1|SSB) + (1|SSBS), 
	family = "binomial", data = PA_11_2014_a_m_b)
anova(m1ti, m2ti)














#PA effectiveness estimates

# benefit for species richness
#
b.sp <- exp(fixef(m1)[2]) -1
b.sp.max <- (exp(fixef(m1)[2]) + 1.96*se.fixef(m1)[2])-1
b.sp.min <- (exp(fixef(m1)[2]) - 1.96*se.fixef(m1)[2])-1

b.a <- exp(fixef(m1a)[2]) -1
b.a.max <- (exp(fixef(m1a)[2]) + 1.96*se.fixef(m1)[2])-1
b.a.min <- (exp(fixef(m1a)[2]) - 1.96*se.fixef(m1)[2])-1

b.vals <- c(b.sp, b.sp.max, b.sp.min, b.a, b.a.max, b.a.min)
all.results <- vector()

b <- b.a.max
b <- b.a


for(b in b.vals){

benefit <- as.numeric(b)		# percentage increase in metric in PAs
PA.pct <- 13 				# percentage of total land area in PAs
global.loss <- 0.116			# global loss of biodiversity (from Newbold et al)

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


all.results <- c(all.results, est)
}

d <- data.frame(estimates= all.results, 
	values = c("b.sp", "b.sp.max", "b.sp.min", "b.a", "b.a.max", "b.a.min"))

write.csv(t(d), "simple.effectiveness.estimates.csv")






# Getting andys numbers - my approach written before I got Andys code. 


# so, if PAs are 7.5% higher in species richness
# Sites outside PAs have 92.5% as many species as sites inside
# we also know 13% of land surface is in Protected Areas
# Globally PREDICTS (nature MS) estimates mean net loss of species from terrestrial sites 
	# across all landuses at 12.9% (relative to pristine)

# what are Protected and Unprotected relevant to pristine
# n = % sp  in unprotected relative to pristine
# y = % sp  in protected relative to pristine

# total loss  = 87% loss in unprotected and 13% loss in protected
# simultaneous equations are
#1 - 0.129 = 0.87*n + 0.13*y
#n * 1.075 = y


# so
#1 - 0.129 = 0.87*n + 0.13*1.075*n


#then
n <- 1 - 0.129/(0.87 + 0.13*1.075) 
y <- n * 1.075

n #0.872
y #0.927

loss of species in unprotected relative to pristine = 0.1277544
loss of species in protected relative to pristine = 

Simultaneous equations suggest:
#PA sites have 93.2% as many species as pristine sites; non-PA sites have 86.2%

# 0% effective is 86.2%, 100% effective is the same as pristine
100 - 86.2 # 13.8
93.2-86.2  #7
7/13.8 = 50.7%
# PAs are therefore 51% effective in retaining site-level species richness

# to vary parameters

total.loss <- 1 - 0.129
p.protect <- 0.13
p.unprotect <- 1- p.protect
x <- as.numeric(exp(fixef(m1)[2]))

n <- total.loss/(p.unprotect + p.protect*x) 
y <- n * x

a <- 1 - n
b <- y - n
b/a


