
# describing data 


setwd("N:/Documents/PREDICTS/WDPA analysis")

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
library(MatchIt)



validate <- function(x) {
  par(mfrow = c(1,2))
  plot(resid(x)~ fitted(x))
  hist(resid(x))
  par(mfrow = c(1,1))
}



setwd("R:/ecocon_d/clg32/GitHub/PREDICTS_WDPA")


# load functions
source("compare_randoms_lmer - with poly.R")
source("model_select.R")


#load data for matched studies
source("WDPA_predicts_prep_matched.landuse_for_analysis.R")
data <- matched.landuse

# or if using LUPA data
source("WDPA_predicts_prep_PA_11_14_for_analysis.R")
data <- PA_11_14


data$SS_PH <- paste(data$SS, data$Predominant_habitat)





# how many sites and studies in each section of the data
nrow(data)
length(unique(data$SS))

nrow(data[which(data$Within_PA == "yes"),])
length(unique(data$SS[which(data$Within_PA == "yes")]))






# how many comparisons in each landuse pairing

LUnames <- unique(data$Predominant_habitat)

LUnames <-  LUnames[c(2,3,6,1,4,5)]
comparisons <- matrix(0,nrow = 6, ncol = 6)
comp <- data.frame(comparisons)
colnames(comp) <-  LUnames
rownames(comp) <- LUnames
comp

studies <- unique(data$SS)
for(s in studies){
	data.subset <- subset(data, SS == s)
	landuses <- unique(data.subset$Predominant_habitat)
	for(l in LUnames){
		if(l %in% landuses){comp[l,l] <- comp[l,l] + 1}
		for(u in LUnames){
			if(l %in% landuses & u %in% landuses){comp[l,u] <- comp[l,u]+1}
	}
   }
}
comp 

write.csv(comp,"studies per landuse comparison 15 09 14.csv")









# how many sites in each land use
 
data$Predominant_habitat <- factor(data$Predominant_habitat)

sites_total <- tapply(data$SSS,data$Predominant_habitat, length)


in_PA <- subset(data, Within_PA == "yes")
sites_in <- tapply(in_PA$SSS, in_PA$Predominant_habitat, length)

s <- data.frame(cbind(sites_total, sites_in))
s$sites_out <- s$sites_total - s$sites_in

write.csv(s, "sites.in.Predominat.Habitats.csv")













# how many studies give these sites

# get table of sites per pred_habitat in each study first


sites_total <- aggregate(SSS ~ SS + Predominant_habitat + SS_PH, data, length)
colnames(sites_total) <- c("SS", "PH", "SS_PH", "sites_total")

in_PA <- subset(data, Within_PA == "yes")
sites_in <- aggregate(SSS ~ SS + Predominant_habitat + SS_PH, in_PA, length)
colnames(sites_in) <- c("SS", "PH", "SS_PH", "sites_in")

View(sites_total)

s <- merge(sites_total, sites_in, by = c("SS", "PH", "SS_PH"))
s$sites_out <- s$sites_total - s$sites_in


# then for each habitat I can get the number of studies contributing data

list <- unique(s$PH)
output <- matrix(nrow = length(list), ncol = 2)

paste(list[1])
df <- subset(s, PH == paste(list[1]))
count <- length(unique(data$SS))
output[1,1] <- paste(list[1])
output[1,2] <- count


for(i in 1:length(list)){
	df <- subset(s, PH == paste(list[i]))
	count <- length(unique(data$SS_PH))
	output[i,1] <- paste(list[i])
	output[i,2] <- count
	print(i)
}
output <- data.frame(output)

write.csv(output, "studies contributing sites in and out of PA matched by Predominant Habitat.csv")











#how many studies have <10 sites remaining? 


sites_per_study<- aggregate(SSS ~ SS, data, length)
View(sites_per_study)

length(sites_per_study$SS)
length(which(sites_per_study$SSS < 10))

length(which(sites_per_study$SSS < 5))






# how many sites per realm


names(in_PA)
sites_total <- tapply(data$SSS, data$Realm, length)

in_PA <- subset(data, Within_PA == "yes")
sites_in <- tapply(in_PA$SSS, in_PA$Realm, length)

per.realm <- data.frame(cbind(sites_total, sites_in))
per.realm$sites_out = per.realm$sites_total - per.realm$sites_in


write.csv(per.realm, "Sites per realm.csv")








# how many sites per biome

sites_total <- tapply(data$SSS, data$Biome, length)

in_PA <- subset(data, Within_PA == "yes")
sites_in <- tapply(in_PA$SSS, in_PA$Biome, length)

per.biome<- data.frame(cbind(sites_total, sites_in))
per.biome$sites_out = per.realm$sites_total - per.realm$sites_in


#write.csv(per.biome, "Sites per biome.csv")










# where are they


# can do this in R
# to get plot of inside and outside points


in_PA <- subset(data, Within_PA == "yes" )
outside_PA <- subset(data, Within_PA == "no" )

long.lat <- in_PA[,c("Longitd", "Latitud")]
spatial_in <- SpatialPoints(long.lat)

long.lat <- outside_PA[,c("Longitd", "Latitud")]
spatial_outside <- SpatialPoints(long.lat)


countries <- readOGR(dsn = "C:/GIS/Data/world_borders", "world_country_admin_boundary_shapefile_with_fips_codes")
str(countries)

par(mfrow = c(2,1))
par(mar=c(0.2,0.1, 1.5,0.2))

plot(countries, lty = NULL, border = "grey", col = "grey")
plot(spatial_in, pch = 21, col =1, bg = "blue", add = TRUE, cex = 0.5)
title(main= "Sites inside PA")

plot(countries, lty = NULL, border = "grey", col = "grey")
plot(spatial_outside, pch = 21, col =1, bg = "red", add = TRUE, cex= 0.5)
title(main= "Sites outside PA")

plot(Latitud ~ Longitd,data) 








# where are the different land use types

list <- c("Primary Forest", 
	"Primary Non-forest", 
	"Mature Secondary Vegetation", 
	"Intermediate Secondary Vegetation", 
	"Young Secondary Vegetation", 
	"Secondary vegetation (indeterminate age)",  
	"Plantation forest", 
	"Cropland", 
	"Pasture", 
	"Urban") 

cols = c(rep("#5B8A3B", 1), 	# primary forest
		rep("#97bb6a", 1), 	# primary non-forest
		rep("#147659", 1), 	# MSV
       	rep("#1B9E77", 1), 	# ISV
		rep("#8ecfbc", 1),	# YSV
		rep("#b4d7db", 1), 	# SV (Ind)
		rep("#7570B3", 1), 	# Plantation
		rep("#E6AB02", 1), 	# cropland
		rep("#D95F02", 1), 	# pasture
		rep("#E7298A", 1))	# urban


i <- 2
for(i in 1:length(list)){
	data <- subset(data, Predominant_habitat == paste(list[i]) )
	long.lat <- data[,c("Longitd", "Latitud")]
	spatial_data <- SpatialPoints(long.lat)
	tiff(paste(list[i],"sites location map.tif"), width = 2500, height = 1200, res = 500, pointsize = 8)
	par(mar=c(0.2,0.1,1.5,0.2))
	plot(countries, lty = NULL, border = "grey", col = "grey", xpd = NULL)
	plot(spatial_data, pch = 21, col = 1, bg = paste(cols[i]), add = TRUE, cex= 0.5)
	title(main= paste(list[i]))
	dev.off()
	print(i)
}








# also, as mentioned by Andy, there may be problems if 
	# any one study has only one site inside or outside


s <- tapply(data$SSS, data$SS, length) # gives how many sites each study has

sites_in <- aggregate(SSS ~ SS, in_PA, length)
names(sites_in)
View(sites_in)
length(which(sites_in$SSS ==1)) #17 studies have only 1 site in a PA

#which land uses are these 
studies <- sites_in$SS[which(sites_in$SSS ==1)]

temp <- subset(data, SS %in% studies == TRUE)
View(temp)
unique(temp$SS_PH)

outside_PA <- subset(data, Within_PA = "no")
s.out <- aggregate(SSS ~ SS, outside_PA, length)
length(which(s.out$SSS ==1)) # 0 studies have only one site outside a PA

length(unique(data$SS))






# how many different land uses in vs outside for each study
# using this dataset
# can we make comparisons between landuses within PA


length(unique(data$SS))

length(unique(data$SS_PH))

studies <- unique(data$SS)
i <- 2
v <- rep(NA, length(studies))
for(i in 1:length(studies)){
	focal <-  studies[i]
	dataset <- subset(data, SS == paste(focal))
	v[i] <- length(unique(dataset$Predominant_habitat))
}
View(v)

length(which(v ==1)) #LUPA = 57 studies only have one habitat
length(which(v >1 )) #LUPA = 144 studies have more than one habitat


studies.overlap <- studies[which(v > 1)]
i <- 2
w <- rep(NA, length(studies.overlap))
for(i in 1:length(studies)){
	focal <-  studies[i]
	dataset <- subset(data, SS == paste(focal))
	w[i] <- paste(unique(dataset$Predominant_habitat))
}

View(w)
#there is pasture there,  a complete lack of overlap between pasture and other landuses is not the problem

length(which(w == "Pasture")) # 14




# looking at pasture
# get the sites in this subset that have pasture in them 
 p <- subset(data, Predominant_habitat == "Pasture")
pasture.SS <- unique(p$SS)

# take these out of the original database
pasture.diversity <- subset(diversity, SS %in% pasture.SS == TRUE)
View(pasture.diversity)

#does sampling effort vary amoung sites
t <- aggregate(Sampling_effort ~ SSS + SS,pasture.diversity, max)
View(t)



#get max and min sampling effort for each study

t.min <- aggregate(Sampling_effort ~ SS ,t, min)
t.max <- aggregate(Sampling_effort ~ SS ,t, max)
t <- merge(t.min, t.max, by = "SS")



#only one study where the sampling effort varies between sites
# AD1_2007__Ockinger 1 

subset(data, SS == "AD1_2007__Ockinger 1")
#this study is a comparison to secondary vegetation






### What is the distance between sites in the same habitat?



earth.dist <- function (long1, lat1, long2, lat2)
#uses haversine formula
{
rad <- pi/180
a1 <- lat1 * rad
a2 <- long1 * rad
b1 <- lat2 * rad
b2 <- long2 * rad
dlon <- b2 - a2
dlat <- b1 - a1
a <- (sin(dlat/2))^2 + cos(a1) * cos(b1) * (sin(dlon/2))^2
c <- 2 * atan2(sqrt(a), sqrt(1 - a))
R <- 6378.145
d <- R * c
return(d)
}




names(data)

results <- data.frame()
SS_PH <- unique(data$SS_PH)
s <- SS_PH[1]

for (p in SS_PH){
	ds <- subset(data, SS_PH == p)

# get names of all comparisons we want

SSS <- ds$SSS
SSS1 <- vector()
SSS2 <- vector()
count <- nrow(ds) -1

for(i in 1:nrow(ds)){
	SSS1 <- c(SSS1, rep(paste(SSS[i]), count))
	count <- count -1
	}

count <- 2
while(count <= nrow(ds)){
	SSS2 <- c(SSS2, paste(SSS[count:nrow(ds)]))
	count <- count + 1
	}

to.compare <- cbind(SSS1, SSS2)

lat1  <- vector()
long1 <- vector()
lat2  <- vector()
long2 <- vector()


# calculate distances
s<- 1
for(s in 1:length(SSS1)){
	lat1[s] <- ds$Latitud[which(SSS1[s] == droplevels(ds$SSS))]
	}
for(s in 1:length(SSS1)){
	long1[s] <- ds$Longitd[which(SSS1[s] == droplevels(ds$SSS))]
	}
for(s in 1:length(SSS2)){
	lat2[s] <- ds$Latitud[which(SSS2[s] == droplevels(ds$SSS))]
	}
for(s in 1:length(SSS2)){
	long2[s] <- ds$Longitd[which(SSS2[s] == droplevels(ds$SSS))]
	}

res <- data.frame(SSS1 = SSS1, SSS2 = SSS2, lat1 = lat1, lat2 = lat2, long1 = long1, long2 = long2)
res$dist_km <- earth.dist(long1, lat1, long2, lat2)

results <- rbind(results, res)
					
print(p)
}

#CHANGE NAME TO REFLECT DATA
write.csv(results, "distances.between.sites.PA_11_2014.csv")



hist(results$dist_km)
quantile(results$dist_km, c(0.05,.1, .12, .15))
# 90% of data is more than 1.6 km apart
# 85% of data is more than 2.6 km apart
# 88% of data more than 2km apart

quantile(results$dist_km)


max(results$dist_km)
min(results$dist_km)
median(results$dist_km)
mean(results$dist_km)

unique(results$SSS1[which(results$dist_km == 0)]) #806 sites are in the same place as another one

unique(results$SSS1) #9028 sites in total









### What is value of key matching variables in protected vs unprotected, max, min, mean, median and the variation? 

hpd.n <- tapply(matched.landuse$hpd, matched.landuse$Within_PA, length)
hpd.mean <- tapply(matched.landuse$hpd, matched.landuse$Within_PA, mean)
hpd.median <- tapply(matched.landuse$hpd, matched.landuse$Within_PA, median)
hpd.max <- tapply(matched.landuse$hpd, matched.landuse$Within_PA, max)
hpd.min <- tapply(matched.landuse$hpd, matched.landuse$Within_PA, min)

hpd.describe <- rbind(hpd.n, hpd.mean, hpd.median, hpd.max, hpd.min)
hist(matched.landuse$hpd[matched.landuse$Within_PA == "yes"])
hist(matched.landuse$hpd[matched.landuse$Within_PA == "no"])


access.n <- tapply(matched.landuse$access, matched.landuse$Within_PA, length)
access.mean <- tapply(matched.landuse$access, matched.landuse$Within_PA, mean)
access.median <- tapply(matched.landuse$access, matched.landuse$Within_PA, median)
access.max <- tapply(matched.landuse$access, matched.landuse$Within_PA, max)
access.min <- tapply(matched.landuse$access, matched.landuse$Within_PA, min)

access.describe <- rbind(access.n, access.mean, access.median, access.max, access.min)
hist(matched.landuse$access[matched.landuse$Within_PA == "yes"])
hist(matched.landuse$access[matched.landuse$Within_PA == "no"])



ag_suit.n <- tapply(matched.landuse$ag_suit, matched.landuse$Within_PA, length)
ag_suit.mean <- tapply(matched.landuse$ag_suit, matched.landuse$Within_PA, mean, na.rm = T)
ag_suit.median <- tapply(matched.landuse$ag_suit, matched.landuse$Within_PA, median, na.rm = T)
ag_suit.max <- tapply(matched.landuse$ag_suit, matched.landuse$Within_PA, max, na.rm = T)
ag_suit.min <- tapply(matched.landuse$ag_suit, matched.landuse$Within_PA, min, na.rm = T)

ag_suit.describe <- rbind(ag_suit.n, ag_suit.mean, ag_suit.median, ag_suit.max, ag_suit.min)
hist(matched.landuse$ag_suit[matched.landuse$Within_PA == "yes"])
hist(matched.landuse$ag_suit[matched.landuse$Within_PA == "no"])


### do protected and unprotected sites differ in elevation/slope/ag suit


plot(log_elevation ~ Within_PA, PA_11_2014)

m <- lmer(log_elevation ~ Within_PA + (Within_PA|SS) + (1|SSB), PA_11_2014)
n <- lmer(log_elevation ~ 1 + (Within_PA|SS) + (1|SSB), PA_11_2014)
summary(m)
anova(m,n)

m <- lmer(log_slope ~ Within_PA + (Within_PA|SS) + (1|SSB), PA_11_2014)
n <- lmer(log_slope ~ 1 + (Within_PA|SS) + (1|SSB), PA_11_2014)
summary(m)
anova(m,n)

m <- lmer(ag_suit ~ Within_PA + (Within_PA|SS) + (1|SSB), PA_11_2014)
n <- lmer(ag_suit ~ 1 + (Within_PA|SS) + (1|SSB), PA_11_2014)
summary(m)
anova(m,n)





### check correlations ###
# matched landuse data

# ag suit and elevation 
plot(ag_suit ~ log_elevation, data)
ag.elev <- lmer(ag_suit ~  log_elevation + (log_elevation|SS) + (1|SSB), data)
ag.elev.null <- lmer(ag_suit ~  1 + (log_elevation|SS) + (1|SSB), data)
anova(ag.elev, ag.elev.null)
#weakly correlated
#get R-squared
cor.test(data$ag_suit, data$log_elevation)$estimate #Pearson's product moment correlation coefficient
# so r2 = 
cor.test(data$ag_suit, data$log_elevation)$estimate^2
# same as estimate from linear model
summary(lm(data$ag_suit ~ data$log_elevation))

# but R2 for GLMM more tricky. one possibility is:
# this is from http://glmm.wikidot.com/faq. this is the regression of the response against fitted values 
r2.corr.mer <- function(m) {
   lmfit <-  lm(model.response(model.frame(m)) ~ fitted(m))
   summary(lmfit)$r.squared
}
m <- ag.elev
r2.corr.mer(ag.elev)

# better version 
source("N:/R notes/R squared for glmm JonLefcheck.R")
rsquared.glmm(ag.elev) #marginal is fixed only, conditional is fixed and random


# slope and elevation 
plot(log_slope ~ log_elevation, data)
slope.elev <- lmer(log_slope ~  log_elevation + (log_elevation|SS) + (1|SSB), data)
slope.elev.null <- lmer(log_slope ~  1 + (log_elevation|SS) + (1|SSB), data)
anova(slope.elev, slope.elev.null)
rsquared.glmm(slope.elev)
# more strongly correlated


# slope and ag suit 
plot(ag_suit ~ log_slope, data)
ag.slope <- lmer(ag_suit ~  log_slope + (log_slope|SS) + (1|SSB), data)
ag.slope.null <- lmer(ag_suit ~  1 + (log_slope|SS) + (1|SSB), data)
anova(ag.slope, ag.slope.null)
rsquared.glmm(ag.slope)



# for PA characteristics, take only points inside PA
PA.data <- subset(data, Within_PA == "yes")

# PA size and elevation 
plot(log_AREA.PA ~ log_elevation, PA.data)
m1 <- lmer(log_AREA.PA ~  log_elevation + (log_elevation|SS) + (1|SSB), PA.data)
m2 <- lmer(log_AREA.PA ~  1 + (log_elevation|SS) + (1|SSB), PA.data)
anova(m1, m2)
rsquared.glmm(m1)
#weakly correlated

# PA size and slope
plot(log_AREA.PA ~ log_slope, PA.data)
m1 <- lmer(log_AREA.PA ~  log_slope + (log_slope|SS) + (1|SSB), PA.data)
m2 <- lmer(log_AREA.PA ~  1 + (log_slope|SS) + (1|SSB), PA.data)
anova(m1, m2)
rsquared.glmm(m1)
# not correlated

# PA size and ag suit 
plot(log_AREA.PA ~ ag_suit, PA.data)
m1 <- lmer(log_AREA.PA ~  ag_suit + (ag_suit|SS) + (1|SSB), PA.data)
m2 <- lmer(log_AREA.PA ~  1 + (ag_suit|SS) + (1|SSB), PA.data)
anova(m1, m2)
rsquared.glmm(m1)
# not correlated

# PA age and elevation 
plot(DoP  ~ log_elevation, PA.data)
m1 <- lmer(DoP ~  log_elevation + (log_elevation|SS) + (1|SSB), PA.data)
m2 <- lmer(DoP  ~  1 + (log_elevation|SS) + (1|SSB), PA.data)
anova(m1, m2)
rsquared.glmm(m1)
#not correlated


# PA age and slope
plot(DoP ~ log_slope, PA.data)
m1 <- lmer(DoP ~  log_slope + (log_slope|SS) + (1|SSB), PA.data)
m2 <- lmer(DoP ~  1 + (log_slope|SS) + (1|SSB), PA.data)
anova(m1, m2)
rsquared.glmm(m1)
# not correlated

# PA age and ag suit 
plot(DoP ~ ag_suit, PA.data)
m1 <- lmer(DoP ~  ag_suit + (ag_suit|SS) + (1|SSB), PA.data)
m2 <- lmer(DoP ~  1 + (ag_suit|SS) + (1|SSB), PA.data)
anova(m1, m2)
rsquared.glmm(m1)
# not correlated


# PA age and size 
plot(DoP ~ log_AREA.PA, PA.data)
m1 <- lmer(DoP ~  log_AREA.PA + (log_AREA.PA|SS) + (1|SSB), PA.data)
m2 <- lmer(DoP ~  1 + (log_AREA.PA|SS) + (1|SSB), PA.data)
anova(m1, m2)
rsquared.glmm(m1)
# weakly correlated - older PAs tend to be larger for this dataset








###


# comparing range and endemicity estimates
# in this data - range can vary between 5.96 and 6.2, for example
# so endemicity measure of 1/5.9 = 0.1694915 = 794328.2  km2 
# and endemicity measure of 1/6.1 = 0.1639344 = 1258925 sq km2 respectively
# so a difference of -0.005557099 in endemicity  = increase in 464597.2 km2 (size of spain is 504,645km2)

# what is a difference of 0.01 in endemicity, so between 0.16 and 0.17
# 0.16 = 1/6.25 = 1778279 km2
# 0.17 = 1/5.882353 = 762698.6
# so endemicity value lower by 0.01 = increase in CWM range of 1015580 km2
# this is over 4 times the size of the UK, and similar to bolivia or ethiopia

# what does this mean in relations to the endemic species? 

# load diversity and get range data

min(diversity$Geographic_range_log10_square_km, na.rm = T)
max(diversity$Geographic_range_log10_square_km, na.rm = T)
hist(diversity$Geographic_range_log10_square_km, na.rm = T)
quantile(diversity$Geographic_range_log10_square_km, na.rm = T)
10^quantile(diversity$Geographic_range_log10_square_km, na.rm = T)
# so more than 50% of the species in the database have a range smaller than 1000000 km2

min.range <- min(diversity$Geographic_range_log10_square_km, na.rm = T)
max.range <- max(diversity$Geographic_range_log10_square_km, na.rm = T)
lower.quart <- quantile(diversity$Geographic_range_log10_square_km, na.rm = T)[2]
10^min.range
10^max.range
10^lower.quart = 208838 km2


#make up some communities
sp.range <- aggregate(Geographic_range_log10_square_km ~ Best_guess_binomial, diversity, unique)
sp.range <- sp.range[order(sp.range$Geographic_range_log10_square_km), ]

mean.range.small <- mean(sp.range$Geographic_range_log10_square_km[1:10])
nrow(sp.range)
mean.range.large <- mean(sp.range$Geographic_range_log10_square_km[16607:16617])
#so a community made up of equal numbers of 10 smallest range species would be
1/mean.range.small ; 10^mean.range.small
#anda community made up of equal numbers of 10 largest range species would be
1/mean.range.large; 10^mean.range.large
# and the difference is
10^mean.range.large - 10^mean.range.small
#31,736,682





### PROPENSITY MATCHING ######
# How many sd of protected area mean is the unprotected data


matched.landuse$protected <- 0
matched.landuse$protected[which(matched.landuse$Within_PA == "yes")] <- 1
names(matched.landuse)
matching.data <- matched.landuse[,c("protected", "ag_suit", "log_access", "log_hpd", 
	"log_elevation", "log_slope",
	"Predominant_habitat", "Country", "SS")]
matching.data <- na.omit(matching.data)

#simplest form
#distance/propensity varies with variables put in to match to
matching <- matchit(protected ~ ag_suit + log_access + log_hpd + log_elevation + log_slope
	+Predominant_habitat + Country + SS, matching.data)
summary(matching$model)
matching.data$distance <- matching$distance
matching.data$discarded
matching.data$weights
matching$nn
plot(distance ~ protected, matching.data)


# this gives a measure of whether we would expect the site to be protected or not
# where 1 = protected


#overall, do unprotected fall within 1sd either side of mean of protected? 

# get the propensity matching score for the protected sites
p.dist <- matching.data$distance[which(matching.data$protected == 1)]
hist(p.dist)
mean.p <- mean(p.dist)
sd.p <- sd(p.dist)
lower <- mean.p - 2*sd.p
upper <- mean.p + 2*sd.p



# what proportion of unprotected sites fall within 2sd of the mean
# get propensity scores for the unprotected data
unp.dist <- matching.data$distance[which(matching.data$protected == 0)]
hist(unp.dist)
mean.unp <- mean(unp.dist)
sd.unp <- sd(unp.dist)
length(which(unp.dist < upper & unp.dist > lower))/length(unp.dist)

x <- seq(0,1,length = 200)
y <- dnorm(x, mean =  mean.unp, sd = sd.unp)
plot(y~x, type = "l", ylim = c(0,5), bty = "l",
	main = "density distribution of propensity scores")


y <- dnorm(x, mean =  mean.p, sd = sd.p)
lines(y~x, lty = 2)

legend("topright", c("unprotected", "protected"), lty = c(1,2))



# 0.25 threshold after carranza et al
lower.25sd <- mean.p - 0.25*sd.p
upper.25sd <- mean.p + 0.25*sd.p
length(which(unp.dist < upper.25sd & unp.dist > lower.25sd))/length(unp.dist)

# 0.5 nelson and chomitz 2011
lower.5sd <- mean.p - 0.5*sd.p
upper.5sd <- mean.p + 0.5*sd.p
length(which(unp.dist < upper.5sd & unp.dist > lower.5sd))/length(unp.dist)


# 1 sd
lower1sd <- mean.p - 1*sd.p
upper1sd <- mean.p + 1*sd.p
length(which(unp.dist < upper1sd & unp.dist > lower1sd))/length(unp.dist)




#tried getting propensity scores within study, doesnt work

# for each study 
# ... harder to do as many studies only have one PA site

df <- data.frame()
for (i in 1:length(studies)){
	data.subset <- subset(matching.data, SS == studies[i])
	
	print(i)
  }

data.subset[order(data.subset$log_access),]

View(df)



# propensity matching in other papers works by country
# could also match by country - the idea is that comparison pixels are always within country
#but here, we have a random factor for study, which is pretty much always within country. so effect of PA
# is estimated first within each study, then overall. 














#############
data only within PAs
show range of values within each study
#############

names(PA.data)

PA.data <- subset(PA.data, DoP > 0)

min.area <- min(PA.data$GIS_AREA) 
max.area <- max(PA.data$GIS_AREA)
mean.area <- mean(PA.data$GIS_AREA)
median.area <- median(PA.data$GIS_AREA, na.rm = T)
hist(PA.data$GIS_AREA)

min.DoP <- min(PA.data$DoP, na.rm = T) 
max.DoP <- max(PA.data$DoP, na.rm = T)
mean.DoP <- mean(PA.data$DoP, na.rm = T)
median.DoP <- median(PA.data$DoP, na.rm = T)
hist(PA.data$DoP)

PA.stats <- data.frame(min.area = min.area, max.area = max.area, mean.area = mean.area, median.area = median.area, 
	min.DoP = min.DoP, max.DoP = max.DoP, mean.DoP = mean.DoP, median.DoP = median.DoP)
t(PA.stats)


table(PA.data$IUCN_CAT)
table(PA.data$IUCN_CAT_number)

tapply(PA.data$GIS_AREA, PA.data$Realm, median)
tapply(PA.data$GIS_AREA, PA.data$Zone, median)
tapply(PA.data$DoP, PA.data$Realm, mean, na.rm = T)
tapply(PA.data$DoP, PA.data$Zone, mean, na.rm = T)
PA.data$IUCN_CAT







# how much overlap between DoP in different studies


tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/xtra - PA DoP values - alphabetical order.tif",
	width = 15, height = 30, units = "cm", pointsize = 12, res = 300)

PA.data <- droplevels(PA.data)


DoP.data <- PA.data[,c("SS","DoP")]
DoP.data <- na.omit(DoP.data)
DoP.data <- droplevels(DoP.data)
s <- split(DoP.data, DoP.data$SS)
min.DoP <- sapply(s, function(rows) min(rows$DoP))
max.DoP <- sapply(s, function(rows) max(rows$DoP))

min.DoP <- data.frame(unlist(min.DoP))
max.DoP <- data.frame(unlist(max.DoP))


study.names <- rownames(min.DoP)

y <- seq(1,length(study.names),1)
xlim <- c(min(DoP.data$DoP, na.rm = T), max(DoP.data$DoP, na.rm = T))

par(mar = c(4.5,10,1,1))
plot(y ~ rep(30,length(y)), type = "n", bty = "l", xlim = xlim, xlab = "Duration of Protection",
		yaxt = "n", ylab = "", ylim = c(1,length(y)),
	main = "Studies listed in alphabetical order")
axis(2, y, study.names, cex.axis = 0.2, las = 1)

 z <- 5

for(z in 1:length(study.names)){
	s.name <- study.names[z]
	arrows(min.DoP$unlist.min.DoP[which(rownames(min.DoP) == s.name)],z,
		max.DoP$unlist.max.DoP[which(rownames(max.DoP) == s.name)],z, length = 0.02, angle = 90, code = 3)
	}


for(z in 1:length(study.names)){
	data.subset <- subset(DoP.data, SS == study.names[z])
	points(data.subset$DoP, rep(z, length(data.subset$DoP)), pch = 16, col = 2, cex = 0.5)
	}

dev.off()





### list in order of min DoP

tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/xtra - PA DoP values - increasing order.tif",
	width = 15, height = 30, units = "cm", pointsize = 12, res = 300)

study.names <- rownames(min.DoP)[order(min.DoP$unlist.min.DoP.)]

y <- seq(1,length(study.names),1)
xlim <- c(min(DoP.data$DoP, na.rm = T), max(DoP.data$DoP, na.rm = T))

par(mar = c(4.5,10,1,1))
plot(y ~ rep(30,length(y)), type = "n", bty = "l", xlim = xlim, xlab = "Duration of Protection",
		yaxt = "n", ylab = "", ylim = c(1,length(y)),
	main = "Studies listed in order of minimum DoP")
axis(2, y, study.names, cex.axis = 0.2, las = 1)

 z <- 5

for(z in 1:length(study.names)){
	s.name <- study.names[z]
	arrows(min.DoP$unlist.min.DoP[which(rownames(min.DoP) == s.name)],z,
		max.DoP$unlist.max.DoP[which(rownames(max.DoP) == s.name)],z, length = 0.01, angle = 90, code = 3)
	}


for(z in 1:length(study.names)){
	data.subset <- subset(DoP.data, SS == study.names[z])
	points(data.subset$DoP, rep(z, length(data.subset$DoP)), pch = 16, col = 2, cex = 0.5)
	}



dev.off()








# how much overlap between GIS AREA in different studies


which(PA.data$log_GIS_AREA == 0)



tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/xtra - PA GIS values - alphabetical order.tif",
	width = 15, height = 30, units = "cm", pointsize = 12, res = 300)


PA.data <- droplevels(PA.data)


GIS.data <- PA.data[,c("SS","log_GIS_AREA")]


GIS.data <- na.omit(GIS.data)
GIS.data <- droplevels(GIS.data)

s <- split(GIS.data, GIS.data$SS)
min.GIS <- sapply(s, function(rows) min(rows$log_GIS_AREA))
max.GIS <- sapply(s, function(rows) max(rows$log_GIS_AREA))

min.GIS <- data.frame(unlist(min.GIS))
max.GIS <- data.frame(unlist(max.GIS))


study.names <- rownames(min.GIS)#[order(min.GIS$unlist.min.GIS.)]

y <- seq(1,length(study.names),1)
xlim <- c(min(GIS.data$log_GIS_AREA, na.rm = T), max(GIS.data$log_GIS_AREA, na.rm = T))

par(mar = c(4.5,10,1,1))
plot(y ~ rep(30,length(y)), type = "n", bty = "l", xlim = xlim, xlab = "Log PA size (km2)",
		yaxt = "n", ylab = "", ylim = c(1,length(y)),
	main = "Studies listed in alphabetical order")
axis(2, y, study.names, cex.axis = 0.2, las = 1)

 z <- 5

for(z in 1:length(study.names)){
	s.name <- study.names[z]
	arrows(min.GIS$unlist.min.GIS[which(rownames(min.GIS) == s.name)],z,
		max.GIS$unlist.max.GIS[which(rownames(max.GIS) == s.name)],z, length = 0.02, angle = 90, code = 3)
	}


for(z in 1:length(study.names)){
	data.subset <- subset(GIS.data, SS == study.names[z])
	points(data.subset$log_GIS_AREA, rep(z, length(data.subset$log_GIS_AREA)), pch = 16, col = 2, cex = 0.5)
	}

dev.off()




### list in order of min GIS AREA
study.names <- rownames(min.GIS)[order(min.GIS$unlist.min.GIS.)]



tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/xtra - PA GIS values - increasing order.tif",
	width = 15, height = 30, units = "cm", pointsize = 12, res = 300)



par(mar = c(4.5,10,1,1))
plot(y ~ rep(30,length(y)), type = "n", bty = "l", xlim = xlim, xlab = "log PA size (km2)",
		yaxt = "n", ylab = "", ylim = c(1,length(y)),
	main = "Studies listed in order of minimum GIS")
axis(2, y, study.names, cex.axis = 0.2, las = 1)

 z <- 1

for(z in 1:length(study.names)){
	s.name <- study.names[z]
	arrows(min.GIS$unlist.min.GIS[which(rownames(min.GIS) == s.name)],z,
		max.GIS$unlist.max.GIS[which(rownames(max.GIS) == s.name)],z, length = 0.01, angle = 90, code = 3)
	}


for(z in 1:length(study.names)){
	data.subset <- subset(GIS.data, SS == study.names[z])
	points(data.subset$log_GIS_AREA, rep(z, length(data.subset$log_GIS_AREA)), pch = 16, col = 2, cex = 0.5)
	}

dev.off()






### list in order of min log dist boundary

tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/xtra - PA multi lbd values - increasing order.tif",
	width = 15, height = 30, units = "cm", pointsize = 12, res = 300)



PA.data <- droplevels(PA.data)


lbd.data <- PA.data[,c("SS","log_bound_dist_km")]
lbd.data <- na.omit(lbd.data)
lbd.data <- droplevels(lbd.data)
s <- split(lbd.data, lbd.data$SS)
min.lbd <- sapply(s, function(rows) min(rows$log_bound_dist_km))
max.lbd <- sapply(s, function(rows) max(rows$log_bound_dist_km))

min.lbd <- data.frame(unlist(min.lbd))
max.lbd <- data.frame(unlist(max.lbd))


study.names <- rownames(min.lbd)[order(min.lbd$unlist.min.lbd)]

y <- seq(1,length(study.names),1)
xlim <- c(min(lbd.data$log_bound_dist_km, na.rm = T), max(lbd.data$log_bound_dist_km, na.rm = T))

par(mar = c(4.5,10,1,1))
plot(y ~ rep(30,length(y)), type = "n", bty = "l", xlim = xlim, xlab = "Log distance to boundary (km)",
		yaxt = "n", ylab = "", ylim = c(1,length(y)),
	main = "Studies listed in order of minimum lbd")
axis(2, y, study.names, cex.axis = 0.2, las = 1)

 z <- 5

for(z in 1:length(study.names)){
	s.name <- study.names[z]
	arrows(min.lbd$unlist.min.lbd[which(rownames(min.lbd) == s.name)],z,
		max.lbd$unlist.max.lbd[which(rownames(max.lbd) == s.name)],z, length = 0.01, angle = 90, code = 3)
	}


for(z in 1:length(study.names)){
	data.subset <- subset(lbd.data, SS == study.names[z])
	points(data.subset$log_bound_dist_km, rep(z, length(data.subset$log_bound_dist_km)), pch = 16, col = 2, cex = 0.5)
	}




dev.off()










### list in order of min hpd

tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/xtra - PA hpd values - increasing order.tif",
	width = 15, height = 30, units = "cm", pointsize = 12, res = 300)



PA.data <- droplevels(PA.data)


hpd.data <- PA.data[,c("SS","log_hpd")]
hpd.data <- na.omit(hpd.data)
hpd.data <- droplevels(hpd.data)
s <- split(hpd.data, hpd.data$SS)
min.hpd <- sapply(s, function(rows) min(rows$log_hpd))
max.hpd <- sapply(s, function(rows) max(rows$log_hpd))

min.hpd <- data.frame(unlist(min.hpd))
max.hpd <- data.frame(unlist(max.hpd))


study.names <- rownames(min.hpd)[order(min.hpd$unlist.min.hpd)]

y <- seq(1,length(study.names),1)
xlim <- c(min(hpd.data$log_hpd, na.rm = T), max(hpd.data$log_hpd, na.rm = T))

par(mar = c(4.5,10,1,1))
plot(y ~ rep(30,length(y)), type = "n", bty = "l", xlim = xlim, xlab = "Duration of Protection",
		yaxt = "n", ylab = "", ylim = c(1,length(y)),
	main = "Studies listed in order of minimum hpd")
axis(2, y, study.names, cex.axis = 0.2, las = 1)

 z <- 5

for(z in 1:length(study.names)){
	s.name <- study.names[z]
	arrows(min.hpd$unlist.min.hpd[which(rownames(min.hpd) == s.name)],z,
		max.hpd$unlist.max.hpd[which(rownames(max.hpd) == s.name)],z, length = 0.01, angle = 90, code = 3)
	}


for(z in 1:length(study.names)){
	data.subset <- subset(hpd.data, SS == study.names[z])
	points(data.subset$log_hpd, rep(z, length(data.subset$log_hpd)), pch = 16, col = 2, cex = 0.5)
	}



dev.off()








### list in order of min hpd

tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/xtra - PA hpd values - increasing order.tif",
	width = 15, height = 30, units = "cm", pointsize = 12, res = 300)



PA.data <- droplevels(PA.data)


hpd.data <- PA.data[,c("SS","log_hpd")]
hpd.data <- na.omit(hpd.data)
hpd.data <- droplevels(hpd.data)
s <- split(hpd.data, hpd.data$SS)
min.hpd <- sapply(s, function(rows) min(rows$log_hpd))
max.hpd <- sapply(s, function(rows) max(rows$log_hpd))

min.hpd <- data.frame(unlist(min.hpd))
max.hpd <- data.frame(unlist(max.hpd))


study.names <- rownames(min.hpd)[order(min.hpd$unlist.min.hpd)]

y <- seq(1,length(study.names),1)
xlim <- c(min(hpd.data$log_hpd, na.rm = T), max(hpd.data$log_hpd, na.rm = T))

par(mar = c(4.5,10,1,1))
plot(y ~ rep(30,length(y)), type = "n", bty = "l", xlim = xlim, xlab = "Log hpd",
		yaxt = "n", ylab = "", ylim = c(1,length(y)),
	main = "Studies listed in order of minimum hpd")
axis(2, y, study.names, cex.axis = 0.2, las = 1)

 z <- 5

for(z in 1:length(study.names)){
	s.name <- study.names[z]
	arrows(min.hpd$unlist.min.hpd[which(rownames(min.hpd) == s.name)],z,
		max.hpd$unlist.max.hpd[which(rownames(max.hpd) == s.name)],z, length = 0.01, angle = 90, code = 3)
	}


for(z in 1:length(study.names)){
	data.subset <- subset(hpd.data, SS == study.names[z])
	points(data.subset$log_hpd, rep(z, length(data.subset$log_hpd)), pch = 16, col = 2, cex = 0.5)
	}



dev.off()







### list in order of min access

tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/xtra - PA access values - increasing order.tif",
	width = 15, height = 30, units = "cm", pointsize = 12, res = 300)



PA.data <- droplevels(PA.data)


access.data <- PA.data[,c("SS","log_access")]
access.data <- na.omit(access.data)
access.data <- droplevels(access.data)
s <- split(access.data, access.data$SS)
min.access <- sapply(s, function(rows) min(rows$log_access))
max.access <- sapply(s, function(rows) max(rows$log_access))

min.access <- data.frame(unlist(min.access))
max.access <- data.frame(unlist(max.access))


study.names <- rownames(min.access)[order(min.access$unlist.min.access)]

y <- seq(1,length(study.names),1)
xlim <- c(min(access.data$log_access, na.rm = T), max(access.data$log_access, na.rm = T))

par(mar = c(4.5,10,1,1))
plot(y ~ rep(30,length(y)), type = "n", bty = "l", xlim = xlim, xlab = "log access",
		yaxt = "n", ylab = "", ylim = c(1,length(y)),
	main = "Studies listed in order of minimum access")
axis(2, y, study.names, cex.axis = 0.2, las = 1)

 z <- 5

for(z in 1:length(study.names)){
	s.name <- study.names[z]
	arrows(min.access$unlist.min.access[which(rownames(min.access) == s.name)],z,
		max.access$unlist.max.access[which(rownames(max.access) == s.name)],z, length = 0.01, angle = 90, code = 3)
	}


for(z in 1:length(study.names)){
	data.subset <- subset(access.data, SS == study.names[z])
	points(data.subset$log_access, rep(z, length(data.subset$log_access)), pch = 16, col = 2, cex = 0.5)
	}



dev.off()







### list in order of min ag_suit

tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/xtra - PA ag_suit values - increasing order.tif",
	width = 15, height = 30, units = "cm", pointsize = 12, res = 300)



PA.data <- droplevels(PA.data)


ag_suit.data <- PA.data[,c("SS","ag_suit")]
ag_suit.data <- na.omit(ag_suit.data)
ag_suit.data <- droplevels(ag_suit.data)
s <- split(ag_suit.data, ag_suit.data$SS)
min.ag_suit <- sapply(s, function(rows) min(rows$ag_suit))
max.ag_suit <- sapply(s, function(rows) max(rows$ag_suit))

min.ag_suit <- data.frame(unlist(min.ag_suit))
max.ag_suit <- data.frame(unlist(max.ag_suit))


study.names <- rownames(min.ag_suit)[order(min.ag_suit$unlist.min.ag_suit)]

y <- seq(1,length(study.names),1)
xlim <- c(min(ag_suit.data$ag_suit, na.rm = T), max(ag_suit.data$ag_suit, na.rm = T))

par(mar = c(4.5,10,1,1))
plot(y ~ rep(30,length(y)), type = "n", bty = "l", xlim = xlim, xlab = "ag_suit",
		yaxt = "n", ylab = "", ylim = c(1,length(y)),
	main = "Studies listed in order of minimum ag_suit")
axis(2, y, study.names, cex.axis = 0.2, las = 1)

 z <- 5

for(z in 1:length(study.names)){
	s.name <- study.names[z]
	arrows(min.ag_suit$unlist.min.ag_suit[which(rownames(min.ag_suit) == s.name)],z,
		max.ag_suit$unlist.max.ag_suit[which(rownames(max.ag_suit) == s.name)],z, length = 0.01, angle = 90, code = 3)
	}


for(z in 1:length(study.names)){
	data.subset <- subset(ag_suit.data, SS == study.names[z])
	points(data.subset$ag_suit, rep(z, length(data.subset$ag_suit)), pch = 16, col = 2, cex = 0.5)
	}



dev.off()








### how many studies have more than one IUCN CAT in this dataset - 


p <- split(PA.multi, PA.multi$SS)
n.IUCN <- lapply(p, function(rows) length(unique(na.omit(rows$IUCN_CAT_number))))
ni <- names(which(n.IUCN > 1)) #36 studies, of 62
length(unique(PA.multi$SS))

x <- subset(PA.multi, SS %in% ni)
nrow(x)
View(x)





### whats the data spread in PA.multi

length(unique(PA.multi$SS))


veg <- subset(PA.multi, taxon_of_interest == "Plants and Fungi")
length(unique(veg$SS))





### are there correlations between the different PA characteristics 


plot(range ~ access, PA.multi, log = "x")
plot(log_GIS_AREA ~ IUCN_CAT_number, PA.multi)
plot(DoP ~ IUCN_CAT_number, PA.multi)
plot(DoP ~ log_GIS_AREA, PA.multi)

M1 <- lmer(log_GIS_AREA ~ IUCN_CAT_number + (IUCN_CAT_number|Country), PA.multi)
M2 <- lmer(log_GIS_AREA ~ 1 + (IUCN_CAT_number|Country), PA.multi)
anova(M1, M2)

M1 <- lmer(DoP ~ IUCN_CAT_number + (IUCN_CAT_number|Country), PA.multi)
M2 <- lmer(DoP ~ 1 + (IUCN_CAT_number|Country), PA.multi)
anova(M1, M2)

M1 <- lmer(DoP ~ log_GIS_AREA + (log_GIS_AREA|Country), PA.multi)
M2 <- lmer(DoP ~ 1 + (log_GIS_AREA|Country), PA.multi)
anova(M1, M2)
summary(M1)

