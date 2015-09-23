
rm(list=ls()) 
library(rgdal)
library(yarg)
library(roquefort)
library(gamm4)

setwd("")
source("compare_randoms.R")
source("model_select.R")
setwd("R:/ecocon_d/clg32/PREDICTS/WDPA analysis")
PREDICTS_WDPA <- read.csv("PREDICTS_WDPA.csv")


matched.landuse <- subset(PREDICTS_WDPA, matched.landuse == "yes")
nrow(matched.landuse) #5015


### sending these to arc for isolating corresponding PAs
# create shapefile for PA_11_14 data
long.lat <- PA_11_14[,c("Longitd", "Latitud")]
spatial <- SpatialPoints(long.lat)
plot(spatial)

spdf <- SpatialPointsDataFrame(spatial, PA_11_14)
#writeOGR(obj = spdf, dsn = "C:/GIS/PA_predicts_mapping", "PA_11_2014_no_ind_sec_veg", driver = "ESRI Shapefile")


# create shapefile for matched.landuse data
long.lat <- matched.landuse[,c("Longitd", "Latitud")]
spatial <- SpatialPoints(long.lat)
plot(spatial)


spdf <- SpatialPointsDataFrame(spatial,matched.landuse)
#writeOGR(obj = spdf, dsn = "C:/GIS/PA_predicts_mapping", "matched.landuse", driver = "ESRI Shapefile")


# create shapefile of sites in PA_11_14 but not in matched landuse
additional <- PA_11_14[which(PA_11_14$SSS %in% matched.landuse$SSS == F),]
nrow(additional)
long.lat <- additional[,c("Longitd", "Latitud")]
spatial <- SpatialPoints(long.lat)

spdf <- SpatialPointsDataFrame(spatial,additional)
#writeOGR(obj = spdf, dsn = "C:/GIS/PA_predicts_mapping", "additional", driver = "ESRI Shapefile")


# equal plot
countries <- readOGR(dsn = "C:/GIS/Data/world_borders", "BORDERS_without_antarctica_equl")
#this is the publicly accessible data layer at http://thematicmapping.org/downloads/world_borders.php
points_ML <- readOGR(dsn = "C:/GIS/PA_predicts_mapping", "matched.landuse")
points_additional <- readOGR(dsn = "C:/GIS/PA_predicts_mapping", "additional")
lines <- readOGR("C:/GIS/Data/world_borders", "equator_tropics_equl")
PAs_all <- readOGR(dsn = "C:/GIS/Data/WDPA/WDPAmollmergeJULY", "WDPAmollmergeJULY_terrestrial_equl")
#this is the publicly accessible world database on protected areas at http://www.protectedplanet.net/



tiff("R:/ecocon_d/clg32/PREDICTS/WDPA analysis/plots/06_15/PA_11_14 map equal.tif",
	width = 25, height =21, units = "cm", pointsize = 12, res = 300)

par(mfrow = c(1,1),mar = c(0,0,0,0), oma=c(0, 0, 0, 0),
    pty = "m",
    xpd = NA)

PA.col <- rgb(0.66,0.86,0.6)

plot(countries, lty = NULL, border = "grey", col = "grey")
plot(PAs_all, lty = NULL, border = PA.col, col = PA.col, add = T, lwd = 0.01)
plot(lines, add = T, lty = 2)
plot(points_additional, pch = 16, col = rgb(0.1,0.1,0.1,0.5), add = T, cex = 1)
plot(points_ML, pch = 16, col = rgb(0.1,0.1,0.1,0.5), add = T, cex = 1)


dev.off()

tiff("R:/ecocon_d/clg32/PREDICTS/WDPA analysis/plots/02_15/PA_11_14 map legend.tif",
	width = 15, height =15, units = "cm", pointsize = 12, res = 300)
plot(1~1, type = "n", axes = F, xlab = "", ylab = "")
legend(0.75,1, cex = 1,
	pch = c(16, 21,15), 
	c("sites matched by land use", "sites not matched by land use", "protected areas"), 
	pt.bg = rgb(0.9,0.9,0.9,0.5), col = c(rgb(0.1,0.1,0.1,0.5),rgb(0.1,0.1,0.1,0.5),PA.col))

dev.off()



# histogram of points

cols <- brewer.pal(8, "Paired")
taxa.cols <- cols[c(4,2,8)]
taxa.cols.ci <- c("#33A02C44", "#1F78B444", "#FF7F0044")
taxa <- c("Plants", "Invertebrates", "Vertebrates")

n.breaks <- 30

min(matched.landuse$Latitud)
max(matched.landuse$Latitud)

#want breaks that fall between 0 and 23.5 neatly
step <- 23.5/4
breaks <- seq(-8*step, 12*step, step)

ML_verts_cut <- hist(matched.landuse$Latitud[which(matched.landuse$taxon_of_interest == "Vertebrates")], breaks)
ML_inverts_cut <- hist(matched.landuse$Latitud[which(matched.landuse$taxon_of_interest == "Invertebrates")], breaks)
ML_plants_cut <- hist(matched.landuse$Latitud[which(matched.landuse$taxon_of_interest == "Plants")], breaks)

LUPA_verts_cut <- hist(PA_11_14$Latitud[which(PA_11_14$taxon_of_interest == "Vertebrates")], breaks)
LUPA_inverts_cut <- hist(PA_11_14$Latitud[which(PA_11_14$taxon_of_interest == "Invertebrates")], breaks)
LUPA_plants_cut <- hist(PA_11_14$Latitud[which(PA_11_14$taxon_of_interest == "Plants")], breaks)

LUPA_verts_cut$mids


tiff( "R:/ecocon_d/clg32/PREDICTS/WDPA analysis/plots/06_15/latitudinal bar plot.tif",
	width = 16, height =14, units = "cm", pointsize = 18, res = 300)

par(mfrow = c(1,2))
#add the zeros at the ends where breaks not made
LUPA_cut <- rbind(LUPA_verts_cut$counts, c(LUPA_inverts_cut$counts), c(LUPA_plants_cut$counts))
ML_cut <- rbind(c(ML_verts_cut$counts), c(ML_inverts_cut$counts), c(ML_plants_cut$counts))

colSums(LUPA_cut)
par(mar = c(4,3,1,0), oma = c(1,1,1,1))
barplot(LUPA_cut, beside = F, horiz = T, space = 0, col = rev(taxa.cols), border = NA,
	#main = "all sites", cex.main = 1.2,
	ylim = c(0,20), 
	xlim = c(1550,0), cex.axis = 1, las = 2, xaxt = "n")
axis(1, c(0,500,1000,1500), c(0,500,1000,1500), las = 2) 
abline(h = (0 - c)/step, col = 1, lty = 2)
abline(h = (23.5 - c)/step, col = 1, lty = 2)
abline(h = (-23.5 - c)/step, col = 1, lty = 2)

# get axis positions for integer values
# intercept is lowest break, gradient = step
c <- -8*step
x <- seq(0,20,1)
#check this gives original breaks
y = step*x + c
y == breaks

# reverse formula
#(y - c)/step = x
# so if
y <- c(-50, -25, 0, 25, 50, 75)
new.x <- (y - c)/step

axis(2, new.x, y, cex.axis = 1, las = 2)

par(mar = c(4,0,1,3))
barplot(ML_cut, beside = F, horiz = T, space = 0, col = rev(taxa.cols), border = NA,
	#main = "sites matched \n by landuse",  cex.main = 1.2, 
	ylim = c(0,20), xlim = c(0,1550), cex.axis = 1, las = 2, xaxt = "n")
axis(1, c(0,500,1000,1500), c(0,500,1000,1500), las = 2) 
abline(h = (0 - c)/step, col = 1, lty = 2)
abline(h = (23.5 - c)/step, col = 1, lty = 2)
abline(h = (-23.5 - c)/step, col = 1, lty = 2)

abline(v = 0, col = 1)

dev.off()

LU1 <- rgb(0.4,0.4,0.4)
LU2 <- rgb(0.7,0.7,0.7)
LU3 <- rgb(0.95,0.95,0.95)

tiff("R:/ecocon_d/clg32/PREDICTS/WDPA analysis/plots/02_15/latitudinal histogram legend.tif",
	width = 15, height =15, units = "cm", pointsize = 12, res = 300)
plot(1~1, type = "n", axes = F, xlab = "", ylab = "")
legend(0.75,1, cex = 1,
	c("Land use X", "Land use Y", "Land use Z",
	"Protected areas", "Sites matched by land use", "Sites not matched by land use", 
	"Plants", "Invertebrates", "Vertebrates"), 
	col = c(LU1,LU2,LU3,PA.col, 1,1, taxa.cols), 
	pch = c(15,15,15,15, 16, 21, 15, 15,15), bty = "n")

dev.off()

