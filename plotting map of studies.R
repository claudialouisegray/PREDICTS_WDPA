
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
source("compare_randoms.R")
source("model_select.R")


#load data for matched studies
setwd("R:/ecocon_d/clg32/GitHub/PREDICTS_WDPA")
source("prep_matched.landuse_for_analysis.R")
data <- matched.landuse

# or if using LUPA data
setwd("R:/ecocon_d/clg32/GitHub/PREDICTS_WDPA")
source("prep_PA_11_14_for_analysis.R")
data <- PA_11_14
nrow(PA_11_14)


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

# map these
countries <- readOGR(dsn = "C:/GIS/Data/world_borders", "BORDERS_without_antarctica_moll")
PAs_all <- readOGR(dsn = "C:/GIS/Data/WDPA/WDPAmollmergeJULY", "WDPAmollmergeJULY_terrestrial")
#PAs_LUPA <- readOGR(dsn = "C:/GIS/PA_predicts_mapping", "PA_11_2014_no_ind_sec_veg_PAs_moll")
#PAs_ML <- readOGR(dsn = "C:/GIS/PA_predicts_mapping", "matched.landuse_PAs_moll")
points_LUPA <- readOGR(dsn = "C:/GIS/PA_predicts_mapping", "PA_11_2014_no_ind_sec_veg_moll")
points_ML <- readOGR(dsn = "C:/GIS/PA_predicts_mapping", "matched.landuse_moll")
points_additional <- readOGR(dsn = "C:/GIS/PA_predicts_mapping", "additional_moll")
lines <- readOGR("C:/GIS/Data/world_borders", "equator_tropics_moll")



tiff("N:/Documents/PREDICTS/WDPA analysis/plots/02_15/PA_11_14 map.tif",
	width = 25, height =20, units = "cm", pointsize = 12, res = 300)
par(mfrow = c(1,1))
par(mar=c(0.2,0.1,1.5,0.2))

plot(countries, lty = NULL, border = "grey", col = "grey")
plot(PAs_all, lty = NULL, border = rgb(0.4,0.7,0.5), col = rgb(0.4,0.7,0.5), add = T, lwd = 0.01)
plot(lines, add = T, lty = 2)
#plot(points_LUPA, pch = 21, col = rgb(0.1,0.1,0.1,0.5), bg = rgb(0.9,0.9,0.9,0.5), add = T, cex = 1)
plot(points_additional, pch = 21, col = rgb(0.1,0.1,0.1,0.5), bg = rgb(0.9,0.9,0.9,0.5), add = T, cex = 1)
plot(points_ML, pch = 16, col = rgb(0.1,0.1,0.1,0.5), add = T, cex = 1)

dev.off()

tiff("N:/Documents/PREDICTS/WDPA analysis/plots/02_15/PA_11_14 map legend.tif",
	width = 15, height =15, units = "cm", pointsize = 12, res = 300)
plot(1~1, type = "n", axes = F, xlab = "", ylab = "")
legend(0.75,1, cex = 1,
	pch = c(16, 21,15), 
	c("sites matched by land use", "sites not matched by land use", "protected areas"), 
	pt.bg = rgb(0.9,0.9,0.9,0.5), col = c(rgb(0.1,0.1,0.1,0.5),rgb(0.1,0.1,0.1,0.5),rgb(0.4,0.7,0.5)))

dev.off()

### make border thickness less



# histogram of points

cols <- brewer.pal(8, "Paired")
taxa.cols <- cols[c(4,2,8)]
taxa.cols.ci <- c("#33A02C44", "#1F78B444", "#FF7F0044")
taxa <- c("Plants", "Invertebrates", "Vertebrates")

n.breaks <- 30

ML_verts_cut <- hist(matched.landuse$Latitud[which(matched.landuse$taxon_of_interest == "Vertebrates")], n.breaks)
ML_inverts_cut <- hist(matched.landuse$Latitud[which(matched.landuse$taxon_of_interest == "Invertebrates")], n.breaks)
ML_plants_cut <- hist(matched.landuse$Latitud[which(matched.landuse$taxon_of_interest == "Plants")], n.breaks)


ML_verts_cut$breaks
ML_inverts_cut$breaks
ML_plants_cut$breaks

LUPA_verts_cut <- hist(PA_11_14$Latitud[which(PA_11_14$taxon_of_interest == "Vertebrates")], n.breaks)
LUPA_inverts_cut <- hist(PA_11_14$Latitud[which(PA_11_14$taxon_of_interest == "Invertebrates")], n.breaks)
LUPA_plants_cut <- hist(PA_11_14$Latitud[which(PA_11_14$taxon_of_interest == "Plants")], n.breaks)

LUPA_verts_cut$breaks
LUPA_inverts_cut$breaks
LUPA_plants_cut$breaks

LUPA_verts_cut$mids


tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/02_15/latitudinal bar plot.tif",
	width = 10, height =20, units = "cm", pointsize = 12, res = 300)

par(mfrow = c(1,2))
#add the zeros at the ends where breaks not made
LUPA_cut <- rbind(LUPA_verts_cut$counts, c(0,LUPA_inverts_cut$counts,0), c(0,LUPA_plants_cut$counts,0,0))
ML_cut <- rbind(c(0,ML_verts_cut$counts), c(0,ML_inverts_cut$counts,0), c(0,ML_plants_cut$counts,0,0))

colSums(LUPA_cut)
par(mar = c(4,2,2,0))
barplot(LUPA_cut, beside = F, horiz = T, space = 0, col = rev(taxa.cols), border = NA,
	main = "all sites", cex.main = 1.2,
	ylim = c(-6,26), xlim = c(1350,0), cex.axis = 1, las = 2)
#all breaks = c(seq(-75,-55,5), LUPA_verts_cut$breaks), at -5:25
axis(2, seq(-4,24,2), seq(-70,70,10), cex.axis = 1)
par(mar = c(4,0,2,2))
barplot(ML_cut, beside = F, horiz = T, space = 0, col = rev(taxa.cols), border = NA,
	main = "sites matched \n by landuse",  cex.main = 1.2, 
	ylim = c(-6,26), xlim = c(0,1350), cex.axis = 1, las = 2)
#axis(2, -6:23, c(seq(-75,-50,5) ML_verts_cut$breaks))
abline(v = 0, col = 1)

dev.off()



tiff("N:/Documents/PREDICTS/WDPA analysis/plots/02_15/latitudinal histogram legend.tif",
	width = 15, height =15, units = "cm", pointsize = 12, res = 300)
plot(1~1, type = "n", axes = F, xlab = "", ylab = "")
legend(0.75,1, cex = 1,
	c("Plants", "Invertebrates", "Vertebrates"), 
	col = taxa.cols, pch = 15)

dev.off()



