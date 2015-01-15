#Distances between sites

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




# load functions

source("N:/Documents/PREDICTS/WDPA analysis/compare_randoms_pairwise.R")
source("N:/Documents/PREDICTS/WDPA analysis/get_pairwise_independent.R")
source("N:/Documents/PREDICTS/WDPA analysis/multiplot_pairwise_all.R")
source("N:/Documents/PREDICTS/WDPA analysis/bray_curtis_dissimilarity.R")
source("N:/Documents/PREDICTS/WDPA analysis/model_select_pairwise.R")


# Load dataset on taxa_split matched all


matched.landuse<- read.csv("matched.landuse_07_2014.csv")




# get distances


pairwise_landuse_abundance <- get_pairwise(matched.landuse, "Total_abundance")

# how many dists = 0 
which(pairwise_landuse_abundance$dist_km == 0)
#none. great


min(pairwise_landuse_abundance$dist_km)
max(pairwise_landuse_abundance$dist_km)
mean(pairwise_landuse_abundance$dist_km)
median(pairwise_landuse_abundance$dist_km)
quantile(pairwise_landuse_abundance$dist_km, c(0.01, 0.05, 0.1, .25, 0.5, 0.75, 0.98))
pairwise_landuse_abundance$SS[which(pairwise_landuse_abundance$dist_km == max(pairwise_landuse_abundance$dist_km))]
View(subset(pairwise_landuse_abundance, SS == "HP1_2008__Grogan 1"))

plot <- qplot(dist_km, data=pairwise_landuse_abundance, binwidth = 5, geom="histogram")
plot + theme_bw()