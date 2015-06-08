


### LOAD DATA on confounding variables



setwd("R:/ecocon_d/clg32/PREDICTS/WDPA analysis/matching data")


access <- read.table("bfer_1km_acc50k_moll_11_14.txt", header = T, sep = ",")
hpd <- read.table("bfer_1km_HPD_moll_11_14.txt", header = T, sep = ",")
elevation <- read.table("bfer_1km_mn30_elevation_moll_11_14.txt", header = T, sep = ",")
slope <- read.table("bfer_1km_mn30_slope_moll_11_14.txt", header = T, sep = ",")


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





### LOAD DIST TO PA BOUNDARY DATA ####




setwd("R:/ecocon_d/clg32/PREDICTS/WDPA analysis")



#PA_dists <- read.csv("matched_landuse_07_2014_dists_to_PA_boundary_same_habitat.csv")


PA_dists <- read.csv("PA_11_2014_dists_to_PA_boundary_same_habitat.csv")


# this is exported from the table of the shapefile in arcMAP, 
# then the blank cells (" ") replaced with NAs and FID (empty) column removed in excel
# -1 values for near dist (when PA in same habitat is >500000 m away) to NA

#get columns needed
PA_dists_data <- PA_dists[,c("SSS", "NEAR_DIST")]










### Load dataset on matched landuse

matched.landuse<- read.csv("matched.landuse_11_2014.csv")
nrow(matched.landuse)

#make plants the reference
matched.landuse$taxon_of_interest <- relevel(matched.landuse$taxon_of_interest, "Plants")


matched.landuse$Predominant_habitat <- gsub("Primary forest", "Primary Vegetation", matched.landuse$Predominant_habitat)
matched.landuse$Predominant_habitat <- gsub("Primary non-forest", "Primary Vegetation", matched.landuse$Predominant_habitat)
#matched.landuse$Predominant_habitat <- gsub("Young secondary vegetation", "Secondary Vegetation",matched.landuse$Predominant_habitat)
#matched.landuse$Predominant_habitat <- gsub("Mature secondary vegetation", "Secondary Vegetation",matched.landuse$Predominant_habitat)
#matched.landuse$Predominant_habitat <- gsub("Intermediate secondary vegetation","Secondary Vegetation", matched.landuse$Predominant_habitat)
#matched.landuse$Predominant_habitat[which(matched.landuse$Predominant_habitat == "Secondary vegetation (indeterminate age)")] <- "Secondary Vegetation"

matched.landuse <- subset(matched.landuse, Predominant_habitat != "Secondary vegetation (indeterminate age)")

matched.landuse$Predominant_habitat <- factor(matched.landuse$Predominant_habitat)

nrow(matched.landuse)#5491 to 5015 after sec veg split and indet. sec veg dropped

m <- merge(matched.landuse, access.1, "SSS", all.x = T)
colnames(m)[which(colnames(m) == "MEAN")] <- "access"
m <- merge(m, hpd.1, "SSS", all.x = T)
colnames(m)[which(colnames(m) == "MEAN")] <- "hpd"
m <- merge(m, elevation.1, "SSS", all.x = T)
colnames(m)[which(colnames(m) == "MEAN")] <- "elevation"
m <- merge(m, slope.1, "SSS", all.x = T)
colnames(m)[which(colnames(m) == "MEAN")] <- "slope"
m <- merge(m, ag.1, "SSS", all.x = T)

nrow(m)

matched.landuse <- m

length(matched.landuse[,1]) 
names(matched.landuse)



#merge distance to boundary data onto matched landuse
nrow(PA_dists_data)
nrow(matched.landuse)

m <- merge(matched.landuse, PA_dists_data, by = "SSS")
nrow(m)
matched.landuse <- m


matched.landuse$bound_dist_km <- matched.landuse$NEAR_DIST/1000
matched.landuse$log_bound_dist_km <- log(matched.landuse$bound_dist_km+1)

#make variables where the sign is negative inside PAs
matched.landuse$bound_dist_km_PA_neg <- matched.landuse$bound_dist_km
matched.landuse$bound_dist_km_PA_neg[which(matched.landuse$Within_PA == "yes")] <- -1*matched.landuse$bound_dist_km_PA_neg[which(matched.landuse$Within_PA == "yes")]
matched.landuse$log_bound_dist_km_PA_neg <- matched.landuse$log_bound_dist_km
matched.landuse$log_bound_dist_km_PA_neg[which(matched.landuse$Within_PA == "yes")] <- -1*matched.landuse$log_bound_dist_km_PA_neg[which(matched.landuse$Within_PA == "yes")]





# sort out explanatory variables

matched.landuse$IUCN_CAT_number <- factor(matched.landuse$IUCN_CAT_number) # they arent really in an order

matched.landuse$log_slope <- log(matched.landuse$slope +1)
matched.landuse$log_elevation <- log(matched.landuse$elevation +1)
matched.landuse$log_hpd<- log(matched.landuse$hpd +1)
matched.landuse$log_access <- log(matched.landuse$access +1)

matched.landuse$DoP.PA <- matched.landuse$DoP
matched.landuse$DoP.PA[which(matched.landuse$Within_PA == "no")] <- 0
# make bins
# based on violin plots try 0 - 10, 11 - 20, 21 - 40, 41 - 80
matched.landuse$DoP.PA.f <- matched.landuse$DoP.PA
matched.landuse$DoP.PA.f[which(matched.landuse$DoP.PA >=0 & matched.landuse$DoP.PA <= 20)] <- "young"
matched.landuse$DoP.PA.f[which(matched.landuse$DoP.PA >= 21 & matched.landuse$DoP.PA <= max(matched.landuse$DoP.PA, na.rm = T))] <- "old"
matched.landuse$DoP.PA.f[which(matched.landuse$Within_PA == "no")] <- "outside"
matched.landuse$DoP.PA.f <- factor(matched.landuse$DoP.PA.f)
matched.landuse$DoP.PA.f <- relevel(matched.landuse$DoP.PA.f, "outside")

matched.landuse$AREA.PA <- matched.landuse$GIS_AREA
matched.landuse$AREA.PA[which(matched.landuse$Within_PA == "no")] <- 0
# make bins
matched.landuse$AREA.PA.f <- matched.landuse$AREA.PA
matched.landuse$AREA.PA.f[which(matched.landuse$AREA.PA >= 0 & matched.landuse$AREA.PA <= 400)] <- "small"
matched.landuse$AREA.PA.f[which(matched.landuse$AREA.PA >= 401 & matched.landuse$AREA.PA <= max(matched.landuse$AREA.PA, na.rm = T))] <- "large"
matched.landuse$AREA.PA.f[which(matched.landuse$Within_PA == "no")] <- "outside"
matched.landuse$AREA.PA.f <- factor(matched.landuse$AREA.PA.f)
matched.landuse$AREA.PA.f <- relevel(matched.landuse$AREA.PA.f, "outside")

#make combined interaction term variable
matched.landuse$AREA_DoP <- paste(matched.landuse$AREA.PA.f, matched.landuse$DoP.PA.f, sep = "_")
matched.landuse$AREA_DoP[which(matched.landuse$Within_PA == "no")] <- "outside"
matched.landuse$AREA_DoP[is.na(matched.landuse$DoP.PA)] <- "outside"
matched.landuse$AREA_DoP <- factor(matched.landuse$AREA_DoP)
matched.landuse$AREA_DoP <- relevel(matched.landuse$AREA_DoP, "outside")


#check outside numbers are the same
#aggregate(SSS~ AREA.PA.f + taxon_of_interest, matched.landuse, length)
#aggregate(SSS~ DoP.PA.f + taxon_of_interest, matched.landuse, length)
table(matched.landuse$AREA_DoP, matched.landuse$taxon_of_interest)
table(matched.landuse$AREA_DoP, matched.landuse$Zone)



matched.landuse$log_GIS_AREA <- log(matched.landuse$GIS_AREA+1)
matched.landuse$log_AREA.PA <- matched.landuse$log_GIS_AREA
matched.landuse$log_AREA.PA[which(matched.landuse$Within_PA == "no")] <- 0


# make IUCN cat variable of I or II vs III to VI vs unknown vs unprotected
matched.landuse$IUCN_CAT_number <- factor(matched.landuse$IUCN_CAT_number) # they arent really in an order
matched.landuse$IUCN_CAT <- matched.landuse$IUCN_CAT_number 
levels(matched.landuse$IUCN_CAT) <- c(levels(matched.landuse$IUCN_CAT), "0")
matched.landuse$IUCN_CAT[which(matched.landuse$Within_PA == "no")] <- 0


matched.landuse$IUCN.PA <- matched.landuse$IUCN_CAT_number
levels(matched.landuse$IUCN.PA) <- c( "1.5", "4.5", "7", "0")
matched.landuse$IUCN.PA[which(matched.landuse$Within_PA == "no")] <- "0"


#make response variables
matched.landuse$log_abundance <- log(matched.landuse$Total_abundance +1)
matched.landuse$range <- matched.landuse$CWM_Geographic_range_log10_square_km
matched.landuse$mass <- matched.landuse$CWM_Mass_log10_g  
matched.landuse$veg <- matched.landuse$CWM_Vegetative_height_log10_m
matched.landuse$vol <- matched.landuse$CWM_Length_derived_volume_3log10_mm

pos <- which(matched.landuse$CWM_Adult_wet_mass_log10_g >0 & matched.landuse$CWM_Maximum_wet_mass_log10_g >0) 
length(pos)






### CREATE MULTIPLE TAXA PER STUDY DATASET 
# create dataset for species richness analysis that without all studies that are only on one taxon




studies.taxa <- read.csv("Number of taxa per study split_taxa_coarse 11_2014.csv")

which(studies.taxa$number.taxa == 1)
length(which(studies.taxa$number.taxa == 1)) #25

more.than.one.taxa <- studies.taxa$SS[which(studies.taxa$number.taxa != 1)]

multiple.taxa.matched.landuse <- subset(matched.landuse, SS %in% more.than.one.taxa)
length(multiple.taxa.matched.landuse[,1]) #5311

multiple.taxa.matched.landuse <- droplevels(multiple.taxa.matched.landuse)












### Load dataset on matched landuse for amphibian, mammal and bird data only



matched.landuse_amp_mam_bir <- read.csv("matched.landuse_amphib_mammal_bird_11_2014.csv")
nrow(matched.landuse_amp_mam_bir)#2074, down to 1983 if 8 landuses



matched.landuse_amp_mam_bir$Predominant_habitat <- gsub("Primary forest", "Primary Vegetation", matched.landuse_amp_mam_bir$Predominant_habitat)
matched.landuse_amp_mam_bir$Predominant_habitat <- gsub("Primary non-forest", "Primary Vegetation", matched.landuse_amp_mam_bir$Predominant_habitat)

matched.landuse_amp_mam_bir <- subset(matched.landuse_amp_mam_bir, Predominant_habitat != "Secondary vegetation (indeterminate age)")

matched.landuse_amp_mam_bir$Predominant_habitat <- factor(matched.landuse_amp_mam_bir$Predominant_habitat)


m <- merge(matched.landuse_amp_mam_bir, access.1, "SSS")
colnames(m)[which(colnames(m) == "MEAN")] <- "access"
m <- merge(m, hpd.1, "SSS")
colnames(m)[which(colnames(m) == "MEAN")] <- "hpd"
m <- merge(m, elevation.1, "SSS")
colnames(m)[which(colnames(m) == "MEAN")] <- "elevation"
m <- merge(m, slope.1, "SSS")
colnames(m)[which(colnames(m) == "MEAN")] <- "slope"
m <- merge(m, ag.1, "SSS")

nrow(m)


matched.landuse_amp_mam_bir <- m

length(matched.landuse_amp_mam_bir[,1]) 
names(matched.landuse_amp_mam_bir)



#merge distance to boundary data onto matched landuse
nrow(PA_dists_data)
nrow(matched.landuse_amp_mam_bir)

m <- merge(matched.landuse_amp_mam_bir, PA_dists_data, by = "SSS")
nrow(m)
matched.landuse_amp_mam_bir <- m


matched.landuse_amp_mam_bir$bound_dist_km <- matched.landuse_amp_mam_bir$NEAR_DIST/1000
matched.landuse_amp_mam_bir$log_bound_dist_km <- log(matched.landuse_amp_mam_bir$bound_dist_km+1)

#make variables where the sign is negative inside PAs
matched.landuse_amp_mam_bir$bound_dist_km_PA_neg <- matched.landuse_amp_mam_bir$bound_dist_km
matched.landuse_amp_mam_bir$bound_dist_km_PA_neg[which(matched.landuse_amp_mam_bir$Within_PA == "yes")] <- -1*matched.landuse_amp_mam_bir$bound_dist_km_PA_neg[which(matched.landuse_amp_mam_bir$Within_PA == "yes")]
matched.landuse_amp_mam_bir$log_bound_dist_km_PA_neg <- matched.landuse_amp_mam_bir$log_bound_dist_km
matched.landuse_amp_mam_bir$log_bound_dist_km_PA_neg[which(matched.landuse_amp_mam_bir$Within_PA == "yes")] <- -1*matched.landuse_amp_mam_bir$log_bound_dist_km_PA_neg[which(matched.landuse_amp_mam_bir$Within_PA == "yes")]





# sort out explanatory variables

matched.landuse_amp_mam_bir$IUCN_CAT_number <- factor(matched.landuse_amp_mam_bir$IUCN_CAT_number) # they arent really in an order

matched.landuse_amp_mam_bir$log_slope <- log(matched.landuse_amp_mam_bir$slope +1)
matched.landuse_amp_mam_bir$log_elevation <- log(matched.landuse_amp_mam_bir$elevation +1)
matched.landuse_amp_mam_bir$log_hpd<- log(matched.landuse_amp_mam_bir$hpd +1)
matched.landuse_amp_mam_bir$log_access <- log(matched.landuse_amp_mam_bir$access +1)

matched.landuse_amp_mam_bir$log_GIS_AREA <- log(matched.landuse_amp_mam_bir$GIS_AREA+1)


#make response variables
matched.landuse_amp_mam_bir$range <- matched.landuse_amp_mam_bir$CWM_Geographic_range_log10_square_km
matched.landuse_amp_mam_bir$mass <- matched.landuse_amp_mam_bir$CWM_Mass_log10_g  ### UPDATE
matched.landuse_amp_mam_bir$veg <- matched.landuse_amp_mam_bir$CWM_Vegetative_height_log10_m
matched.landuse_amp_mam_bir$vol <- matched.landuse_amp_mam_bir$CWM_Length_derived_volume_3log10_mm

pos <- which(matched.landuse_amp_mam_bir$CWM_Adult_wet_mass_log10_g >0 & matched.landuse_amp_mam_bir$CWM_Maximum_wet_mass_log10_g >0) 
length(pos)

matched.landuse_amp_mam_bir$log_abundance <- log(matched.landuse_amp_mam_bir$Total_abundance +1)
matched.landuse_amp_mam_bir$log_sp_rich <- log(matched.landuse_amp_mam_bir$Species_richness+1)



matched.landuse_amp_mam_bir$DoP.PA <- matched.landuse_amp_mam_bir$DoP
matched.landuse_amp_mam_bir$DoP.PA[which(matched.landuse_amp_mam_bir$Within_PA == "no")] <- 0

matched.landuse_amp_mam_bir$AREA.PA <- matched.landuse_amp_mam_bir$GIS_AREA
matched.landuse_amp_mam_bir$AREA.PA[which(matched.landuse_amp_mam_bir$Within_PA == "no")] <- 0

matched.landuse_amp_mam_bir$log_AREA.PA <- matched.landuse_amp_mam_bir$log_GIS_AREA
matched.landuse_amp_mam_bir$log_AREA.PA[which(matched.landuse_amp_mam_bir$Within_PA == "no")] <- 0

matched.landuse_amp_mam_bir$IUCN.PA <- matched.landuse_amp_mam_bir$IUCN_CAT_number
levels(matched.landuse_amp_mam_bir$IUCN.PA) <- c( "1.5", "4.5", "7", "0")
matched.landuse_amp_mam_bir$IUCN.PA[which(matched.landuse_amp_mam_bir$Within_PA == "no")] <- "0"





### load dataset on endangered species only
# no longer using this - Jan 2015
#matched.landuse_VU_EN_CR <- read.csv("matched.landuse_VU_EN_CR_11_2014.csv")

#remove problem studies from checking abstracts

#nrow(matched.landuse_VU_EN_CR) #796


#m <- merge(matched.landuse_VU_EN_CR, access.1, "SSS")
#colnames(m)[which(colnames(m) == "MEAN")] <- "access"
#m <- merge(m, hpd.1, "SSS")
#colnames(m)[which(colnames(m) == "MEAN")] <- "hpd"
#m <- merge(m, elevation.1, "SSS")
#colnames(m)[which(colnames(m) == "MEAN")] <- "elevation"
#m <- merge(m, slope.1, "SSS")
#colnames(m)[which(colnames(m) == "MEAN")] <- "slope"
#m <- merge(m, ag.1, "SSS")
#nrow(m)

#matched.landuse_VU_EN_CR <- m

#merge distance to boundary data onto matched landuse
#nrow(PA_dists_data)
#nrow(matched.landuse_VU_EN_CR)

#m <- merge(matched.landuse_VU_EN_CR, PA_dists_data, by = "SSS")
#nrow(m)
#matched.landuse_VU_EN_CR <- m


#matched.landuse_VU_EN_CR$bound_dist_km <- matched.landuse_VU_EN_CR$NEAR_DIST/1000
#matched.landuse_VU_EN_CR$log_bound_dist_km <- log(matched.landuse_VU_EN_CR$bound_dist_km+1)

#make variables where the sign is negative inside PAs
#matched.landuse_VU_EN_CR$bound_dist_km_PA_neg <- matched.landuse_VU_EN_CR$bound_dist_km
#matched.landuse_VU_EN_CR$bound_dist_km_PA_neg[which(matched.landuse_VU_EN_CR$Within_PA == "yes")] <- -1*matched.landuse_VU_EN_CR$bound_dist_km_PA_neg[which(matched.landuse_VU_EN_CR$Within_PA == "yes")]
#matched.landuse_VU_EN_CR$log_bound_dist_km_PA_neg <- matched.landuse_VU_EN_CR$log_bound_dist_km
#matched.landuse_VU_EN_CR$log_bound_dist_km_PA_neg[which(matched.landuse_VU_EN_CR$Within_PA == "yes")] <- -1*matched.landuse_VU_EN_CR$log_bound_dist_km_PA_neg[which(matched.landuse_VU_EN_CR$Within_PA == "yes")]

# sort out explanatory variables

#matched.landuse_VU_EN_CR$IUCN_CAT_number <- factor(matched.landuse_VU_EN_CR$IUCN_CAT_number) # they arent really in an order

#matched.landuse_VU_EN_CR$log_slope <- log(matched.landuse_VU_EN_CR$slope +1)
#matched.landuse_VU_EN_CR$log_elevation <- log(matched.landuse_VU_EN_CR$elevation +1)
#matched.landuse_VU_EN_CR$log_hpd<- log(matched.landuse_VU_EN_CR$hpd +1)
#matched.landuse_VU_EN_CR$log_access <- log(matched.landuse_VU_EN_CR$access +1)
#matched.landuse_VU_EN_CR$log_GIS_AREA <- log(matched.landuse_VU_EN_CR$GIS_AREA+1)

#matched.landuse_VU_EN_CR$DoP.PA <- matched.landuse_VU_EN_CR$DoP
#matched.landuse_VU_EN_CR$DoP.PA[which(matched.landuse_VU_EN_CR$Within_PA == "no")] <- 0

#matched.landuse_VU_EN_CR$AREA.PA <- matched.landuse_VU_EN_CR$GIS_AREA
#matched.landuse_VU_EN_CR$AREA.PA[which(matched.landuse_VU_EN_CR$Within_PA == "no")] <- 0

#matched.landuse_VU_EN_CR$log_AREA.PA <- matched.landuse_VU_EN_CR$log_GIS_AREA
#matched.landuse_VU_EN_CR$log_AREA.PA[which(matched.landuse_VU_EN_CR$Within_PA == "no")] <- 0

#matched.landuse_VU_EN_CR$IUCN.PA <- matched.landuse_VU_EN_CR$IUCN_CAT_number
#levels(matched.landuse_VU_EN_CR$IUCN.PA) <- c( "1.5", "4.5", "7", "0")
#matched.landuse_VU_EN_CR$IUCN.PA[which(matched.landuse_VU_EN_CR$Within_PA == "no")] <- "0"

#make response variables
#matched.landuse_VU_EN_CR$log_abundance <- log(matched.landuse_VU_EN_CR$Total_abundance +1)

### make version with only species that studies more than one species
#multiple.taxa.matched.landuse_VU_EN_CR <- subset(matched.landuse_VU_EN_CR, SS %in% more.than.one.taxa)
#multiple.taxa.matched.landuse_VU_EN_CR <- droplevels(multiple.taxa.matched.landuse_VU_EN_CR)





