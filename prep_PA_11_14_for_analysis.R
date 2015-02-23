


### LOAD DATA on confounding variables


setwd("N:/Documents/PREDICTS/WDPA analysis/matching data")


access <- read.table("bfer_1km_acc50k_moll_11_14.txt", header = T, sep = ",")
hpd <- read.table("bfer_1km_HPD_moll_11_14.txt", header = T, sep = ",")
elevation <- read.table("bfer_1km_mn30_elevation_moll_11_14.txt", header = T, sep = ",")
slope <- read.table("bfer_1km_mn30_slope_moll_11_14.txt", header = T, sep = ",")


ag_suit <- read.csv("ag_suitability_11_2014_moll.csv") 

# this is taken from extract values to points function, saved as text file, then opened in excel
# FID column deleted and -9999 values in RASTERVALU switched to NA

ag.1 <- ag_suit[,c("SSS", "RASTERVALU")] 
colnames(ag.1) <- c("SSS", "ag_suit")


#values of 9 are water, sites have fallen on the edge of lakes/rivers
ag.1$ag_suit[which(ag.1$ag_suit==9)] <- NA

# switch scale to be more intuitive - higher values more agriculturally suitable
ag.1$ag_suit <- 9 - ag.1$ag_suit

access.1 <- access[,c("SSS", "MEAN")]
hpd.1 <- hpd[,c("SSS", "MEAN")]
elevation.1 <- elevation[,c("SSS", "MEAN")]
slope.1 <- slope[,c("SSS", "MEAN")]




### Load PA_11_14 dataset 

setwd("N:/Documents/PREDICTS/WDPA analysis")

PA_11_14<- read.csv("PA_11_2014.csv")
nrow(PA_11_14) #7077

# group primary landuses and remove secondary indeterminate
PA_11_14$Predominant_habitat <- gsub("Primary forest", "Primary Vegetation", PA_11_14$Predominant_habitat)
PA_11_14$Predominant_habitat <- gsub("Primary non-forest", "Primary Vegetation", PA_11_14$Predominant_habitat)
PA_11_14 <- subset(PA_11_14, Predominant_habitat != "Secondary vegetation (indeterminate age)")
PA_11_14$Predominant_habitat <- factor(PA_11_14$Predominant_habitat)
nrow(PA_11_14)

#PA_11_14$Predominant_habitat <- relevel(PA_11_14$Predominant_habitat, "Primary Vegetation")

# merge on the confounding variables
m <- merge(PA_11_14, access.1, "SSS", all.x = T)
colnames(m)[which(colnames(m) == "MEAN")] <- "access"
m <- merge(m, hpd.1, "SSS", all.x = T)
colnames(m)[which(colnames(m) == "MEAN")] <- "hpd"
m <- merge(m, elevation.1, "SSS", all.x = T)
colnames(m)[which(colnames(m) == "MEAN")] <- "elevation"
m <- merge(m, slope.1, "SSS", all.x = T)
colnames(m)[which(colnames(m) == "MEAN")] <- "slope"
m <- merge(m, ag.1, "SSS", all.x = T)

nrow(m) #6631
PA_11_14 <- m

# sort out explanatory variables
#PA_11_14$IUCN_CAT_number <- factor(PA_11_14$IUCN_CAT_number) # they arent really in an order
# make IUCN cat variable of I or II vs III to VI vs unknown vs unprotected
PA_11_14$IUCN_CAT <- PA_11_14$IUCN_CAT_number 
levels(PA_11_14$IUCN_CAT) <- c(levels(PA_11_14$IUCN_CAT), "0")
PA_11_14$IUCN_CAT[which(PA_11_14$Within_PA == "no")] <- 0
PA_11_14$IUCN_CAT  <- factor(PA_11_14$IUCN_CAT)

PA_11_14$log_slope <- log(PA_11_14$slope +1)
PA_11_14$log_elevation <- log(PA_11_14$elevation +1)
PA_11_14$log_hpd<- log(PA_11_14$hpd +1)
PA_11_14$log_access <- log(PA_11_14$access +1)
PA_11_14$log_AREA.PA <- log(PA_11_14$GIS_AREA+1)

PA_11_14$LU_3 <- PA_11_14$Predominant_habitat
PA_11_14$LU_3 <- gsub("Cropland", "Human_dominated", PA_11_14$LU_3)
PA_11_14$LU_3 <- gsub("Plantation forest", "Human_dominated", PA_11_14$LU_3)
PA_11_14$LU_3 <- gsub("Pasture", "Human_dominated", PA_11_14$LU_3)
PA_11_14$LU_3 <- gsub("Urban", "Human_dominated", PA_11_14$LU_3)
PA_11_14$LU_3 <- gsub("Mature secondary vegetation", "Secondary", PA_11_14$LU_3)
PA_11_14$LU_3 <- gsub("Intermediate secondary vegetation", "Secondary", PA_11_14$LU_3)
PA_11_14$LU_3 <- gsub("Young secondary vegetation", "Secondary", PA_11_14$LU_3)
unique(PA_11_14$LU_3)

PA_11_14$PA <- PA_11_14$Within_PA
levels(PA_11_14$PA) <- c(levels(PA_11_14$PA), "IN", "OUT")
PA_11_14$PA[which(PA_11_14$Within_PA == "yes")] <- "IN"
PA_11_14$PA[which(PA_11_14$Within_PA == "no")] <- "OUT"
PA_11_14$LUPA <- factor(paste(PA_11_14$PA, PA_11_14$Predominant_habitat))

#make response variables
PA_11_14$log_abundance <- log(PA_11_14$Total_abundance +1)
PA_11_14$range <- PA_11_14$CWM_Geographic_range_log10_square_km

#set reference levels
PA_11_14 <- droplevels(PA_11_14)
PA_11_14$Predominant_habitat <- relevel(PA_11_14$Predominant_habitat, "Primary Vegetation")
PA_11_14$Within_PA <- relevel(PA_11_14$Within_PA, "no")
PA_11_14$LUPA <- relevel(PA_11_14$LUPA, "OUT Primary Vegetation")



### CREATE MULTIPLE TAXA PER STUDY DATASET 
# create dataset for species richness analysis that without all studies that are only on one taxon

studies.taxa <- read.csv("Number of taxa per study split_taxa_coarse 11_2014.csv")

which(studies.taxa$number.taxa == 1)
length(which(studies.taxa$number.taxa == 1)) #25

more.than.one.taxa <- studies.taxa$SS[which(studies.taxa$number.taxa != 1)]

multiple.taxa.PA_11_14 <- subset(PA_11_14, SS %in% more.than.one.taxa)
length(multiple.taxa.PA_11_14[,1]) #6461

# this is the dataset to use for species richness and rarefied richness analysis
multiple.taxa.PA_11_14 <- droplevels(multiple.taxa.PA_11_14)







### prep the data for the proportion threatened analysis
### Load dataset on matched landuse for amphibian, mammal and bird data only

PA_11_14_amp_mam_bir <- read.csv("PA_11_2014_amph_mamm_bird.csv")

# group primary landuses and remove secondary indeterminate
PA_11_14_amp_mam_bir$Predominant_habitat <- gsub("Primary forest", "Primary Vegetation", PA_11_14_amp_mam_bir$Predominant_habitat)
PA_11_14_amp_mam_bir$Predominant_habitat <- gsub("Primary non-forest", "Primary Vegetation", PA_11_14_amp_mam_bir$Predominant_habitat)
PA_11_14_amp_mam_bir <- subset(PA_11_14_amp_mam_bir, Predominant_habitat != "Secondary vegetation (indeterminate age)")
PA_11_14_amp_mam_bir$Predominant_habitat <- factor(PA_11_14_amp_mam_bir$Predominant_habitat)
nrow(PA_11_14_amp_mam_bir)#2617

m <- merge(PA_11_14_amp_mam_bir, access.1, "SSS")
colnames(m)[which(colnames(m) == "MEAN")] <- "access"
m <- merge(m, hpd.1, "SSS")
colnames(m)[which(colnames(m) == "MEAN")] <- "hpd"
m <- merge(m, elevation.1, "SSS")
colnames(m)[which(colnames(m) == "MEAN")] <- "elevation"
m <- merge(m, slope.1, "SSS")
colnames(m)[which(colnames(m) == "MEAN")] <- "slope"
m <- merge(m, ag.1, "SSS")

nrow(m) # check still 2617
PA_11_14_amp_mam_bir <- m

# sort out explanatory variables

PA_11_14_amp_mam_bir$IUCN_CAT_number <- factor(PA_11_14_amp_mam_bir$IUCN_CAT_number) # they arent really in an order
PA_11_14_amp_mam_bir$log_slope <- log(PA_11_14_amp_mam_bir$slope +1)
PA_11_14_amp_mam_bir$log_elevation <- log(PA_11_14_amp_mam_bir$elevation +1)
PA_11_14_amp_mam_bir$log_hpd<- log(PA_11_14_amp_mam_bir$hpd +1)
PA_11_14_amp_mam_bir$log_access <- log(PA_11_14_amp_mam_bir$access +1)
PA_11_14_amp_mam_bir$log_GIS_AREA <- log(PA_11_14_amp_mam_bir$GIS_AREA+1)

#make response variables
PA_11_14_amp_mam_bir$log_abundance <- log(PA_11_14_amp_mam_bir$Total_abundance +1)





