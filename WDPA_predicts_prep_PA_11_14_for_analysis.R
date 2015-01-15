


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

m <- merge(PA_11_14, access.1, "SSS", all.x = T)
colnames(m)[which(colnames(m) == "MEAN")] <- "access"
m <- merge(m, hpd.1, "SSS", all.x = T)
colnames(m)[which(colnames(m) == "MEAN")] <- "hpd"
m <- merge(m, elevation.1, "SSS", all.x = T)
colnames(m)[which(colnames(m) == "MEAN")] <- "elevation"
m <- merge(m, slope.1, "SSS", all.x = T)
colnames(m)[which(colnames(m) == "MEAN")] <- "slope"
m <- merge(m, ag.1, "SSS", all.x = T)

nrow(m) #check still 7077

PA_11_14 <- m



# sort out explanatory variables
PA_11_14$IUCN_CAT_number <- factor(PA_11_14$IUCN_CAT_number) # they arent really in an order
PA_11_14$log_slope <- log(PA_11_14$slope +1)
PA_11_14$log_elevation <- log(PA_11_14$elevation +1)
PA_11_14$log_hpd<- log(PA_11_14$hpd +1)
PA_11_14$log_access <- log(PA_11_14$access +1)
PA_11_14$log_AREA.PA <- log(PA_11_14$GIS_AREA+1)


#make response variables
# SAM - are you also now doing this? (remember us talking about how you were dropping 0 values before - cant remember what we decided in the end?)
PA_11_14$log_abundance <- log(PA_11_14$Total_abundance +1)





### CREATE MULTIPLE TAXA PER STUDY DATASET 
# create dataset for species richness analysis that without all studies that are only on one taxon
# SAM - possibly another thing we need to discuss here. I've dropped all studies that are only on one taxon from 
# the species richness and the rarefied richness analyses, as they cant tell the difference in species richness


studies.taxa <- read.csv("Number of taxa per study split_taxa_coarse 11_2014.csv")

which(studies.taxa$number.taxa == 1)
length(which(studies.taxa$number.taxa == 1)) #25

more.than.one.taxa <- studies.taxa$SS[which(studies.taxa$number.taxa != 1)]

multiple.taxa.PA_11_14 <- subset(PA_11_14, SS %in% more.than.one.taxa)
length(multiple.taxa.PA_11_14[,1]) #6874

# this is the dataset to use for species richness and rarefied richness analysis
multiple.taxa.PA_11_14 <- droplevels(multiple.taxa.PA_11_14)











### prep the data for the proportion threatened analysis
### Load dataset on matched landuse for amphibian, mammal and bird data only



PA_11_14_amp_mam_bir <- read.csv("PA_11_2014_amph_mamm_bird.csv")

nrow(PA_11_14_amp_mam_bir)#2695


m <- merge(PA_11_14_amp_mam_bir, access.1, "SSS")
colnames(m)[which(colnames(m) == "MEAN")] <- "access"
m <- merge(m, hpd.1, "SSS")
colnames(m)[which(colnames(m) == "MEAN")] <- "hpd"
m <- merge(m, elevation.1, "SSS")
colnames(m)[which(colnames(m) == "MEAN")] <- "elevation"
m <- merge(m, slope.1, "SSS")
colnames(m)[which(colnames(m) == "MEAN")] <- "slope"
m <- merge(m, ag.1, "SSS")

nrow(m) # check still 2695
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





