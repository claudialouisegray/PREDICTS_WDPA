


# load file from spatial join
# give PA attributes
# split to matched all and matched landuse datasets


#### read in ESRI file with WDPA data on each site, with site metrics already calculated

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


setwd("R:/ecocon_d/clg32/PREDICTS/WDPA analysis")

# in Arc, I ran a spatial join one to many so that can select PA characteristics below. 
# All PAs included
# so now I have a dataset where the WDPA cat is on each site

d.metrics.PA <- readOGR(dsn = "C:/GIS/PA_predicts_mapping", "predicts_wdpa_join_taxa_split_coarse_WDPA_nobiosphere_inscribed_adopted_11_2014")


# get just the data not the spatial info
all_PA_info <- d.metrics.PA@data

names(all_PA_info)
head(all_PA_info)
length(all_PA_info[,1]) #16533


# make binary in out variable
all_PA_info$Within_PA[all_PA_info$WDPAID == 0] <- "no" 
all_PA_info$Within_PA[all_PA_info$WDPAID > 0]<-"yes"


# get number of PAs per site:
# take within_PA sites first
in_PA <- subset(all_PA_info, Within_PA == "yes")
out_PA <- subset(all_PA_info, Within_PA == "no")


#get number of records per site (this is equal to the number of PAs it intercepts)
PAs_per_site <- table(droplevels(in_PA$SSS))
t = tabulate(PAs_per_site)
hist(PAs_per_site)



# rename fields concatenated by ARC

new.names <- gsub("Sorc_ID", "Source_ID", names(all_PA_info))
new.names <- gsub("Stdy_nmb", "Study_number", new.names)
new.names <- gsub("Study_nm", "Study_name", new.names)
new.names <- gsub("St_nmbr", "Site_number", new.names)
new.names <- gsub("Site_nm", "Site_name", new.names)
new.names <- gsub("Prdmnn_", "Predominant_habitat", new.names)
new.names <- gsub("Us_ntns", "UseIntensity", new.names)
new.names <- gsub("Smpl_s_", "Sample_start_earliest", new.names)
new.names <- gsub("Smpl_n_", "Sample_end_latest", new.names)
new.names <- gsub("Dvrst__", "Diversity_metric_type", new.names)


#new.names <- gsub("Ttl_bnd", "Total_abundance", new.names)
#new.names <- gsub("Spcs_rc", "Species_richness", new.names)
#new.names <- gsub("Smpsn_d", "Simpson_diversity", new.names)  
#new.names <- gsub("Rchnss_", "Richness_rarefied", new.names)
#new.names <- gsub("CWM_A__",  "CWM_Adult_wet_mass_log10_g", new.names)
#new.names <- gsub("PWD_A__", "PropWithData_Adult_wet_mass_log10_g" , new.names)
#new.names <- gsub("CWM_G__",  "CWM_Geographic_range_log10_square_km" , new.names)
#new.names <- gsub("PWD_G__",  "PropWithData_Geographic_range_log10_square_km", new.names)
#new.names <- gsub("IPR_A__",   "IPR_Adult_wet_mass_log10_g", new.names)
#new.names <- gsub("IPR_G__", "IPR_Geographic_range_log10_square_km" , new.names)


colnames(all_PA_info) <- new.names



# for each site, give the PA characteristics an order corresponding to the selection we want to make 
# Max area of all PAs, Lowest IUCN number (highest protection)

#set orders of category so that min and max can run on these variables
levels(all_PA_info$IUCN_CAT)
all_PA_info$IUCN_CAT_ordered <- factor(all_PA_info$IUCN_CAT, ordered = T, 
	levels = rev(c("Ia","Ib","II","III","IV","V","VI","Not Applicable","Not Reported")))

levels(all_PA_info$STATUS)
all_PA_info$STATUS_ordered <- factor(all_PA_info$STATUS, ordered = T, 
	levels = rev(c("Designated","Inscribed","Proposed","Not Reported")))

levels(all_PA_info$DESIG_TYPE)
all_PA_info$DESIG_TYPE_ordered <- factor(all_PA_info$DESIG_TYPE, ordered = T, 
	levels = (c("International","Regional","National")))









### JOIN ON DATA FOR SITES REQUIRED

# join on site data to restore it
# this is needed as site_metrics code has been updated to change simp. diversity calculation
# and because mass data is updated with all data (may need to be changed if permissons not obtained)
# also as ARC replaces NA with 0 and NaN with -1

d.metrics <- read.csv("d.metrics_11_2014_all_trait_data.csv")

names(d.metrics)

traits <- c("CWM_Adult_wet_mass_log10_g", "PropWithData_Adult_wet_mass_log10_g",
	"CWM_Geographic_range_log10_square_km", "PropWithData_Geographic_range_log10_square_km",
	"CWM_Maximum_wet_mass_log10_g", "PropWithData_Maximum_wet_mass_log10_g",
	"CWM_Mass_log10_g", "PropWithData_Mass_log10_g",
	"CWM_Vegetative_height_log10_m", "PropWithData_Vegetative_height_log10_m",
	"CWM_Length_derived_volume_3log10_mm", "PropWithData_Length_derived_volume_3log10_mm",
	"IPR_Adult_wet_mass_log10_g", "IPR_Geographic_range_log10_square_km",
	"IPR_Maximum_wet_mass_log10_g", "IPR_Vegetative_height_log10_m", 
	"IPR_Length_derived_volume_3log10_mm", "IPR_Mass_log10_g")

site.data <- d.metrics[,c("SSS", "taxon_of_interest",
	"Total_abundance", "Species_richness", "Simpson_diversity",  "Richness_rarefied", traits)]
length(site.data[,1]) #14792

	
# merge onto data for each site
all_PA_info_merge <- merge(all_PA_info, site.data, "SSS", all.y = T)
length(all_PA_info_merge[,1]) #16533


# for buffered points, take rep_area and fill in GIS_AREA with this
# buffered points have field "radius" as this is what was used to calculate buffer
names(all_PA_info_merge)
pos <- which(all_PA_info_merge$Radius > 0)

all_PA_info_merge$GIS_AREA[pos] <- all_PA_info_merge$REP_AREA[pos]


# are there values for all GIS_AREA inside PAs? 
in_PA <- subset(all_PA_info_merge, Within_PA == "yes")
which(in_PA$GIS_AREA == 0) #yes


# select data per site

sites <- unique(all_PA_info_merge$SSS)
PA.predicts <- data.frame()

#s <- 1

all_PA_info_merge$SSB <- paste(all_PA_info_merge$SS, all_PA_info_merge$Block)
all_PA_info_merge$SSBS <- paste(all_PA_info_merge$SS, all_PA_info_merge$Block, all_PA_info_merge$Site_number)

### Allocate PA attributes to each site

for (s in 1:length(sites)){
	data.subset <- subset(all_PA_info_merge, SSS == sites[s])

	Source_ID <- unique(data.subset$Source_ID)
	SS <- unique(data.subset$SS)
	SSS <- unique(data.subset$SSS)
	SSB <- unique(data.subset$SSB)
	SSBS <- unique(data.subset$SSBS)
	Country <- unique(data.subset$Country)
	Study_number <- unique(data.subset$Study_number)
	Site_number <- unique(data.subset$Site_number)
	Long <- unique(data.subset$Longitd)
	Lat <- unique(data.subset$Latitud)
	Diversity_metric_type <- unique(data.subset$Diversity_metric_type)
	Predominant_habitat <- unique(data.subset$Predominant_habitat)
	UseIntensity <- unique(data.subset$UseIntensity)
	Sample_start_earliest <- unique(data.subset$Sample_start_earliest)
	Realm <- unique(data.subset$Realm)
	Zone <- unique(data.subset$Zone)
	taxon_of_interest <- unique(data.subset$taxon_of_interest)
	Total_abundance <- unique(data.subset$Total_abundance)
	Species_richness <- unique(data.subset$Species_richness)
	ChaoR <- unique(data.subset$ChaoR)
	Richness_rarefied <- unique(data.subset$Richness_rarefied)
	Simpson_diversity <- unique(data.subset$Simpson_diversity)

	CWM_Adult_wet_mass_log10_g <- unique(data.subset$CWM_Adult_wet_mass_log10_g)         
	PropWithData_Adult_wet_mass_log10_g <- unique(data.subset$PropWithData_Adult_wet_mass_log10_g)         
	CWM_Geographic_range_log10_square_km <- unique(data.subset$CWM_Geographic_range_log10_square_km)
	PropWithData_Geographic_range_log10_square_km <- unique(data.subset$PropWithData_Geographic_range_log10_square_km)
	CWM_Maximum_wet_mass_log10_g <- unique(data.subset$CWM_Maximum_wet_mass_log10_g)
	PropWithData_Maximum_wet_mass_log10_g <- unique(data.subset$PropWithData_Maximum_wet_mass_log10_g)
	CWM_Vegetative_height_log10_m <- unique(data.subset$CWM_Vegetative_height_log10_m)
	PropWithData_Vegetative_height_log10_m <- unique(data.subset$PropWithData_Vegetative_height_log10_m)
	CWM_Length_derived_volume_3log10_mm <- unique(data.subset$CWM_Length_derived_volume_3log10_mm)
	PropWithData_Length_derived_volume_3log10_mm <- unique(data.subset$PropWithData_Length_derived_volume_3log10_mm)
	CWM_Mass_log10_g <- unique(data.subset$CWM_Mass_log10_g)
	PropWithData_Mass_log10_g <- unique(data.subset$PropWithData_Mass_log10_g)

	IPR_Adult_wet_mass_log10_g <- unique(data.subset$IPR_Adult_wet_mass_log10_g)
	IPR_Geographic_range_log10_square_km <- unique(data.subset$IPR_Geographic_range_log10_square_km)
	IPR_Maximum_wet_mass_log10_g <- unique(data.subset$IPR_Maximum_wet_mass_log10_g)
	IPR_Vegetative_height_log10_m <- unique(data.subset$IPR_Vegetative_height_log10_m)
	IPR_Length_derived_volume_3log10_mm <- unique(data.subset$IPR_Length_derived_volume_3log10_mm)
	IPR_IUCN_Red_List_Score <- unique(data.subset$IPR_IUCN_Red_List_Score)

	Within_PA <- unique(data.subset$Within_PA)

	IUCN_CAT <- max(data.subset$IUCN_CAT_ordered)
	STATUS <- max(data.subset$STATUS_ordered)
	DESIG_TYPE <- max(data.subset$DESIG_TYPE_ordered)
	GIS_AREA <- max(data.subset$GIS_AREA)
	STATUS_YR <- min(data.subset$STATUS_YR)

	
	df <- data.frame(Source_ID = Source_ID,
		Study_number = Study_number, Site_number = Site_number,
		SS = SS, SSS = SSS, SSB = SSB, SSBS = SSBS,
		Country = Country,
		Longitd = Long,
		Latitud = Lat,
		Diversity_metric_type = Diversity_metric_type, 
		Predominant_habitat = Predominant_habitat, 
		UseIntensity  = UseIntensity, 
		Sample_start_earliest = Sample_start_earliest, 
		Realm = Realm, 
		Zone = Zone, 
		taxon_of_interest = taxon_of_interest, 
		Total_abundance = Total_abundance, 
		Species_richness = Species_richness, 
		ChaoR = ChaoR, 
		Richness_rarefied = Richness_rarefied, 
		Simpson_diversity = Simpson_diversity, 
		CWM_Adult_wet_mass_log10_g = CWM_Adult_wet_mass_log10_g,
		PropWithData_Adult_wet_mass_log10_g = PropWithData_Adult_wet_mass_log10_g,
		CWM_Geographic_range_log10_square_km = CWM_Geographic_range_log10_square_km,
		PropWithData_Geographic_range_log10_square_km = PropWithData_Geographic_range_log10_square_km,
		CWM_Maximum_wet_mass_log10_g = CWM_Maximum_wet_mass_log10_g,
		PropWithData_Maximum_wet_mass_log10_g = PropWithData_Maximum_wet_mass_log10_g,
		CWM_Vegetative_height_log10_m = CWM_Vegetative_height_log10_m,
		PropWithData_Vegetative_height_log10_m = PropWithData_Vegetative_height_log10_m,
		CWM_Length_derived_volume_3log10_mm = CWM_Length_derived_volume_3log10_mm,
		PropWithData_Length_derived_volume_3log10_mm = PropWithData_Length_derived_volume_3log10_mm,
		CWM_Mass_log10_g = CWM_Mass_log10_g,
		PropWithData_Mass_log10_g = PropWithData_Mass_log10_g,


		IPR_Adult_wet_mass_log10_g = IPR_Adult_wet_mass_log10_g, 
		IPR_Geographic_range_log10_square_km = IPR_Geographic_range_log10_square_km,
		IPR_Maximum_wet_mass_log10_g = IPR_Maximum_wet_mass_log10_g,
		IPR_Vegetative_height_log10_m = IPR_Vegetative_height_log10_m,
		IPR_Length_derived_volume_3log10_mm = IPR_Length_derived_volume_3log10_mm,
		IUCN_CAT = IUCN_CAT, 
		STATUS = STATUS, 
		DESIG_TYPE  = DESIG_TYPE, 
		GIS_AREA = GIS_AREA, 
		STATUS_YR = STATUS_YR, 
		Within_PA = Within_PA)

	PA.predicts<- rbind(PA.predicts, df)
	#print(s)
	}

length(unique(PA.predicts$SS)) # 468

# take out sites that are Predominant_habitat "Cannot decide"
unique(PA.predicts$Predominant_habitat)
PA.predicts.subset <- subset(PA.predicts, Predominant_habitat != "Cannot decide")
unique(PA.predicts.subset$Predominant_habitat)

# can leave in the ones with intensity "Cannot decide"
length(PA.predicts.subset[,1]) # 14331

# change the classifications to new groups based on previous summary of how many sites in each group
# exact numbers of sites in each class may differ from totals of the previous summary tables
# as reclassifying here may allow more studies to be included
# e.g. previously MSV inside and ISV outside would not have been included. Now it can be. 

#PA.predicts.subset$Predominant_habitat <- gsub("Primary forest", "Primary Vegetation", PA.predicts.subset$Predominant_habitat)
#PA.predicts.subset$Predominant_habitat <- gsub("Primary non-forest", "Primary Vegetation", PA.predicts.subset$Predominant_habitat)
#PA.predicts.subset$Predominant_habitat <- gsub("Young secondary vegetation", "Secondary Vegetation",PA.predicts.subset$Predominant_habitat)
#PA.predicts.subset$Predominant_habitat <- gsub("Mature secondary vegetation", "Secondary Vegetation",PA.predicts.subset$Predominant_habitat)
#PA.predicts.subset$Predominant_habitat <- gsub("Intermediate secondary vegetation","Secondary Vegetation", PA.predicts.subset$Predominant_habitat)
#PA.predicts.subset$Predominant_habitat[which(PA.predicts.subset$Predominant_habitat == "Secondary vegetation (indeterminate age)")] <- "Secondary Vegetation"

# MSV = 629, ISV = 901, YSC = 1133, IndYS = 1051  - i.e. a large chunk is indeterminate sec veg. 
unique(PA.predicts.subset$Predominant_habitat)

### calculate duration of protection (DoP)

earliest_sample_date <- as.Date(PA.predicts.subset$Sample_start_earliest)
sample_year <- format(earliest_sample_date, "%Y")
PA.predicts.subset$DoP <- as.numeric(sample_year) - PA.predicts.subset$STATUS_YR
PA.predicts.subset$DoP[which(PA.predicts.subset$STATUS_YR == 0)] <- NA
PA.predicts.subset$STATUS_YR[which(PA.predicts.subset$STATUS_YR == 0)] <- NA

# if any sites are DoP <0 , make unprotected
PA.predicts.subset$Within_PA[which(PA.predicts.subset$DoP < 0)] <- "no"



### get distance to other sites in each study

earth.dist <- function (long1, lat1, long2, lat2)
#uses haversine formula, gives distance in km
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


studies <- unique(PA.predicts.subset$SS)
s <- studies[2]
n <- 3
i <- sites[1]
j <- sites[1]
df <- data.frame()
results <- list()

for(s in studies){
	ds <- subset(PA.predicts.subset, SS == s)
	sites <- unique(ds$SSS)
	distances <- data.frame(matrix(nrow = length(sites), ncol = length(sites)))
	colnames(distances) <- sites ; rownames(distances) <- sites
		for(i in sites){
			for(j in sites){
		 		dist <- earth.dist(ds$Longitd[which(ds$SSS == i)],ds$Latitud[which(ds$SSS == i)], ds$Longitd[which(ds$SSS == j)],ds$Latitud[which(ds$SSS == j)])
				distances[i,j] <- dist
			}
		}
	for(n in 1:length(sites)){
		row.data <- distances[n, -n] # get the distances without the distance to self
		mean.dist.to.other.sites <- mean(as.numeric(row.data))
		results[[paste(sites[n])]] <- mean.dist.to.other.sites
	}
   print(s)
}
			
too.far <- names(which(results > 150))
hist(as.numeric(results))
nrow(PA.predicts.subset) # 14331
PA.predicts.subset<- subset(PA.predicts.subset, SSS %in% too.far == F)
nrow(PA.predicts.subset) #13908
# studies with only one site are those that were merged together in merge sites

# are there values for all GIS_AREA inside PAs? 
in_PA <- subset(PA.predicts.subset, Within_PA == "yes")
which(in_PA$GIS_AREA == 0) #yes

long.lat <- PA.predicts[,c("Longitd", "Latitud")]
spatial <- SpatialPoints(long.lat)
plot(spatial)
PA.predicts$IUCN_CAT <- factor(PA.predicts$IUCN_CAT, ordered = F)
PA.predicts$STATUS <- factor(PA.predicts$STATUS, ordered = F)
PA.predicts$DESIG_TYPE <- factor(PA.predicts$DESIG_TYPE, ordered = F)
spdf <- SpatialPointsDataFrame(spatial, PA.predicts)
writeOGR(obj = spdf, dsn = "C:/GIS/PA_predicts_mapping", "PA.predicts_11_2014", driver = "ESRI Shapefile")

setwd("R:/ecocon_d/clg32/PREDICTS/WDPA analysis")
######
SUBSET TO MATCHED PAIRS ALL
all sites within matched studies
###### 

# get subset of only sites within PAs
within_PA <- subset(PA.predicts.subset, Within_PA == "yes")

#calculate total and number inside PA for each study
sites_total <-  aggregate(Within_PA ~  SS, PA.predicts.subset,length)
#aggregate(PA.predicts.subset$Within_PA, by = list( SS = PA.predicts.subset$SS), length)

colnames(sites_total) <- c("SS", "sites_total")
length(sites_total[,1]) #459

sites_in <- aggregate(Within_PA ~  SS, within_PA,length)
colnames(sites_in) <- c( "SS","sites_in")

#combine into dataframe
sites_in_out <- merge(sites_total, sites_in, by = c("SS"), all.x = TRUE)
length(sites_in_out[,1]) #still 459

# replace NAs with 0 as there are no sites in PAs in these studies
sites_in_out$sites_in[is.na(sites_in_out$sites_in)] <- 0

# calculate number of sites outside of PAs
sites_in_out$sites_out <- sites_in_out$sites_total - sites_in_out$sites_in
#View(sites_in_out)

#make index of which studies meet inclusion criteria (i.e. > 1 site in and out)
include <- unique(sites_in_out$SS[which(sites_in_out$sites_in > 0 
		& sites_in_out$sites_out > 0 )])
length(include) # check = 187 studies that have sites in and out

#use this index to select site level data that matches inclusion criteria
PA.predicts.subset1 <- subset(PA.predicts.subset, SS %in% include)

unique(PA.predicts.subset1$Predominant_habitat) #check no "Cannot decide" or indeterminate
unique(PA.predicts.subset1$UseIntensity) #doesnt matter if cannot decide as not including as fixed factor

# change IUCN cats to numbers so can express mean of all categories succintly
# and combine so that 
# 1.5 = Ia, Ib and II
# 4.5 = III and IV, V, VI

PA.predicts.subset1$IUCN_CAT_number <- NA

# NB CHECK ALL uniquePA.predicts.subset1$IUCN_CAT_number are accounted for

PA.predicts.subset1$IUCN_CAT_number[which(PA.predicts.subset1$IUCN_CAT == "Ia")] <- 1.5
PA.predicts.subset1$IUCN_CAT_number[which(PA.predicts.subset1$IUCN_CAT == "Ib")] <- 1.5
PA.predicts.subset1$IUCN_CAT_number[which(PA.predicts.subset1$IUCN_CAT == "II")] <- 1.5
PA.predicts.subset1$IUCN_CAT_number[which(PA.predicts.subset1$IUCN_CAT == "III")] <- 4.5
PA.predicts.subset1$IUCN_CAT_number[which(PA.predicts.subset1$IUCN_CAT == "IV")] <- 4.5
PA.predicts.subset1$IUCN_CAT_number[which(PA.predicts.subset1$IUCN_CAT == "V")] <- 4.5
PA.predicts.subset1$IUCN_CAT_number[which(PA.predicts.subset1$IUCN_CAT == "VI")] <- 4.5
PA.predicts.subset1$IUCN_CAT_number[which(PA.predicts.subset1$IUCN_CAT == "Not Reported")] <- 7
PA.predicts.subset1$IUCN_CAT_number[which(PA.predicts.subset1$IUCN_CAT == "Not Applicable")] <- 7



PA.predicts.subset1$Within_PA <- factor(PA.predicts.subset1$Within_PA)
PA.predicts.subset1$Within_PA <- relevel(PA.predicts.subset1$Within_PA, "yes")

PA.predicts.subset1$Predominant_habitat <- factor(PA.predicts.subset1$Predominant_habitat)
PA.predicts.subset1$Predominant_habitat <- relevel(PA.predicts.subset1$Predominant_habitat, "Primary Vegetation")

matched.all <- PA.predicts.subset1
length(matched.all[,1]) #7077
matched.all$Richness_rarefied <- round(matched.all$Richness_rarefied)


#save - previously matched.all
write.csv(matched.all, "PA_11_2014.csv")

#when data updates have occurred, check differences
#PA_11_2014 <- read.csv("PA_11_2014.csv")
#sites.old <-  unique(PA_11_2014$SSS)
#sites.new <- unique(matched.all$SSS)
#sites.new[which(sites.new %in% sites.old == F)]

# save as shapefile
long.lat <- matched.all[,c("Longitd", "Latitud")]
spatial <- SpatialPoints(long.lat)
plot(spatial)
matched.all$IUCN_CAT <- factor(matched.all$IUCN_CAT, ordered = F)
matched.all$STATUS <- factor(matched.all$STATUS, ordered = F)
matched.all$DESIG_TYPE <- factor(matched.all$DESIG_TYPE, ordered = F)

spdf <- SpatialPointsDataFrame(spatial, matched.all)
writeOGR(obj = spdf, dsn = "C:/GIS/PA_predicts_mapping", "PA_11_2014", driver = "ESRI Shapefile")

######
SUBSET TO MATCHED LANDUSE PAIRS
######

# get subset of only sites within PAs
within_PA <- subset(PA.predicts.subset, Within_PA == "yes")

#calculate total and number inside PA for each landuse for each study
sites_PH_total <- aggregate(Within_PA ~ Predominant_habitat + SS, PA.predicts.subset,length)
colnames(sites_PH_total) <- c("Predominant_habitat","SS", "sites_total")
length(sites_PH_total[,1]) #836 diff levels of SS and PH
sites_PH_in <- aggregate(Within_PA ~  Predominant_habitat + SS, within_PA,length)
colnames(sites_PH_in) <- c( "Predominant_habitat", "SS","sites_in")


#combine into dataframe
sites_PH_in_out <- merge(sites_PH_total, sites_PH_in, by = c("SS", "Predominant_habitat"), all.x = TRUE)
length(sites_PH_in_out[,1]) 

# replace NAs with 0 as there are no sites in PAs in these studies
sites_PH_in_out$sites_in[is.na(sites_PH_in_out$sites_in)] <- 0

# calculate number of sites outside of PAs
sites_PH_in_out$sites_out <- sites_PH_in_out$sites_total - sites_PH_in_out$sites_in
#View(sites_PH_in_out)

#make index of which habitats in which studies meet inclusion criteria (i.e. > 1 site in and out)
sites_PH_in_out$SS_PH <- paste(sites_PH_in_out$SS, sites_PH_in_out$Predominant_habitat)
include <- unique(sites_PH_in_out$SS_PH[which(sites_PH_in_out$sites_in > 0 
		& sites_PH_in_out$sites_out > 0 )])
length(include) 

#use this index to select site level data that matches inclusion criteria
PA.predicts.subset$SS_PH <- paste(PA.predicts.subset$SS, PA.predicts.subset$Predominant_habitat)
PA.predicts.subset1 <- subset(PA.predicts.subset, SS_PH %in% include)

# change IUCN cats to numbers so can take the highest one per study
# and combine so that 
# 1.5 = Ia, Ib and II
# 4.5 = III,  IV, V and VI

PA.predicts.subset1$IUCN_CAT_number <- NA
PA.predicts.subset1$IUCN_CAT_number[which(PA.predicts.subset1$IUCN_CAT == "Ia")] <- 1.5
PA.predicts.subset1$IUCN_CAT_number[which(PA.predicts.subset1$IUCN_CAT == "Ib")] <- 1.5
PA.predicts.subset1$IUCN_CAT_number[which(PA.predicts.subset1$IUCN_CAT == "II")] <- 1.5
PA.predicts.subset1$IUCN_CAT_number[which(PA.predicts.subset1$IUCN_CAT == "III")] <- 4.5
PA.predicts.subset1$IUCN_CAT_number[which(PA.predicts.subset1$IUCN_CAT == "IV")] <- 4.5
PA.predicts.subset1$IUCN_CAT_number[which(PA.predicts.subset1$IUCN_CAT == "V")] <- 4.5
PA.predicts.subset1$IUCN_CAT_number[which(PA.predicts.subset1$IUCN_CAT == "VI")] <- 4.5
PA.predicts.subset1$IUCN_CAT_number[which(PA.predicts.subset1$IUCN_CAT == "Not Reported")] <- 7
PA.predicts.subset1$IUCN_CAT_number[which(PA.predicts.subset1$IUCN_CAT == "Not Applicable")] <- 7

PA.predicts.subset1$Within_PA <- factor(PA.predicts.subset1$Within_PA)
PA.predicts.subset1$Within_PA <- relevel(PA.predicts.subset1$Within_PA, "yes")

PA.predicts.subset1$Predominant_habitat <- factor(PA.predicts.subset1$Predominant_habitat)
PA.predicts.subset1$Predominant_habitat <- relevel(PA.predicts.subset1$Predominant_habitat, "Primary Vegetation")

unique(PA.predicts.subset1$Predominant_habitat) #check no "Cannot decide" or indeterminate
matched.landuse <- PA.predicts.subset1 #incase overwriting during coding
names(matched.landuse)
length(matched.landuse[,1]) # 5366

#save
matched.landuse$Richness_rarefied <- round(matched.landuse$Richness_rarefied)
write.csv(matched.landuse, "matched.landuse_11_2014.csv")

#when data updates have occurred, check differences
#old <- read.csv("matched.landuse_11_2014.csv")
#sites.old <- unique(old$SSS)
#sites.new <- unique(matched.landuse$SSS)
#sites.new[which(sites.new %in% sites.old == F)]


# save as shapefile
long.lat <- matched.landuse[,c("Longitd", "Latitud")]
spatial <- SpatialPoints(long.lat)
plot(spatial)
matched.landuse$IUCN_CAT <- factor(matched.landuse$IUCN_CAT, ordered = F)
matched.landuse$STATUS <- factor(matched.landuse$STATUS, ordered = F)
matched.landuse$DESIG_TYPE <- factor(matched.landuse$DESIG_TYPE, ordered = F)
spdf <- SpatialPointsDataFrame(spatial, matched.landuse)
writeOGR(obj = spdf, dsn = "C:/GIS/PA_predicts_mapping", "matched.landuse", driver = "ESRI Shapefile")



