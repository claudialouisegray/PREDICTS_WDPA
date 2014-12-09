

setwd("N:/Documents/PREDICTS/WDPA analysis")

rm(list=ls()) 

library(yarg)
library(roquefort)
library(rgdal)
library(sp)
library(gplots)
library(ggplot2)
library(scales)
library(gridExtra)
library(optimx)
library(SDMTools)
library(data.table)



#set up useful functions

validate <- function(x) {
	par(mfrow = c(1,2))
	plot(resid(x)~ fitted(x))
	hist(resid(x))
	par(mfrow = c(1,1))
	}




source("N:/Documents/PREDICTS/WDPA analysis/site_metrics_per_taxon.R")





### load diversity file with PREDICTS data


setwd("N:/Documents/PREDICTS/data from Lawrence")
 
diversity <- readRDS('diversity-2014-11-13-03-40-23.rds')


# Drop studies that we do not have permission to use the data from
diversity<-diversity[((diversity$Insightly_category=="Data from paper") | (
  diversity$Insightly_category=="Data received/permission obtained") | (
    diversity$Insightly_category=="Email address invalid") | (
      diversity$Insightly_category=="Email sent - no response received")),]

# Drop a study with uncertain permission
diversity<-diversity[(diversity$Source_ID!="VK1_2007__StLaurent"),]


# Drop some studies where the land-use classification is problematic
problem.studies<-c("AD1_2011__Hanley","AD1_2011__Yoon","VK1_2007__StLaurent",
                   "VK1_2013__ABMIboreal","SE1_2012__Rosselli","HW1_2008__Pillsbury",
                   "MG1_2008__Boutin","MG1_2011__Schon","MG1_2012__Schmidt",
                   "DI1_2013__Munyuli")
diversity<-diversity[!(diversity$Source_ID %in% problem.studies),]


#Drop studies where all sites should be either inside or outside a PA, or which are unlikely to be entered correctly from 2014 students (GP)
to.remove <- c("AD1_2007__Winfree", "AD1_2008__Kohler", 
			"DL1_2007__Verdu", "DL1_2010__Carvalho", "DL1_2011__Mallari", "DL1_2013__Bartolommei", 
			"GP1_2009__Devineau", "GP1_2011__Schumann",
			"HB1_2008__Sedlock", "HP1_2005__Kessler",
			"HP1_2005__Kumar", "HP1_2006__Lachat",  "HP1_2007__Shahabuddin", "HP1_2009__Kessler",
			"HW1_2008__Cagle", "LK1_2012__Wiafe", "MJ1_2013__Ndanganga", "SC1_2010__Marsh")

diversity <- subset(diversity, Source_ID %in% to.remove == F)



# correct for sampling effort
# function divides abundances by the sampling effort so that all "Measurement" values are standardised to one sampling effort unit
cat('Correcting for sampling effort\n')
diversity <- CorrectSamplingEffort(diversity)


# merge sites that are very close
# i.e. merge that have exactly the same values of
#  match.cols <- c('Longitude','Latitude','Predominant_habitat', 
#    'Use_intensity', 'Sample_start_earliest', 'Sample_end_latest',
#    'Sample_date_resolution', 'Block')
# NB note in this code to improve by using coordinate precision
# this must be done before any extra columns are added

cat('Merging sites\n')

#remove column that is not included in MergeSites column
diversity <- diversity[, - which(names(diversity) == "Insightly_restrictions")]
diversity <- MergeSites(diversity)







# Load traits
#traits.dir <- 'N:/Documents/PREDICTS/data from Lawrence'
#pantheria <- PantheriaMammalianBodyMass(file.path(traits.dir, 'PANTHERIAMammalianBodyMass/RawData/PanTHERIA_1-0_WR05_Aug2008.txt'))
#gbif <- GBIFGeographicRange(file.path(traits.dir, 'GBIF/1_query_gbif'))



# Load traits extract
load('traits.RData')

ls()
# c.values - side project of Rebeccas on genome size, not to use here                   
# pantheria - mammals, gives "Adult_wet_mass_log10_g"
# su - ants, gives "Length_derived_volume_3log10_mm", "Head_width_log10_mm", "Head_length_log10_mm"
# edgar - beetles, gives "Length_derived_volume_3log10_mm"
# horton - spiders, "Diet_type_Category", "Hunting_method_Category", "Length_derived_volume_3log10_mm", "Ecological_specialism_Category"   
# middleton welling - butterflies, gives "Flight_duration_months", "Wingspan_log10_mm", "Forewing_length_log10_mm", 
	#"Host_diversity_log10_N_hosts", "Host_specificity_Category", "Univoltine_Logical","Found_in_artificial_habitat_Logical"
# butchart - birds, gives "Adult_wet_mass_log10_g"
# gilbert - hoverflies, gives "Thorax_volume_log10_mm.3"
# cooper - amphibs  gives  "Length_derived_volume_3log10_mm"
# meiri.and.feldman - reptiles .... , gives "Maximum_wet_mass_log10_g"
# senior - amphibs.... gives "Length_derived_volume_3log10_mm"
# try - plants, gives vegetative height


cat('Adding traits to diversity data\n')

#using this to check what each dataset adds
#diversity1 <- AddTraits(diversity, iucn) 
#names(diversity1) 


# in IUCN, some version 2.3 categories remain see http://www.iucnredlist.org/static/categories_criteria_2_3#categories
# LR/cd Lower Risk/conservation dependent
# LR/nt = Lower Risk/near threatened
# Lower Risk/least concern

diversity <- AddTraits(diversity, butchart, cooper,  gilbert, horton, edgar, su,
	meiri.and.feldman, middleton, pantheria, senior,  try, try.inferred, gbif, iucn)


cat('Adding genus-averaged traits to diversity data\n')
diversity <- AddGenusAveragedTraits(diversity, pantheria)


names(diversity)





# do any species have both adult wet mass and maximum wet mass
which(diversity$Adult_wet_mass_log10_g>0 & diversity$Maximum_wet_mass_log10_g >0)


# make new mass variable that combines both
diversity$Mass_log10_g <- NA
diversity$Mass_log10_g[which(diversity$Adult_wet_mass_log10_g>0)] <- diversity$Adult_wet_mass_log10_g[which(diversity$Adult_wet_mass_log10_g>0)]
diversity$Mass_log10_g[which(diversity$Maximum_wet_mass_log10_g>0)] <- diversity$Maximum_wet_mass_log10_g[which(diversity$Maximum_wet_mass_log10_g>0)]


traits <- c('Adult_wet_mass_log10_g',
            'Geographic_range_log10_square_km',
		'Maximum_wet_mass_log10_g',
		'Vegetative_height_log10_m',
		'Length_derived_volume_3log10_mm',
		'Mass_log10_g')


# check that there is now more info than for just mammals

#with.mass <- subset(diversity, Adult_wet_mass_log10_g > 0)
#unique(with.mass$Class) # just birds and mammmals

#with.max.mass <- subset(diversity, Maximum_wet_mass_log10_g > 0)
#unique(with.max.mass$Class) # reptiles 

#with.thorax <- subset(diversity, Thorax_volume_log10_mm.3 > 0)
#unique(with.thorax$Order) # flies 

#with.lengthvol <- subset(diversity, Length_derived_volume_3log10_mm > 0)
#unique(with.lengthvol$Order) # Coleoptera Araneae Sarcoptiformes (mites) Hymenoptera        amphibs: Anura, Caudata, Gymnophiona(caecilians)
#unique(with.lengthvol$Class)

#with.height <- subset(diversity, Vegetative_height_log10_m > 0)
#unique(with.height$Class) # try.inferred




# get taxon info

unique(diversity$Kingdom)
which(diversity$Kingdom == "Protozoa") 

# drop protozoa, as not enough to work with
diversity <- subset(diversity, Kingdom != "Protozoa")
diversity$taxon_of_interest <- factor(diversity$Kingdom, levels = c("Animalia", "Fungi", "Plantae", "Plants and Fungi", "Invertebrates", "Vertebrates"))

# drop study where only kingdom Animalia is known
index <- rep(1, length(diversity[,1]))
index[which(diversity$Phylum == "" & diversity$Kingdom == "Animalia")]  <- 0
diversity <- diversity[which(index == 1),]


# drop Fungi, as not enough to work with after matching has been done
diversity <- subset(diversity, Kingdom != "Fungi")
diversity$taxon_of_interest <- factor(diversity$Kingdom, levels = c("Animalia", "Plantae", "Plants", "Invertebrates", "Vertebrates"))




# split verts and inverts

diversity$taxon_of_interest[which(diversity$Phylum == "Arthropoda")] <- "Invertebrates"
diversity$taxon_of_interest[which(diversity$Phylum == "Mollusca")] <- "Invertebrates"
diversity$taxon_of_interest[which(diversity$Phylum == "Platyhelminthes")] <- "Invertebrates"
diversity$taxon_of_interest[which(diversity$Phylum == "Nematoda")] <- "Invertebrates"
diversity$taxon_of_interest[which(diversity$Phylum == "Annelida")] <- "Invertebrates"
diversity$taxon_of_interest[which(diversity$Phylum == "Onychophora")] <- "Invertebrates"
diversity$taxon_of_interest[which(diversity$Phylum == "Chordata")] <- "Vertebrates"
diversity$taxon_of_interest[which(diversity$Kingdom == "Plantae")] <- "Plants"



# check that none have a blank kingdom or unassigned kingdom
which(diversity$Kingdom == "")
which(diversity$Kingdom == "Not assigned")

# check that all the cases where phylum is not known can be allocated based on kingdom
unique(diversity$Kingdom[which(diversity$Phylum == "")]) # yes all are plants





# check how many studies have more than one taxon
# need each study to be on a different taxon

study.data <- split(diversity, diversity$SS)
numbers <- sapply(study.data, function(rows) length(unique(rows$taxon_of_interest)))
more.than.one.taxon <- names(which(numbers >1))
more.than.one.taxon
#"AD1_2002__Vazquez 1"    "SC1_2005__Richardson 1" - only two studies


# look at these - sanity check
#View(subset(diversity, SS %in% more.than.one.taxon)) 
#yes, reasonable. e.g. pollinators inc hummingbirds


# split studies into sub_studies based on these groups

# now create new diversity file where study has taxon names
# ideally would have new study numbers for each level of SS_taxa, with corresponding site numbers
# but this will work 

diversity.original <- diversity
#diversity <- diversity.original

diversity$SS_original <- diversity$SS
diversity$SSS_original <- diversity$SSS
diversity$SS <- paste(diversity$SS_original, diversity$taxon_of_interest)
diversity$SSS <- paste(diversity$SSS_original, diversity$taxon_of_interest)





# Make Zone - tropical versus temperate

diversity$Zone <- rep("Temperate", length(diversity[,1]))
diversity$Zone[which(diversity$Latitude > -23.4378 & diversity$Latitude <  23.4378)] <- "Tropical"



### would any study now be both tropical and temperate?

study.data <- split(diversity, diversity$SS)
climate.zone <- sapply(study.data, function(rows) length(unique(rows$Zone)))
which(climate.zone == 2)
# great, each study only has one level of zone








#####  get SiteMetrics


diversity <- droplevels(diversity)

# make SS and SSS factors so that site metrics code works
diversity$SS <- factor(diversity$SS)
diversity$SSS <- factor(diversity$SSS)



### now get SiteMetrics!

d.metrics <- SiteMetrics_taxa_split(diversity, extra.cols = c("Realm", "Zone", "taxon_of_interest", "Country"),  traits = traits)

length(d.metrics[,1]) #14792

# now that there are more sites, send to ARC to get WDPA data on each site. 

long.lat <- d.metrics[,c("Longitude", "Latitude")]
spatial <- SpatialPoints(long.lat)
plot(spatial)

spdf <- SpatialPointsDataFrame(spatial, d.metrics)

writeOGR(obj = spdf, dsn = "C:/GIS/PA_predicts_mapping", "all_sites_taxa_split_coarse_11_2014_new", driver = "ESRI Shapefile")





# write d.metrics file to add back on to site level data after ArcMap processing
setwd("N:/Documents/PREDICTS/WDPA analysis")

write.csv(d.metrics, "d.metrics_11_2014_all_trait_data.csv")
	

#d.metrics <- read.csv("d.metrics_all_mass_data.csv")






### CREATE MULTIPLE TAXA PER STUDY DATASET 
# create dataset for species richness analysis that without all studies that are only on one taxon



#take out the studies that only compare one species
# calculate how many species per study 
#can also do this as highest taxa level = species




study.data <- split(diversity, diversity$SS)
numbers <- sapply(study.data, function(rows) length(unique(rows$Taxon_number)))
names(study.data)



studies.taxa <- data.frame(SS = names(study.data), number.taxa = numbers)
nrow(studies.taxa)

setwd("N:/Documents/PREDICTS/WDPA analysis")
write.csv(studies.taxa, "Number of taxa per study split_taxa_coarse 11_2014.csv")



