

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




source("R:/ecocon_d/clg32/GitHub/PREDICTS_WDPA/site_metrics_per_taxon.R")





### load diversity file with PREDICTS data


setwd("N:/Documents/PREDICTS/data from Lawrence")

diversity <- readRDS('diversity-2014-11-13-03-40-23.rds')



# View(diversity)
# diversity is a dataset with each row = one species, therefore multiple rows per site




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


#Drop studies where all sites should be either inside or outside a PA, or which are unlikely to be entered correctly (GP)
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


diversity1 <- AddTraits(diversity, iucn) 
names(diversity1) #using this to check what each dataset adds

str(diversity1$IUCN_Red_List_Status) #ordered factor 


#in IUCN, some version 2.3 categories remain see http://www.iucnredlist.org/static/categories_criteria_2_3#categories
# LR/cd Lower Risk/conservation dependent
# LR/nt = Lower Risk/near threatened
# Lower Risk/least concern

diversity <- AddTraits(diversity, butchart, cooper,  gilbert, horton, edgar, su,
	meiri.and.feldman, middleton, pantheria, senior,  try, try.inferred, gbif, iucn)


cat('Adding genus-averaged traits to diversity data\n')
diversity <- AddGenusAveragedTraits(diversity, pantheria)


names(diversity)




# do any species have both adult wet mass and maxiumum wet mass
which(diversity$Adult_wet_mass_log10_g>0 & diversity$Maximum_wet_mass_log10_g >0)


# make new mass variable that combines both
diversity$Mass_log10_g <- NA
diversity$Mass_log10_g[which(diversity$Adult_wet_mass_log10_g>0)] <- diversity$Adult_wet_mass_log10_g[which(diversity$Adult_wet_mass_log10_g>0)]
diversity$Mass_log10_g[which(diversity$Maximum_wet_mass_log10_g>0)] <- diversity$Maximum_wet_mass_log10_g[which(diversity$Maximum_wet_mass_log10_g>0)]

length(which(diversity$Adult_wet_mass_log10_g>0))
length(which(diversity$Maximum_wet_mass_log10_g >0))
length(which(diversity$Mass_log10_g>0))


# make RedListScore

diversity$IUCN_Red_List_Score <- NA
which(diversity$IUCN_Red_List_Status == "EX") # none
which(diversity$IUCN_Red_List_Status == "EW") # none

#diversity$IUCN_Red_List_Score[which(diversity$IUCN_Red_List_Status == "CR")] <- 4
#diversity$IUCN_Red_List_Score[which(diversity$IUCN_Red_List_Status == "EN")] <- 3
#diversity$IUCN_Red_List_Score[which(diversity$IUCN_Red_List_Status == "VU")] <- 2
#diversity$IUCN_Red_List_Score[which(diversity$IUCN_Red_List_Status == "NT")] <- 1
#diversity$IUCN_Red_List_Score[which(diversity$IUCN_Red_List_Status == "LC")] <- 0
#diversity$IUCN_Red_List_Score[which(diversity$IUCN_Red_List_Status == "LR/cd")] <- 1
#diversity$IUCN_Red_List_Score[which(diversity$IUCN_Red_List_Status == "LR/nt")] <- 1
#diversity$IUCN_Red_List_Score[which(diversity$IUCN_Red_List_Status == "LR/lc")] <- 0
#diversity$IUCN_Red_List_Score[which(diversity$IUCN_Red_List_Status == "DD")] <- NA
#diversity$IUCN_Red_List_Score[which(diversity$IUCN_Red_List_Status == "NE")] <- NA


diversity$IUCN_Red_List_Score[which(diversity$IUCN_Red_List_Status == "CR")] <- 1
diversity$IUCN_Red_List_Score[which(diversity$IUCN_Red_List_Status == "EN")] <- 2
diversity$IUCN_Red_List_Score[which(diversity$IUCN_Red_List_Status == "VU")] <- 3
diversity$IUCN_Red_List_Score[which(diversity$IUCN_Red_List_Status == "NT")] <- 4
diversity$IUCN_Red_List_Score[which(diversity$IUCN_Red_List_Status == "LC")] <- 5
diversity$IUCN_Red_List_Score[which(diversity$IUCN_Red_List_Status == "LR/cd")] <- 4
diversity$IUCN_Red_List_Score[which(diversity$IUCN_Red_List_Status == "LR/nt")] <- 4
diversity$IUCN_Red_List_Score[which(diversity$IUCN_Red_List_Status == "LR/lc")] <- 5
diversity$IUCN_Red_List_Score[which(diversity$IUCN_Red_List_Status == "DD")] <- NA
diversity$IUCN_Red_List_Score[which(diversity$IUCN_Red_List_Status == "NE")] <- NA


# make binary endangered index
diversity$endangered <- NA

diversity$endangered[which(diversity$IUCN_Red_List_Status == "CR")] <- 1
diversity$endangered[which(diversity$IUCN_Red_List_Status == "EN")] <- 1
diversity$endangered[which(diversity$IUCN_Red_List_Status == "VU")] <- 1
diversity$endangered[which(diversity$IUCN_Red_List_Status == "NT")] <- 0
diversity$endangered[which(diversity$IUCN_Red_List_Status == "LC")] <- 0
diversity$endangered[which(diversity$IUCN_Red_List_Status == "LR/cd")] <- 0
diversity$endangered[which(diversity$IUCN_Red_List_Status == "LR/nt")] <- 0
diversity$endangered[which(diversity$IUCN_Red_List_Status == "LR/lc")] <- 0
diversity$endangered[which(diversity$IUCN_Red_List_Status == "DD")] <- NA
diversity$endangered[which(diversity$IUCN_Red_List_Status == "NE")] <- NA




traits <- c('Adult_wet_mass_log10_g',
            'Geographic_range_log10_square_km',
		'Maximum_wet_mass_log10_g',
		'Vegetative_height_log10_m',
		'Length_derived_volume_3log10_mm',
		'Mass_log10_g',
		'IUCN_Red_List_Score',
		'endangered')





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


# split studies into taxonomically unique studies based on these groups

# now create new diversity file where study has taxon names
# ideally would have new study numbers for each level of SS_taxa, with corresponding site numbers
# but this will work 

diversity.original <- diversity
#diversity <- diversity.original
diversity$SS_original <- diversity$SS
diversity$SSS_original <- diversity$SSS
diversity$SS <- as.factor(paste(diversity$SS_original, diversity$taxon_of_interest))
diversity$SSS <- as.factor(paste(diversity$SSS_original, diversity$taxon_of_interest))





# Make Zone - tropical versus temperate

names(diversity)


diversity$Zone <- rep("Temperate", length(diversity[,1]))
diversity$Zone[which(diversity$Latitude > -23.4378 & diversity$Latitude <  23.4378)] <- "Tropical"

#View(diversity[1:5000,c("Country", "Zone")])
#View(diversity[which(diversity$Country == "Australia"),c("Country", "Zone")])


### would any study now be both tropical and temperate?

study.data <- split(diversity, diversity$SS)
climate.zone <- sapply(study.data, function(rows) length(unique(rows$Zone)))
which(climate.zone == 2)
# great, each study only has one level of zone





### this is where code starts being different from main taxa split code ###


### get CWM IUCN red list score for all birds, mammals, amphibs ###

comprehensively.assessed <- c("Aves", "Mammalia", "Amphibia")
Red.List.assessed <- subset(diversity, Class %in% comprehensively.assessed)



Red.List.assessed <- droplevels(Red.List.assessed)




# how many species have redlist status

s <- split(Red.List.assessed, Red.List.assessed$Parsed_name)
length(s)

test <- sapply(s, function(rows) unique(rows$IUCN_Red_List_Status))
test <- unlist(test)

names(which(is.na(test) == T)) # some of these are on the IUCN Red List but have no information here???

length(which(test == "DD"))
length(which(test == "NE"))
length(which(test == "CR"))
length(which(test == "EN"))
length(which(test == "VU"))
length(which(test == "NT"))
length(which(test == "LC"))
length(which(test == "LR/nt"))
length(which(test == "LR/cd"))
length(which(test == "LR/lc"))




#diversity.original <- diversity
# diversity <- Red.List.assessed 




### now get SiteMetrics!

d.metrics_all_amphib_mammal_bird <- SiteMetrics_taxa_split(Red.List.assessed, 
	extra.cols = c("Realm", "Zone", "taxon_of_interest", "Country"),  traits = traits)


names(d.metrics_all_amphib_mammal_bird)
length(d.metrics_all_amphib_mammal_bird [,1]) #5387 sites in total 


# write d.metrics file to add back on to site level data after ArcMap processing
setwd("N:/Documents/PREDICTS/WDPA analysis")

write.csv(d.metrics_all_amphib_mammal_bird, "d.metrics_11_14_all_amphib_mammal_bird.csv")
	









##### only endangered species ####
# d.metrics_all_amphib_mammal_bird includes all the data on amphibs, birds, mammals
# to compare species richness of only endangered, we need to get site level stats for only endangered species
# but again use only these taxonomic groups, so that if needed, we can compare to whats going on for all species in the group
# the comparison then easier to interpret as there are not missing endangered species that could completely alter results
# So keep only data on species that are endangered, i.e. in one of the endangered categories: CR EN VU

unique(diversity.no.na$IUCN_Red_List_Status)

endangered.cats <- c("CR", "EN", "VU")
endangered.sp <- subset(Red.List.assessed, IUCN_Red_List_Status %in% endangered.cats)

length(unique(endangered.sp$SSS)) #1747


#get site metrics for this data only


traits <- c('Adult_wet_mass_log10_g',
            'Geographic_range_log10_square_km',
		'Maximum_wet_mass_log10_g',
		'Vegetative_height_log10_m',
		'Length_derived_volume_3log10_mm',
		'Mass_log10_g')



endangered.sp <- droplevels(endangered.sp)
names(endangered.sp)
d.metrics_VU_EN_CR <- SiteMetrics_taxa_split(endangered.sp, 
	extra.cols = c("Realm", "Zone", "taxon_of_interest", "Country"),  traits = traits)


# write this to file
setwd("N:/Documents/PREDICTS/WDPA analysis")
write.csv(d.metrics_VU_EN_CR, "d.metrics_11_14_VU_EN_CR.csv")
nrow(d.metrics_VU_EN_CR) #1747







### CREATE MULTIPLE TAXA PER STUDY DATASET 

#take out the studies that only compare one species
# calculate how many species per study 
#can also do this as highest taxa level = species




study.data <- split(endangered.sp, endangered.sp$SS)
numbers <- sapply(study.data, function(rows) length(unique(rows$Taxon_number)))
names(study.data)


studies.taxa <- data.frame(SS = names(study.data), number.taxa = numbers)

setwd("N:/Documents/PREDICTS/WDPA analysis")
write.csv(studies.taxa, "Number of taxa per study VU_EN_CR split_taxa_coarse 11_2014.csv")











