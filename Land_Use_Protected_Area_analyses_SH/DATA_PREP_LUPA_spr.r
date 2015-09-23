library(foreign)

PA_11_2014<-read.csv("PA_11_2014.csv")
data<-PA_11_2014

### CREATE MULTIPLE TAXA PER STUDY DATASET 
# create dataset for species richness analysis that without all studies that are only on one taxon

studies.taxa <- read.csv("Number of taxa per study split_taxa_coarse 11_2014.csv")
which(studies.taxa$number.taxa == 1)
length(which(studies.taxa$number.taxa == 1)) #25
more.than.one.taxa <- studies.taxa$SS[which(studies.taxa$number.taxa != 1)]
multiple.taxa.PA_11_14 <- subset(data, SS %in% more.than.one.taxa)
length(multiple.taxa.PA_11_14[,1]) #6874

data.spr.rr<-multiple.taxa.PA_11_14

data.spr.rr$Site_ID<-factor(paste(data.spr.rr$Source_ID, data.spr.rr$Study_number, data.spr.rr$Site_number))

#remove "secondary vegetation indeterminate" and merge primary forest and non-forest
data.spr.rr$Predominant_habitat<-as.character(data.spr.rr$Predominant_habitat)
data.spr.rr$Predominant_habitat[which(data.spr.rr$Predominant_habitat=="Primary forest")]<-"Primary vegetation"
data.spr.rr$Predominant_habitat[which(data.spr.rr$Predominant_habitat=="Primary non-forest")]<-"Primary vegetation"
data.spr.rr$Predominant_habitat[which(data.spr.rr$Predominant_habitat=="Secondary vegetation (indeterminate age)")]<-NA
data.spr.rr<-data.spr.rr[!is.na(data.spr.rr$Predominant_habitat),]
data.spr.rr<-droplevels(data.spr.rr)
data.spr.rr$Predominant_habitat<-as.factor(data.spr.rr$Predominant_habitat)
data.spr.rr$Predominant_habitat<-relevel(data.spr.rr$Predominant_habitat,ref="Primary vegetation")

#make LUPA 
data.spr.rr$Within_PA<-as.character(data.spr.rr$Within_PA)
data.spr.rr$Within_PA[which(data.spr.rr$Within_PA=="yes")]<-"In"
data.spr.rr$Within_PA[which(data.spr.rr$Within_PA=="no")]<-"Out"
data.spr.rr$Within_PA<-as.factor(data.spr.rr$Within_PA)
data.spr.rr$Within_PA<-relevel(data.spr.rr$Within_PA, ref="Out")
data.spr.rr$LUPA<-factor(paste(data.spr.rr$Predominant_habitat, data.spr.rr$Within_PA))
data.spr.rr$LUPA<-relevel(data.spr.rr$LUPA,ref="Primary vegetation Out")

#LogAbund
data.spr.rr$LogAbund<-log(data.spr.rr$Total_abundance+1)

# LOAD DATA on confounding variables
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
match.sl.site<-match(data.spr.rr$SSS, slope.1$SSS)
data.spr.rr$slope<-slope.1$MEAN[match.sl.site]
match.el.site<-match(data.spr.rr$SSS, elevation.1$SSS)
data.spr.rr$elevation<-elevation.1$MEAN[match.el.site]
match.ag.site<-match(data.spr.rr$SSS, ag.1$SSS)
data.spr.rr$ag_suit<-ag.1$ag_suit[match.ag.site]

# sort out explanatory variables
data.spr.rr$IUCN_CAT_number <- factor(data.spr.rr$IUCN_CAT_number) # they arent really in an order
data.spr.rr$log_slope <- log(data.spr.rr$slope +1)
data.spr.rr$log_elevation <- log(data.spr.rr$elevation +1)

#combined dataset with 3 levels####
data.spr.rr.com<-data.spr.rr
data.spr.rr.com$Predominant_habitat<-as.character(data.spr.rr.com$Predominant_habitat)
data.spr.rr.com$Predominant_habitat[which(data.spr.rr.com$Predominant_habitat == "Cropland")] <- "Human dominated"
data.spr.rr.com$Predominant_habitat[which(data.spr.rr.com$Predominant_habitat == "Pasture")] <- "Human dominated"
data.spr.rr.com$Predominant_habitat[which(data.spr.rr.com$Predominant_habitat == "Plantation forest")] <- "Human dominated"
data.spr.rr.com$Predominant_habitat[which(data.spr.rr.com$Predominant_habitat == "Urban")]<-"Human dominated"
secondary<-grep(pattern="econdary",x=data.spr.rr.com$Predominant_habitat,fixed=F)
data.spr.rr.com$Predominant_habitat[secondary]<-"Secondary vegetation"
data.spr.rr.com$Predominant_habitat<-as.factor(data.spr.rr.com$Predominant_habitat)
data.spr.rr.com$Predominant_habitat<-relevel(data.spr.rr.com$Predominant_habitat,ref="Primary vegetation")
#make LUPA, LUUI, LUUIPA, LUPAtax, LUPAzone
data.spr.rr.com<-subset(data.spr.rr.com,UseIntensity=="Intense use" |UseIntensity=="Light use"|UseIntensity== "Minimal use")
data.spr.rr.com<-droplevels(data.spr.rr.com)
data.spr.rr.com$LUUI<-factor(paste(data.spr.rr.com$Predominant_habitat, data.spr.rr.com$UseIntensity))
data.spr.rr.com$LUUI<-relevel(data.spr.rr.com$LUUI,ref="Primary vegetation Minimal use")

data.spr.rr.com$LUPA<-factor(paste(data.spr.rr.com$Predominant_habitat, data.spr.rr.com$Within_PA))
data.spr.rr.com$LUPA<-relevel(data.spr.rr.com$LUPA,ref="Primary vegetation Out")

data.spr.rr.com$LUUIPA<-factor(paste(data.spr.rr.com$LUUI, data.spr.rr.com$Within_PA))
data.spr.rr.com$LUUIPA<-relevel(data.spr.rr.com$LUUIPA,ref="Primary vegetation Minimal use Out")

data.spr.rr.com$LUPAtax<-factor(paste(data.spr.rr.com$LUPA, data.spr.rr.com$taxon_of_interest))
data.spr.rr.com$LUPAzone<-factor(paste(data.spr.rr.com$LUPA, data.spr.rr.com$Zone))
