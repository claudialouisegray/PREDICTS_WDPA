library(foreign)

#8 level LUPA dataset####
#read in WDPA/PREDICTS dataset
data<-read.csv("PA_11_2014.csv")

data$Site_ID<-factor(paste(data$Source_ID, data$Study_number, data$Site_number))

#remove "secondary vegetation indeterminate" and merge primary forest and non-forest
data$Predominant_habitat<-as.character(data$Predominant_habitat)
data$Predominant_habitat[which(data$Predominant_habitat=="Primary forest")]<-"Primary vegetation"
data$Predominant_habitat[which(data$Predominant_habitat=="Primary non-forest")]<-"Primary vegetation"
data$Predominant_habitat[which(data$Predominant_habitat=="Secondary vegetation (indeterminate age)")]<-NA
data<-data[!is.na(data$Predominant_habitat),]
data<-droplevels(data)
data$Predominant_habitat<-as.factor(data$Predominant_habitat)
data$Predominant_habitat<-relevel(data$Predominant_habitat,ref="Primary vegetation")

#make LUPA 
data$Within_PA<-as.character(data$Within_PA)
data$Within_PA[which(data$Within_PA=="yes")]<-"In"
data$Within_PA[which(data$Within_PA=="no")]<-"Out"
data$Within_PA<-as.factor(data$Within_PA)
data$Within_PA<-relevel(data$Within_PA, ref="Out")
data$LUPA<-factor(paste(data$Predominant_habitat, data$Within_PA))
data$LUPA<-relevel(data$LUPA,ref="Primary vegetation Out")

#LogAbund
data$LogAbund<-log(data$Total_abundance+1)

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
match.sl.site<-match(data$SSS, slope.1$SSS)
data$slope<-slope.1$MEAN[match.sl.site]
match.el.site<-match(data$SSS, elevation.1$SSS)
data$elevation<-elevation.1$MEAN[match.el.site]
match.ag.site<-match(data$SSS, ag.1$SSS)
data$ag_suit<-ag.1$ag_suit[match.ag.site]

# sort out explanatory variables
data$IUCN_CAT_number <- factor(data$IUCN_CAT_number) # they arent really in an order
data$log_slope <- log(data$slope +1)
data$log_elevation <- log(data$elevation +1)

#combined dataset with 3 levels####
data.com<-data
data.com$Predominant_habitat<-as.character(data.com$Predominant_habitat)
data.com$Predominant_habitat[which(data.com$Predominant_habitat == "Cropland")] <- "Human dominated"
data.com$Predominant_habitat[which(data.com$Predominant_habitat == "Pasture")] <- "Human dominated"
data.com$Predominant_habitat[which(data.com$Predominant_habitat == "Plantation forest")] <- "Human dominated"
data.com$Predominant_habitat[which(data.com$Predominant_habitat == "Urban")]<-"Human dominated"
secondary<-grep(pattern="econdary",x=data.com$Predominant_habitat,fixed=F)
data.com$Predominant_habitat[secondary]<-"Secondary vegetation"
data.com$Predominant_habitat<-as.factor(data.com$Predominant_habitat)
data.com$Predominant_habitat<-relevel(data.com$Predominant_habitat,ref="Primary vegetation")
#make LUPA, LUUI, LUUIPA, LUPAtax, LUPAzone
data.com<-subset(data.com,UseIntensity=="Intense use" |UseIntensity=="Light use"|UseIntensity== "Minimal use")
data.com<-droplevels(data.com)
data.com$LUUI<-factor(paste(data.com$Predominant_habitat, data.com$UseIntensity))
data.com$LUUI<-relevel(data.com$LUUI,ref="Primary vegetation Minimal use")

data.com$LUPA<-factor(paste(data.com$Predominant_habitat, data.com$Within_PA))
data.com$LUPA<-relevel(data.com$LUPA,ref="Primary vegetation Out")

data.com$LUUIPA<-factor(paste(data.com$LUUI, data.com$Within_PA))
data.com$LUUIPA<-relevel(data.com$LUUIPA,ref="Primary vegetation Minimal use Out")

data.com$LUPAtax<-factor(paste(data.com$LUPA, data.com$taxon_of_interest))
data.com$LUPAzone<-factor(paste(data.com$LUPA, data.com$Zone))

