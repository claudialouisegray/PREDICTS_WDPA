library(foreign)

#8 level LUPA dataset####
#read in the dataframe from NHM data portal
setwd("R:/ecocon_d/shared/PREDICTS PA paper/NHM data portal")
data<-read.csv("PREDICTS_WDPA.csv")

#change column names to match analyses scripts
colnames(data) <- gsub("log_abundance", "LogAbund", colnames(data))
data$Site_ID<-factor(paste(data$Source_ID, data$Study_number, data$Site_number))

#change Predominant_habitat strings to match those used in analyses scripts
data$Predominant_habitat <-as.factor(gsub("Primary Vegetation", "Primary vegetation", data$Predominant_habitat))
data$Predominant_habitat<-relevel(data$Predominant_habitat,ref="Primary vegetation")
data$LU_3 <-as.factor(gsub("Primary Vegetation", "Primary vegetation", data$Predominant_habitat))

#change Within_PA and remake LUPA to match strings used in analyses scripts
data$Within_PA<-as.character(data$Within_PA)
data$Within_PA[which(data$Within_PA=="yes")]<-"In"
data$Within_PA[which(data$Within_PA=="no")]<-"Out"
data$Within_PA<-as.factor(data$Within_PA)
data$Within_PA<-relevel(data$Within_PA, ref="Out")
data$LUPA<-factor(paste(data$Predominant_habitat, data$Within_PA))
data$LUPA<-relevel(data$LUPA,ref="Primary vegetation Out")

#8 level LUPA dataset####
#combined dataset with 3 levels
data.com<-data
data.com$LU_3 <- gsub("Secondary", "Secondary vegetation", data.com$LU_3)
data.com$LU_3 <- gsub("Human_dominated", "Human dominated", data.com$LU_3)
data.com$Predominant_habitat<-as.factor(data.com$LU_3)
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

