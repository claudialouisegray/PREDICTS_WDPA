
rm(list=ls()) 
library(yarg)
lirary(roquefort)

### load diversity file with PREDICTS data from other code
setwd("R:/ecocon_d/clg32/GitHub/PREDICTS_WDPA")

# get diversity file for LUPA analysis

source("prep_PA_11_14_for_analysis.R")

#make original SSS
PA_11_14$SSS_original <- paste(PA_11_14$Source_ID, 
						PA_11_14$Study_number,
						PA_11_14$Site_number)

diversity.LUPA <- subset(diversity, SSS %in% PA_11_14$SSS)
diversity.LUPA <- subset(diversity, SSS %in% PA_11_14$SSS_original)
nrow(diversity.LUPA)
saveRDS(diversity.LUPA, file = "R:/ecocon_d/clg32/PREDICTS/WDPA analysis/diversity.PA_11_14.rds")

length(unique(diversity.LUPA$SSS))
length(unique(diversity.LUPA$SS))
length(unique(diversity.LUPA$Parsed_name))
length(unique(diversity.LUPA$Country))
length(unique(diversity.LUPA$Ecoregion))
length(unique(diversity.LUPA$Hotspot))

#check number of contributers
names(diversity)

sort(table(diversity.LUPA$Source_ID)/nrow(diversity.LUPA))
#no source more than 3% of data


# do the same for matched landuse

setwd("R:/ecocon_d/clg32/GitHub/PREDICTS_WDPA")
source("prep_matched.landuse_for_analysis.R")
names(matched.landuse)


#make original SSS
matched.landuse$SSS_original <- paste(matched.landuse$Source_ID, 
						matched.landuse$Study_number,
						matched.landuse$Site_number)


diversity.matched.landuse <- subset(diversity, SSS %in% matched.landuse$SSS_original)
nrow(diversity.matched.landuse)

saveRDS(diversity.matched.landuse, file = "R:/ecocon_d/clg32/PREDICTS/WDPA analysis/diversity.matched.landuse_11_14.rds")


length(unique(diversity.matched.landuse$SSS))
length(unique(diversity.matched.landuse$SS))
length(unique(diversity.matched.landuse$Parsed_name))
length(unique(diversity.matched.landuse$Country))
length(unique(diversity.matched.landuse$Ecoregion))
length(unique(diversity.matched.landuse$Hotspot))


### get number of records per site

s <- split(diversity, diversity$SSS)
s <- split(diversity, diversity$SSS_original)

length(s)
l <- lapply(s, function(x) nrow(x))
length(l)
max(as.numeric(l))
min(as.numeric(l))

pos <- which(as.numeric(l) == 1)
to.check <- names(s)[pos
View(subset(diversity, SSS_original %in% to.check))


# some sites only have one record as they are the split ones, e.g. the hummingbirds out of pollinator study
# other sites only have one record as the study was on one particular species

mean(as.numeric(l))
sd(as.numeric(l))


