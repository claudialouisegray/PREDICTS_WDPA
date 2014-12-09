


rm(list=ls()) 

library(yarg)
library(roquefort)



### load diversity file with PREDICTS data


setwd("N:/Documents/PREDICTS/data from Lawrence")

diversity <- readRDS('diversity-2014-11-13-03-40-23.rds')




setwd("N:/Documents/PREDICTS/WDPA analysis")





# get diversity file for LUPA analysis

PA_11_14 <- read.csv("PA_11_2014.csv")

#make original SSS
PA_11_14$SSS_original <- paste(PA_11_14$Source_ID, 
						PA_11_14$Study_number,
						PA_11_14$Site_number)

diversity.LUPA <- subset(diversity, SSS %in% PA_11_14$SSS_original)
nrow(diversity.LUPA)
saveRDS(diversity.LUPA, file = "diversity.PA_11_14.rds")







# do the same for matched landuse


matched.landuse <- read.csv("matched.landuse_11_2014.csv")
names(matched.landuse)


#make original SSS

matched.landuse$SSS_original <- paste(matched.landuse$Source_ID, 
						matched.landuse$Study_number,
						matched.landuse$Site_number)


diversity.matched.landuse <- subset(diversity, SSS %in% matched.landuse$SSS_original)
nrow(diversity.matched.landuse)

saveRDS(diversity.matched.landuse, file = "diversity.matched.landuse_11_14.rds")



