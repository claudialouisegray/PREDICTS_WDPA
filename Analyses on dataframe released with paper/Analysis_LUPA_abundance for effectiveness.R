rm(list=ls()) 

library(yarg)
library(roquefort)
library(gamm4)


setwd("R:/ecocon_d/clg32/GitHub/PREDICTS_WDPA")
source("compare_randoms.R")
source("model_select.R")
setwd("R:/ecocon_d/clg32/PREDICTS/WDPA analysis")
PA_11_14 <- read.csv("PREDICTS_WDPA.csv")



###combine secondary for land use estimate

#make new datasets
PA_11_14_sec <- PA_11_14
PA_11_14_sec$Predominant_habitat <- gsub("Young secondary vegetation", "Secondary vegetation", PA_11_14_sec$Predominant_habitat)
PA_11_14_sec$Predominant_habitat <- gsub("Intermediate secondary vegetation", "Secondary vegetation", PA_11_14_sec$Predominant_habitat)
PA_11_14_sec$Predominant_habitat <- gsub("Mature secondary vegetation", "Secondary vegetation", PA_11_14_sec$Predominant_habitat)

multiple.taxa.PA_11_14_sec <- multiple.taxa.PA_11_14
multiple.taxa.PA_11_14_sec$Predominant_habitat <- gsub("Young secondary vegetation", "Secondary vegetation", multiple.taxa.PA_11_14_sec$Predominant_habitat)
multiple.taxa.PA_11_14_sec$Predominant_habitat <- gsub("Intermediate secondary vegetation", "Secondary vegetation", multiple.taxa.PA_11_14_sec$Predominant_habitat)
multiple.taxa.PA_11_14_sec$Predominant_habitat <- gsub("Mature secondary vegetation", "Secondary vegetation", multiple.taxa.PA_11_14_sec$Predominant_habitat)

#make a factor and set primary as reference
PA_11_14_sec$Predominant_habitat <- factor(PA_11_14_sec$Predominant_habitat)
PA_11_14_sec$Predominant_habitat <- relevel(PA_11_14_sec$Predominant_habitat, "Primary Vegetation")

#also need to make new LUPA
PA_11_14_sec$LUPA <- factor(paste(PA_11_14_sec$PA, PA_11_14_sec$Predominant_habitat))


#check non linear relationships

fF <- c("Zone", "taxon_of_interest", "Within_PA", "Predominant_habitat") 
fT <- list("ag_suit" = "3", "log_slope" = "3", "log_elevation" = "3")
keepVars <- character(0)
fI <- character(0)
RS <-  c("Within_PA")
#"log_abundance~Predominant_habitat+taxon_of_interest+Zone+(1+Within_PA|SS)+(1|SSB)"

log_abundance.best.random <- compare_randoms(PA_11_14_sec, "log_abundance",
                                             fixedFactors=fF,
                                             fixedTerms=fT,
                                             keepVars = keepVars,
                                             fixedInteractions=fI,
                                             otherRandoms=c("Predominant_habitat"),
                                             fixed_RandomSlopes = RS,
                                             fitInteractions=FALSE,
                                             verbose=TRUE)


log_abundance.best.random$best.random #


# model select

log_abundance.poly <- model_select(all.data  = PA_11_14_sec, 
                                   responseVar = "log_abundance",
                                   alpha = 0.05,
                                   fixedFactors= fF,
                                   fixedTerms= fT,
                                   keepVars = keepVars,
                                   fixedInteractions=fI,
                                   randomStruct = log_abundance.best.random$best.random,
                                   otherRandoms=c("Predominant_habitat"),
                                   verbose=TRUE)
log_abundance.poly$stats

# add interactions
fF <- c("Zone", "taxon_of_interest", "Within_PA", "Predominant_habitat") 
fT <- character(0)
keepVars <- list("ag_suit" = "1", "log_elevation" = "1", "log_slope" = "1") 
fI <- c("Within_PA:Predominant_habitat")
RS <-  c("Within_PA")
#"log_abundance~Predominant_habitat+taxon_of_interest+Zone+Within_PA:Predominant_habitat+Within_PA+poly(ag_suit,1)+poly(log_elevation,1)+poly(log_slope,1)+(1+Within_PA|SS)+(1|SSB)"


log_abundance.model.sec <- model_select(all.data  = PA_11_14_sec, 
                                        responseVar = "log_abundance",
                                        alpha = 0.05,
                                        fixedFactors= fF,
                                        fixedTerms= fT,
                                        keepVars = keepVars,
                                        fixedInteractions=fI,
                                        randomStruct = log_abundance.best.random$best.random,
                                        otherRandoms=c("Predominant_habitat"),
                                        verbose=TRUE)
summary(log_abundance.model.sec$model)




### get effectiveness estimate using proportion land use in each land use type

crop.PAs <- read.csv("R:/ecocon_d/clg32/PREDICTS/WDPA analysis/effectiveness estimates/proportion land use/crop_WDPA.txt", header = T)
pasture.PAs <- read.csv("R:/ecocon_d/clg32/PREDICTS/WDPA analysis/effectiveness estimates/proportion land use/pasture_WDPA.txt", header = T)
primary.PAs <- read.csv("R:/ecocon_d/clg32/PREDICTS/WDPA analysis/effectiveness estimates/proportion land use/primary_WDPA.txt", header = T)
secondary.PAs <- read.csv("R:/ecocon_d/clg32/PREDICTS/WDPA analysis/effectiveness estimates/proportion land use/secondary_WDPA.txt", header = T)
urban.PAs <- read.csv("R:/ecocon_d/clg32/PREDICTS/WDPA analysis/effectiveness estimates/proportion land use/urban_WDPA.txt", header = T)

crop.out <- read.csv("R:/ecocon_d/clg32/PREDICTS/WDPA analysis/effectiveness estimates/proportion land use/crop_outside.txt", header = T)
pasture.out <- read.csv("R:/ecocon_d/clg32/PREDICTS/WDPA analysis/effectiveness estimates/proportion land use/pasture_outside.txt", header = T)
primary.out <- read.csv("R:/ecocon_d/clg32/PREDICTS/WDPA analysis/effectiveness estimates/proportion land use/primary_outside.txt", header = T)
secondary.out <- read.csv("R:/ecocon_d/clg32/PREDICTS/WDPA analysis/effectiveness estimates/proportion land use/secondary_outside.txt", header = T)
urban.out <- read.csv("R:/ecocon_d/clg32/PREDICTS/WDPA analysis/effectiveness estimates/proportion land use/urban_outside.txt", header = T)



landuses <- c("Primary", "Secondary", "Pasture", "Cropland", "Urban")
LU.props <- expand.grid(landuse = landuses, outside = NA, PAs = NA)

#get proportion land use
#this equals sum of number of cells with each proportion divided by total number of cells
# total number of cells multiplied by 1000 as proportions out of 1000, not 1. 

LU.props[which(landuses == "Cropland"), which(names(LU.props) == "PAs")] <- 
  sum(crop.PAs$COUNT_*crop.PAs$VALUE_)/(sum(crop.PAs$COUNT_)*1000)
LU.props[which(landuses == "Primary"), which(names(LU.props) == "PAs")] <- 
  sum(as.numeric(primary.PAs$COUNT_*primary.PAs$VALUE_))/(sum(primary.PAs$COUNT_)*1000)
LU.props[which(landuses == "Secondary"), which(names(LU.props) == "PAs")] <- 
  sum(as.numeric(secondary.PAs$COUNT_*secondary.PAs$VALUE_))/(sum(secondary.PAs$COUNT_)*1000)
LU.props[which(landuses == "Pasture"), which(names(LU.props) == "PAs")] <- 
  sum(as.numeric(pasture.PAs$COUNT_*pasture.PAs$VALUE_))/(sum(pasture.PAs$COUNT_)*1000)
LU.props[which(landuses == "Urban"), which(names(LU.props) == "PAs")] <- 
  sum(urban.PAs$COUNT_*urban.PAs$VALUE_)/(sum(urban.PAs$COUNT_)*1000)

LU.props[which(landuses == "Cropland"), which(names(LU.props) == "outside")] <- 
  sum(as.numeric(crop.out$COUNT_*crop.out$VALUE_))/(sum(crop.out$COUNT_)*1000)
LU.props[which(landuses == "Primary"), which(names(LU.props) == "outside")] <- 
  sum(as.numeric(primary.out$COUNT_*primary.out$VALUE_))/(sum(primary.out$COUNT_)*1000)
LU.props[which(landuses == "Secondary"), which(names(LU.props) == "outside")] <- 
  sum(as.numeric(secondary.out$COUNT_*secondary.out$VALUE_))/(sum(secondary.out$COUNT_)*1000)
LU.props[which(landuses == "Pasture"), which(names(LU.props) == "outside")] <- 
  sum(as.numeric(pasture.out$COUNT_*pasture.out$VALUE_))/(sum(pasture.out$COUNT_)*1000)
LU.props[which(landuses == "Urban"), which(names(LU.props) == "outside")] <- 
  sum(urban.out$COUNT_*urban.out$VALUE_)/(sum(urban.out$COUNT_)*1000)

LU.props

sec_out <- fixef(log_abundance.model.sec$model)[
  which(names(fixef(log_abundance.model.sec$model)) == "Predominant_habitatSecondary vegetation")]
crop_out <- fixef(log_abundance.model.sec$model)[
  which(names(fixef(log_abundance.model.sec$model)) == "Predominant_habitatCropland")]
pas_out <- fixef(log_abundance.model.sec$model)[
  which(names(fixef(log_abundance.model.sec$model)) == "Predominant_habitatPasture")]
urb_out <- fixef(log_abundance.model.sec$model)[
  which(names(fixef(log_abundance.model.sec$model)) == "Predominant_habitatUrban")]

pri_in <- fixef(log_abundance.model.sec$model)[
  which(names(fixef(log_abundance.model.sec$model)) == "Within_PAyes")]
sec_in <- fixef(log_abundance.model.sec$model)[
  which(names(fixef(log_abundance.model.sec$model)) == "Within_PAyes")] +
  fixef(log_abundance.model.sec$model)[
    which(names(fixef(log_abundance.model.sec$model)) == "Predominant_habitatSecondary vegetation")] +
  fixef(log_abundance.model.sec$model)[
    which(names(fixef(log_abundance.model.sec$model)) == "Predominant_habitatSecondary vegetation:Within_PAyes")]
crop_in <- fixef(log_abundance.model.sec$model)[
  which(names(fixef(log_abundance.model.sec$model)) == "Within_PAyes")] +
  fixef(log_abundance.model.sec$model)[
    which(names(fixef(log_abundance.model.sec$model)) == "Predominant_habitatCropland")] +
  fixef(log_abundance.model.sec$model)[
    which(names(fixef(log_abundance.model.sec$model)) == "Predominant_habitatCropland:Within_PAyes")]
pas_in <- fixef(log_abundance.model.sec$model)[
  which(names(fixef(log_abundance.model.sec$model)) == "Within_PAyes")] +
  fixef(log_abundance.model.sec$model)[
    which(names(fixef(log_abundance.model.sec$model)) == "Predominant_habitatPasture")] +
  fixef(log_abundance.model.sec$model)[
    which(names(fixef(log_abundance.model.sec$model)) == "Predominant_habitatPasture:Within_PAyes")]
urb_in <- fixef(log_abundance.model.sec$model)[
  which(names(fixef(log_abundance.model.sec$model)) == "Within_PAyes")] +
  fixef(log_abundance.model.sec$model)[
    which(names(fixef(log_abundance.model.sec$model)) == "Predominant_habitatUrban")] +
  fixef(log_abundance.model.sec$model)[
    which(names(fixef(log_abundance.model.sec$model)) == "Predominant_habitatUrban:Within_PAyes")]

var.cov <- vcov(log_abundance.model.sec$model)

#max values for CI

sec_out_max <- sec_out +
  2*se.fixef(log_abundance.model.sec$model)[
    which(names(fixef(log_abundance.model.sec$model)) == "Predominant_habitatSecondary vegetation")]
crop_out_max <- crop_out+
  2*se.fixef(log_abundance.model.sec$model)[
    which(names(fixef(log_abundance.model.sec$model)) == "Predominant_habitatCropland")]
pas_out_max <- pas_out+
  2*se.fixef(log_abundance.model.sec$model)[
    which(names(fixef(log_abundance.model.sec$model)) == "Predominant_habitatPasture")]
urb_out_max <- urb_out+
  2*se.fixef(log_abundance.model.sec$model)[
    which(names(fixef(log_abundance.model.sec$model)) == "Predominant_habitatUrban")]

pri_in_max <- pri_in+
  2*se.fixef(log_abundance.model.sec$model)[
    which(names(fixef(log_abundance.model.sec$model)) == "Within_PAyes")]

cov<-var.cov[which(row.names(var.cov)=="Within_PAyes"),
             which(names(var.cov[1,])=="Predominant_habitatSecondary vegetation:Within_PAyes")]
se1<-se.fixef(log_abundance.model.sec$model)[which(names(fixef(log_abundance.model.sec$model)) == "Within_PAyes")]
se2<-se.fixef(log_abundance.model.sec$model)[which(names(fixef(log_abundance.model.sec$model)) == "Predominant_habitatSecondary vegetation:Within_PAyes")]
se.sec_in<-sqrt(se1^2+se2^2+2*cov)
sec_in_max <- sec_in + 2*se.sec_in

cov<-var.cov[which(row.names(var.cov)=="Within_PAyes"),
             which(names(var.cov[1,])=="Predominant_habitatCropland:Within_PAyes")]
se1<-se.fixef(log_abundance.model.sec$model)[which(names(fixef(log_abundance.model.sec$model)) == "Within_PAyes")]
se2<-se.fixef(log_abundance.model.sec$model)[which(names(fixef(log_abundance.model.sec$model)) == "Predominant_habitatCropland:Within_PAyes")]
se.crop_in<-sqrt(se1^2+se2^2+2*cov)
crop_in_max <- crop_in + 2*se.crop_in

cov<-var.cov[which(row.names(var.cov)=="Within_PAyes"),
             which(names(var.cov[1,])=="Predominant_habitatPasture:Within_PAyes")]
se1<-se.fixef(log_abundance.model.sec$model)[which(names(fixef(log_abundance.model.sec$model)) == "Within_PAyes")]
se2<-se.fixef(log_abundance.model.sec$model)[which(names(fixef(log_abundance.model.sec$model)) == "Predominant_habitatPasture:Within_PAyes")]
se.pas_in<-sqrt(se1^2+se2^2+2*cov)
pas_in_max <- pas_in + 2*se.pas_in

cov<-var.cov[which(row.names(var.cov)=="Within_PAyes"),
             which(names(var.cov[1,])=="Predominant_habitatUrban:Within_PAyes")]
se1<-se.fixef(log_abundance.model.sec$model)[which(names(fixef(log_abundance.model.sec$model)) == "Within_PAyes")]
se2<-se.fixef(log_abundance.model.sec$model)[which(names(fixef(log_abundance.model.sec$model)) == "Predominant_habitatUrban:Within_PAyes")]
se.urb_in<-sqrt(se1^2+se2^2+2*cov)
urb_in_max <- urb_in + 2*se.urb_in




# min vals for CI

sec_out_min <- sec_out-
  2*se.fixef(log_abundance.model.sec$model)[
    which(names(fixef(log_abundance.model.sec$model)) == "Predominant_habitatSecondary vegetation")]
crop_out_min <- crop_out -
  2*se.fixef(log_abundance.model.sec$model)[
    which(names(fixef(log_abundance.model.sec$model)) == "Predominant_habitatCropland")]
pas_out_min <- pas_out -
  2*se.fixef(log_abundance.model.sec$model)[
    which(names(fixef(log_abundance.model.sec$model)) == "Predominant_habitatPasture")]
urb_out_min <- urb_out -
  2*se.fixef(log_abundance.model.sec$model)[
    which(names(fixef(log_abundance.model.sec$model)) == "Predominant_habitatUrban")]

pri_in_min <- pri_in - 
  2*se.fixef(log_abundance.model.sec$model)[
    which(names(fixef(log_abundance.model.sec$model)) == "Within_PAyes")]
sec_in_min <- sec_in - 2*se.sec_in
crop_in_min <- crop_in - 2*se.crop_in
pas_in_min <- pas_in - 2*se.pas_in
urb_in_min <- urb_in - 2*se.urb_in


#weight the biodiversity response by proportion of area outside and inside PAs in that land use
response_out <- exp(0 + 
                      sec_out*LU.props[which(landuses == "Secondary"), which(names(LU.props) == "outside")] +
                      crop_out*LU.props[which(landuses == "Cropland"), which(names(LU.props) == "outside")] +
                      pas_out*LU.props[which(landuses == "Pasture"), which(names(LU.props) == "outside")]  + 
                      urb_out*LU.props[which(landuses == "Urban"), which(names(LU.props) == "outside")])
response_out_max <- exp(0 + 
                          sec_out_max*LU.props[which(landuses == "Secondary"), which(names(LU.props) == "outside")] +
                          crop_out_max*LU.props[which(landuses == "Cropland"), which(names(LU.props) == "outside")] +
                          pas_out_max*LU.props[which(landuses == "Pasture"), which(names(LU.props) == "outside")]  + 
                          urb_out_max*LU.props[which(landuses == "Urban"), which(names(LU.props) == "outside")])
response_out_min <- exp(0 + 
                          sec_out_min*LU.props[which(landuses == "Secondary"), which(names(LU.props) == "outside")] +
                          crop_out_min*LU.props[which(landuses == "Cropland"), which(names(LU.props) == "outside")] +
                          pas_out_min*LU.props[which(landuses == "Pasture"), which(names(LU.props) == "outside")]  + 
                          urb_out_min*LU.props[which(landuses == "Urban"), which(names(LU.props) == "outside")])

response_in <- exp(pri_in*LU.props[which(landuses == "Primary"), which(names(LU.props) == "PAs")] + 
                     sec_in*LU.props[which(landuses == "Secondary"), which(names(LU.props) == "PAs")] +
                     crop_in*LU.props[which(landuses == "Cropland"), which(names(LU.props) == "PAs")] +
                     pas_in*LU.props[which(landuses == "Pasture"), which(names(LU.props) == "PAs")]  + 
                     urb_in*LU.props[which(landuses == "Urban"), which(names(LU.props) == "PAs")])
response_in_max <- exp(pri_in_max*LU.props[which(landuses == "Primary"), which(names(LU.props) == "PAs")] + 
                         sec_in_max*LU.props[which(landuses == "Secondary"), which(names(LU.props) == "PAs")] +
                         crop_in_max*LU.props[which(landuses == "Cropland"), which(names(LU.props) == "PAs")] +
                         pas_in_max*LU.props[which(landuses == "Pasture"), which(names(LU.props) == "PAs")]  + 
                         urb_in_max*LU.props[which(landuses == "Urban"), which(names(LU.props) == "PAs")])
response_in_min <- exp(pri_in_min*LU.props[which(landuses == "Primary"), which(names(LU.props) == "PAs")] + 
                         sec_in_min*LU.props[which(landuses == "Secondary"), which(names(LU.props) == "PAs")] +
                         crop_in_min*LU.props[which(landuses == "Cropland"), which(names(LU.props) == "PAs")] +
                         pas_in_min*LU.props[which(landuses == "Pasture"), which(names(LU.props) == "PAs")]  + 
                         urb_in_min*LU.props[which(landuses == "Urban"), which(names(LU.props) == "PAs")])


b.a <-  response_in/response_out -1
b.a_max <-  response_in_max/response_out_max -1
b.a_min <-  response_in_min/response_out_min -1

benefit <- as.numeric(c(b.a, b.a_max, b.a_min))		# percentage increase in metric in PAs
PA.pct <- 15.4				# percentage of total land area in PAs
# global loss of biodiversity (from Newbold et al)
global.loss <-  0.126 


NPA.rel <- 1-benefit			# relative biodiversity in unprotected sites
PA.rel <- 1					# biodiversity in protected sites
NPA.pct <- 100-PA.pct			# land area unprotected 
global.int <- 1-global.loss		# overall status of biodiversity relative to pristine


# we want NPA.abs and PA.abs - where these are the biodiv metrics in unprotected and protected relative to pristine
# simultaneous equations are

# global.int = NPA.pct/100*NPA.abs + PA.pct/100*PA.abs #overall loss is loss in PAs and nonPAs relative to pristine
# NPA.abs = PA.abs*NPA.rel			     


# so
#global.int = (NPA.pct/100)*PA.abs*NPA.rel + (PA.pct/100)*PA.abs
# which is the same as
#global.int = PA.abs*(NPA.pct/100*(NPA.rel) + PA.pct/100)

# then
PA.abs <- global.int/(PA.pct/100 + (NPA.pct/100 * (NPA.rel)))
NPA.abs <- PA.abs*NPA.rel

# if pristine is 1, then where between 1 and NPA.abs does PA.abs fall?
# get difference between PA.abs and NPA.abs as a percentage of NPA.abs
#((1-NPA.abs)-(1-PA.abs))/(1-NPA.abs)

est <- 1-(1-PA.abs)/(1-NPA.abs)

# ie
#est <- ((1-NPA.abs)-(1-PA.abs))/(1-NPA.abs)

est



