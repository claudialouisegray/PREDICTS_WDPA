# PLOTS
# Within_PA

load("\\\\smbhome.uscs.susx.ac.uk\\clg32\\Documents\\PREDICTS\\WDPA analysis\\RData files\\rarefied.richness_with_block_and_keeping_confounding_vars.RData")

library(reshape)
library(lme4)
library(RColorBrewer)

setwd("R:\\ecocon_d\\clg32\\GitHub\\PREDICTS_WDPA")
source("addHistogram.R")
source("addDataDistribution.R")

model.data <- Richness_rarefied.model2$data

data <- multiple.taxa.matched.landuse[,c("Richness_rarefied", "Zone", "taxon_of_interest", "ag_suit", "log_slope", "slope", 
	"DoP.PA", "AREA.PA", "log_AREA.PA", "IUCN.PA",
	 "log_elevation", "elevation", "log_bound_dist_km_PA_neg", "bound_dist_km_PA_neg", "Within_PA", "Predominant_habitat", "SS", "SSBS")]

data <- na.omit(data) 

names(matched.landuse)

nrow(data)
length(unique(data$SS))

nrow(model.data)
length(unique(model.data$SS))





L = 100

#make colours

display.brewer.all()
cols <- brewer.pal(8, "Paired")
display.brewer.pal(8, "Paired")
#drop the red to avoid red-green colorblindness issues

taxa.cols <- cols[c(4,2,8)]
taxa.cols.ci <- c("#33A02C44", "#1F78B444", "#FF7F0044")
taxa <- c("Plants", "Invertebrates", "Vertebrates")

cols <- brewer.pal(8, "Paired")
#drop the red to avoid red-green colorblindness issues

zone.cols <- cols[c(6,2)]
zone.cols.ci <- c("#E31A1C44", "#1F78B444")

inside.col <- cols[4]
inside.col.ci <- "#33A02C44"
outside.col <- 1
outside.col.ci <- "#33333344"


#land use colours for 6 landuses
#lu <- c("Primary Vegetation", "Secondary Vegetation", "Plantation forest", "Cropland", "Pasture", "Urban")
#lu.cols = c("#5B8A3B", "#1B9E77", "#7570B3", "#E6AB02", "#D95F02", "#E7298A")
#lu.cols2 = c("#66A61E", "#8ecfbc", "#7570B3","#E6AB02","#D95F02", "#E7298A")
#lu.cols2.ci <- c("#66A61E90","#8ecfbc90","#7570B390","#E6AB0290","#D95F0290","#E7298A90")

#landuse colours for 8 land uses
lu <- c("Primary Vegetation", "Mature secondary vegetation", "Intermediate secondary vegetation", "Young secondary vegetation",
	"Plantation forest", "Cropland", "Pasture", "Urban")
lu.cols = c("#5B8A3B","#147659", "#1B9E77","#8ecfbc", "#7570B3", "#E6AB02", "#D95F02", "#E7298A")
lu.cols2.ci <- c("#5B8A3B90","#14765990", "#1B9E7790","#8ecfbc90", "#7570B390", "#E6AB0290", "#D95F0290", "#E7298A90")


ylims <- c(0.15,0.22)
#slope.lim <- log(c(0,25)+1)
#elev.lim <- c()
#size.lim <- log(c(0,10000)+1)
#age.lim <- c(0,85)

#new limits for panel plots with data distribution
slope.lim <- c(-0.3, log(100))
elev.lim <- c(-0.3,log(300000))
size.lim <- c(-0.3,log(300000))
age.lim <- c(-5,110)
ag.lim <- c(0.3,131737-500-700.5)


#get model objects
model.data$Zone<- relevel(model.data$Zone , "Tropical")
model.data$Within_PA<- relevel(model.data$Within_PA , "yes")
mam.in <- glmer(Richness_rarefied.model2$final.call, model.data, family = "poisson")

model.data$Within_PA<- relevel(model.data$Within_PA , "no")
mam.out <- glmer(Richness_rarefied.model2$final.call, model.data, family = "poisson")


#### ag suit ###########


tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/01_15/Richness_rarefied vs ag suit.tif",
	width = 18, height = 20, units = "cm", pointsize = 12, res = 300)
  par(mfrow = c(3,1))
  par(mar=c(4,4.5,1,2))
  par(mgp=c(2.5,1,0))
  
  ag <-seq(from=min(model.data$ag_suit[which(model.data$Within_PA == "yes")]),
	to=max(model.data$ag_suit[which(model.data$Within_PA == "yes")]),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  
 for (x in 1:L)
  {
    Richness_rarefied <-0
    log_elevation <- mean(model.data$log_elevation)
    log_AREA.PA <- mean(model.data$log_AREA.PA)
    DoP.PA <- mean(model.data$DoP.PA)
    IUCN.PA <- "1.5"
    Zone <-"Tropical"
    taxon_of_interest<- "Vertebrates"
    Predominant_habitat <-"Primary Vegetation"
    ag_suit <- ag[x]
    log_slope <- mean(model.data$log_slope)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(Richness_rarefied, log_elevation,log_AREA.PA, DoP.PA, log_slope, ag_suit))
    newdat.f<-data.frame(cbind( Zone,taxon_of_interest, Predominant_habitat, Within_PA))
    levels(newdat.f$Zone)<-levels(model.data$Zone)
    levels(newdat.f$taxon_of_interest)<-levels(model.data$taxon_of_interest)
    levels(newdat.f$Predominant_habitat)<-levels(model.data$Predominant_habitat)
    levels(newdat.f$ag_suit) <- levels(model.data$ag_suit)
    levels(newdat.f$Within_PA) <- levels(model.data$Within_PA)
    newdat<-data.frame(newdat.f,newdat.n)
    mm<-model.matrix(terms(mam.in),newdat)
    pvar1 <- diag(mm %*% tcrossprod(vcov(mam.in),mm))
    z[x]<-mm %*% fixef(mam.in)
    zu[x]<-z[x]+sqrt(pvar1)
    zl[x]<-z[x]-sqrt(pvar1) 
  }
  

  z <- exp(z)
  zu <- exp(zu)
  zl <- exp(zl)


  plot(ag,z, ylim=c(0,20), xlim = ag.lim, axes = F,
		col =outside.col, lwd = 1, 
		bty = "l", 
		type = "l",ylab = "Rarefied species richness per site ± s.e", xlab="Agricultural suitability (higher = more suitable)")
  polygon(c(ag,rev(ag)),c(zu, rev(zl)),lty=0, col = outside.col.ci)
  axis(1,at = seq(1,8,1), seq(1,8,1))
  axis(2,at = seq(0,20,5), seq(0,20,5))

ag2 <- rep(seq(0.5, 8.5, 1),each = 2)
	end <- length(ag2)-1
ag2 <- ag2[2:end]

addHistogram(data = model.data,
			var = "Predominant_habitat",
			x =   "ag_suit",
			xlim = ag.lim,
			levels = lu,
			levels.col = lu.cols2.ci,
			bar.breaks = ag2)


addHistogram(data = model.data,
			var = "taxon_of_interest",
			x =   "ag_suit",
			xlim = ag.lim,
			levels = taxa,
			levels.col = taxa.cols,
			bar.breaks = ag2)


dev.off()





###### PA area ##########

#tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/01_15/Richness_rarefied vs PA size.tif",
#	width = 12, height = 10, units = "cm", pointsize = 12, res = 300)

tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/01_15/Richness_rarefied vs PA size EXTRA.tif",
	width = 20, height = 20, units = "cm", pointsize = 12, res = 300)

  par(mfrow = c(3,1))
  par(mar=c(4,4.5,4,1.5))
  par(mgp=c(2.5,1,0))
  
  gis <-seq(from=min(model.data$log_AREA.PA),to=max(model.data$log_AREA.PA),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)


 for (x in 1:L)
  {
    Richness_rarefied <-0
    log_elevation <- mean(model.data$log_elevation)
    log_AREA.PA <- gis[x]
    DoP.PA <- mean(model.data$DoP.PA)
    IUCN.PA <- "1.5"
    Zone <-"Tropical"
    taxon_of_interest<- "Invertebrates"
    Predominant_habitat <-"Primary Vegetation"
    ag_suit <- median(model.data$ag_suit)
    log_slope <- mean(model.data$log_slope)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(Richness_rarefied, log_elevation,log_AREA.PA, DoP.PA, log_slope, ag_suit))
    newdat.f<-data.frame(cbind( Zone,taxon_of_interest, Predominant_habitat, Within_PA))
    levels(newdat.f$Zone)<-levels(model.data$Zone)
    levels(newdat.f$taxon_of_interest)<-levels(model.data$taxon_of_interest)
    levels(newdat.f$Predominant_habitat)<-levels(model.data$Predominant_habitat)
    levels(newdat.f$ag_suit) <- levels(model.data$ag_suit)
   # levels(newdat.f$IUCN.PA) <- levels(model.data$IUCN.PA)
    levels(newdat.f$Within_PA) <- levels(model.data$Within_PA)
    newdat<-data.frame(newdat.f,newdat.n)
    mm<-model.matrix(terms(mam.in),newdat)
    pvar1 <- diag(mm %*% tcrossprod(vcov(mam.in),mm))
    z[x]<-mm %*% fixef(mam.in)
    zu[x]<-z[x]+sqrt(pvar1)
    zl[x]<-z[x]-sqrt(pvar1) 
  }
  
#  gis <- exp(gis) - 1

  z <- exp(z)
  zu <- exp(zu)
  zl <- exp(zl)



  plot(gis,z, ylim=c(0,20), xlim = size.lim, col = 1, lwd = 1,
		bty = "l", axes = F, #log = "x",
		type = "l",ylab = "Rarefied species richness per site ± s.e", xlab="PA size (km2)")
  rug(data$log_AREA.PA, ticksize = 0.03, side = 1, lwd = 1, col = 8)
  axis(1,at = log(c(0,10,100,1000,10000)+1), c(0,10,100,1000,10000))
  axis(2,at = seq(0,20,5), seq(0,20,5))
#  points(gis,zu,type="l",lty=2,  lwd = 1, col =  1)
#  points(gis,zl,type="l",lty=2,  lwd = 1, col =  1)
  polygon(c(gis,rev(gis)),c(zu, rev(zl)),lty=0, col = outside.col.ci)


addDataDistribution(data = model.data, 
			b = 50,
			x = "log_AREA.PA",
			xlim = size.lim,
			var = "Predominant_habitat",
			include.lowest = F,
			levels = lu,
			legend.spacing = 0.1,
			levels.col = lu.cols2.ci,
			axis.text.pos = log(c(0,5,10,100,1000,10000)+1),
			axis.text = c(0,5,10,100,1000,10000))



addDataDistribution(data = model.data,
			b = 50, 
			x = "log_AREA.PA",
			xlim = size.lim,
			var = "taxon_of_interest",
			include.lowest = F,
			levels = taxa,
			legend.spacing = 0.1,
			levels.col = taxa.cols,
			axis.text.pos = log(c(0,5,10,100,1000,10000)+1),
			axis.text = c(0,5,10,100,1000,10000))


dev.off()




###### elevation ##########

tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/01_15/Richness_rarefied vs elevation.tif",
	width = 20, height = 20, units = "cm", pointsize = 12, res = 300)

  par(mfrow = c(3,1))
  par(mar=c(4,4.5,4,1.5))
  par(mgp=c(2.5,1,0))

  elevation <-seq(from=min(model.data$log_elevation),to=max(model.data$log_elevation),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)

 for (x in 1:L)
  {
    Richness_rarefied <-0
    log_AREA.PA <- mean(model.data$log_AREA.PA)
    log_elevation <- elevation[x]
    DoP.PA <- mean(model.data$DoP.PA)
    IUCN.PA <- "1.5"
    Zone <-"Tropical"
    taxon_of_interest<- "Invertebrates"
    Predominant_habitat <-"Primary Vegetation"
    ag_suit <- median(model.data$ag_suit)
    log_slope <- mean(model.data$log_slope)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(Richness_rarefied,log_AREA.PA,log_elevation, DoP.PA, log_slope, ag_suit))
    newdat.f<-data.frame(cbind( Zone,taxon_of_interest, Predominant_habitat, Within_PA))
    levels(newdat.f$Zone)<-levels(model.data$Zone)
    levels(newdat.f$taxon_of_interest)<-levels(model.data$taxon_of_interest)
    levels(newdat.f$Predominant_habitat)<-levels(model.data$Predominant_habitat)
    levels(newdat.f$ag_suit) <- levels(model.data$ag_suit)
   # levels(newdat.f$IUCN.PA) <- levels(model.data$IUCN.PA)
    levels(newdat.f$Within_PA) <- levels(model.data$Within_PA)
    newdat<-data.frame(newdat.f,newdat.n)
    mm<-model.matrix(terms(mam.in),newdat)
    pvar1 <- diag(mm %*% tcrossprod(vcov(mam.in),mm))
    z[x]<-mm %*% fixef(mam.in)
    zu[x]<-z[x]+sqrt(pvar1)
    zl[x]<-z[x]-sqrt(pvar1) 
  }
  
  #elevation <- exp(elevation) - 1

  z <- exp(z)
  zu <- exp(zu)
  zl <- exp(zl)



  plot(elevation,z, ylim=c(0,20), xlim = elev.lim, col = 1, lwd = 1,
		bty = "l", axes = F,
		type = "l",ylab = "Rarefied species richness per site ± s.e", xlab="Elevation (m)")
  rug(data$log_elevation, ticksize = 0.03, side = 1, lwd = 1, col = 8)
  polygon(c(elevation,rev(elevation)),c(zu, rev(zl)),lty=0, col = outside.col.ci)
  axis(2,at = seq(0,20,5), seq(0,20,5))
  axis(1, log(c(0,1,5,50,500,2000, 5000)+1), c(0,1, 5, 50 , 500, 2000, 5000))


addDataDistribution(data = model.data, 
			b = 50,
			x = "log_elevation",
			xlim = elev.lim,
			var = "Predominant_habitat",
			include.lowest = F,
			levels = lu,
			legend.spacing = 0.1,
			levels.col = lu.cols2.ci,
			axis.text.pos = log(c(0,1,5,50,500,2000, 5000)+1),
			axis.text = c(0,1, 5, 50 , 500, 2000, 5000))



addDataDistribution(data = model.data,
			b = 50, 
			x = "log_elevation",
			xlim = elev.lim,
			var = "taxon_of_interest",
			include.lowest = F,
			levels = taxa,
			legend.spacing = 0.1,
			levels.col = taxa.cols,
			axis.text.pos = log(c(0,1,5,50,500,2000, 5000)+1),
			axis.text = c(0,1, 5, 50 , 500, 2000, 5000))


dev.off()



### slope ###

tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/01_15/Richness_rarefied vs slope.tif",
	width = 20, height = 20, units = "cm", pointsize = 12, res = 300)

  par(mfrow = c(3,1))
  par(mar=c(4,4.5,1,1.5))
  par(mgp=c(2.5,1,0))

  slope <-seq(from=min(model.data$log_slope),to=max(model.data$log_slope),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  
 for (x in 1:L)
  {
    Richness_rarefied <-0
    log_AREA.PA <- mean(model.data$log_AREA.PA)
    log_slope <- slope[x]
    DoP.PA <- mean(model.data$DoP.PA)
    IUCN.PA <- "1.5"
    Zone <-"Tropical"
    taxon_of_interest<- "Invertebrates"
    Predominant_habitat <-"Primary Vegetation"
    ag_suit <- median(model.data$ag_suit)
    log_elevation <- mean(model.data$log_elevation)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(Richness_rarefied,log_AREA.PA,log_elevation, DoP.PA, log_slope, ag_suit))
    newdat.f<-data.frame(cbind( Zone,taxon_of_interest, Predominant_habitat, Within_PA))
    levels(newdat.f$Zone)<-levels(model.data$Zone)
    levels(newdat.f$taxon_of_interest)<-levels(model.data$taxon_of_interest)
    levels(newdat.f$Predominant_habitat)<-levels(model.data$Predominant_habitat)
    levels(newdat.f$ag_suit) <- levels(model.data$ag_suit)
   # levels(newdat.f$IUCN.PA) <- levels(model.data$IUCN.PA)
    levels(newdat.f$Within_PA) <- levels(model.data$Within_PA)
    newdat<-data.frame(newdat.f,newdat.n)
    mm<-model.matrix(terms(mam.in),newdat)
    pvar1 <- diag(mm %*% tcrossprod(vcov(mam.in),mm))
    z[x]<-mm %*% fixef(mam.in)
    zu[x]<-z[x]+sqrt(pvar1)
    zl[x]<-z[x]-sqrt(pvar1) 
  }
  
  #slope <- exp(slope) - 1

  z <- exp(z)
  zu <- exp(zu)
  zl <- exp(zl)



  plot(slope,z, ylim=c(0,20), xlim = slope.lim, col = 1, lwd = 1,
		bty = "l",axes = F,
		type = "l",ylab = "Rarefied species richness per site ± s.e", xlab="Slope (degrees)")
  rug(data$log_slope, ticksize = 0.03, side = 1, lwd = 1, col = 8)
  polygon(c(slope,rev(slope)),c(zu, rev(zl)),lty=0, col = outside.col.ci)
  axis(1, log(c(0,2.5,5,10,20,30)+1), c(0,2.5,5,10,20,30))
  axis(2,at = seq(0,20,5), seq(0,20,5))

addDataDistribution(data = model.data, 
			b = 50,
			x = "log_slope",
			xlim = slope.lim,
			var = "Predominant_habitat",
			include.lowest = F,
			levels = lu,
			legend.spacing = 0.1,
			levels.col = lu.cols2.ci,
			axis.text.pos =  log(c(0,2.5,5,10,20,30)+1),
			axis.text = c(0,2.5,5,10,20,30))



addDataDistribution(data = model.data,
			b = 50, 
			x = "log_slope",
			xlim = slope.lim,
			var = "taxon_of_interest",
			include.lowest = F,
			levels = taxa,
			legend.spacing = 0.1,
			levels.col = taxa.cols,
			axis.text.pos =  log(c(0,2.5,5,10,20,30)+1),
			axis.text = c(0,2.5,5,10,20,30))




dev.off()












##### PA and non PA ######
### different taxa #### 



tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/12_14/Richness_rarefied vs PA size.tif",
	width = 12, height = 10, units = "cm", pointsize = 12, res = 300)


PA <- vector()
taxon <- vector()
estimate <- vector()
se <- vector()

model.data$taxon_of_interest <- relevel(model.data$taxon_of_interest , "Invertebrates")
model.data$Within_PA <- relevel(model.data$Within_PA, "yes")
mam <- glmer(Richness_rarefied.model2$final.call, model.data, family = "poisson")
estimate[1] <- fixef(mam)[1]
se[1] <- se.fixef(mam)[1]
PA[1] <- "Protected"
taxon[1] <- "Invertebrates"


model.data$Within_PA <- relevel(model.data$Within_PA, "no")
mam <- glmer(Richness_rarefied.model2$final.call, model.data, family = "poisson")
estimate[2] <- fixef(mam)[1]
se[2] <- se.fixef(mam)[1]
PA[2] <- "Unprotected"
taxon[2] <- "Invertebrates"

model.data$taxon_of_interest <- relevel(model.data$taxon_of_interest , "Vertebrates")
model.data$Within_PA <- relevel(model.data$Within_PA, "yes")
mam <- glmer(Richness_rarefied.model2$final.call, model.data, family = "poisson")
estimate[3] <- fixef(mam)[1]
se[3] <- se.fixef(mam)[1]
PA[3] <- "Protected"
taxon[3] <- "Vertebrates"


model.data$Within_PA <- relevel(model.data$Within_PA, "no")
mam <- glmer(Richness_rarefied.model2$final.call, model.data, family = "poisson")
estimate[4] <- fixef(mam)[1]
se[4] <- se.fixef(mam)[1]
PA[4] <- "Unprotected"
taxon[4] <- "Vertebrates"


model.data$taxon_of_interest <- relevel(model.data$taxon_of_interest , "Plants")
model.data$Within_PA <- relevel(model.data$Within_PA, "yes")
mam <- glmer(Richness_rarefied.model2$final.call, model.data, family = "poisson")
estimate[5] <- fixef(mam)[1]
se[5] <- se.fixef(mam)[1]
PA[5] <- "Protected"
taxon[5] <- "Plants and Fungi"


multiple.taxa.matched.landuse$Within_PA <- relevel(multiple.taxa.matched.landuse$Within_PA, "no")
mam <- glmer(Richness_rarefied.model2$final.call, model.data, family = "poisson")
estimate[6] <- fixef(mam)[1]
se[6] <- se.fixef(mam)[1]
PA[6] <- "Unprotected"
taxon[6] <- "Plants and Fungi"



plot.data <- data.frame(PA = PA, taxon = taxon, se = se, estimate = estimate)
plot.data$upper <- plot.data$estimate + plot.data$se
plot.data$lower <- plot.data$estimate - plot.data$se
pos <- c(1:6)

plot(plot.data$estimate ~ pos, pch = c(16,21,16,21), cex = 1.5,
	ylim = c(0, 6), 
	bty = "l", xaxt = "n", xlab = "", ylab = "Rarefied species richness per site ±se")
arrows(pos, plot.data$lower, pos, plot.data$upper, angle = 90, code = 3, length = 0.05)
axis(1, c(1.5, 3.5, 5.5), c("Invertebrates", "Vertebrates", "Plants"))
legend("topleft", pch = c(16, 21), pt.cex = 1.5, c("Protected", "Unprotected"), bty = "n")

dev.off()












### older plots incase useful ###


###### DoP vs taxon of interest ##########


tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/11_14/Richness_rarefied vs DoP vs taxon .tif",
	width = 12, height = 25, units = "cm", pointsize = 12, res = 300)


taxa <- levels(model.data$taxon_of_interest)
 t <- taxa[2]


par(mfrow = c(3,1))

for(t in taxa){ 


  par(mar=c(4,4.5,4,1.5))
  par(mgp=c(2.5,1,0))
  
 L = 30

  dop <-seq(from=min(model.data$DoP.PA[which(model.data$taxon_of_interest == t)])
		,to=max(model.data$DoP.PA[which(model.data$taxon_of_interest == t)]),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  


model.data$Zone<- relevel(model.data$Zone , "Tropical")
model.data$taxon_of_interest <- relevel(model.data$taxon_of_interest , t)
mam <- glmer(Richness_rarefied.model2$final.call, model.data, family = "poisson")

 for (x in 1:L)
  {
    log_bound_dist_km_PA_neg <-  mean(model.data$log_bound_dist_km_PA_neg)
    Richness_rarefied <-0
    log_elevation <- mean(model.data$log_elevation)
    log_AREA.PA <- mean(model.data$log_AREA.PA)
    DoP.PA <- dop[x]
    IUCN.PA <- "1.5"
    Zone <-"Tropical"
    taxon_of_interest<- "Invertebrates"
    Predominant_habitat <-"Primary Vegetation"
    ag_suit <- median(model.data$ag_suit)
    log_slope <- mean(model.data$log_slope)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, Richness_rarefied, log_elevation,log_AREA.PA, DoP.PA, log_slope, ag_suit))
    newdat.f<-data.frame(cbind( Zone,taxon_of_interest, Predominant_habitat, Within_PA))
    levels(newdat.f$Zone)<-levels(model.data$Zone)
    levels(newdat.f$taxon_of_interest)<-levels(model.data$taxon_of_interest)
    levels(newdat.f$Predominant_habitat)<-levels(model.data$Predominant_habitat)
    levels(newdat.f$ag_suit) <- levels(model.data$ag_suit)
   # levels(newdat.f$IUCN.PA) <- levels(model.data$IUCN.PA)
    levels(newdat.f$Within_PA) <- levels(model.data$Within_PA)
    newdat<-data.frame(newdat.f,newdat.n)
    mm<-model.matrix(terms(mam),newdat)
    pvar1 <- diag(mm %*% tcrossprod(vcov(mam),mm))
    z[x]<-mm %*% fixef(mam)
    zu[x]<-z[x]+sqrt(pvar1)
    zl[x]<-z[x]-sqrt(pvar1) 
  }
  

  z <- exp(z)
  zu <- exp(zu)
  zl <- exp(zl)



  plot(dop,z, ylim=c(0,17), xlim = c(0,81),
		col = taxa.cols[which(taxa == t)], lwd = 2, main =t,
		bty = "l",  
		type = "l",ylab = "Rarefied species richness per site ± s.e", xlab="Duration of protection (yr)")
  rug(model.data$DoP.PA[which(model.data$taxon_of_interest == t)], ticksize = 0.03, side = 1, lwd = 1.5, col = taxa.cols[which(taxa == t)])
  points(dop,zu,type="l",lty=2,  lwd = 2, col =  taxa.cols[which(taxa == t)])
  points(dop,zl,type="l",lty=2,  lwd = 2.5, col =  taxa.cols[which(taxa == t)])


}




dev.off()













###### DoP ##########


tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/Richness_rarefied vs DoP.tif",
	width = 12, height = 10, units = "cm", pointsize = 12, res = 300)



  par(mar=c(4,4.5,4,1.5))
  par(mgp=c(2.5,1,0))
  
 L = 30

  dop <-seq(from=min(model.data$DoP.PA),to=max(model.data$DoP.PA),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  


model.data$Zone<- relevel(model.data$Zone , "Tropical")
mam <- glmer(Richness_rarefied.model2$final.call, model.data, family = "poisson")

 for (x in 1:L)
  {
    log_bound_dist_km_PA_neg <-  mean(model.data$log_bound_dist_km_PA_neg)
    Richness_rarefied <-0
    log_elevation <- mean(model.data$log_elevation)
    log_AREA.PA <- mean(model.data$log_AREA.PA)
    DoP.PA <- dop[x]
    IUCN.PA <- "1.5"
    Zone <-"Tropical"
    taxon_of_interest<- "Invertebrates"
    Predominant_habitat <-"Primary Vegetation"
    ag_suit <- median(model.data$ag_suit)
    log_slope <- mean(model.data$log_slope)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, Richness_rarefied, log_elevation,log_AREA.PA, DoP.PA, log_slope, ag_suit))
    newdat.f<-data.frame(cbind( Zone,taxon_of_interest, Predominant_habitat, Within_PA))
    levels(newdat.f$Zone)<-levels(model.data$Zone)
    levels(newdat.f$taxon_of_interest)<-levels(model.data$taxon_of_interest)
    levels(newdat.f$Predominant_habitat)<-levels(model.data$Predominant_habitat)
    levels(newdat.f$ag_suit) <- levels(model.data$ag_suit)
   # levels(newdat.f$IUCN.PA) <- levels(model.data$IUCN.PA)
    levels(newdat.f$Within_PA) <- levels(model.data$Within_PA)
    newdat<-data.frame(newdat.f,newdat.n)
    mm<-model.matrix(terms(mam),newdat)
    pvar1 <- diag(mm %*% tcrossprod(vcov(mam),mm))
    z[x]<-mm %*% fixef(mam)
    zu[x]<-z[x]+sqrt(pvar1)
    zl[x]<-z[x]-sqrt(pvar1) 
  }
  

  z <- exp(z)
  zu <- exp(zu)
  zl <- exp(zl)



  plot(dop,z, ylim=c(0,17), col = 1, lwd = 2,
		bty = "l",  
		type = "l",ylab = "Rarefied species richness per site ± s.e", xlab="Duration of protection (yr)")
  rug(data$DoP.PA, ticksize = 0.03, side = 1, lwd = 1.5, col = 1)
  points(dop,zu,type="l",lty=2,  lwd = 2, col =  1)
  points(dop,zl,type="l",lty=2,  lwd = 2.5, col =  1)



dev.off()













###### PA size in different zones  #############
###### for different taxa ##################

tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/Richness_rarefied vs zone vs PA size vs taxon.tif",
	width = 12, height = 10, units = "cm", pointsize = 12, res = 300)



  par(mar=c(4,4.5,4,1.5))
  par(mgp=c(2.5,1,0))
  
 L = 30

  gis <-seq(from=min(model.data$log_AREA.PA),to=max(model.data$log_AREA.PA),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  

par(mfrow = c(3,1))
taxa <- levels(matched.landuse$taxon_of_interest)
i <- 0
t <- taxa[1]

for(t in taxa){

model.data$taxon_of_interest <- relevel(model.data$taxon_of_interest, t)
model.data$Zone<- relevel(model.data$Zone , "Tropical")
mam <- glmer(Richness_rarefied.model$final.call, model.data, family = "poisson")

 for (x in 1:L)
  {
    log_bound_dist_km_PA_neg <-  mean(model.data$log_bound_dist_km_PA_neg)
    Richness_rarefied <-0
    log_elevation <- mean(model.data$log_elevation)
    log_AREA.PA <-gis[x]
    DoP.PA <- mean(model.data$DoP.PA)
    IUCN.PA <- "1.5"
    Zone <-"Tropical"
    taxon_of_interest<- t
    Predominant_habitat <-"Primary Vegetation"
    ag_suit <- median(model.data$ag_suit)
    log_slope <- mean(model.data$log_slope)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, Richness_rarefied, log_elevation,log_AREA.PA, DoP.PA, log_slope, ag_suit))
    newdat.f<-data.frame(cbind( Zone,taxon_of_interest, Predominant_habitat, IUCN.PA))
    levels(newdat.f$Zone)<-levels(model.data$Zone)
    levels(newdat.f$taxon_of_interest)<-levels(model.data$taxon_of_interest)
    levels(newdat.f$Predominant_habitat)<-levels(model.data$Predominant_habitat)
    levels(newdat.f$ag_suit) <- levels(model.data$ag_suit)
    levels(newdat.f$IUCN.PA) <- levels(model.data$IUCN.PA)
   # levels(newdat.f$Within_PA) <- levels(model.data$Within_PA)
    newdat<-data.frame(newdat.f,newdat.n)
    mm<-model.matrix(terms(mam),newdat)
    pvar1 <- diag(mm %*% tcrossprod(vcov(mam),mm))
    z[x]<-mm %*% fixef(mam)
    zu[x]<-z[x]+sqrt(pvar1)
    zl[x]<-z[x]-sqrt(pvar1) 
  }
  

  z <- exp(z)
  zu <- exp(zu)
  zl <- exp(zl)

 gis_plot <- exp(gis) -1


  plot(gis_plot,z, ylim=c(0,25), col = zone.cols[1], lwd = 2, main = t,
		bty = "l", log = "x", #yaxt = "n", 
		type = "l",ylab = "Richness_rarefied size Log10(sq km) ± s.e", xlab="PA size (km2)")
  rug(data$AREA.PA[which(data$Zone == "Tropical" & data$taxon_of_interest == t)]
	, ticksize = 0.03, side = 1, lwd = 1.5, col = zone.cols[1])
  points(gis_plot,zu,type="l",lty=2,  lwd = 2, col =  zone.cols[1])
  points(gis_plot,zl,type="l",lty=2,  lwd = 2.5, col =  zone.cols[1])



  y<-vector(mode="numeric",length=L)
  yu<-vector(mode="numeric",length=L)
  yl<-vector(mode="numeric",length=L)
  

model.data$taxon_of_interest <- relevel(model.data$taxon_of_interest, t)
model.data$Zone<- relevel(model.data$Zone , "Temperate")
mam <- glmer(Richness_rarefied.model$final.call, model.data, family = "poisson")

 for (x in 1:L)
  {
    log_bound_dist_km_PA_neg <-  mean(model.data$log_bound_dist_km_PA_neg)
    Richness_rarefied <-0
    log_elevation <- mean(model.data$log_elevation)
    log_AREA.PA <-gis[x]
    DoP.PA <- mean(model.data$DoP.PA)
    IUCN.PA <- "1.5"
    Zone <-"Tropical"
    taxon_of_interest<- t
    Predominant_habitat <-"Primary Vegetation"
    ag_suit <- median(model.data$ag_suit)
    log_slope <- mean(model.data$log_slope)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, Richness_rarefied, log_elevation,log_AREA.PA, DoP.PA, log_slope, ag_suit))
    newdat.f<-data.frame(cbind( Zone,taxon_of_interest, Predominant_habitat, IUCN.PA))
    levels(newdat.f$Zone)<-levels(model.data$Zone)
    levels(newdat.f$taxon_of_interest)<-levels(model.data$taxon_of_interest)
    levels(newdat.f$Predominant_habitat)<-levels(model.data$Predominant_habitat)
    levels(newdat.f$ag_suit) <- levels(model.data$ag_suit)
    levels(newdat.f$IUCN.PA) <- levels(model.data$IUCN.PA)
   # levels(newdat.f$Within_PA) <- levels(model.data$Within_PA)
    newdat<-data.frame(newdat.f,newdat.n)
    mm<-model.matrix(terms(mam),newdat)
    pvar1 <- diag(mm %*% tcrossprod(vcov(mam),mm))
    y[x]<-mm %*% fixef(mam)
    yu[x]<-y[x]+sqrt(pvar1)
    yl[x]<-y[x]-sqrt(pvar1) 
  }
  

  y <- exp(y)
  yu <- exp(yu)
  yl <- exp(yl)


  points(gis_plot,y,type="l",lty=1, col =  zone.cols[2])
  points(gis_plot,yu,type="l",lty=2, col =  zone.cols[2])
  points(gis_plot,yl,type="l",lty=2, col =  zone.cols[2])
  rug(data$AREA.PA[which(data$Zone == "Temperate" & data$taxon_of_interest == t)],
		col = zone.cols[2], lwd = 0.8, pos = 0)

}

legend("topright", c("Tropical", "Temperate") , col = zone.cols, lty = 1, lwd = c(2,1))

dev.off()




### ag suit ####






tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/Richness_rarefied vs ag suit.tif",
	width = 12, height = 10, units = "cm", pointsize = 12, res = 300)



  par(mar=c(4,4.5,4,1.5))
  par(mgp=c(2.5,1,0))
  
 L = 30

  ag <-seq(from=min(model.data$ag_suit),to=max(model.data$ag_suit),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  




model.data$Zone<- relevel(model.data$Zone , "Tropical")
mam <- glmer(Richness_rarefied.model2$final.call, model.data, family = "poisson")

 for (x in 1:L)
  {
    log_bound_dist_km_PA_neg <-  mean(model.data$log_bound_dist_km_PA_neg)
    Richness_rarefied <-0
    log_elevation <- mean(model.data$log_elevation)
    log_AREA.PA <- mean(model.data$log_AREA.PA)
    DoP.PA <- mean(model.data$DoP.PA)
    IUCN.PA <- "1.5"
    Zone <-"Tropical"
    taxon_of_interest<- t
    Predominant_habitat <-"Primary Vegetation"
    ag_suit <- ag[x]
    log_slope <- mean(model.data$log_slope)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, Richness_rarefied, log_elevation,log_AREA.PA, DoP.PA, log_slope, ag_suit))
    newdat.f<-data.frame(cbind( Zone,taxon_of_interest, Predominant_habitat, Within_PA))
    levels(newdat.f$Zone)<-levels(model.data$Zone)
    levels(newdat.f$taxon_of_interest)<-levels(model.data$taxon_of_interest)
    levels(newdat.f$Predominant_habitat)<-levels(model.data$Predominant_habitat)
    levels(newdat.f$ag_suit) <- levels(model.data$ag_suit)
    levels(newdat.f$Within_PA) <- levels(model.data$Within_PA)
    newdat<-data.frame(newdat.f,newdat.n)
    mm<-model.matrix(terms(mam),newdat)
    pvar1 <- diag(mm %*% tcrossprod(vcov(mam),mm))
    z[x]<-mm %*% fixef(mam)
    zu[x]<-z[x]+sqrt(pvar1)
    zl[x]<-z[x]-sqrt(pvar1) 
  }
  

  z <- exp(z)
  zu <- exp(zu)
  zl <- exp(zl)


  plot(ag,z, ylim=c(0,25), col =1, lwd = 1, 
		bty = "l", 
		type = "l",ylab = "Rarefied species richness per site ± s.e", xlab="Agricultural suitability (higher = more suitable)")
  rug(data$ag_suit
	, ticksize = 0.03, side = 1, lwd = 1.5, col = 1)
  points(ag,zu,type="l",lty=2,  lwd = 1, col =  1)
  points(ag,zl,type="l",lty=2,  lwd = 1, col =  1)


dev.off()

