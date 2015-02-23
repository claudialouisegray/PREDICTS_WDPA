# PLOTS
# distance to boundary model

load("N:\\Documents\\PREDICTS\\WDPA analysis\\RData_files\\species.richness_boundary_distance.RData")

#double check: Species_richness.model$final.call should match that in excel file of all analyses


library(reshape)
library(lme4)
library(RColorBrewer)

setwd("R:\\ecocon_d\\clg32\\GitHub\\PREDICTS_WDPA")
source("addHistogram.R")
source("addDataDistribution.R")




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


# set plot limits
ylims <- c(0,50)


slope.lim <- log(c(0,25)+1)
elev.lim <- c(0,log(5000))
size.lim <- log(c(0,10000)+1)
age.lim <- c(0,85)


#new limits for panel plots with data distribution
slope.lim <- c(-0.3, log(100))
elev.lim <- c(-0.3,log(300000))
size.lim <- c(-0.3,log(300000))
age.lim <- c(-5,110)
ag.lim <- c(0.3,13.5)




###### dist to boundary  #############
###### for different taxa ##################




model.data <- Species_richness.model$data

data <- multiple.taxa.matched.landuse[,c("Species_richness", "Zone", "taxon_of_interest", "ag_suit", "log_elevation", "elevation",
	 "log_AREA.PA", "DoP.PA", "AREA.PA",
	 "log_slope", "slope", "log_bound_dist_km_PA_neg", "bound_dist_km_PA_neg", "Within_PA", "Predominant_habitat", 
	 "SS", "SSBS", "SSB")]

data <- na.omit(data) #have to get rid of NA to try poly on ag_suit


nrow(data)
length(unique(data$SS))

nrow(model.data)
length(unique(model.data$SS))




tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/01_15/Species_richness vs taxon vs dist to boundary.tif",
	width = 25, height = 10, units = "cm", pointsize = 12, res = 300)

par(mfrow = c(1,3))

for(t in taxa){

i <- which(taxa == t)

  par(mar=c(4,4.5,4,1.5))
  par(mgp=c(2.5,1,0))
  
 L = 100

  lbd <-seq(from=min(model.data$log_bound_dist_km_PA_neg[which(model.data$taxon_of_interest == t)]),
		to=max(model.data$log_bound_dist_km_PA_neg[which(model.data$taxon_of_interest == t)]),
		length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  



model.data$Zone<- relevel(model.data$Zone , "Tropical")
model.data$taxon_of_interest <- relevel(model.data$taxon_of_interest , t)
mam <- glmer(Species_richness.model$final.call, model.data, family = "poisson")

 for (x in 1:L)
  {
    log_bound_dist_km_PA_neg <- lbd[x]
    Species_richness <-0
    log_slope <- mean(model.data$log_slope)
    log_AREA.PA <- mean(model.data$log_AREA.PA )
    DoP.PA <- mean(model.data$DoP.PA)
    IUCN.PA <- "1.5"
    Zone <-"Tropical"
    taxon_of_interest<- t
    Predominant_habitat <-"Primary Vegetation"
    ag_suit <- median(model.data$ag_suit)
    log_elevation <- mean(model.data$log_elevation)
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, Species_richness, log_slope,log_AREA.PA, DoP.PA, log_elevation, ag_suit))
    newdat.f<-data.frame(cbind( Zone,taxon_of_interest, Predominant_habitat))
    levels(newdat.f$Zone)<-levels(model.data$Zone)
    levels(newdat.f$taxon_of_interest)<-levels(model.data$taxon_of_interest)
    levels(newdat.f$Predominant_habitat)<-levels(model.data$Predominant_habitat)
    levels(newdat.f$ag_suit) <- levels(model.data$ag_suit)
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

 lbd1 <- exp(abs(lbd)) -1
 inside <- which(lbd < 0)
 lbd1[inside] <- lbd1[inside]*-1

#  plot(lbd1,z, ylim=ylims, xlim = c(-80, 200), col = cols[i],
#		bty = "l", #log = "x", #yaxt = "n", 
#		type = "l",ylab = "Species richness per site ± s.e", xlab="Distance to PA boundary (km)", main = t)
#  rug(data$bound_dist_km_PA_neg[which(data$taxon_of_interest == t)]
#	, ticksize = 0.03, side = 1, lwd = 0.5, pos = 0, col = cols[i])
#  points(lbd1,zu,type="l",lty=2, col = cols[i])
#  points(lbd1,zl,type="l",lty=2, col = cols[i])


 plot(lbd,z, ylim=ylims, xlim = c(-1*log(50+1), log(500+1)), col = taxa.cols[i], main = t,
		bty = "l", xaxt = "n", #log = "x",
		type = "l",ylab = "Species richness ± s.e", xlab=" Distance to PA boundary (km)")
 axis(1,at = log(c(0,1,10,50,200)+1), c(0,1,10,50,200)) 
 axis(1,at = -1*(log(c(1,10,50,200)+1)), c(-1,-10,-50,-200)) 
 rug(model.data$log_bound_dist_km_PA_neg[which(data$taxon_of_interest == t)], ticksize = 0.03, side = 1, lwd = 0.5, col = taxa.cols[i]) 
 polygon(c(lbd,rev(lbd)),c(zu, rev(zl)),lty=0, col = taxa.cols.ci[i])



 abline(v = 0, lty = 2, col = 8)

}



dev.off()













#PLOTs

#Within_PA

load("\\\\smbhome.uscs.susx.ac.uk\\clg32\\Documents\\PREDICTS\\WDPA analysis\\RData files\\species.richness_with_block_and_keeping_confounding_vars.RData")


model.data <- Species_richness.model2$data

data <- multiple.taxa.matched.landuse[,c("Species_richness", "Zone", "taxon_of_interest", "ag_suit", "log_elevation", "elevation", "DoP.PA", "AREA.PA", "log_AREA.PA",
	 "log_slope", "slope", "log_bound_dist_km_PA_neg", "bound_dist_km_PA_neg", "Within_PA", "Predominant_habitat", "SS", "SSBS")]

data <- na.omit(data) #have to get rid of NA to try poly on ag_suit


nrow(data)
length(unique(data$SS))

nrow(model.data)
length(unique(model.data$SS))




model.data$Within_PA <- relevel(model.data$Within_PA, "yes")
model.data$Zone<- relevel(model.data$Zone , "Tropical")
model.data$taxon_of_interest <- relevel(model.data$taxon_of_interest, "Invertebrates")
mam.in <- glmer(Species_richness.model2$final.call, model.data,  family = "poisson", control= glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))


model.data$Within_PA <- relevel(model.data$Within_PA, "no")
mam.out <- glmer(Species_richness.model2$final.call, model.data, family = "poisson")






#### within PA and ag suit ##########

tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/01_15/Species_richness within pa vs  ag_suit EXTRA.tif",
	width = 20, height = 20, units = "cm", pointsize = 12, res = 300)

#tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/01_15/Species_richness within pa vs  ag_suit.tif",
#	width = 12, height = 10, units = "cm", pointsize = 12, res = 300)


  par(mfrow = c(3,1))
  par(mar=c(4,4.5,1,1.5))
  par(mgp=c(2.5,1,0))
  
  ag <-seq(from=min(model.data$ag_suit[which(model.data$Within_PA == "no")]),
		to=max(model.data$ag_suit[which(model.data$Within_PA == "no")]),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  

 for (x in 1:L)
  {
    ag_suit <- ag[x]
    Species_richness <-0
    log_slope <- mean(model.data$log_slope)
    log_AREA.PA <-mean(model.data$log_AREA.PA)
    DoP.PA <- mean(model.data$DoP.PA)
    IUCN.PA <- "1.5"
    Zone <-"Tropical"
    taxon_of_interest<- "Vertebrates"
    Predominant_habitat <-"Primary Vegetation"
    log_elevation <- mean(model.data$log_elevation)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(Species_richness, log_slope,log_AREA.PA, DoP.PA, log_elevation, ag_suit))
    newdat.f<-data.frame(cbind( Zone,taxon_of_interest, Predominant_habitat, Within_PA))
    levels(newdat.f$Zone)<-levels(model.data$Zone)
    levels(newdat.f$taxon_of_interest)<-levels(model.data$taxon_of_interest)
    levels(newdat.f$Predominant_habitat)<-levels(model.data$Predominant_habitat)
    levels(newdat.f$ag_suit) <- levels(model.data$ag_suit)
    levels(newdat.f$Within_PA) <- levels(model.data$Within_PA)
    newdat<-data.frame(newdat.f,newdat.n)
    mm<-model.matrix(terms(mam.out),newdat)
    pvar1 <- diag(mm %*% tcrossprod(vcov(mam.out),mm))
    z[x]<-mm %*% fixef(mam.out)
    zu[x]<-z[x]+sqrt(pvar1)
    zl[x]<-z[x]-sqrt(pvar1) 
  }
  
  #backtransform from poisson log-link
  z <- exp(z)
  zu <- exp(zu)
  zl <- exp(zl)


  plot(ag,z, ylim=c(0,50), xlim = ag.lim, col = outside.col,
		bty = "l", axes = F,
		type = "l",ylab = "Species richness per site ± s.e", 
		xlab="Agricultural suitability (higher = more suitable)")
  axis(1,at = seq(1,8,1), seq(1,8,1))
  axis(2,at = seq(0,50,10), seq(0,50,10))
  #points(ag,zu,type="l",lty=2, col = outside.col)
  #points(ag,zl,type="l",lty=2, col = outside.col)

  polygon(c(ag,rev(ag)),c(zu, rev(zl)),lty=0, col = outside.col.ci)






  ag <-seq(from=min(model.data$ag_suit[which(model.data$Within_PA == "yes")]),
		to=max(model.data$ag_suit[which(model.data$Within_PA == "yes")]),length=L)
  y <-vector(mode="numeric",length=length(slope))
  yu <-vector(mode="numeric",length=length(slope))
  yl <-vector(mode="numeric",length=length(slope))


 for (x in 1:L)
  {
    ag_suit <- ag[x]
    Species_richness <-0
    log_slope <- mean(model.data$log_slope)
    log_AREA.PA <-mean(model.data$log_AREA.PA)
    DoP.PA <- mean(model.data$DoP.PA)
    IUCN.PA <- "1.5"
    Zone <-"Tropical"
    taxon_of_interest<- "Vertebrates"
    Predominant_habitat <-"Primary Vegetation"
    log_elevation <- mean(model.data$log_elevation)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(Species_richness, log_slope,log_AREA.PA, DoP.PA, log_elevation, ag_suit))
    newdat.f<-data.frame(cbind( Zone,taxon_of_interest, Predominant_habitat, Within_PA))
    levels(newdat.f$Zone)<-levels(model.data$Zone)
    levels(newdat.f$taxon_of_interest)<-levels(model.data$taxon_of_interest)
    levels(newdat.f$Predominant_habitat)<-levels(model.data$Predominant_habitat)
    levels(newdat.f$ag_suit) <- levels(model.data$ag_suit)
    levels(newdat.f$Within_PA) <- levels(model.data$Within_PA)
    newdat<-data.frame(newdat.f,newdat.n)
    mm<-model.matrix(terms(mam.in),newdat)
    pvar1 <- diag(mm %*% tcrossprod(vcov(mam.in),mm))
    y[x]<-mm %*% fixef(mam.in)
    yu[x]<-y[x]+sqrt(pvar1)
    yl[x]<-y[x]-sqrt(pvar1) 
  }

  #backtransform from poisson log-link
  y <- exp(y)
  yu <- exp(yu)
  yl <- exp(yl)

  
  points(ag,y,type="l",lty=1 , col = inside.col, lwd = 1)
 # points(ag,yu,type="l",lty=3, lwd = 2, col = inside.col)
 # points(ag,yl,type="l",lty=3, lwd = 2, col = inside.col)
  polygon(c(ag,rev(ag)),c(yu, rev(yl)),lty=0, col = inside.col.ci)


legend(x = 8.7, y = 25, c("Protected", "Unprotected") , cex = 1,
	col = c(inside.col, outside.col), bty = "n", lty = c(1,1), lwd = c(1,1))

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








#### within PA and slope ##########


#tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/01_15/Species_richness within pa vs slope.tif",
#	width = 12, height = 10, units = "cm", pointsize = 12, res = 300)

tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/01_15/Species_richness within pa vs slope EXTRA.tif",
	width = 20, height = 20, units = "cm", pointsize = 12, res = 300)

  par(mfrow = c(3,1))
  par(mar=c(4,4.5,4,1.5))
  par(mgp=c(2.5,1,0))
  
  slope <-seq(from=min(model.data$log_slope[which(model.data$Within_PA == "no")]),
		to=max(model.data$log_slope[which(model.data$Within_PA == "no")]),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  


# it doesnt matter what the values for taxon and zone are, as there is no interaction with these



 for (x in 1:L)
  {
    ag_suit <- median(model.data$ag_suit)
    Species_richness <-0
    log_slope <- slope[x]
    log_AREA.PA <-mean(model.data$log_AREA.PA)
    DoP.PA <- mean(model.data$DoP.PA)
    IUCN.PA <- "1.5"
    Zone <-"Tropical"
    taxon_of_interest<- "Vertebrates"
    Predominant_habitat <-"Primary Vegetation"
    log_elevation <- mean(model.data$log_elevation)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(Species_richness, log_slope,log_AREA.PA, DoP.PA, log_elevation, ag_suit))
    newdat.f<-data.frame(cbind( Zone,taxon_of_interest, Predominant_habitat, Within_PA))
    levels(newdat.f$Zone)<-levels(model.data$Zone)
    levels(newdat.f$taxon_of_interest)<-levels(model.data$taxon_of_interest)
    levels(newdat.f$Predominant_habitat)<-levels(model.data$Predominant_habitat)
    levels(newdat.f$slope) <- levels(model.data$slope)
    levels(newdat.f$Within_PA) <- levels(model.data$Within_PA)
    newdat<-data.frame(newdat.f,newdat.n)
    mm<-model.matrix(terms(mam.out),newdat)
    pvar1 <- diag(mm %*% tcrossprod(vcov(mam.out),mm))
    z[x]<-mm %*% fixef(mam.out)
    zu[x]<-z[x]+sqrt(pvar1)
    zl[x]<-z[x]-sqrt(pvar1) 
  }
  

 # slope <- exp(slope) -1


  #backtransform from poisson log-link
  z <- exp(z)
  zu <- exp(zu)
  zl <- exp(zl)


  plot(slope,z, ylim=ylims, xlim = slope.lim, col = outside.col,
		bty = "l", axes = F, #log = "x",
		type = "l",ylab = "Species richness per site ± s.e", 
		xlab="Slope (degrees)")
  axis(1, log(c(0,2.5,5,10,20)+1), c(0,2.5,5,10,20))
  axis(2,at = seq(0,50,10), seq(0,50,10))
  polygon(c(slope,rev(slope)),c(zu, rev(zl)),lty=0, col = outside.col.ci)
 # points(slope,zu,type="l",lty=2, col = outside.col)
 # points(slope,zl,type="l",lty=2, col = outside.col)
  rug(data$log_slope[which(data$Within_PA == "no")],
	ticksize = 0.03, side = 1,   lwd = 1, col = 8, pos = min(ylims)) 




  slope <-seq(from=min(model.data$log_slope[which(model.data$Within_PA == "yes")]),
		to=max(model.data$log_slope[which(model.data$Within_PA == "yes")]),length=L)
  y <-vector(mode="numeric",length=length(slope))
  yu <-vector(mode="numeric",length=length(slope))
  yl <-vector(mode="numeric",length=length(slope))


 for (x in 1:L)
  {
    ag_suit <- median(model.data$ag_suit)
    Species_richness <-0
    log_slope <- slope[x]
    log_AREA.PA <-mean(model.data$log_AREA.PA)
    DoP.PA <- mean(model.data$DoP.PA)
    IUCN.PA <- "1.5"
    Zone <-"Tropical"
    taxon_of_interest<- "Vertebrates"
    Predominant_habitat <-"Primary Vegetation"
    log_elevation <- mean(model.data$log_elevation)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(Species_richness, log_slope,log_AREA.PA, DoP.PA, log_elevation, ag_suit))
    newdat.f<-data.frame(cbind( Zone,taxon_of_interest, Predominant_habitat, Within_PA))
    levels(newdat.f$Zone)<-levels(model.data$Zone)
    levels(newdat.f$taxon_of_interest)<-levels(model.data$taxon_of_interest)
    levels(newdat.f$Predominant_habitat)<-levels(model.data$Predominant_habitat)
    levels(newdat.f$slope) <- levels(model.data$slope)
    levels(newdat.f$Within_PA) <- levels(model.data$Within_PA)
    newdat<-data.frame(newdat.f,newdat.n)
    mm<-model.matrix(terms(mam.in),newdat)
    pvar1 <- diag(mm %*% tcrossprod(vcov(mam.in),mm))
    y[x]<-mm %*% fixef(mam.in)
    yu[x]<-y[x]+sqrt(pvar1)
    yl[x]<-y[x]-sqrt(pvar1) 
  }


 # slope <- exp(slope) -1


  #backtransform from poisson log-link
  y <- exp(y)
  yu <- exp(yu)
  yl <- exp(yl)

  
  points(slope,y,type="l",lty=1 , col = inside.col, lwd = 1)
  polygon(c(slope,rev(slope)),c(yu, rev(yl)),lty=0, col = inside.col.ci)
#  points(slope,yu,type="l",lty=3, lwd = 2, col = inside.col)
#  points(slope,yl,type="l",lty=3, lwd = 2, col = inside.col)
  rug(data$log_slope[which(data$Within_PA == "yes")],
	ticksize = 0.03, side = 1,   lwd = 1, col = inside.col) 

legend(x = log(30), y = 40,
	 c("Protected", "Unprotected") , col = c(inside.col, outside.col), bty = "n", lty = c(1,1), lwd = c(1,1))

addDataDistribution(data = model.data, 
			b = 50,
			x = "log_slope",
			xlim = slope.lim,
			var = "Predominant_habitat",
			levels = lu,
			legend.spacing = 0.25,
			levels.col = lu.cols2.ci,
			axis.text.pos = log(c(0,2.5,5, 10, 20)+1),
			axis.text = c(0,2.5,5,10,20) )



addDataDistribution(data = model.data,
			b = 50, 
			x = "log_slope",
			xlim = slope.lim,
			var = "taxon_of_interest",
			levels = taxa,
			legend.spacing = 0.25,
			levels.col = taxa.cols,
			axis.text.pos = log(c(0,2.5,5, 10, 20)+1),
			axis.text = c(0,2.5,5,10,20) )

dev.off()







#### PA size #####

#tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/01_15/Species richness vs PA size.tif",
#	width = 15, height = 12, units = "cm", pointsize = 12, res = 300)


tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/01_15/Species richness vs PA size EXTRA.tif",
	width = 20, height = 20, units = "cm", pointsize = 12, res = 300)

  par(mfrow = c(3,1))
  par(mar=c(4,4.5,4,1.5))
  par(mgp=c(2.5,1,0))
  
  gis <-seq(from=min(model.data$log_AREA.PA),
		to=max(model.data$log_AREA.PA),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  

model.data$Zone<- relevel(model.data$Zone , "Tropical")
model.data$taxon_of_interest <- relevel(model.data$taxon_of_interest , "Invertebrates")

 for (x in 1:L)
  {
    Species_richness <-0
    log_slope <- mean(model.data$log_slope)
    log_elevation <- mean(model.data$log_elevation)
    log_AREA.PA <-gis[x]
    DoP.PA <- mean(model.data$DoP.PA)
    IUCN.PA <- "1.5"
    Zone <-"Tropical"
    taxon_of_interest<- "Vertebrates"
    Predominant_habitat <-"Primary Vegetation"
    ag_suit <- median(model.data$ag_suit)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(Species_richness, log_slope,log_AREA.PA, DoP.PA, log_elevation, ag_suit))
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

# gis_plot <- exp(gis) -1
  gis_plot <- gis

  plot(gis_plot,z, ylim=ylims, xlim = size.lim,  col = 1, lwd = 1,
		bty = "l",  axes = F, #log = "x",
		type = "l",ylab = "Species richness per site ± s.e", xlab="PA size (km2)")
  rug(data$log_AREA.PA, ticksize = 0.03, side = 1, lwd = 1, col = 8)
  axis(1,at = log(c(0,10,100,1000,10000)+1), c(0,10,100,1000,10000))
  axis(2,at = seq(0,50,10), seq(0,50,10))
  polygon(c(gis_plot,rev(gis_plot)),c(zu, rev(zl)),lty=0, col = outside.col.ci)
#  points(gis_plot,zu,type="l",lty=2,  lwd = 1, col = 1)
#  points(gis_plot,zl,type="l",lty=2,  lwd = 1, col =  1)


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




#### elevation #####

#tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/01_15/Species richness vs elevation.tif",
#	width = 15, height = 12, units = "cm", pointsize = 12, res = 300)


tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/01_15/Species richness vs elevation EXTRA.tif",
	width = 20, height = 20, units = "cm", pointsize = 12, res = 300)

  par(mfrow = c(3,1))
  par(mar=c(4,4.5,4,1.5))
  par(mgp=c(2.5,1,0))
  
  elevation <-seq(from=min(model.data$log_elevation),
		to=max(model.data$log_elevation),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  


 for (x in 1:L)
  {
    Species_richness <-0
    log_slope <- mean(model.data$log_slope)
    log_AREA.PA<- mean(model.data$log_AREA.PA)
    log_elevation <-elevation[x]
    DoP.PA <- mean(model.data$DoP.PA)
    IUCN.PA <- "1.5"
    Zone <-"Tropical"
    taxon_of_interest<- "Vertebrates"
    Predominant_habitat <-"Primary Vegetation"
    ag_suit <- median(model.data$ag_suit)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(Species_richness, log_slope,log_elevation, DoP.PA, log_AREA.PA, ag_suit))
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

# elevation_plot <- exp(elevation) -1
  elevation_plot <- elevation

  plot(elevation_plot,z, ylim=ylims, xlim = elev.lim,  col = 1, lwd = 1,
		bty = "l",  axes = F, #log = "x",
		type = "l",ylab = "Species richness per site ± s.e", xlab="Elevation (m)")
  rug(data$log_elevation, ticksize = 0.03, side = 1, lwd = 1, col = 8)
  axis(1, log(c(0,1,5,50,500,2000, 5000)+1), c(0,1, 5, 50 , 500, 2000, 5000))
  axis(2,at = seq(0,50,10), seq(0,50,10))
  polygon(c(elevation_plot,rev(elevation_plot)),c(zu, rev(zl)),lty=0, col = outside.col.ci)
#  points(elevation_plot,zu,type="l",lty=2,  lwd = 1, col = 1)
#  points(elevation_plot,zl,type="l",lty=2,  lwd = 1, col =  1)


addDataDistribution(data = model.data, 
			b = 50,
			x = "log_elevation",
			xlim = elev.lim,
			var = "Predominant_habitat",
			include.lowest = F,
			levels = lu,
			legend.spacing = 0.1,
			levels.col = lu.cols2.ci,
			axis.text.pos = log(c(0,1,5,50,500,2000,5000)+1),
			axis.text =c(0,1, 5, 50 , 500, 2000, 5000))


addDataDistribution(data = model.data, 
			b = 50,
			x = "log_elevation",
			xlim = elev.lim,
			var = "taxon_of_interest",
			include.lowest = F,
			levels = taxa,
			legend.spacing = 0.1,
			levels.col = taxa.cols,
			axis.text.pos = log(c(0,1,5,50,500,2000,5000)+1),
			axis.text =c(0,1, 5, 50 , 500, 2000, 5000))

dev.off()






### 07_14 plots ###



#### zone vs PA size #####

tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/Species richness vs zone vs PA size.tif",
	width = 15, height = 12, units = "cm", pointsize = 12, res = 300)





  par(mar=c(4,4.5,4,1.5))
  par(mgp=c(2.5,1,0))
  
 L = 30

  gis <-seq(from=min(model.data$log_AREA.PA[which(model.data$Zone == "Tropical")]),
		to=max(model.data$log_AREA.PA[which(model.data$Zone == "Tropical")]),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  





model.data$Zone<- relevel(model.data$Zone , "Tropical")
mam <- glmer(Species_richness.model$final.call, model.data, family = "poisson")

 for (x in 1:L)
  {
    log_bound_dist_km_PA_neg <-  mean(model.data$log_bound_dist_km_PA_neg)
    Species_richness <-0
    log_slope <- mean(model.data$log_slope)
    log_AREA.PA <-gis[x]
    DoP.PA <- mean(model.data$DoP.PA)
    IUCN.PA <- "1.5"
    Zone <-"Tropical"
    taxon_of_interest<- t
    Predominant_habitat <-"Primary Vegetation"
    ag_suit <- median(model.data$ag_suit)
    log_elevation <- mean(model.data$log_elevation)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, Species_richness, log_slope,log_AREA.PA, DoP.PA, log_elevation, ag_suit))
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

 gis_plot <- exp(gis) -1


  plot(gis_plot,z, ylim=c(0,35),  col = zone.cols[1], lwd = 2,
		bty = "l", log = "x", #yaxt = "n", 
		type = "l",ylab = "Species richness per site ± s.e", xlab="PA size (km2)")
  rug(data$AREA.PA[which(data$Zone == "Tropical")]
	, ticksize = 0.03, side = 1, lwd = 1.5, col = zone.cols[1])
  points(gis_plot,zu,type="l",lty=2,  lwd = 2, col =  zone.cols[1])
  points(gis_plot,zl,type="l",lty=2,  lwd = 2.5, col =  zone.cols[1])


  gis <-seq(from=min(model.data$log_AREA.PA[which(model.data$Zone == "Temperate")]),
		to=max(model.data$log_AREA.PA[which(model.data$Zone == "Temperate")]),length=L)
  y<-vector(mode="numeric",length=L)
  yu<-vector(mode="numeric",length=L)
  yl<-vector(mode="numeric",length=L)
  


model.data$Zone<- relevel(model.data$Zone , "Temperate")
mam <- glmer(Species_richness.model$final.call, model.data, family = "poisson")

 for (x in 1:L)
  {
    log_bound_dist_km_PA_neg <- lbd[x]
    Species_richness <-0
    log_slope <- mean(model.data$log_slope)
    log_AREA.PA <- mean(model.data$log_AREA.PA )
    DoP.PA <- mean(model.data$DoP.PA )
    IUCN.PA <- "1.5"
    Zone <-"Tropical"
    taxon_of_interest<- t
    Predominant_habitat <-"Primary Vegetation"
    ag_suit <- median(model.data$ag_suit)
    log_elevation <- mean(model.data$log_elevation)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, Species_richness, log_slope,log_AREA.PA, DoP.PA, log_elevation, ag_suit))
    newdat.f<-data.frame(cbind( Zone,taxon_of_interest, Predominant_habitat, Within_PA))
    levels(newdat.f$Zone)<-levels(model.data$Zone)
    levels(newdat.f$taxon_of_interest)<-levels(model.data$taxon_of_interest)
    levels(newdat.f$Within_PA) <- levels(model.data$Within_PA)
    levels(newdat.f$Predominant_habitat)<-levels(model.data$Predominant_habitat)
    levels(newdat.f$ag_suit) <- levels(model.data$ag_suit)
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

  gis_plot <- exp(gis) -1


  points(gis_plot,y,type="l",lty=2, col =  zone.cols[2])
  points(gis_plot,yu,type="l",lty=2, col =  zone.cols[2])
  points(gis_plot,yl,type="l",lty=2, col =  zone.cols[2])
  rug(data$AREA.PA[which(data$Zone == "Temperate")],
		col = zone.cols[2], lwd = 0.8, pos = 0)



legend("topright", c("Tropical", "Temperate") , col = zone.cols, lty = 1, lwd = c(2,1))


dev.off()




#### slope ##########



tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/Species_richness vs slope.tif",
	width = 12, height = 10, units = "cm", pointsize = 12, res = 300)



par(mfrow = c(1,1))


  par(mar=c(4,4.5,4,1.5))
  par(mgp=c(2.5,1,0))
  
 L = 30

  slope <-seq(from=min(model.data$log_slope),
		to=max(model.data$log_slope),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  


# it doesnt matter what the values for taxon and zone are, as there is no interaction with these



model.data$Within_PA <- relevel(model.data$Within_PA, "no")
model.data$Zone<- relevel(model.data$Zone , "Tropical")
model.data$taxon_of_interest <- relevel(model.data$taxon_of_interest , "Invertebrates")
mam <- glmer(Species_richness.model$final.call, model.data, family = "poisson")

 for (x in 1:L)
  {
    ag_suit <- median(model.data$ag_suit)
    log_bound_dist_km_PA_neg <-  mean(model.data$log_bound_dist_km_PA_neg)
    Species_richness <-0
    log_slope <- slope[x]
    log_AREA.PA <-mean(model.data$log_AREA.PA)
    DoP.PA <- mean(model.data$DoP.PA)
    IUCN.PA <- "1.5"
    Zone <-"Tropical"
    taxon_of_interest<- t
    Predominant_habitat <-"Primary Vegetation"
    log_elevation <- mean(model.data$log_elevation)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, Species_richness, log_slope,log_AREA.PA, DoP.PA, log_elevation, ag_suit))
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
  

  slope <- exp(slope) -1

  #backtransform from poisson log-link
  z <- exp(z)
  zu <- exp(zu)
  zl <- exp(zl)


  plot(slope,z, ylim=c(0, 30), col = 1,
		bty = "l", log = "x",
		type = "l",ylab = "Species richness per site ± s.e", 
		xlab="Human population density (per km2)")
  points(slope,zu,type="l",lty=2, col = 1)
  points(slope,zl,type="l",lty=2, col = 1)
  rug(data$slope,
		col = 1, lwd = 0.8, pos = 0)


dev.off()









#### PA age, for different zones #####

tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/Species richness vs zone vs DoP.tif",
	width = 12, height = 10, units = "cm", pointsize = 12, res = 300)



  par(mar=c(4,4.5,4,1.5))
  par(mgp=c(2.5,1,0))
  
 L = 30

  dop <-seq(from=min(model.data$DoP.PA),to=max(model.data$DoP.PA),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  




model.data$Zone<- relevel(model.data$Zone , "Tropical")
mam <- glmer(Species_richness.model.int$final.call, model.data, family = "poisson")


 for (x in 1:L)
  {
    log_bound_dist_km_PA_neg <- mean(model.data$log_bound_dist_km_PA_neg)
    Species_richness <-0
    log_slope <- mean(model.data$log_slope)
    log_AREA.PA <- mean(model.data$log_AREA.PA )
    DoP.PA <- dop[x]
    IUCN.PA <- "1.5"
    Zone <-"Tropical"
    taxon_of_interest<- "Vertebrates"
    Predominant_habitat <-"Primary Vegetation"
    ag_suit <- median(model.data$ag_suit)
    log_elevation <- mean(model.data$log_elevation)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, Species_richness, log_slope,log_AREA.PA, DoP.PA, log_elevation, ag_suit))
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



  plot(dop,z, ylim=c(0,35), col = zone.cols[1], lwd = 2,
		bty = "l",
		type = "l",ylab = "Species richness per site ± s.e", xlab="Duration of protection (yr)")
  rug(data$DoP.PA[which(data$Zone == "Tropical")]
	, ticksize = 0.03, side = 1, lwd = 1.5, col = zone.cols[1])
  points(dop,zu,type="l",lty=2,  lwd = 2, col =  zone.cols[1])
  points(dop,zl,type="l",lty=2,  lwd = 2, col =  zone.cols[1])



  y<-vector(mode="numeric",length=L)
  yu<-vector(mode="numeric",length=L)
  yl<-vector(mode="numeric",length=L)
  



model.data$Zone<- relevel(model.data$Zone , "Temperate")
mam <- glmer(Species_richness.model.int$final.call, model.data, family = "poisson")

 for (x in 1:L)
  {
    log_bound_dist_km_PA_neg <- lbd[x]
    Species_richness <-0
    log_slope <- mean(model.data$log_slope)
    log_AREA.PA <- mean(model.data$log_AREA.PA )
    DoP.PA <- mean(model.data$DoP.PA )
    IUCN.PA <- "1.5"
    Zone <-"Temperate"
    taxon_of_interest<- "Vertebrates"
    Predominant_habitat <-"Primary Vegetation"
    ag_suit <- median(model.data$ag_suit)
    log_elevation <- mean(model.data$log_elevation)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, Species_richness, log_slope,log_AREA.PA,  DoP.PA, log_elevation, ag_suit))
    newdat.f<-data.frame(cbind( Zone,taxon_of_interest, Predominant_habitat, IUCN.PA))
    levels(newdat.f$Zone)<-levels(model.data$Zone)
    levels(newdat.f$taxon_of_interest)<-levels(model.data$taxon_of_interest)
    levels(newdat.f$Predominant_habitat)<-levels(model.data$Predominant_habitat)
    levels(newdat.f$ag_suit) <- levels(model.data$ag_suit)
    levels(newdat.f$IUCN.PA) <- levels(model.data$IUCN.PA)
  #  levels(newdat.f$Within_PA) <- levels(model.data$Within_PA)
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

  points(dop,y,type="l",lty=2, col =  zone.cols[2])
  points(dop,yu,type="l",lty=2, col =  zone.cols[2])
  points(dop,yl,type="l",lty=2, col =  zone.cols[2])
  rug(data$DoP.PA[which(data$Zone == "Temperate")],
		col = zone.cols[2], lwd = 1, pos = 0)



legend("topright", c("Tropical", "Temperate") , col = zone.cols, lty = 1, lwd = c(2,1))

dev.off()







