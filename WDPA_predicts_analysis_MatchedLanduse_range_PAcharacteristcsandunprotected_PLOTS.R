
# PLOTS

library(reshape)

model.data <- range.model2$data

data <- matched.landuse[,c("range", "Source_ID", "Zone", "taxon_of_interest", "ag_suit", "log_elevation", "log_slope", 
	"DoP.PA", "log_AREA.PA", "AREA.PA", "IUCN.PA",
	 "slope", "elevation", "log_bound_dist_km_PA_neg", "bound_dist_km_PA_neg", "Within_PA", "Predominant_habitat", "SS", "SSBS")]

data <- na.omit(data) #have to get rid of NA to try poly on ag_suit

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


cols <- brewer.pal(8, "Paired")
#drop the red to avoid red-green colorblindness issues

zone.cols <- cols[c(6,2)]
zone.cols.ci <- c("#E31A1C44", "#1F78B444")

inside.col <- cols[4]
inside.col.ci <- "#33A02C44"
outside.col <- 1
outside.col.ci <- "#33333344"


lu <- c("Primary Vegetation", "Secondary Vegetation", "Plantation forest", "Cropland", "Pasture", "Urban")
lu.cols = c("#5B8A3B", "#1B9E77", "#7570B3", "#E6AB02", "#D95F02", "#E7298A")


ylims <- c(0.1,0.25)
slope.lim <- log(c(0,25)+1)
elev.lim <- c()
size.lim <- log(c(0,10000)+1)
age.lim <- c(0,85)


# is data dominated by invertebrates
# yes majority of studies are inverts, though more sites are verts 
 
inv <- subset(model.data, taxon_of_interest == "Invertebrates")
nrow(inv)# 1432
length(unique(inv$SS)) #66

vert <- subset(model.data, taxon_of_interest == "Vertebrates")
nrow(vert)# 1697
length(unique(vert$SS)) #45

plant <- subset(model.data, taxon_of_interest == "Plants")
nrow(plant)# 1214
length(unique(plant$SS)) #31






#### PA size, for different zones #####

tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/01_15/endemicity vs zone vs PA size.tif",
	width = 12, height = 10, units = "cm", pointsize = 12, res = 300)



  par(mar=c(4,4.5,4,1.5))
  par(mgp=c(2.5,1,0))
  
 L = 30

  gis <-seq(from=min(model.data$log_AREA.PA[which(model.data$Zone == "Tropical")]),
		to=max(model.data$log_AREA.PA[which(model.data$Zone == "Tropical")]),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  





model.data$Zone<- relevel(model.data$Zone , "Tropical")
mam <- lmer(range.model2$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

 for (x in 1:L)
  {
    range <-0
    log_elevation <- mean(model.data$log_elevation)
    log_AREA.PA <-gis[x]
    DoP.PA <- mean(model.data$DoP.PA)
    IUCN.PA <- "1.5"
    Zone <-"Tropical"
    taxon_of_interest<- "Vertebrates"
    Predominant_habitat <-"Primary Vegetation"
    ag_suit <- median(model.data$ag_suit)
    log_elevation <- mean(model.data$log_elevation)
    log_slope <- mean(model.data$log_slope)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(range, log_elevation,log_AREA.PA, DoP.PA, log_slope, ag_suit))
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
  
  z <- 1/z
  zu <- 1/zu
  zl <- 1/zl

# gis_plot <- exp(gis) -1
  gis_plot <- gis

  plot(gis_plot,z, ylim = ylims, xlim = size.lim, col = zone.cols[1], lwd = 1,
		bty = "l",  xaxt = "n", #yaxt = "n", log = "x",
		type = "l",ylab = "Endemicity (1/10^(CWM log10 range)) ± s.e", xlab="PA size (km2)")
  axis(1,at = log(c(0,10,100,1000,10000)+1), c(0,10,100,1000,10000))
  rug(data$log_AREA.PA[which(data$Zone == "Tropical")]
	, ticksize = 0.03, side = 1, lwd = 1, col = zone.cols[1])
  polygon(c(gis_plot,rev(gis_plot)),c(zu, rev(zl)),lty=0, col = zone.cols.ci[1])
#  points(gis_plot,zu,type="l",lty=2,  lwd = 1, col =  zone.cols[1])
#  points(gis_plot,zl,type="l",lty=2,  lwd = 1, col =  zone.cols[1])


 gis <-seq(from=min(model.data$log_AREA.PA[which(model.data$Zone == "Temperate")]),
		to=max(model.data$log_AREA.PA[which(model.data$Zone == "Temperate")]),length=L)

  y<-vector(mode="numeric",length=L)
  yu<-vector(mode="numeric",length=L)
  yl<-vector(mode="numeric",length=L)
  


model.data$Zone<- relevel(model.data$Zone , "Temperate")
mam <- lmer(range.model2$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

 for (x in 1:L)
  {
    range <-0
    log_elevation <- mean(model.data$log_elevation)
    log_slope<- mean(model.data$log_slope)
    log_AREA.PA <-gis[x]
    DoP.PA <- mean(model.data$DoP.PA)
    IUCN.PA <- "1.5"
    Zone <-"Tropical"
    taxon_of_interest<- "Vertebrates"
    Predominant_habitat <-"Primary Vegetation"
    ag_suit <- median(model.data$ag_suit)
    log_elevation <- mean(model.data$log_elevation)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind( range, log_elevation,log_AREA.PA, DoP.PA, log_slope, ag_suit))
    newdat.f<-data.frame(cbind( Zone,taxon_of_interest, Predominant_habitat, Within_PA))
    levels(newdat.f$Zone)<-levels(model.data$Zone)
    levels(newdat.f$taxon_of_interest)<-levels(model.data$taxon_of_interest)
    levels(newdat.f$Predominant_habitat)<-levels(model.data$Predominant_habitat)
    levels(newdat.f$ag_suit) <- levels(model.data$ag_suit)
    levels(newdat.f$Within_PA) <- levels(model.data$Within_PA)
    newdat<-data.frame(newdat.f,newdat.n)
    mm<-model.matrix(terms(mam),newdat)
    pvar1 <- diag(mm %*% tcrossprod(vcov(mam),mm))
    y[x]<-mm %*% fixef(mam)
    yu[x]<-y[x]+sqrt(pvar1)
    yl[x]<-y[x]-sqrt(pvar1) 
  }
  
  y <- 1/y
  yu <- 1/yu
  yl <- 1/yl

 # gis_plot <- exp(gis) -1
  gis_plot <- gis

  points(gis_plot,y,type="l",lty=1, col =  zone.cols[2])
#  points(gis_plot,yu,type="l",lty=2, col =  zone.cols[2])
#  points(gis_plot,yl,type="l",lty=2, col =  zone.cols[2])
  polygon(c(gis_plot,rev(gis_plot)),c(yu, rev(yl)),lty=0, col = zone.cols.ci[2])
  rug(data$log_AREA.PA[which(data$Zone == "Temperate")],
		col = zone.cols[2], lwd = 1, pos = min(ylims))



legend("topright", c("Tropical", "Temperate"), col = zone.cols, lty = 1, lwd = c(1,1))

dev.off()






#### PA age, for different zones #####

tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/01_15/endemicity vs zone vs PA age.tif",
	width = 12, height = 10, units = "cm", pointsize = 12, res = 300)



  par(mar=c(4,4.5,4,1.5))
  par(mgp=c(2.5,1,0))
  
 L = 30

  dop <-seq(from=min(model.data$DoP.PA[which(model.data$Zone == "Tropical")]),
		to=max(model.data$DoP.PA[which(model.data$Zone == "Tropical")]),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  





model.data$Zone<- relevel(model.data$Zone , "Tropical")
mam <- lmer(range.model2$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

 for (x in 1:L)
  {
    range <-0
    log_elevation <- mean(model.data$log_elevation)
    log_AREA.PA <- mean(model.data$log_AREA.PA)
    DoP.PA <- dop[x]
    IUCN.PA <- "1.5"
    Zone <-"Tropical"
    taxon_of_interest<- "Vertebrates"
    Predominant_habitat <-"Primary Vegetation"
    ag_suit <- median(model.data$ag_suit)
    log_elevation <- mean(model.data$log_elevation)
    log_slope <- mean(model.data$log_slope)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(range, log_elevation,log_AREA.PA, DoP.PA, log_slope, ag_suit))
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
  

  z <- 1/z
  zu <- 1/zu
  zl <- 1/zl


  plot(dop,z, ylim=ylims, xlim = age.lim, col = zone.cols[1], lwd = 1,
		bty = "l",  #yaxt = "n", log = "x",
		type = "l",ylab = "Endemicity (1/CWM range) ± s.e", xlab="PA age (years)")
  rug(data$DoP.PA[which(data$Zone == "Tropical")]
	, ticksize = 0.03, side = 1, lwd = 1, col = zone.cols[1])
  polygon(c(dop,rev(dop)),c(zu, rev(zl)),lty=0, col = zone.cols.ci[1])
#  points(dop,zu,type="l",lty=2,  lwd = 1, col =  zone.cols[1])
#  points(dop,zl,type="l",lty=2,  lwd = 1, col =  zone.cols[1])


 dop <-seq(from=min(model.data$DoP.PA[which(model.data$Zone == "Temperate")]),
		to=max(model.data$DoP.PA[which(model.data$Zone == "Temperate")]),length=L)

  y<-vector(mode="numeric",length=L)
  yu<-vector(mode="numeric",length=L)
  yl<-vector(mode="numeric",length=L)
  


model.data$Zone<- relevel(model.data$Zone , "Temperate")
mam <- lmer(range.model2$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

 for (x in 1:L)
  {
    range <-0
    log_elevation <- mean(model.data$log_elevation)
    log_slope<- mean(model.data$log_slope)
    log_AREA.PA <- mean(model.data$log_AREA.PA)
    DoP.PA <- dop[x]
    IUCN.PA <- "1.5"
    Zone <-"Tropical"
    taxon_of_interest<- "Vertebrates"
    Predominant_habitat <-"Primary Vegetation"
    ag_suit <- median(model.data$ag_suit)
    log_elevation <- mean(model.data$log_elevation)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind( range, log_elevation,log_AREA.PA, DoP.PA, log_slope, ag_suit))
    newdat.f<-data.frame(cbind( Zone,taxon_of_interest, Predominant_habitat, Within_PA))
    levels(newdat.f$Zone)<-levels(model.data$Zone)
    levels(newdat.f$taxon_of_interest)<-levels(model.data$taxon_of_interest)
    levels(newdat.f$Predominant_habitat)<-levels(model.data$Predominant_habitat)
    levels(newdat.f$ag_suit) <- levels(model.data$ag_suit)
    levels(newdat.f$Within_PA) <- levels(model.data$Within_PA)
    newdat<-data.frame(newdat.f,newdat.n)
    mm<-model.matrix(terms(mam),newdat)
    pvar1 <- diag(mm %*% tcrossprod(vcov(mam),mm))
    y[x]<-mm %*% fixef(mam)
    yu[x]<-y[x]+sqrt(pvar1)
    yl[x]<-y[x]-sqrt(pvar1) 
  }
  
  y <- 1/y
  yu <- 1/yu
  yl <- 1/yl


  points(dop,y,type="l",lty=1, col =  zone.cols[2])
#  points(dop,yu,type="l",lty=2, col =  zone.cols[2])
#  points(dop,yl,type="l",lty=2, col =  zone.cols[2])
  polygon(c(dop,rev(dop)),c(yu, rev(yl)),lty=0, col = zone.cols.ci[2])
  rug(data$DoP.PA[which(data$Zone == "Temperate")],
		col = zone.cols[2], lwd = 1, pos = min(ylims))



legend("topright", c("Tropical", "Temperate"), col = zone.cols, lty = 1, lwd = c(1,1))

dev.off()











####  Within PA and slope ##########



tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/01_15/endemicity within pa vs slope.tif",
	width = 12, height = 10, units = "cm", pointsize = 12, res = 300)



par(mfrow = c(1,1))


  par(mar=c(4,4.5,1,1.5))
  par(mgp=c(2.5,1,0))
  
 L = 100

 slope <-seq(from=min(model.data$log_slope[which(model.data$Within_PA == "no")]),
		to=max(model.data$log_slope[which(model.data$Within_PA == "no")]),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  


model.data$Within_PA <- relevel(model.data$Within_PA, "no")
model.data$Zone<- relevel(model.data$Zone , "Tropical")
model.data$taxon_of_interest <- relevel(model.data$taxon_of_interest , "Invertebrates")
mam <- lmer(range.model2$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

 for (x in 1:L)
  {
    log_slope <- slope[x]
    range <-0
    log_elevation <- mean(model.data$log_elevation)
    ag_suit <- median(model.data$ag_suit)
    log_AREA.PA <- mean(model.data$log_AREA.PA )
    DoP.PA <- mean(model.data$DoP.PA)
    Zone <-"Temperate"
    taxon_of_interest<- "Vertebrates"
    Predominant_habitat <-"Primary Vegetation"
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_elevation, log_AREA.PA, DoP.PA, range, log_slope, ag_suit))
    newdat.f<-data.frame(cbind( Zone,taxon_of_interest, Predominant_habitat, Within_PA))
    levels(newdat.f$Zone)<-levels(model.data$Zone)
    levels(newdat.f$taxon_of_interest)<-levels(model.data$taxon_of_interest)
    levels(newdat.f$Predominant_habitat)<-levels(model.data$Predominant_habitat)
    levels(newdat.f$Within_PA) <- levels(model.data$Within_PA)
    newdat<-data.frame(newdat.f,newdat.n)
    mm<-model.matrix(terms(mam),newdat)
    pvar1 <- diag(mm %*% tcrossprod(vcov(mam),mm))
    z[x]<-mm %*% fixef(mam)
    zu[x]<-z[x]+sqrt(pvar1)
    zl[x]<-z[x]-sqrt(pvar1) 
  }
  

  z <- 1/z
  zu <- 1/zu
  zl <- 1/zl

  plot(slope,z, ylim=ylims, xlim = slope.lim,
		col = outside.col, #xlim = c(0.001,1500),
		bty = "l",  xaxt = "n", #log = "x",
		type = "l",ylab = "Endemicity (1/CWM range) ± s.e", 
		xlab="Slope (degrees)")
   rug(data$log_slope[which(data$Within_PA == "no")]
	, ticksize = 0.03, side = 1,   lwd = 1, col = 8, pos = min(ylims)) 
  polygon(c(slope,rev(slope)),c(zu, rev(zl)),lty=0, col = outside.col.ci)
  axis(1, log(c(0,2.5,5,10,20)+1), c(0,2.5,5,10,20))


model.data$Within_PA <- relevel(model.data$Within_PA, "yes")

mam<- lmer(range.model2$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

  slope <-seq(from=min(model.data$log_slope[which(model.data$Within_PA == "yes")]),
		to=max(model.data$log_slope[which(model.data$Within_PA == "yes")]),length=L)
  y<-vector(mode="numeric",length=L)
  yu<-vector(mode="numeric",length=L)
  yl<-vector(mode="numeric",length=L)


 for (x in 1:L)
  {
    log_slope <- slope[x]
    range <-0
    log_elevation <- mean(model.data$log_elevation)
    ag_suit <- median(model.data$ag_suit)
    log_AREA.PA <- mean(model.data$log_AREA.PA )
    DoP.PA <- mean(model.data$DoP.PA)
    Zone <-"Temperate"
    taxon_of_interest<- "Vertebrates"
    Predominant_habitat <-"Primary Vegetation"
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_elevation, log_AREA.PA, DoP.PA, range, log_slope, ag_suit))
    newdat.f<-data.frame(cbind( Zone,taxon_of_interest, Predominant_habitat, Within_PA))
    levels(newdat.f$Zone)<-levels(model.data$Zone)
    levels(newdat.f$taxon_of_interest)<-levels(model.data$taxon_of_interest)
    levels(newdat.f$Predominant_habitat)<-levels(model.data$Predominant_habitat)
    levels(newdat.f$Within_PA) <- levels(model.data$Within_PA)
    newdat<-data.frame(newdat.f,newdat.n)
    mm<-model.matrix(terms(mam),newdat)
    pvar1 <- diag(mm %*% tcrossprod(vcov(mam),mm))
    y[x]<-mm %*% fixef(mam)
    yu[x]<-y[x]+sqrt(pvar1)
    yl[x]<-y[x]-sqrt(pvar1) 
  }


  y <- 1/y
  yu <- 1/yu
  yl <- 1/yl

  points(slope,y,type="l",lty=1 , col = inside.col, lwd = 1)
  polygon(c(slope,rev(slope)),c(yu, rev(yl)),lty=0, col = inside.col.ci)
   rug(data$log_slope[which(data$Within_PA == "yes")]
	, ticksize = 0.03, side = 1,   lwd = 1, col = inside.col) 





legend("topleft", c("Protected", "Unprotected") , col = c(inside.col,outside.col), lty = c(1,1))


dev.off()





#### within PA and ag suit ##########



tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/01_15/endemicity within pa vs  ag_suit.tif",
	width = 15, height = 12, units = "cm", pointsize = 12, res = 300)

par(mfrow = c(1,1))


  par(mar=c(4,4.5,4,1.5))
  par(mgp=c(2.5,1,0))
  
 L = 30

  ag <-seq(from=min(model.data$ag_suit[which(model.data$Within_PA == "no")]),
		to=max(model.data$ag_suit[which(model.data$Within_PA == "no")]),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  


# it doesnt matter what the values for taxon and zone are, as there is no interaction with these



model.data$Within_PA <- relevel(model.data$Within_PA, "no")
model.data$Zone<- relevel(model.data$Zone, "Tropical")
model.data$taxon_of_interest <- relevel(model.data$taxon_of_interest , "Invertebrates")
mam <- lmer(range.model2$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

 for (x in 1:L)
  {
    ag_suit <- ag[x]
    range <-0
    log_elevation <- mean(model.data$log_elevation)
    log_slope <- mean(model.data$log_slope)
    log_AREA.PA <- mean(model.data$log_AREA.PA )
    DoP.PA <- mean(model.data$DoP.PA )
    Zone <-"Temperate"
    taxon_of_interest<- "Invertebrates"
    Predominant_habitat <-"Primary Vegetation"
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_slope, log_AREA.PA, DoP.PA, range, log_elevation, ag_suit))
    newdat.f<-data.frame(cbind( Zone,taxon_of_interest, Predominant_habitat, Within_PA))
    levels(newdat.f$Zone)<-levels(model.data$Zone)
    levels(newdat.f$taxon_of_interest)<-levels(model.data$taxon_of_interest)
    levels(newdat.f$Predominant_habitat)<-levels(model.data$Predominant_habitat)
    levels(newdat.f$Within_PA) <- levels(model.data$Within_PA)
    newdat<-data.frame(newdat.f,newdat.n)
    mm<-model.matrix(terms(mam),newdat)
    pvar1 <- diag(mm %*% tcrossprod(vcov(mam),mm))
    z[x]<-mm %*% fixef(mam)
    zu[x]<-z[x]+sqrt(pvar1)
    zl[x]<-z[x]-sqrt(pvar1) 
  }
  
  z <- 1/z
  zu <- 1/zu
  zl <- 1/zl

  plot(ag,z, ylim=ylims, col = outside.col,
		bty = "l", #log = "x",
		type = "l",ylab = "Endemicity (1/CWM range) ± s.e", xlab="Agricultural suitability (higher = more suitable)")
  polygon(c(ag,rev(ag)),c(zu, rev(zl)),lty=0, col = outside.col.ci)
#  points(ag,zu,type="l",lty=2, col = 8)
#  points(ag,zl,type="l",lty=2, col = 8)




model.data$Within_PA <- relevel(model.data$Within_PA, "yes")

mam<- lmer(range.model2$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

  ag <-seq(from=min(model.data$ag_suit[which(model.data$Within_PA == "yes")]),
		to=max(model.data$ag_suit[which(model.data$Within_PA == "yes")]),length=L)
  y<-vector(mode="numeric",length=L)
  yu<-vector(mode="numeric",length=L)
  yl<-vector(mode="numeric",length=L)


 for (x in 1:L)
  {
    ag_suit <- ag[x]
    range <-0
    log_elevation <- mean(model.data$log_elevation)
    log_slope <- mean(model.data$log_slope)
    log_AREA.PA <- mean(model.data$log_AREA.PA )
    DoP.PA <- mean(model.data$DoP.PA )
    Zone <-"Temperate"
    taxon_of_interest<- "Invertebrates"
    Predominant_habitat <-"Primary Vegetation"
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_slope, log_AREA.PA, DoP.PA, range, log_elevation, ag_suit))
    newdat.f<-data.frame(cbind( Zone,taxon_of_interest, Predominant_habitat, Within_PA))
    levels(newdat.f$Zone)<-levels(model.data$Zone)
    levels(newdat.f$taxon_of_interest)<-levels(model.data$taxon_of_interest)
    levels(newdat.f$Predominant_habitat)<-levels(model.data$Predominant_habitat)
    levels(newdat.f$Within_PA) <- levels(model.data$Within_PA)
    newdat<-data.frame(newdat.f,newdat.n)
    mm<-model.matrix(terms(mam),newdat)
    pvar1 <- diag(mm %*% tcrossprod(vcov(mam),mm))
    y[x]<-mm %*% fixef(mam)
    yu[x]<-y[x]+sqrt(pvar1)
    yl[x]<-y[x]-sqrt(pvar1) 
  }


  y <- 1/y
  yu <- 1/yu
  yl <- 1/yl
  
  points(ag,y,type="l",lty=1 , col = inside.col, lwd = 1)
  polygon(c(ag,rev(ag)),c(yu, rev(yl)),lty=0, col = inside.col.ci)
#  points(ag,yu,type="l",lty=3, lwd = 2, col = 1)
#  points(ag,yl,type="l",lty=3, lwd = 2, col = 1)


legend("topright", c("Protected", "Unprotected") , col = c(inside.col,outside.col), lty = c(1,1), lwd = c(1,1))


# add half violin plots to show data spread of different landuses
t1 <- data.frame()
y <- 0
ag <- seq(1,8,1)
for (i in 1:length(lu)){
	t2 <- aggregate(Predominant_habitat ~ ag_suit, data = model.data[which(model.data$Predominant_habitat == lu[i]),], length)
	names(t2) <- c("ag_suit", lu[i])
	#combine into table so missing values can be replaced with 0 where appropriate
		if(i == 1){
		t1 <- t2
		}else{
		t1 <- merge(t1,t2, by = "ag_suit", all.x = T)
	}
}

t1[is.na(t1)] <- 0

y <- 0.05
t3 <- t1
for(i in 2:ncol(t1)){
	t3[,i] <- t3[,i]/40000 + min(ylims) + y 
	polygon(c(ag,rev(ag)), c(t3[,i],rep(min(ylims)+y,length(ag))),lty=1, col = lu.cols[i-1])
	y <- y-0.01
	}

scale.bar <- 500/40000
start <- 0.165
arrows(1,start,1,start + scale.bar, length = 0.06, code = 3, angle = 90)
text(x = 1.5, y = start + scale.bar/2, pos = 4, cex = 0.7,
	" polygons below show agricultural suitability values for different landuses, \n scale bar height equals 500 sites")

legend("topleft",lu, cex = 0.7, col = lu.cols, pch = 16, bty = "n")
#x = 6.5, y = 0.185

dev.off()























### not needed ####


#### PA age #####

tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/12_14/endemicity vs DoP.tif",
	width = 12, height = 10, units = "cm", pointsize = 12, res = 300)



  par(mar=c(4,4.5,4,1.5))
  par(mgp=c(2.5,1,0))
  
 L = 100

  dop <-seq(from=min(model.data$DoP.PA),
		to=max(model.data$DoP.PA),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  

mam <- lmer(range.model$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

 for (x in 1:L)
  {
    range <-0
    log_elevation <- mean(model.data$log_elevation)
    log_AREA.PA <- mean(model.data$log_AREA.PA )
    DoP.PA <- dop[x]
    IUCN.PA <- "1.5"
    Zone <-"Tropical"
    taxon_of_interest<- "Vertebrates"
    Predominant_habitat <-"Primary Vegetation"
    ag_suit <- median(model.data$ag_suit)
    log_slope <- mean(model.data$log_slope)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind( range, log_elevation,log_AREA.PA, DoP.PA, log_slope, ag_suit))
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
  
  z <- 1/z
  zu <- 1/zu
  zl <- 1/zl



  plot(dop,z, ylim=c(0.15,0.22), lwd = 2,
		bty = "l",
		type = "l",ylab = " Endemicity (1/CWM range) ± s.e", xlab="Duration of protection (yr)")
  rug(data$DoP.PA[which(data$Zone == "Tropical")]
	, ticksize = 0.03, side = 1, lwd = 2, )
  points(dop,zu,type="l",lty=2,  lwd = 2)
  points(dop,zl,type="l",lty=2,  lwd = 2)




#### slope ##########



tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/12_14/endemicity vs slope.tif",
	width = 12, height = 10, units = "cm", pointsize = 12, res = 300)

par(mfrow = c(1,1))


  par(mar=c(4,4.5,4,1.5))
  par(mgp=c(2.5,1,0))
  
 L = 30

  slope <-seq(from=min(model.data$log_slope),to=max(model.data$log_slope),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  


# it doesnt matter what the values for taxon and zone are, as there is no interaction with these



model.data$Within_PA <- relevel(model.data$Within_PA, "no")
model.data$Zone<- relevel(model.data$Zone , "Tropical")
model.data$taxon_of_interest <- relevel(model.data$taxon_of_interest , "Invertebrates")
mam <- lmer(range.model$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

 for (x in 1:L)
  {
    log_slope <- slope[x]
    range <-0
    log_elevation <- mean(model.data$log_elevation)
    ag_suit <- median(model.data$ag_suit)
    log_AREA.PA <- mean(model.data$log_AREA.PA )
    DoP.PA <- mean(model.data$DoP.PA )
    Zone <-"Temperate"
    taxon_of_interest<- "Invertebrates"
    Predominant_habitat <-"Primary Vegetation"
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_slope, range, log_AREA.PA, DoP.PA, log_elevation, ag_suit))
    newdat.f<-data.frame(cbind( Zone,taxon_of_interest, Predominant_habitat, Within_PA))
    levels(newdat.f$Zone)<-levels(model.data$Zone)
    levels(newdat.f$taxon_of_interest)<-levels(model.data$taxon_of_interest)
    levels(newdat.f$Predominant_habitat)<-levels(model.data$Predominant_habitat)
    levels(newdat.f$Within_PA) <- levels(model.data$Within_PA)
    newdat<-data.frame(newdat.f,newdat.n)
    mm<-model.matrix(terms(mam),newdat)
    pvar1 <- diag(mm %*% tcrossprod(vcov(mam),mm))
    z[x]<-mm %*% fixef(mam)
    zu[x]<-z[x]+sqrt(pvar1)
    zl[x]<-z[x]-sqrt(pvar1) 
  }
  
  z <- 1/z
  zu <- 1/zu
  zl <- 1/zl

  slope <- exp(slope) - 1


  plot(slope,z, ylim=c(0.16,0.21), col = 1 ,
		bty = "l", log = "x",
		type = "l",ylab = "Endemicity (1/CWM range) ± s.e", xlab="Slope")
   rug(data$slope, ticksize = 0.03, side = 1,  lwd = 0.5, col = 1 ) 


  points(slope,zu,type="l",lty=2, col = 1 )
  points(slope,zl,type="l",lty=2, col = 1 )







dev.off()




####





####  Within PA and elevation ##########



tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/12_14/endemicity within pa vs elevation.tif",
	width = 12, height = 10, units = "cm", pointsize = 12, res = 300)



par(mfrow = c(1,1))


  par(mar=c(4,4.5,1,1.5))
  par(mgp=c(2.5,1,0))
  
 L = 100

  elevation <-seq(from=min(model.data$log_elevation[which(model.data$Within_PA == "no")]),
		to=max(model.data$log_elevation[which(model.data$Within_PA == "no")]),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  


# it doesnt matter what the values for taxon and zone are, as there is no interaction with these



model.data$Within_PA <- relevel(model.data$Within_PA, "no")
model.data$Zone<- relevel(model.data$Zone , "Tropical")
model.data$taxon_of_interest <- relevel(model.data$taxon_of_interest , "Invertebrates")
mam <- lmer(range.model$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

 for (x in 1:L)
  {
    log_elevation <- elevation[x]
    range <-0
    log_slope <- mean(model.data$log_slope)
    ag_suit <- median(model.data$ag_suit)
    log_AREA.PA <- mean(model.data$log_AREA.PA )
    DoP.PA <- mean(model.data$DoP.PA)
    Zone <-"Temperate"
    taxon_of_interest<- "Vertebrates"
    Predominant_habitat <-"Primary Vegetation"
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_elevation, log_AREA.PA, DoP.PA, range, log_slope, ag_suit))
    newdat.f<-data.frame(cbind( Zone,taxon_of_interest, Predominant_habitat, Within_PA))
    levels(newdat.f$Zone)<-levels(model.data$Zone)
    levels(newdat.f$taxon_of_interest)<-levels(model.data$taxon_of_interest)
    levels(newdat.f$Predominant_habitat)<-levels(model.data$Predominant_habitat)
    levels(newdat.f$Within_PA) <- levels(model.data$Within_PA)
    newdat<-data.frame(newdat.f,newdat.n)
    mm<-model.matrix(terms(mam),newdat)
    pvar1 <- diag(mm %*% tcrossprod(vcov(mam),mm))
    z[x]<-mm %*% fixef(mam)
    zu[x]<-z[x]+sqrt(pvar1)
    zl[x]<-z[x]-sqrt(pvar1) 
  }
  
  elev_plot <- (exp(elevation) -1)

  z <- 1/z
  zu <- 1/zu
  zl <- 1/zl

  plot(elev_plot,z, ylim=c(0.17,0.22),  col = outside.col, #xlim = c(0.001,1500),
		bty = "l",  xaxt = "n", log = "x",
		type = "l",ylab = "Endemicity (1/CWM range) ± s.e", 
		xlab="Elevation (m)")
   rug(data$elevation[which(data$Within_PA == "no")]
	, ticksize = 0.03, side = 1,   lwd = 1, col = outside.col) 
  axis(1, c(0.1, 1, 2, 5, 20, 100, 1000), c(0.1,1, 2, 5, 20,  100,  1000))
  points(elev_plot,zu,type="l",lty=2, col = outside.col)
  points(elev_plot,zl,type="l",lty=2, col = outside.col)


min(data$elevation[which(data$Within_PA == "no")])


model.data$Within_PA <- relevel(model.data$Within_PA, "yes")

mam<- lmer(range.model$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

  elev <-seq(from=min(model.data$log_elevation[which(model.data$Within_PA == "yes")]),
		to=max(model.data$log_elevation[which(model.data$Within_PA == "yes")]),length=L)
  y<-vector(mode="numeric",length=length(elevation))
  yu<-vector(mode="numeric",length=length(elevation))
  yl<-vector(mode="numeric",length=length(elevation))


 for (x in 1:L)
  {
    log_elevation <- elevation[x]
    range <-0
    log_slope <- mean(model.data$log_slope)
    ag_suit <- median(model.data$ag_suit)
    log_AREA.PA <- mean(model.data$log_AREA.PA )
    DoP.PA <- mean(model.data$DoP.PA)
    Zone <-"Temperate"
    taxon_of_interest<- "Vertebrates"
    Predominant_habitat <-"Primary Vegetation"
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_elevation, log_AREA.PA, DoP.PA, range, log_slope, ag_suit))
    newdat.f<-data.frame(cbind( Zone,taxon_of_interest, Predominant_habitat, Within_PA))
    levels(newdat.f$Zone)<-levels(model.data$Zone)
    levels(newdat.f$taxon_of_interest)<-levels(model.data$taxon_of_interest)
    levels(newdat.f$Predominant_habitat)<-levels(model.data$Predominant_habitat)
    levels(newdat.f$Within_PA) <- levels(model.data$Within_PA)
    newdat<-data.frame(newdat.f,newdat.n)
    mm<-model.matrix(terms(mam),newdat)
    pvar1 <- diag(mm %*% tcrossprod(vcov(mam),mm))
    y[x]<-mm %*% fixef(mam)
    yu[x]<-y[x]+sqrt(pvar1)
    yl[x]<-y[x]-sqrt(pvar1) 
  }

  elev_plot <- (exp(elevation) -1)

  y <- 1/y
  yu <- 1/yu
  yl <- 1/yl

  points(elev_plot,y,type="l",lty=1 , col = inside.col, lwd = 1)
  points(elev_plot,yu,type="l",lty=3, lwd = 1, col = inside.col)
  points(elev_plot,yl,type="l",lty=3, lwd = 1, col = inside.col)
   rug(data$elevation[which(data$Within_PA == "yes")]
	, ticksize = 0.03, side = 1,   lwd = 2, col = inside.col, pos = 0.17) 





legend("topleft", c("Protected", "Unprotected") , col = c(inside.col,outside.col), lty = c(1,1), lwd = c(2,1))


dev.off()



