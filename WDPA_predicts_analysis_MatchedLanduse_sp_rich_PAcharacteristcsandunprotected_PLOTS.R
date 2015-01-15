# PLOTS
# distance to boundary model

load("N:\\Documents\\PREDICTS\\WDPA analysis\\RData_files\\species.richness_boundary_distance.RData")

#double check: Species_richness.model$final.call should match that in excel file of all analyses

library(RColorBrewer)
library(lme4)

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


ylims <- c(0,50)
slope.lim <- log(c(0,25)+1)
elev.lim <- c()
size.lim <- log(c(0,10000)+1)
age.lim <- c(0,85)






###### dist to boundary  #############
###### for different taxa ##################



tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/12_14/Species_richness vs taxon vs dist to boundary.tif",
	width = 25, height = 10, units = "cm", pointsize = 12, res = 300)

par(mfrow = c(1,3))
taxa <- levels(matched.landuse$taxon_of_interest)
i <- 0
t <- taxa[1]

for(t in taxa){

i <- i +1

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

  plot(lbd1,z, ylim=ylims, xlim = c(-80, 200), col = cols[i],
		bty = "l", #log = "x", #yaxt = "n", 
		type = "l",ylab = "Species richness per site ± s.e", xlab="Distance to PA boundary (km)", main = t)
  rug(data$bound_dist_km_PA_neg[which(data$taxon_of_interest == t)]
	, ticksize = 0.03, side = 1, lwd = 0.5, pos = 0, col = cols[i])
  points(lbd1,zu,type="l",lty=2, col = cols[i])
  points(lbd1,zl,type="l",lty=2, col = cols[i])

 abline(v = 0, lty = 2, col = 8)

}



dev.off()













#PLOTs

#Within_PA

load("N:\\Documents\\PREDICTS\\WDPA analysis\\RData files\\species.richness.RData")

double check:  Species_richness.model2$final.call should match that in excel file of all analyses



model.data <- Species_richness.model2$data

data <- multiple.taxa.matched.landuse[,c("Species_richness", "Zone", "taxon_of_interest", "ag_suit", "log_elevation", "elevation", "DoP.PA", "AREA.PA", "log_AREA.PA",
	 "log_slope", "slope", "log_bound_dist_km_PA_neg", "bound_dist_km_PA_neg", "Within_PA", "Predominant_habitat", "SS", "SSBS")]

data <- na.omit(data) #have to get rid of NA to try poly on ag_suit


nrow(data)
length(unique(data$SS))

nrow(model.data)
length(unique(model.data$SS))



#### within PA and ag suit ##########



tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/01_15/Species_richness within pa vs  ag_suit.tif",
	width = 12, height = 10, units = "cm", pointsize = 12, res = 300)

par(mfrow = c(1,1))


  par(mar=c(4,4.5,4,1.5))
  par(mgp=c(2.5,1,0))
  
 L = 30

  ag <-seq(from=min(model.data$ag_suit[which(model.data$Within_PA == "no")]),
		to=max(model.data$ag_suit[which(model.data$Within_PA == "no")]),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  

model.data$Within_PA <- relevel(model.data$Within_PA, "no")
model.data$Zone<- relevel(model.data$Zone , "Tropical")
model.data$taxon_of_interest <- relevel(model.data$taxon_of_interest , "Invertebrates")
mam <- glmer(Species_richness.model2$final.call, model.data, family = "poisson")

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
    mm<-model.matrix(terms(mam),newdat)
    pvar1 <- diag(mm %*% tcrossprod(vcov(mam),mm))
    z[x]<-mm %*% fixef(mam)
    zu[x]<-z[x]+sqrt(pvar1)
    zl[x]<-z[x]-sqrt(pvar1) 
  }
  
  #backtransform from poisson log-link
  z <- exp(z)
  zu <- exp(zu)
  zl <- exp(zl)


  plot(ag,z, ylim=c(0,50), col = outside.col,
		bty = "l", 
		type = "l",ylab = "Species richness per site ± s.e", 
		xlab="Agricultural suitability (higher = more suitable)")
  #points(ag,zu,type="l",lty=2, col = outside.col)
  #points(ag,zl,type="l",lty=2, col = outside.col)

  polygon(c(ag,rev(ag)),c(zu, rev(zl)),lty=0, col = outside.col.ci)





model.data$Within_PA <- relevel(model.data$Within_PA, "yes")
model.data$Zone<- relevel(model.data$Zone , "Tropical")
model.data$taxon_of_interest <- relevel(model.data$taxon_of_interest, "Invertebrates")
mam<- glmer(Species_richness.model2$final.call, model.data,  family = "poisson", control= glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

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
    mm<-model.matrix(terms(mam),newdat)
    pvar1 <- diag(mm %*% tcrossprod(vcov(mam),mm))
    y[x]<-mm %*% fixef(mam)
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


legend("topright", c("Protected", "Unprotected") , cex = 1,
	col = c(inside.col, outside.col), lty = c(1,1), lwd = c(1,1))



dev.off()


### extra stuff on ag suit values - replace with violin plot style


# add on distribution of ag_suit values for landuses
#ph <- coef(Species_richness.model$model)$Predominant_habitat[,1]
#exp(ph)

#from data 

#lu.means <- aggregate(ag_suit ~ Predominant_habitat, model.data, mean)
#lu.sd <- aggregate(ag_suit ~ Predominant_habitat, model.data, sd)
#lu.n <- aggregate(ag_suit ~ Predominant_habitat, model.data, length)

#lu.se <- data.frame(Predominant_habitat = lu.sd$Predominant_habitat, 
		ag_suit = lu.sd$ag_suit/sqrt(lu.n$ag_suit))

#for(y in 1:6){
#	ylim.min <- 5
#	xpoint <- lu.means$ag_suit[which(lu.means$Predominant_habitat == lu[y])]
#	xse <- lu.sd$ag_suit[which(lu.sd$Predominant_habitat == lu[y])]
#	points(y = ylim.min + y, x = xpoint, pch = 16, col = lu.cols[y])
#	xmax <- xpoint + xse
#	xmin <- xpoint - xse
#	arrows(xmax,ylim.min+y,xmin,ylim.min+y, 
#	cex = 0.7, code = 3, angle = 90, length = 0.04, col = lu.cols[y])
#}

#legend("topleft", rev(lu), cex = 0.7,
#	col = rev(lu.cols), pch = 16)







#### within PA and slope ##########



tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/01_15/Species_richness within pa vs slope.tif",
	width = 12, height = 10, units = "cm", pointsize = 12, res = 300)

par(mfrow = c(1,1))


  par(mar=c(4,4.5,4,1.5))
  par(mgp=c(2.5,1,0))
  
 L = 30

  slope <-seq(from=min(model.data$log_slope[which(model.data$Within_PA == "no")]),
		to=max(model.data$log_slope[which(model.data$Within_PA == "no")]),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  


# it doesnt matter what the values for taxon and zone are, as there is no interaction with these



model.data$Within_PA <- relevel(model.data$Within_PA, "no")
model.data$Zone<- relevel(model.data$Zone , "Tropical")
model.data$taxon_of_interest <- relevel(model.data$taxon_of_interest , "Invertebrates")
mam <- glmer(Species_richness.model2$final.call, model.data, family = "poisson")

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
    mm<-model.matrix(terms(mam),newdat)
    pvar1 <- diag(mm %*% tcrossprod(vcov(mam),mm))
    z[x]<-mm %*% fixef(mam)
    zu[x]<-z[x]+sqrt(pvar1)
    zl[x]<-z[x]-sqrt(pvar1) 
  }
  

 # slope <- exp(slope) -1


  #backtransform from poisson log-link
  z <- exp(z)
  zu <- exp(zu)
  zl <- exp(zl)


  plot(slope,z, ylim=ylims, xlim = slope.lim, col = outside.col,
		bty = "l", xaxt = "n", #log = "x",
		type = "l",ylab = "Species richness per site ± s.e", 
		xlab="Slope (degrees)")
  axis(1, log(c(0,2.5,5,10,20)+1), c(0,2.5,5,10,20))
  polygon(c(slope,rev(slope)),c(zu, rev(zl)),lty=0, col = outside.col.ci)
 # points(slope,zu,type="l",lty=2, col = outside.col)
 # points(slope,zl,type="l",lty=2, col = outside.col)
  rug(data$log_slope[which(data$Within_PA == "no")],
	ticksize = 0.03, side = 1,   lwd = 1, col = 8, pos = min(ylims)) 





model.data$Within_PA <- relevel(model.data$Within_PA, "yes")

mam<- glmer(Species_richness.model2$final.call, model.data,  family = "poisson", control= glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

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
    mm<-model.matrix(terms(mam),newdat)
    pvar1 <- diag(mm %*% tcrossprod(vcov(mam),mm))

    y[x]<-mm %*% fixef(mam)
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

legend("topright", c("Protected", "Unprotected") , col = c(inside.col, outside.col), lty = c(1,1), lwd = c(1,1))


dev.off()







#### PA size #####

tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/01_15/Species richness vs PA size.tif",
	width = 15, height = 12, units = "cm", pointsize = 12, res = 300)



  par(mar=c(4,4.5,4,1.5))
  par(mgp=c(2.5,1,0))
  
 L = 30

  gis <-seq(from=min(model.data$log_AREA.PA),
		to=max(model.data$log_AREA.PA),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  

model.data$Zone<- relevel(model.data$Zone , "Tropical")
model.data$taxon_of_interest <- relevel(model.data$taxon_of_interest , "Invertebrates")
mam <- glmer(Species_richness.model2$final.call, model.data, family = "poisson")

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
    mm<-model.matrix(terms(mam),newdat)
    pvar1 <- diag(mm %*% tcrossprod(vcov(mam),mm))
    z[x]<-mm %*% fixef(mam)
    zu[x]<-z[x]+sqrt(pvar1)
    zl[x]<-z[x]-sqrt(pvar1) 
  }
  
  z <- exp(z)
  zu <- exp(zu)
  zl <- exp(zl)

# gis_plot <- exp(gis) -1
  gis_plot <- gis

  plot(gis_plot,z, ylim=ylims, xlim = size.lim,  col = 1, lwd = 1,
		bty = "l",  xaxt = "n", #log = "x",
		type = "l",ylab = "Species richness per site ± s.e", xlab="PA size (km2)")
  rug(data$log_AREA.PA, ticksize = 0.03, side = 1, lwd = 1, col = 1)
  axis(1,at = log(c(0,10,100,1000,10000)+1), c(0,10,100,1000,10000))
  polygon(c(gis_plot,rev(gis_plot)),c(zu, rev(zl)),lty=0, col = outside.col.ci)
#  points(gis_plot,zu,type="l",lty=2,  lwd = 1, col = 1)
#  points(gis_plot,zl,type="l",lty=2,  lwd = 1, col =  1)



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







