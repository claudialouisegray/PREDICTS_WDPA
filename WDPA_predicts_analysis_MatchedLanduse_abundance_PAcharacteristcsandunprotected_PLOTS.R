library(lme4)
library(RColorBrewer)

setwd("R:\\ecocon_d\\clg32\\GitHub\\PREDICTS_WDPA")
source("addHistogram.R")
source("addDataDistribution.R")


### dist to boundary PLOTS

load("N:\\Documents\\PREDICTS\\WDPA analysis\\RData files\\abundance_with_block_and_keeping_confounding_vars.RData")

model.data <- log_abundance.model$data

data <- matched.landuse[,c("log_abundance", "Zone", "taxon_of_interest", "ag_suit", "log_elevation", "elevation", 
	"log_AREA.PA", "AREA.PA", "DoP.PA", "IUCN.PA",
	 "log_slope", "slope", "bound_dist_km_PA_neg",  "Predominant_habitat", "SS", "SSBS")]

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
lu.cols2 = c("#66A61E", "#8ecfbc", "#7570B3","#E6AB02","#D95F02", "#E7298A")
lu.cols2.ci <- c("#66A61E90","#8ecfbc90","#7570B390","#E6AB0290","#D95F0290","#E7298A90")

ylims <- c(0,400) #c(0,500)
slope.lim <- c(-0.3, log(60))
elev.lim <- c(-0.3,log(300000))
size.lim <- c(-0.3,log(500000))
age.lim <- c(-5,100)
ag.lim <- c(0.3,10.5)







###### dist to boundary vs taxon  #############



tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/01_15/abundance vs dist to boundary vs taxon.tif",
	width = 25, height = 10, units = "cm", pointsize = 12, res = 300)


par(mar=c(4,4.5,4,1.5))
par(mgp=c(2.5,1,0))
  
par(mfrow = c(1,3))


taxa <- c("Plants", "Invertebrates", "Vertebrates")

t <- taxa[1]
i <- 0


for(t in taxa) {


i <- i + 1

 L = 30

  lbd <-seq(from=min(model.data$log_bound_dist_km_PA_neg[which(model.data$taxon_of_interest == t)]),
	to=max(model.data$log_bound_dist_km_PA_neg[which(model.data$taxon_of_interest == t)]),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  




model.data$Zone<- relevel(model.data$Zone , "Tropical")
model.data$taxon_of_interest <- relevel(model.data$taxon_of_interest , t)
mam <- lmer(log_abundance.model$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

 for (x in 1:L)
  {
    log_bound_dist_km_PA_neg <- lbd[x]
    log_abundance <-0
#    log_slope <- mean(model.data$log_slope)
    log_elevation<- mean(model.data$log_elevation)
    log_AREA.PA <- mean(model.data$log_AREA.PA)
    IUCN.PA <- "1.5"
    DoP.PA <- mean(model.data$DoP.PA)
    Zone <-"Tropical"
    taxon_of_interest<- "Invertebrates"
    Predominant_habitat <-"Primary Vegetation"
    ag_suit <- median(model.data$ag_suit)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, log_abundance, log_AREA.PA, DoP.PA, log_elevation,  ag_suit)) #log_slope,
    newdat.f<-data.frame(cbind(Zone,taxon_of_interest, Predominant_habitat, Within_PA))
    levels(newdat.f$Zone)<-levels(model.data$Zone)
    levels(newdat.f$taxon_of_interest)<-levels(model.data$taxon_of_interest)
    levels(newdat.f$Predominant_habitat)<-levels(model.data$Predominant_habitat)
    levels(newdat.f$ag_suit) <- levels(model.data$ag_suit)
   # levels(newdat.f$IUCN.PA) <- levels(model.data$IUCN.PA)
    newdat<-data.frame(newdat.f,newdat.n)
    mm<-model.matrix(terms(mam),newdat)
    pvar1 <- diag(mm %*% tcrossprod(vcov(mam),mm))
    z[x]<-mm %*% fixef(mam)
    zu[x]<-z[x]+sqrt(pvar1)
    zl[x]<-z[x]-sqrt(pvar1) 
  }
  
 


 # transforming back - take - off first and put back
 inside <- which(lbd <0)
 lbd1 <- exp(abs(lbd)) -1
 lbd1[inside] <- -1* lbd1[inside] 

 	
  z  <- exp(z) 
  zu <- exp(zu)
  zl <- exp(zl)



#  plot(lbd1,z, ylim=ylims, xlim = c(-50,200), col = taxa.cols[i], main = t,lwd = 1,
#		bty = "l", #yaxt = "n", 
#		type = "l",ylab = "Abundance per site ± s.e", xlab="Distance to PA boundary (km)")
#  rug(data$bound_dist_km_PA_neg[which(data$taxon_of_interest == t)]
#	, ticksize = 0.03, side = 1, lwd = 0.5, col = taxa.cols[i]) 
#  axis(2,at = c(100,200,500,1000), c(100,200,500,1000))
#  points(lbd1,zu,type="l",lty=2, col = taxa.cols[i])
#  points(lbd1,zl,type="l",lty=2, col = taxa.cols[i])
#  polygon(c(lbd1,rev(lbd1)),c(zu, rev(zl)),lty=0, col = taxa.cols.ci[i])



  plot(lbd,z, ylim=ylims, xlim = c(-1*log(50+1), log(500+1)), col = taxa.cols[i], main = t,
		bty = "l", xaxt = "n", #log = "x",
		type = "l",ylab = "Abundance per site ± s.e", xlab="(Log) Distance to PA boundary (km)")
 axis(1,at = log(c(0,1,10,50,200)+1), c(0,1,10,50,200)) 
 axis(1,at = -1*(log(c(1,10,50,200)+1)), c(-1,-10,-50,-200)) 
 rug(model.data$log_bound_dist_km_PA_neg[which(data$taxon_of_interest == t)], ticksize = 0.03, side = 1, lwd = 0.5, col = taxa.cols[i]) 
# points(lbd,zu,type="l",lty=2, col =8)
# points(lbd,zl,type="l",lty=2, col = 8)
 polygon(c(lbd,rev(lbd)),c(zu, rev(zl)),lty=0, col = taxa.cols.ci[i])

 abline(v = 0, lty = 2, col = 8)

}

max(data$bound_dist_km_PA_neg[which(data$taxon_of_interest == t)])


dev.off()









### Within_PA PLOTS


model.data <- log_abundance.model2$data



data <- matched.landuse[,c("log_abundance", "Zone", "taxon_of_interest", "ag_suit", "log_elevation", "elevation", 
	"log_AREA.PA", "AREA.PA", "DoP.PA", "IUCN.PA",
	 "log_slope", "slope", "bound_dist_km_PA_neg", "Within_PA", "Predominant_habitat", "SS", "SSBS")]

data <- na.omit(data) #have to get rid of NA to try poly on ag_suit

nrow(data)
length(unique(data$SS))

nrow(model.data)
length(unique(model.data$SS))



model.data$taxon_of_interest <- relevel(model.data$taxon_of_interest, "Invertebrates")




#### PA Size in different zones #####


#tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/01_15/abundance vs PA.AREA vs zone.tif",
#	width = 12, height = 10, units = "cm", pointsize = 12, res = 300)

tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/01_15/abundance vs PA.AREA vs zone EXTRA.tif",
	width = 20, height = 20, units = "cm", pointsize = 12, res = 300)


  par(mfrow = c(3,1))
  par(mar=c(4,4.5,1,1.5))
  par(mgp=c(2.5,1,0))
  
  L = 30

  gis <-seq(from=min(model.data$log_AREA.PA[which(model.data$Zone == "Tropical")]),
		to=max(model.data$log_AREA.PA[which(model.data$Zone == "Tropical")]),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  



model.data$Zone<- relevel(model.data$Zone , "Tropical")
mam <- lmer(log_abundance.model2$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

 for (x in 1:L)
  {
    log_abundance <-0
    log_slope <- mean(model.data$log_slope)
    log_elevation <- mean(model.data$log_elevation)
    log_AREA.PA <- gis[x]
    DoP.PA <- mean(model.data$DoP.PA)
    IUCN.PA <- "1.5"
    Zone <-"Tropical"
    taxon_of_interest<- "Vertebrates"
    Predominant_habitat <-"Primary Vegetation"
    ag_suit <- median(model.data$ag_suit)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind( log_abundance, log_slope,log_AREA.PA, DoP.PA, log_elevation, ag_suit))
    newdat.f<-data.frame(cbind( Zone,taxon_of_interest, Predominant_habitat,  Within_PA))
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
  
 #gis <- exp(gis) -1

 z <- exp(z)
 zu <- exp(zu)
 zl <- exp(zl)

  plot(gis,z, ylim=ylims, xlim = size.lim, col = zone.cols[1], lwd = 1,
		bty = "l", axes = F, #log = "x", yaxt = "n", 
		type = "l",ylab = "Abundance per site ± s.e", xlab="PA size (km2)")
  axis(1, log(c(0,5,10,100,1000,10000)+1), c(0,5,10,100,1000,10000))
  axis(2, c(0,100,200,300,400),c(0,100,200,300,400))
  rug(data$log_AREA.PA[which(data$Zone == "Tropical")]
	, ticksize = 0.03, side = 1, lwd = 1, col = zone.cols[1])
#  points(gis,zu,type="l",lty=2,  lwd = 1, col =  1)
#  points(gis,zl,type="l",lty=2,  lwd = 1, col =  1)
  polygon(c(gis,rev(gis)),c(zu, rev(zl)),lty=0, col = zone.cols.ci[1])



  gis <-seq(from=min(model.data$log_AREA.PA[which(model.data$Zone == "Temperate")]),
		to=max(model.data$log_AREA.PA[which(model.data$Zone == "Temperate")]),length=L)
  y<-vector(mode="numeric",length=L)
  yu<-vector(mode="numeric",length=L)
  yl<-vector(mode="numeric",length=L)
  



model.data$Zone<- relevel(model.data$Zone , "Temperate")
mam <- lmer(log_abundance.model2$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

 for (x in 1:L)
  {
    log_abundance <-0
    log_slope <- mean(model.data$log_slope)
    log_elevation <- mean(model.data$log_elevation)
    log_AREA.PA <- gis[x]
    DoP.PA <- mean(model.data$DoP.PA)
    IUCN.PA <- "1.5"
    Zone <-"Temperate"
    taxon_of_interest<- "Vertebrates"
    Predominant_habitat <-"Primary Vegetation"
    ag_suit <- median(model.data$ag_suit)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind( log_abundance, log_slope,log_AREA.PA, DoP.PA, log_elevation, ag_suit))
    newdat.f<-data.frame(cbind( Zone,taxon_of_interest, Predominant_habitat,  Within_PA))
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
  
 #gis <- exp(gis) -1

 y <- exp(y)
 yu <- exp(yu)
 yl <- exp(yl)

  points(gis,y, type = "l", col = zone.cols[2])
  rug(data$log_AREA.PA[which(data$Zone == "Temperate")], pos= min(ylims), 
	, ticksize = 0.03, side = 1, lwd = 1, col = zone.cols[2])
#  points(gis,zu,type="l",lty=2,  lwd = 1, col =  1)
#  points(gis,zl,type="l",lty=2,  lwd = 1, col =  1)
  polygon(c(gis,rev(gis)),c(yu, rev(yl)),lty=0, col = zone.cols.ci[2])

legend(x = log(10000 +1) + 0.1, y = mean(ylims),
	c("Tropical", "Temperate") , col = zone.cols, bty = "n", lty = 1, lwd = c(2,1))




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







#### PA age, for different zones #####


#tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/01_15/abundance vs zone vs DoP.tif",
#	width = 12, height = 10, units = "cm", pointsize = 12, res = 300)

tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/01_15/abundance vs zone vs DoP EXTRA.tif",
	width = 20, height = 20, units = "cm", pointsize = 12, res = 300)


par(mfrow = c(3,1))


  par(mar=c(4,4.5,1,1.5))
  par(mgp=c(2.5,1,0))
  
 L = 30

  dop <-seq(from=min(model.data$DoP.PA[which(model.data$Zone == "Tropical")]),
		to=max(model.data$DoP.PA[which(model.data$Zone == "Tropical")]),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  



model.data$Zone<- relevel(model.data$Zone , "Tropical")
mam <- lmer(log_abundance.model2$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

 for (x in 1:L)
  {
    log_abundance <-0
    log_slope <- mean(model.data$log_slope)
    log_AREA.PA <- mean(model.data$log_AREA.PA)
    DoP.PA <- dop[x]
    IUCN.PA <- "1.5"
    Zone <-"Tropical"
    taxon_of_interest<- "Vertebrates"
    Predominant_habitat <-"Primary Vegetation"
    ag_suit <- median(model.data$ag_suit)
    log_elevation <- mean(model.data$log_elevation)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_abundance, log_AREA.PA, DoP.PA, log_elevation, log_slope, ag_suit))
    newdat.f<-data.frame(cbind( Zone,taxon_of_interest, Predominant_habitat,  Within_PA))
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

  plot(dop,z, ylim=ylims, xlim = age.lim, col = zone.cols[1], lwd = 1,
		bty = "l", axes = F, #log = "x",  
		type = "l",ylab = "Abundance per site ± s.e", xlab="Duration of protection (yr)")
  rug(data$DoP.PA[which(data$Zone == "Tropical")]
	, ticksize = 0.03, side = 1, lwd = 1, col = zone.cols[1])
  axis(1, c(0,20,40,60,80),c(0,20,40,60,80))
  axis(2, c(0,100,200,300,400),c(0,100,200,300,400))

#  points(dop,zu,type="l",lty=2,  lwd = 2, col =  zone.cols[1])
#  points(dop,zl,type="l",lty=2,  lwd = 2, col =  zone.cols[1])
  polygon(c(dop,rev(dop)),c(zu, rev(zl)),lty=0, col = zone.cols.ci[1])

  dop <-seq(from=min(model.data$DoP.PA[which(model.data$Zone == "Temperate")]),
		to=max(model.data$DoP.PA[which(model.data$Zone == "Temperate")]),length=L)
  y<-vector(mode="numeric",length=L)
  yu<-vector(mode="numeric",length=L)
  yl<-vector(mode="numeric",length=L)
  


model.data$Zone<- relevel(model.data$Zone , "Temperate")
mam <- lmer(log_abundance.model2$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

 for (x in 1:L)
  {
    log_abundance <-0
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
    newdat.n<-data.frame(cbind(log_abundance, log_slope,log_AREA.PA, DoP.PA, log_elevation, ag_suit))
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
  
 y <- exp(y)
 yu <- exp(yu)
 yl <- exp(yl)

  points(dop,y,type="l",lty=1, col =  zone.cols[2])
#  points(dop,yu,type="l",lty=2, col =  zone.cols[2])
#  points(dop,yl,type="l",lty=2, col =  zone.cols[2])
  rug(data$DoP.PA[which(data$Zone == "Temperate")],
		col = zone.cols[2], lwd = 1, pos = min(ylims))
  polygon(c(dop,rev(dop)),c(yu, rev(yl)),lty=0, col = zone.cols.ci[2])


legend("right", c("Tropical", "Temperate") , col = zone.cols, lty = 1, lwd = c(2,1), bty = "n")



addDataDistribution(data = model.data, 
			b = 50,
			x = "DoP.PA",
			xlim = age.lim,
			var = "Predominant_habitat",
			include.lowest = F,
			levels = lu,
			legend.spacing = 3,
			levels.col = lu.cols2.ci,
			axis.text.pos = c(0,20,40,60,80),
			axis.text = c(0,20,40,60,80) )


addDataDistribution(data = model.data, 
			b = 50,
			x = "DoP.PA",
			xlim = age.lim,
			var = "taxon_of_interest",
			include.lowest = F,
			levels = taxa,
			legend.spacing = 3,
			levels.col = taxa.cols,
			axis.text.pos = c(0,20,40,60,80),
			axis.text = c(0,20,40,60,80))


dev.off()






#### within PA and slope ##########


tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/01_15/abundance vs within pa vs slope EXTRA.tif",
	width = 20, height = 20, units = "cm", pointsize = 12, res = 300)

#tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/01_15/abundance vs within pa vs slope.tif",
#	width = 12, height = 10, units = "cm", pointsize = 12, res = 300)

par(mfrow = c(3,1))

  par(mar=c(4,4.5,1,1.5))
  par(mgp=c(2.5,1,0))
  

  slope <-seq(from=min(model.data$log_slope[which(data$Within_PA == "yes")])
		,to=max(model.data$log_slope[which(data$Within_PA == "yes")]),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  


model.data$Zone<- relevel(model.data$Zone , "Temperate")
model.data$Within_PA <- relevel(model.data$Within_PA, "yes")
mam <- lmer(log_abundance.model2$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

 for (x in 1:L)
  {
    ag_suit <- median(model.data$ag_suit)
    log_abundance <-0
    log_slope <- slope[x]
    log_elevation <- mean(model.data$log_elevation)
    log_AREA.PA <- mean(model.data$log_AREA.PA)
    DoP.PA <- mean(model.data$DoP.PA)
    Zone <- "Tropical"
    taxon_of_interest<- "Invertebrates"
    Predominant_habitat <-"Primary Vegetation"
    Within_PA <- "yes"
    newdat.n<-data.frame(cbind(log_abundance, log_AREA.PA, DoP.PA, log_elevation, log_slope, ag_suit))
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
  
 # original backtransform, for occurence(binomial) data
 # z<-1/(1+exp(-(z)))
 # zu<-1/(1+exp(-(zu)))
 # zl<-1/(1+exp(-(zl)))


 # transforming back 

 # slope_plot <- exp(slope) -1
 slope_plot <- slope

  z  <- exp(z)
  zu <- exp(zu)
  zl <- exp(zl)


  plot(slope_plot,z, ylim=ylims, xlim = slope.lim,
		 col =inside.col, 
		bty = "l",  lwd =1, axes = F, #log = "x",
		type = "l",ylab = "Abundance per site ± s.e", xlab="Slope (degrees)")
#  points(slope_plot,zu,type="l",lty=2,lwd = 1, col = inside.col)
#  points(slope_plot,zl,type="l",lty=2, lwd = 1,col = inside.col)
  polygon(c(slope_plot,rev(slope_plot)),c(zu, rev(zl)),lty=0, col = inside.col.ci)
  axis(1, log(c(0,2.5,5, 10, 20)+1), c(0,2.5,5,10,20))
  axis(2, c(0,100,200,300,400),c(0,100,200,300,400))
  rug(data$log_slope[which(data$Within_PA == "yes")],
		col = inside.col, lwd = 1)




  slope <-seq(from=min(model.data$log_slope[which(data$Within_PA == "no")])
		,to=max(model.data$log_slope[which(data$Within_PA == "no")]),length=L)
  y<-vector(mode="numeric",length=L)
  yu<-vector(mode="numeric",length=L)
  yl<-vector(mode="numeric",length=L)
  





model.data$Within_PA <- relevel(model.data$Within_PA, "no")
mam <- lmer(log_abundance.model2$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

 for (x in 1:L)
  {
    ag_suit <- median(model.data$ag_suit)
    log_abundance <-0
    log_slope <- slope[x]
    log_elevation <- mean(model.data$log_elevation)
    log_AREA.PA <- mean(model.data$log_AREA.PA)
    DoP.PA <- mean(model.data$DoP.PA)
    Zone <- "Tropical"
    taxon_of_interest<- "Invertebrates"
    Predominant_habitat <-"Primary Vegetation"
    Within_PA <- "yes"
    newdat.n<-data.frame(cbind(log_abundance, log_AREA.PA, DoP.PA, log_elevation, log_slope, ag_suit))
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
  
 y <- exp(y)
 yu <- exp(yu)
 yl <- exp(yl)

 # slope_plot <- exp(slope) -1
 slope_plot <- slope

  points(slope_plot,y,type="l",lty=1, col =  outside.col)
#  points(slope_plot,yu,type="l",lty=2, col =  outside.col)
#  points(slope_plot,yl,type="l",lty=2, col =  outside.col)
   polygon(c(slope_plot,rev(slope_plot)),c(yu, rev(yl)),lty=0, col = outside.col.ci)
   rug(data$log_slope[which(data$Within_PA == "no")],
		col = 8, lwd = 1, pos = min(ylims))


legend("right", c("Protected", "Unprotected") , col = c(inside.col,outside.col), bty = "n", lty = 1, lwd = c(2))


# add half violin plots




par(mar=c(2,4.5,1,3))


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









#### within PA and ag suit ##########



tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/01_15/abundance within pa vs ag_suit EXTRA.tif",
	width = 20, height = 20, units = "cm", pointsize = 12, res = 300)


par(mfrow = c(3,1))


  par(mar=c(2,4.5,1,3))
  par(mgp=c(2.5,1,0))
  
 L = 30

  ag <-seq(from=min(model.data$ag_suit),to=max(model.data$ag_suit),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  


model.data$Within_PA <- relevel(model.data$Within_PA, "no")
model.data$Zone<- relevel(model.data$Zone , "Tropical")
model.data$taxon_of_interest <- relevel(model.data$taxon_of_interest , "Invertebrates")
mam <- lmer(log_abundance.model2$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

 for (x in 1:L)
  {
    ag_suit <- ag[x]
    log_abundance <-0
    log_slope <- mean(model.data$log_slope)
    log_elevation<- mean(model.data$log_elevation)
    log_AREA.PA <- mean(model.data$log_AREA.PA)
    DoP.PA <- mean(model.data$DoP.PA)
    Zone <-"Tropical"
    taxon_of_interest<- "Invertebrates"
    Predominant_habitat <-"Primary Vegetation"
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_abundance, log_AREA.PA, DoP.PA, log_elevation, log_slope, ag_suit)) 
    newdat.f<-data.frame(cbind(Zone,taxon_of_interest, Predominant_habitat, Within_PA))
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
  
 # original backtransform, for occurence(binomial) data
 # z<-1/(1+exp(-(z)))
 # zu<-1/(1+exp(-(zu)))
 # zl<-1/(1+exp(-(zl)))


 # transforming back 

  z  <- exp(z)
  zu <- exp(zu)
  zl <- exp(zl)

  plot(ag,z, ylim=ylims, xlim = ag.lim, col = outside.col,
		bty = "l", axes = F,
		type = "l",ylab = "Abundance per site ± s.e", xlab="Agricultural suitability (higher = more suitable)")
  axis(1, seq(1,8,1), seq(1,8,1))
  axis(2, seq(min(ylims),max(ylims),50), seq(min(ylims),max(ylims),50))
#  points(ag,zu,type="l",lty=2, col = outside.col)
#  points(ag,zl,type="l",lty=2, col = outside.col)
   polygon(c(ag,rev(ag)),c(zu, rev(zl)),lty=0, col = outside.col.ci)






model.data$Within_PA <- relevel(model.data$Within_PA, "yes")
mam<- lmer(log_abundance.model2$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))


  y<-vector(mode="numeric",length=L)
  yu<-vector(mode="numeric",length=L)
  yl<-vector(mode="numeric",length=L)


 for (x in 1:L)
  {
    ag_suit <- ag[x]
    log_abundance <-0
    log_slope <- mean(model.data$log_slope)
    log_elevation<- mean(model.data$log_elevation)
    log_AREA.PA <- mean(model.data$log_AREA.PA)
    IUCN.PA <- "1.5"
    DoP.PA <- mean(model.data$DoP.PA)
    Zone <-"Tropical"
    taxon_of_interest<- "Invertebrates"
    Predominant_habitat <-"Primary Vegetation"
    Within_PA <- "yes"
    newdat.n<-data.frame(cbind(log_abundance, log_AREA.PA, DoP.PA, log_elevation, log_slope, ag_suit))
    newdat.f<-data.frame(cbind(Zone,taxon_of_interest, Predominant_habitat, Within_PA))
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


 # backtransform abundance values
   y  <- exp(y)
  yu <- exp(yu)
  yl <- exp(yl)


  
  points(ag,y,type="l",lty=1 , col = inside.col, lwd = 1)
#  points(ag,yu,type="l",lty=3, lwd = 2, col = inside.col)
#  points(ag,yl,type="l",lty=3, lwd = 2, col = inside.col)
   polygon(c(ag,rev(ag)),c(yu, rev(yl)),lty=0, col = inside.col.ci)

legend(x = 8.5, y = 175 , c("Protected", "Unprotected") , col = c(inside.col, outside.col), lty = c(1,1), lwd = c(1,1), bty = "n")


# add half violin plots to show data spread of different landuses and taxa 

ag2 <- rep(seq(0.5, 8.5, 1),each = 2)
	end <- length(ag2)-1
ag2 <- ag2[2:end]

par(mar=c(2,4.5,1,3))

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









#### elevation ##########


#tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/01_15/abundance vs elevation.tif",
#	width = 12, height = 10, units = "cm", pointsize = 12, res = 300)

tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/01_15/abundance vs elevation EXTRA.tif",
	width = 20, height = 20, units = "cm", pointsize = 12, res = 300)

  par(mfrow = c(3,1))
  par(mar=c(4,4.5,4,1.5))
  par(mgp=c(2.5,1,0))

  elevation <-seq(from=min(model.data$log_elevation),to=max(model.data$log_elevation),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  

model.data$Within_PA <- relevel(model.data$Within_PA, "no")
model.data$Zone<- relevel(model.data$Zone , "Tropical")
model.data$taxon_of_interest <- relevel(model.data$taxon_of_interest , "Invertebrates")
mam <- lmer(log_abundance.model$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

 for (x in 1:L)
  {
    ag_suit <- median(model.data$ag_suit)
    log_abundance <-0
    log_slope <- mean(model.data$log_slope)
    log_elevation<- elevation[x]
    log_AREA.PA <- mean(model.data$log_AREA.PA)
    IUCN.PA <- "1.5"
    DoP.PA <- mean(model.data$DoP.PA)
    Zone <-"Tropical"
    taxon_of_interest<- "Invertebrates"
    Predominant_habitat <-"Primary Vegetation"
    Within_PA <- "yes"
    newdat.n<-data.frame(cbind(log_abundance, log_AREA.PA, DoP.PA, log_elevation, log_slope, ag_suit))
    newdat.f<-data.frame(cbind(Zone,taxon_of_interest, Predominant_habitat, Within_PA))
    levels(newdat.f$Zone)<-levels(model.data$Zone)
    levels(newdat.f$taxon_of_interest)<-levels(model.data$taxon_of_interest)
    levels(newdat.f$Predominant_habitat)<-levels(model.data$Predominant_habitat)
    levels(newdat.f$ag_suit) <- levels(model.data$ag_suit)
    levels(newdat.f$Within_PA) <- levels(model.data$Within_PA)
    newdat<-data.frame(newdat.f,newdat.n)
    mm<-model.matrix(terms(mam),newdat)
    pvar1 <- diag(mm %*% tcrossprod(vcov(mam),mm))
    pvar1 <- diag(mm %*% tcrossprod(vcov(mam),mm))
    z[x]<-mm %*% fixef(mam)
    zu[x]<-z[x]+sqrt(pvar1)
    zl[x]<-z[x]-sqrt(pvar1) 
  }
  
 # original backtransform, for occurence(binomial) data
 # z<-1/(1+exp(-(z)))
 # zu<-1/(1+exp(-(zu)))
 # zl<-1/(1+exp(-(zl)))


 # transforming back 

#  elevation <- exp(elevation) -1

  z  <- exp(z)
  zu <- exp(zu)
  zl <- exp(zl)

  plot(elevation,z, ylim=ylims, col = 1,
		bty = "l",  xaxt = "n", #log = "x",
		type = "l",ylab = "Abundance per site ± s.e", xlab="Elevation (m)")
#  points(elevation,zu,type="l",lty=2, col = 1)
#  points(elevation,zl,type="l",lty=2, col = 1)
   polygon(c(elevation,rev(elevation)),c(zu, rev(zl)),lty=0, col = outside.col.ci)
  rug(model.data$log_elevation, ticksize = 0.03, side = 1,  lwd = 0.5, col = 1) 
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























#### older code retained incase it is needed







#### within PA and slope ##########
### and diff zones ####



tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/abundance vs within pa vs slope vs zone.tif",
	width = 25, height = 12, units = "cm", pointsize = 12, res = 300)


zones <- c("Tropical", "Temperate")
t <- zones[1]

par(mfrow = c(1,2))

for(t in zones) { 


  par(mar=c(4,4.5,4,1.5))
  par(mgp=c(2.5,1,0))
  
 L = 100

  slope <-seq(from=min(model.data$log_slope[which(data$Within_PA == "yes" & data$Zone == t)])
		,to=max(model.data$log_slope[which(data$Within_PA == "yes" & data$Zone == t)]),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  



model.data$Within_PA <- relevel(model.data$Within_PA, "yes")
model.data$Zone<- relevel(model.data$Zone , t)
model.data$taxon_of_interest <- relevel(model.data$taxon_of_interest , "Invertebrates")
mam <- lmer(log_abundance.model$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

 for (x in 1:L)
  {
    ag_suit <- median(model.data$ag_suit)
    log_abundance <-0
    log_slope <- slope[x]
    log_elevation <- mean(model.data$log_elevation)
    log_AREA.PA <- mean(model.data$log_AREA.PA)
    DoP.PA <- mean(model.data$DoP.PA)
    Zone <-t
    taxon_of_interest<- "Invertebrates"
    Predominant_habitat <-"Primary Vegetation"
    log_bound_dist_km_PA_neg <- mean(model.data$log_bound_dist_km_PA_neg)
    Within_PA <- "yes"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, log_abundance, log_AREA.PA, DoP.PA, log_elevation, log_slope, ag_suit))
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
  
 # original backtransform, for occurence(binomial) data
 # z<-1/(1+exp(-(z)))
 # zu<-1/(1+exp(-(zu)))
 # zl<-1/(1+exp(-(zl)))


 # transforming back 

  slope_plot <- exp(slope) -1

  z  <- exp(z)
  zu <- exp(zu)
  zl <- exp(zl)


  plot(slope_plot,z, ylim=c(0, 500), xlim = c(0.05,10000),
		 col =1, main = zones[which(zones == t)],
		bty = "l", log = "x", lwd =2, xaxt = "n",
		type = "l",ylab = "Abundance per site ± s.e", xlab="Human population density (per km2)")
  points(slope_plot,zu,type="l",lty=2,lwd = 2, col = 1)
  points(slope_plot,zl,type="l",lty=2, lwd = 2,col = 1)
  axis(1, c(0.1,1,10,100,1000,10000), c(0.1,1,10,100,1000,10000) )
  rug(data$slope[which(data$Within_PA == "yes" & data$Zone == t)],
		col = 1, lwd = 1)



  slope <-seq(from=min(model.data$log_slope[which(data$Within_PA == "no" & data$Zone == t)])
		,to=max(model.data$log_slope[which(data$Within_PA == "no" & data$Zone == t)]),length=L)
  y<-vector(mode="numeric",length=L)
  yu<-vector(mode="numeric",length=L)
  yl<-vector(mode="numeric",length=L)
  





model.data$Within_PA <- relevel(model.data$Within_PA, "no")
model.data$Zone<- relevel(model.data$Zone , t)
mam <- lmer(log_abundance.model$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

 for (x in 1:L)
  {
    ag_suit <- median(model.data$ag_suit)
    log_abundance <-0
    log_slope <- slope[x]
    log_elevation <- mean(model.data$log_elevation)
    log_AREA.PA <- mean(model.data$log_AREA.PA)
    DoP.PA <- mean(model.data$DoP.PA)
    Zone <-"Temperate"
    taxon_of_interest<- "Invertebrates"
    Predominant_habitat <-"Primary Vegetation"
    log_bound_dist_km_PA_neg <- mean(model.data$log_bound_dist_km_PA_neg)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, log_abundance, log_AREA.PA, DoP.PA, log_elevation, log_slope, ag_suit))
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
  
 y <- exp(y)
 yu <- exp(yu)
 yl <- exp(yl)

  slope_plot <- exp(slope) -1

  points(slope_plot,y,type="l",lty=1, col =  8)
  points(slope_plot,yu,type="l",lty=2, col =  8)
  points(slope_plot,yl,type="l",lty=2, col =  8)
  rug(data$slope[which(data$Within_PA == "no" & data$Zone == t)],
		col = 8, lwd = 1, pos = 0)

}



legend("topright", c("Protected", "Unprotected") , col = c(1,8), lty = 1, lwd = c(2))


dev.off()





#### within PA and slope in diff IUCN cats ##########



tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/abundance vs IUCN cat vs slope.tif",
	width = 20, height = 20, units = "cm", pointsize = 12, res = 300)

par(mfrow = c(2,2))


titles <- c("I and II", "III to VI", "NA or Not Reported", "Unprotected")


  par(mar=c(4,4.5,4,1.5))
  par(mgp=c(2.5,1,0))
  
 L = 30

  slope <-seq(from=0,to=max(model.data$log_slope),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  


cats <- levels(matched.landuse$IUCN.PA)
par(mfrow = c(2,2))

for(c in cats){ 


#model.data$Within_PA <- relevel(model.data$Within_PA, "no")
model.data$Zone<- relevel(model.data$Zone , "Tropical")
model.data$IUCN.PA <- relevel(model.data$IUCN.PA , c)
model.data$taxon_of_interest <- relevel(model.data$taxon_of_interest , "Invertebrates")
mam <- lmer(log_abundance.model$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

 for (x in 1:L)
  {
    ag_suit <- median(model.data$ag_suit)
    log_abundance <-0
    log_slope <- slope[x]
    log_elevation <- mean(model.data$log_elevation)
    log_AREA.PA <- mean(model.data$log_AREA.PA)
    IUCN.PA <- c
    DoP.PA <- mean(model.data$DoP.PA)
    Zone <-"Temperate"
    taxon_of_interest<- "Invertebrates"
    Predominant_habitat <-"Primary Vegetation"
    log_bound_dist_km_PA_neg <- mean(model.data$log_bound_dist_km_PA_neg)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, log_abundance, log_AREA.PA, DoP.PA, log_elevation, log_slope, ag_suit))
    newdat.f<-data.frame(cbind( Zone,taxon_of_interest, Predominant_habitat, IUCN.PA))
    levels(newdat.f$Zone)<-levels(model.data$Zone)
    levels(newdat.f$taxon_of_interest)<-levels(model.data$taxon_of_interest)
    levels(newdat.f$Predominant_habitat)<-levels(model.data$Predominant_habitat)
    levels(newdat.f$ag_suit) <- levels(model.data$ag_suit)
 #   levels(newdat.f$Within_PA) <- levels(model.data$Within_PA)
    levels(newdat.f$IUCN.PA) <- levels(model.data$IUCN.PA)
    newdat<-data.frame(newdat.f,newdat.n)
    mm<-model.matrix(terms(mam),newdat)
    pvar1 <- diag(mm %*% tcrossprod(vcov(mam),mm))
    z[x]<-mm %*% fixef(mam)
    zu[x]<-z[x]+sqrt(pvar1)
    zl[x]<-z[x]-sqrt(pvar1) 
  }
  
 # original backtransform, for occurence(binomial) data
 # z<-1/(1+exp(-(z)))
 # zu<-1/(1+exp(-(zu)))
 # zl<-1/(1+exp(-(zl)))


 # transforming back 

  slope_plot <- exp(slope) -1

  z  <- exp(z)
  zu <- exp(zu)
  zl <- exp(zl)

  plot(slope_plot,z, ylim=c(0, 600), col = 8, main = titles[which(cats == c)],
		bty = "l", log = "x", lwd =2,
		type = "l",ylab = "Abundance per site ± s.e", xlab="Human population density (per km2)")
  points(slope_plot,zu,type="l",lty=2,lwd = 2, col = 8)
  points(slope_plot,zl,type="l",lty=2, lwd = 2,col = 8)
  rug(data$slope[which(data$IUCN.PA == c)], pos = 0,
	 ticksize = 0.03, side = 1,  lwd = 0.5, col = 8)


}


dev.off()






#### within PA and slope in diff IUCN cats ##########
### and diff zones ####



tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/abundance vs IUCN cat vs slope vs zone.tif",
	width = 20, height = 20, units = "cm", pointsize = 12, res = 300)

par(mfrow = c(2,2))


titles <- c("I and II", "III to VI", "NA or Not Reported", "Unprotected")


  par(mar=c(4,4.5,4,1.5))
  par(mgp=c(2.5,1,0))
  
 L = 30

  slope <-seq(from=0,to=max(model.data$log_slope),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  


cats <- levels(matched.landuse$IUCN.PA)
par(mfrow = c(2,2))

for(c in cats){ 


#model.data$Within_PA <- relevel(model.data$Within_PA, "no")
model.data$Zone<- relevel(model.data$Zone , "Tropical")
model.data$IUCN.PA <- relevel(model.data$IUCN.PA , c)
model.data$taxon_of_interest <- relevel(model.data$taxon_of_interest , "Invertebrates")
mam <- lmer(log_abundance.model$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

 for (x in 1:L)
  {
    ag_suit <- median(model.data$ag_suit)
    log_abundance <-0
    log_slope <- slope[x]
    log_elevation <- mean(model.data$log_elevation)
    log_AREA.PA <- mean(model.data$log_AREA.PA)
    IUCN.PA <- c
    DoP.PA <- mean(model.data$DoP.PA)
    Zone <-"Temperate"
    taxon_of_interest<- "Invertebrates"
    Predominant_habitat <-"Primary Vegetation"
    log_bound_dist_km_PA_neg <- mean(model.data$log_bound_dist_km_PA_neg)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, log_abundance, log_AREA.PA, DoP.PA, log_elevation, log_slope, ag_suit))
    newdat.f<-data.frame(cbind( Zone,taxon_of_interest, Predominant_habitat, IUCN.PA))
    levels(newdat.f$Zone)<-levels(model.data$Zone)
    levels(newdat.f$taxon_of_interest)<-levels(model.data$taxon_of_interest)
    levels(newdat.f$Predominant_habitat)<-levels(model.data$Predominant_habitat)
    levels(newdat.f$ag_suit) <- levels(model.data$ag_suit)
 #   levels(newdat.f$Within_PA) <- levels(model.data$Within_PA)
    levels(newdat.f$IUCN.PA) <- levels(model.data$IUCN.PA)
    newdat<-data.frame(newdat.f,newdat.n)
    mm<-model.matrix(terms(mam),newdat)
    pvar1 <- diag(mm %*% tcrossprod(vcov(mam),mm))
    z[x]<-mm %*% fixef(mam)
    zu[x]<-z[x]+sqrt(pvar1)
    zl[x]<-z[x]-sqrt(pvar1) 
  }
  
 # original backtransform, for occurence(binomial) data
 # z<-1/(1+exp(-(z)))
 # zu<-1/(1+exp(-(zu)))
 # zl<-1/(1+exp(-(zl)))


 # transforming back 

  slope_plot <- exp(slope) -1

  z  <- exp(z)
  zu <- exp(zu)
  zl <- exp(zl)

  plot(slope_plot,z, ylim=c(-50, 600), col = zone.cols[1], main = titles[which(cats == c)],
		bty = "l", log = "x", lwd =2,
		type = "l",ylab = "Abundance per site ± s.e", xlab="Human population density (per km2)")
  points(slope_plot,zu,type="l",lty=2,lwd = 2, col = zone.cols[1])
  points(slope_plot,zl,type="l",lty=2, lwd = 2,col = zone.cols[1])
  rug(data$slope[which(data$IUCN.PA == c & data$Zone == "Tropical")],
	 ticksize = 0.03, side = 1,  lwd = 0.5, col = zone.cols[1])




model.data$Zone<- relevel(model.data$Zone , "Temperate")
mam <- lmer(log_abundance.model$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

 for (x in 1:L)
  {
    ag_suit <- median(model.data$ag_suit)
    log_abundance <-0
    log_slope <- slope[x]
    log_elevation <- mean(model.data$log_elevation)
    log_AREA.PA <- mean(model.data$log_AREA.PA)
    IUCN.PA <- c
    DoP.PA <- mean(model.data$DoP.PA)
    Zone <-"Temperate"
    taxon_of_interest<- "Invertebrates"
    Predominant_habitat <-"Primary Vegetation"
    log_bound_dist_km_PA_neg <- mean(model.data$log_bound_dist_km_PA_neg)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, log_abundance, log_AREA.PA, DoP.PA, log_elevation, log_slope, ag_suit))
    newdat.f<-data.frame(cbind( Zone,taxon_of_interest, Predominant_habitat, IUCN.PA))
    levels(newdat.f$Zone)<-levels(model.data$Zone)
    levels(newdat.f$taxon_of_interest)<-levels(model.data$taxon_of_interest)
    levels(newdat.f$Predominant_habitat)<-levels(model.data$Predominant_habitat)
    levels(newdat.f$ag_suit) <- levels(model.data$ag_suit)
 #   levels(newdat.f$Within_PA) <- levels(model.data$Within_PA)
    levels(newdat.f$IUCN.PA) <- levels(model.data$IUCN.PA)
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

  points(slope_plot,y,type="l",lty=1, col =  zone.cols[2])
  points(slope_plot,yu,type="l",lty=2, col =  zone.cols[2])
  points(slope_plot,yl,type="l",lty=2, col =  zone.cols[2])
  rug(data$slope[which(data$IUCN.PA == c & data$Zone == "Temperate")],
		col = zone.cols[2], lwd = 1, pos = -50)


}



legend("topright", c("Tropical", "Temperate") , col = zone.cols, lty = 1, lwd = c(2,1))


dev.off()




### PA size vs taxon ####


tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/abundance vs PA size vs taxon.tif",
	width = 12, height = 25, units = "cm", pointsize = 12, res = 300)



  par(mar=c(4,4.5,4,1.5))
  par(mgp=c(2.5,1,0))
  
  par(mfrow = c(3,1))

taxa <- levels(matched.landuse$taxon_of_interest)
t <- taxa[1]
i <- 0

for(t in taxa) {


i <- i + 1

 L = 30

  gis <-seq(from=min(model.data$log_AREA.PA),to=max(model.data$log_AREA.PA),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  




model.data$Zone<- relevel(model.data$Zone , "Tropical")
model.data$taxon_of_interest <- relevel(model.data$taxon_of_interest , t)
mam <- lmer(log_abundance.model$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

 for (x in 1:L)
  {
    log_bound_dist_km_PA_neg <- mean (model.data$log_bound_dist_km_PA_neg)
    log_abundance <-0
    log_slope <- mean(model.data$log_slope)
    log_elevation<- mean(model.data$log_elevation)
    log_AREA.PA <- gis[x]
    IUCN.PA <- "1.5"
    DoP.PA <- mean(model.data$DoP.PA)
    Zone <-"Tropical"
    taxon_of_interest<- "Invertebrates"
    Predominant_habitat <-"Primary Vegetation"
    ag_suit <- median(model.data$ag_suit)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, log_abundance, log_AREA.PA, DoP.PA, log_elevation, log_slope, ag_suit))
    newdat.f<-data.frame(cbind(Zone,taxon_of_interest, Predominant_habitat, IUCN.PA))
    levels(newdat.f$Zone)<-levels(model.data$Zone)
    levels(newdat.f$taxon_of_interest)<-levels(model.data$taxon_of_interest)
    levels(newdat.f$Predominant_habitat)<-levels(model.data$Predominant_habitat)
    levels(newdat.f$ag_suit) <- levels(model.data$ag_suit)
   # levels(newdat.f$Within_PA) <- levels(model.data$Within_PA)
    levels(newdat.f$IUCN.PA) <- levels(model.data$IUCN.PA)
    newdat<-data.frame(newdat.f,newdat.n)
    mm<-model.matrix(terms(mam),newdat)
    pvar1 <- diag(mm %*% tcrossprod(vcov(mam),mm))
    z[x]<-mm %*% fixef(mam)
    zu[x]<-z[x]+sqrt(pvar1)
    zl[x]<-z[x]-sqrt(pvar1) 
  }
  
 


 # transforming back - take - off first and put back
 gis <- exp(gis) - 1

 	
  z  <- exp(z) 
  zu <- exp(zu)
  zl <- exp(zl)

  plot(gis,z, ylim=c(0, 900), col = 1, main = t,
		bty = "l", log = "x", 
		type = "l",ylab = "Abundance per site ± s.e", xlab="Protected area size (km2)")
  rug(data$AREA.PA[which(data$taxon_of_interest == t)]
	, ticksize = 0.03, side = 1, lwd = 0.5, col = 1) 
  axis(2,at = c(100,200,500,1000), c(100,200,500,1000))
  points(gis,zu,type="l",lty=2, col = 1)
  points(gis,zl,type="l",lty=2, col = 1)


}





dev.off()













