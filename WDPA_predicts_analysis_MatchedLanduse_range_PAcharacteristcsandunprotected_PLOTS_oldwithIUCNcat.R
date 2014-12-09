# PLOTS


model.data <- range.model$data

#need dataset with variables before log transform to plot rugs
data <- matched.landuse[,c("range", "Zone", "taxon_of_interest", "ag_suit", "log_access", "access", "DoP.PA", "AREA.PA", "IUCN.PA",
	 "log_hpd", "hpd", "log_bound_dist_km_PA_neg", "bound_dist_km_PA_neg", "Within_PA", "Predominant_habitat", "SS", "SSBS")]

data <- na.omit(data) 

names(matched.landuse)

#check data for rug has same number of sites/studies as model data
nrow(data)
length(unique(data$SS))

nrow(model.data)
length(unique(model.data$SS))



L = 100


### dist to boundary for different taxa 

display.brewer.all()
cols <- brewer.pal(8, "Paired")
display.brewer.pal(8, "Paired")
#drop the red to avoid red-green colorblindness issues

cols <- cols[c(2,4,8)]


tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/range vs taxon vs dist to boundary.tif",
	width = 10, height = 25, units = "cm", pointsize = 12, res = 300)

par(mfrow = c(3,1))
taxa <- levels(matched.landuse$taxon_of_interest)
i <- 0
t <- taxa[1]

for(t in taxa){

i <- i +1

  par(mar=c(4,4.5,4,1.5))
  par(mgp=c(2.5,1,0))
  
 L = 30

  lbd <-seq(from=min(model.data$log_bound_dist_km_PA_neg),
		to=max(model.data$log_bound_dist_km_PA_neg),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  



model.data$Zone<- relevel(model.data$Zone , "Tropical")
model.data$taxon_of_interest <- relevel(model.data$taxon_of_interest , t)
mam <- lmer(range.model$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

 for (x in 1:L)
  {
    log_bound_dist_km_PA_neg <- lbd[x]
    range <-0
    log_hpd <- mean(model.data$log_hpd)
    log_AREA.PA <- mean(model.data$log_AREA.PA )
    DoP.PA <- mean(model.data$DoP.PA)
    IUCN.PA <- "1.5"
    Zone <-"Tropical"
    taxon_of_interest<- t
    Predominant_habitat <-"Primary Vegetation"
    ag_suit <- median(model.data$ag_suit)
    log_access <- mean(model.data$log_access)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, range, log_hpd,log_AREA.PA, DoP.PA, log_access, ag_suit))
    newdat.f<-data.frame(cbind( Zone,taxon_of_interest, Predominant_habitat,  IUCN.PA))
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
  


 lbd1 <- exp(abs(lbd)) -1
 inside <- which(lbd < 0)
 lbd1[inside] <- lbd1[inside]*-1

  plot(lbd1,z, ylim=c(4.7,6.2), col = cols[i],
		bty = "l", #log = "x", #yaxt = "n", 
		type = "l",ylab = "range size Log10(sq km) ± s.e", xlab="Distance to PA boundary (km)", main = t)
  rug(data$bound_dist_km_PA_neg[which(data$taxon_of_interest == t)]
	, ticksize = 0.03, side = 1, lwd = 0.5, pos = 4.74, col = cols[i])
  points(lbd1,zu,type="l",lty=2, col = cols[i])
  points(lbd1,zl,type="l",lty=2, col = cols[i])

 abline(v = 0, lty = 2, col = 8)

}



dev.off()






###### dist to boundary  #############
###### for different zones and taxa ##################

display.brewer.all()
cols <- brewer.pal(8, "Paired")
display.brewer.pal(8, "Paired")
#drop the red to avoid red-green colorblindness issues

zone.cols <- cols[c(6,5)]


tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/endemicity vs zone vs taxon vs dist to boundary LOG SCALE.tif",
	width = 10, height = 25, units = "cm", pointsize = 12, res = 300)



par(mfrow = c(3,1))
taxa <- levels(matched.landuse$taxon_of_interest)
i <- 0
t <- taxa[1]

for(t in taxa){


  par(mar=c(4,4.5,4,1.5))
  par(mgp=c(2.5,1,0))
  
 L = 30

 lbd <-seq(from=min(model.data$log_bound_dist_km_PA_neg[which(data$Zone == "Tropical" & data$taxon_of_interest == t)]),
		to=max(model.data$log_bound_dist_km_PA_neg[which(data$Zone == "Tropical" & data$taxon_of_interest == t)]),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  





model.data$Zone<- relevel(model.data$Zone , "Tropical")
model.data$taxon_of_interest <- relevel(model.data$taxon_of_interest , t)
mam <- lmer(range.model$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

 for (x in 1:L)
  {
    log_bound_dist_km_PA_neg <- lbd[x]
    range <-0
    log_hpd <- mean(model.data$log_hpd)
    log_AREA.PA <- mean(model.data$log_AREA.PA )
    DoP.PA <- mean(model.data$DoP.PA )
    IUCN.PA <- "1.5"
    Zone <-"Tropical"
    taxon_of_interest<- t
    Predominant_habitat <-"Primary Vegetation"
    ag_suit <- median(model.data$ag_suit)
    log_access <- mean(model.data$log_access)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, range, log_hpd,log_AREA.PA, DoP.PA, log_access, ag_suit))
    newdat.f<-data.frame(cbind( Zone,taxon_of_interest, Predominant_habitat, Within_PA))
    levels(newdat.f$Zone)<-levels(model.data$Zone)
    levels(newdat.f$taxon_of_interest)<-levels(model.data$taxon_of_interest)
    levels(newdat.f$Predominant_habitat)<-levels(model.data$Predominant_habitat)
    levels(newdat.f$ag_suit) <- levels(model.data$ag_suit)
#    levels(newdat.f$IUCN.PA) <- levels(model.data$IUCN.PA)
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

 lbd1 <- exp(abs(lbd)) -1
 inside <- which(lbd < 0)
 lbd1[inside] <- lbd1[inside]*-1

  plot(lbd1,z, ylim=c(0.1,0.21), xlim = c(-50, 200), 
		col = zone.cols[1], lwd = 2,
		bty = "l", #log = "x", #yaxt = "n", 
		type = "l",ylab = "range size Log10(sq km) ± s.e", xlab="Distance to PA boundary (km)", main = t)
 rug(data$bound_dist_km_PA_neg[which(data$taxon_of_interest == t & data$Zone == "Tropical")]
	, ticksize = 0.03, side = 1, lwd = 0.5, col = zone.cols[1])
  points(lbd1,zu,type="l",lty=2,  lwd = 2, col =  zone.cols[1])
  points(lbd1,zl,type="l",lty=2,  lwd = 2, col =  zone.cols[1])



#  plot(lbd,z, ylim=c(4.5,8), xlim = c(-1*log(50+1), log(200+1)), 
#		col = zone.cols[1], lwd = 2,
#		bty = "l", #log = "x", 
#		xaxt = "n", 
#		type = "l",ylab = "range size Log10(sq km) ± s.e", xlab="Distance to PA boundary (km)", main = t)
 # rug(data$log_bound_dist_km_PA_neg[which(data$taxon_of_interest == t & data$Zone == "Tropical")]
#	, ticksize = 0.03, side = 1, lwd = 0.5, col = zone.cols[1])
#  points(lbd,zu,type="l",lty=2,  lwd = 2, col =  zone.cols[1])
#  points(lbd,zl,type="l",lty=2,  lwd = 2, col =  zone.cols[1])
#  axis(1, c(-1*log(10+1), 0, log(10+1),log(50+1), log(100+1), log(200+1)), c(-10, 0, 10, 50, 100, 200))

 

  lbd <-seq(from=min(model.data$log_bound_dist_km_PA_neg[which(data$Zone == "Temperate" & data$taxon_of_interest == t)]),
		to=max(model.data$log_bound_dist_km_PA_neg[which(data$Zone == "Temperate" & data$taxon_of_interest == t)]),length=L)
  y<-vector(mode="numeric",length=L)
  yu<-vector(mode="numeric",length=L)
  yl<-vector(mode="numeric",length=L)
  



model.data$Zone<- relevel(model.data$Zone , "Temperate")
mam <- lmer(range.model$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

 for (x in 1:L)
  {
    log_bound_dist_km_PA_neg <- lbd[x]
    range <-0
    log_hpd <- mean(model.data$log_hpd)
    log_AREA.PA <- mean(model.data$log_AREA.PA )
    DoP.PA <- mean(model.data$DoP.PA )
    IUCN.PA <- "1.5"
    Zone <-"Tropical"
    taxon_of_interest<- t
    Predominant_habitat <-"Primary Vegetation"
    ag_suit <- median(model.data$ag_suit)
    log_access <- mean(model.data$log_access)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, range, log_hpd,log_AREA.PA, DoP.PA, log_access, ag_suit))
    newdat.f<-data.frame(cbind( Zone,taxon_of_interest, Predominant_habitat, Within_PA))
    levels(newdat.f$Zone)<-levels(model.data$Zone)
    levels(newdat.f$taxon_of_interest)<-levels(model.data$taxon_of_interest)
    levels(newdat.f$Predominant_habitat)<-levels(model.data$Predominant_habitat)
    levels(newdat.f$ag_suit) <- levels(model.data$ag_suit)
#    levels(newdat.f$IUCN.PA) <- levels(model.data$IUCN.PA)
    levels(newdat.f$Within_PA) <- levels(model.data$Within_PA)
    newdat<-data.frame(newdat.f,newdat.n)
    mm<-model.matrix(terms(mam),newdat)
    pvar1 <- diag(mm %*% tcrossprod(vcov(mam),mm))
    y[x]<-mm %*% fixef(mam)
    yu[x]<-y[x]+sqrt(pvar1)
    yl[x]<-y[x]-sqrt(pvar1) 
  }
  

 lbd1 <- exp(abs(lbd)) -1
 inside <- which(lbd < 0)
 lbd1[inside] <- lbd1[inside]*-1

  y <- 1/y
  yu <- 1/yu
  yl <- 1/yl


  points(lbd1,y,type="l",lty=2, col =  zone.cols[2])
  points(lbd1,yu,type="l",lty=2, col =  zone.cols[2])
  points(lbd1,yl,type="l",lty=2, col =  zone.cols[2])
  rug(data$bound_dist_km_PA_neg[which(data$taxon_of_interest == t & data$Zone == "Temperate")],
		col = zone.cols[2], lwd = 1, pos = 0.1)



#  points(lbd,y,type="l",lty=2, col =  zone.cols[2])
#  points(lbd,yu,type="l",lty=2, col =  zone.cols[2])
#  points(lbd,yl,type="l",lty=2, col =  zone.cols[2])
#  rug(data$log_bound_dist_km_PA_neg[which(data$taxon_of_interest == t & data$Zone == "Temperate")],
#		col = zone.cols[2], lwd = 1, pos = 0.1)


 abline(v = 0, lty = 2, col = 8)

} 

legend("topright", c("Tropical", "Temperate") , col = zone.cols, lty = c(2,1), lwd = c(1,2))





dev.off()








#### PA size, for different zones #####

tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/endemicity vs zone vs PA size.tif",
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
mam <- lmer(range.model$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

 for (x in 1:L)
  {
    log_bound_dist_km_PA_neg <-  mean(model.data$log_bound_dist_km_PA_neg)
    range <-0
    log_hpd <- mean(model.data$log_hpd)
    log_AREA.PA <-gis[x]
    DoP.PA <- mean(model.data$DoP.PA)
    IUCN.PA <- "1.5"
    Zone <-"Tropical"
    taxon_of_interest<- t
    Predominant_habitat <-"Primary Vegetation"
    ag_suit <- median(model.data$ag_suit)
    log_access <- mean(model.data$log_access)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, range, log_hpd,log_AREA.PA, DoP.PA, log_access, ag_suit))
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

 gis_plot <- exp(gis) -1


  plot(gis_plot,z, ylim=c(0.16,0.21), col = zone.cols[1], lwd = 2,
		bty = "l", log = "x", #yaxt = "n", 
		type = "l",ylab = "Endemicity (1/CWM range) ± s.e", xlab="PA size (km2)")
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
mam <- lmer(range.model$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

 for (x in 1:L)
  {
    log_bound_dist_km_PA_neg <-  mean(model.data$log_bound_dist_km_PA_neg)
    range <-0
    log_hpd <- mean(model.data$log_hpd)
    log_AREA.PA <-gis[x]
    DoP.PA <- mean(model.data$DoP.PA)
    IUCN.PA <- "1.5"
    Zone <-"Tropical"
    taxon_of_interest<- t
    Predominant_habitat <-"Primary Vegetation"
    ag_suit <- median(model.data$ag_suit)
    log_access <- mean(model.data$log_access)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, range, log_hpd,log_AREA.PA, DoP.PA, log_access, ag_suit))
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

  gis_plot <- exp(gis) -1

  points(gis_plot,y,type="l",lty=2, col =  zone.cols[2])
  points(gis_plot,yu,type="l",lty=2, col =  zone.cols[2])
  points(gis_plot,yl,type="l",lty=2, col =  zone.cols[2])
  rug(data$AREA.PA[which(data$Zone == "Temperate")],
		col = zone.cols[2], lwd = 0.8, pos = 0.16)



legend("topright", c("Tropical", "Temperate") , col = zone.cols, lty = 1, lwd = c(2,1))

dev.off()







#### PA age, for different zones #####

tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/endemicity vs zone vs DoP.tif",
	width = 12, height = 10, units = "cm", pointsize = 12, res = 300)



  par(mar=c(4,4.5,4,1.5))
  par(mgp=c(2.5,1,0))
  
 L = 100

  dop <-seq(from=min(model.data$DoP.PA[which(model.data$Zone == "Tropical")]),
		to=max(model.data$DoP.PA[which(model.data$Zone == "Tropical")]),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  




model.data$Zone<- relevel(model.data$Zone , "Tropical")
mam <- lmer(range.model$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

 for (x in 1:L)
  {
    log_bound_dist_km_PA_neg <- mean(model.data$log_bound_dist_km_PA_neg)
    range <-0
    log_hpd <- mean(model.data$log_hpd)
    log_AREA.PA <- mean(model.data$log_AREA.PA )
    DoP.PA <- dop[x]
    IUCN.PA <- "1.5"
    Zone <-"Tropical"
    taxon_of_interest<- "Vertebrates"
    Predominant_habitat <-"Primary Vegetation"
    ag_suit <- median(model.data$ag_suit)
    log_access <- mean(model.data$log_access)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, range, log_hpd,log_AREA.PA, DoP.PA, log_access, ag_suit))
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



  plot(dop,z, ylim=c(0.16,0.21), col = zone.cols[1], lwd = 2,
		bty = "l",
		type = "l",ylab = " Endemicity (1/CWM range) ± s.e", xlab="Duration of protection (yr)")
  rug(data$DoP.PA[which(data$Zone == "Tropical")]
	, ticksize = 0.03, side = 1, lwd = 2, col = zone.cols[1])
  points(dop,zu,type="l",lty=2,  lwd = 2, col =  zone.cols[1])
  points(dop,zl,type="l",lty=2,  lwd = 2, col =  zone.cols[1])


 dop <-seq(from=min(model.data$DoP.PA[which(model.data$Zone == "Temperate")]),
		to=max(model.data$DoP.PA[which(model.data$Zone == "Temperate")]),length=L)
   y<-vector(mode="numeric",length=L)
  yu<-vector(mode="numeric",length=L)
  yl<-vector(mode="numeric",length=L)
  



model.data$Zone<- relevel(model.data$Zone , "Temperate")
mam <- lmer(range.model$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

 for (x in 1:L)
  {
    log_bound_dist_km_PA_neg <- mean(model.data$log_bound_dist_km_PA_neg)
    range <-0
    log_hpd <- mean(model.data$log_hpd)
    log_AREA.PA <- mean(model.data$log_AREA.PA )
    DoP.PA <- dop[x]
    IUCN.PA <- "1.5"
    Zone <-"Tropical"
    taxon_of_interest<- "Vertebrates"
    Predominant_habitat <-"Primary Vegetation"
    ag_suit <- median(model.data$ag_suit)
    log_access <- mean(model.data$log_access)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, range, log_hpd,log_AREA.PA, DoP.PA, log_access, ag_suit))
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
  points(dop,yu,type="l",lty=2, col =  zone.cols[2])
  points(dop,yl,type="l",lty=2, col =  zone.cols[2])
  rug(data$DoP.PA[which(data$Zone == "Temperate")],
		col = zone.cols[2], lwd = 1, pos = 0.16)



legend("topright", c("Tropical", "Temperate") , col = zone.cols, lty = 1, lwd = c(2,1))

dev.off()







###### PA age  #############
###### for different taxa ##################

display.brewer.all()
cols <- brewer.pal(8, "Paired")
display.brewer.pal(8, "Paired")
#drop the red to avoid red-green colorblindness issues

cols <- cols[c(2,4,8)]


tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/range vs taxon vs DoP.tif",
	width = 10, height = 25, units = "cm", pointsize = 12, res = 300)

par(mfrow = c(3,1))
taxa <- levels(matched.landuse$taxon_of_interest)
i <- 0
t <- taxa[1]

for(t in taxa){

i <- i +1

  par(mar=c(4,4.5,4,1.5))
  par(mgp=c(2.5,1,0))
  
 L = 30

  dop <-seq(from=min(model.data$DoP.PA),to=max(model.data$DoP.PA),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  




model.data$Within_PA <- relevel(model.data$Within_PA, "no")
model.data$Zone<- relevel(model.data$Zone , "Tropical")
model.data$taxon_of_interest <- relevel(model.data$taxon_of_interest , t)
mam <- lmer(range.model$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

 for (x in 1:L)
  {
    log_bound_dist_km_PA_neg <- mean(model.data$log_bound_dist_km_PA_neg )
    range <-0
    log_hpd <- mean(model.data$log_hpd)
    log_AREA.PA <- mean(model.data$log_AREA.PA)
    DoP.PA <- dop[x]
    Zone <-"Tropical"
    taxon_of_interest<- t
    Predominant_habitat <-"Primary Vegetation"
    ag_suit <- median(model.data$ag_suit)
    log_access <- mean(model.data$log_access)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, range, log_hpd,log_AREA.PA, DoP.PA, log_access, ag_suit))
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
  


  plot(dop,z, ylim=c(4.7,8), col = cols[i],
		bty = "l", #log = "x", #yaxt = "n", 
		type = "l",ylab = "range size Log10(sq km) ± s.e", xlab="Duration of protection (yr)", main = t)
  rug(data$DoP.PA[which(matched.landuse$taxon_of_interest == t)]
	, ticksize = 0.03, side = 1, lwd = 0.5, pos = 4.74, col = cols[i])
  points(dop,zu,type="l",lty=2, col = cols[i])
  points(dop,zl,type="l",lty=2, col = cols[i])



}



dev.off()






####  Within PA and access ##########



tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/endemicity within pa vs accessibility.tif",
	width = 12, height = 10, units = "cm", pointsize = 12, res = 300)



par(mfrow = c(1,1))


  par(mar=c(4,4.5,1,1.5))
  par(mgp=c(2.5,1,0))
  
 L = 100

  acc <-seq(from=min(model.data$log_access[which(model.data$Within_PA == "no")]),
		to=max(model.data$log_access[which(model.data$Within_PA == "no")]),length=L)
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
    log_access <- acc[x]
    range <-0
    log_hpd <- mean(model.data$log_hpd)
    ag_suit <- median(model.data$ag_suit)
    log_AREA.PA <- mean(model.data$log_AREA.PA )
    DoP.PA <- mean(model.data$DoP.PA )
    Zone <-"Temperate"
    taxon_of_interest<- "Invertebrates"
    Predominant_habitat <-"Primary Vegetation"
    log_bound_dist_km_PA_neg <- mean(model.data$log_bound_dist_km_PA_neg)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, log_access, log_AREA.PA, DoP.PA, range, log_hpd, ag_suit))
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
  
  acc_plot <- (exp(acc) -1)/60

  z <- 1/z
  zu <- 1/zu
  zl <- 1/zl

  plot(acc_plot,z, ylim=c(0.16,0.21), xlim = c(0.02,60), col = 8,
		bty = "l", log = "x", xaxt = "n",
		type = "l",ylab = "Endemicity (1/CWM range) ± s.e", 
		xlab="Accessibility (hours to city >50000)")
   rug((data$access[which(data$Within_PA == "no")])/60
	, ticksize = 0.03, side = 1,   lwd = 1, col = 8) 
  axis(1, c(0.05, 0.1, 0.2, 0.5, 1, 2, 5, 20, 50), c(0.05, 0.1, 0.2, 0.5, 1, 2, 5, 20, 50))
  points(acc_plot,zu,type="l",lty=2, col = 8)
  points(acc_plot,zl,type="l",lty=2, col = 8)





model.data$Within_PA <- relevel(model.data$Within_PA, "yes")

mam<- lmer(range.model$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

acc <-seq(from=min(model.data$log_access[which(model.data$Within_PA == "yes")]),
		to=max(model.data$log_access[which(model.data$Within_PA == "yes")]),length=L)
  y<-vector(mode="numeric",length=length(hpd))
  yu<-vector(mode="numeric",length=length(hpd))
  yl<-vector(mode="numeric",length=length(hpd))


 for (x in 1:L)
  {
    log_access <- acc[x]
    range <-0
    log_hpd <- mean(model.data$log_hpd)
    ag_suit <- median(model.data$ag_suit)
    log_AREA.PA <- mean(model.data$log_AREA.PA )
    DoP.PA <- mean(model.data$DoP.PA )
    Zone <-"Temperate"
    taxon_of_interest<- "Invertebrates"
    Predominant_habitat <-"Primary Vegetation"
    log_bound_dist_km_PA_neg <- mean(model.data$log_bound_dist_km_PA_neg)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, log_access, log_AREA.PA, DoP.PA, range, log_hpd, ag_suit))
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

  acc_plot <- (exp(acc) -1)/60

  y <- 1/y
  yu <- 1/yu
  yl <- 1/yl

  points(acc_plot,y,type="l",lty=1 , col = 1, lwd = 2)
  points(acc_plot,yu,type="l",lty=3, lwd = 2, col = 1)
  points(acc_plot,yl,type="l",lty=3, lwd = 2, col = 1)
   rug((data$access[which(data$Within_PA == "yes")])/60
	, ticksize = 0.03, side = 1,   lwd = 1, col = 1, pos = 0.16) 





legend("topright", c("Protected", "Unprotected") , col = c(1,8), lty = c(1,1), lwd = c(1,2))


dev.off()






### range vs accessibility ###



tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/range vs accessibility.tif",
	width = 12, height = 10, units = "cm", pointsize = 12, res = 300)



par(mfrow = c(1,1))

#Zones <- c("Tropical", "Temperate")

#for(o in Zones){ 

  par(mar=c(4,4.5,1,1.5))
  par(mgp=c(2.5,1,0))
  
 L = 30

  acc <-seq(from=min(model.data$log_access),to=max(model.data$log_access),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  


# it doesnt matter what the values for taxon and zone are, as there is no interaction with these




model.data$taxon_of_interest <- relevel(model.data$taxon_of_interest , "Invertebrates")
mam <- lmer(range.model$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

 for (x in 1:L)
  {
    log_access <- acc[x]
    range <-0
    log_hpd <- mean(model.data$log_hpd)
    ag_suit <- median(model.data$ag_suit)
    log_AREA.PA <- mean(model.data$log_AREA.PA )
    DoP.PA <- mean(model.data$DoP.PA )
    IUCN.PA <- "1.5"
    Zone <-"Temperate"
    taxon_of_interest<- "Invertebrates"
    Predominant_habitat <-"Primary Vegetation"
    log_bound_dist_km_PA_neg <- mean(model.data$log_bound_dist_km_PA_neg)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, log_access, log_AREA.PA, DoP.PA, range, log_hpd, ag_suit))
    newdat.f<-data.frame(cbind( Zone,taxon_of_interest, Predominant_habitat, IUCN.PA))
    levels(newdat.f$Zone)<-levels(model.data$Zone)
    levels(newdat.f$taxon_of_interest)<-levels(model.data$taxon_of_interest)
    levels(newdat.f$Predominant_habitat)<-levels(model.data$Predominant_habitat)
   # levels(newdat.f$Within_PA) <- levels(model.data$Within_PA)
    levels(newdat.f$IUCN.PA) <- levels(model.data$IUCN.PA)
    newdat<-data.frame(newdat.f,newdat.n)
    mm<-model.matrix(terms(mam),newdat)
    pvar1 <- diag(mm %*% tcrossprod(vcov(mam),mm))
    z[x]<-mm %*% fixef(mam)
    zu[x]<-z[x]+sqrt(pvar1)
    zl[x]<-z[x]-sqrt(pvar1) 
  }
  
  acc_plot <- (exp(acc) -1)/60

  plot(acc_plot,z, ylim=c(4.7,6.5), col = 1,
		bty = "l", log = "x",
		type = "l",ylab = "range Log10(sq km) ± s.e", xlab="Accessibility (hours to city >50000)")
   rug((data$access)/60
	, ticksize = 0.03, side = 1,   lwd = 1, col = 8) 

  points(acc_plot,zu,type="l",lty=2, col = 1)
  points(acc_plot,zl,type="l",lty=2, col = 1)




dev.off()







#### within PA and ag suit ##########



tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/range within pa vs  ag_suit.tif",
	width = 12, height = 10, units = "cm", pointsize = 12, res = 300)

par(mfrow = c(1,1))


  par(mar=c(4,4.5,4,1.5))
  par(mgp=c(2.5,1,0))
  
 L = 30

  ag <-seq(from=min(model.data$ag_suit),to=max(model.data$ag_suit),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  


# it doesnt matter what the values for taxon and zone are, as there is no interaction with these



model.data$Within_PA <- relevel(model.data$Within_PA, "no")
model.data$Zone<- relevel(model.data$Zone, "Tropical")
model.data$taxon_of_interest <- relevel(model.data$taxon_of_interest , "Invertebrates")
mam <- lmer(range.model$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

 for (x in 1:L)
  {
    ag_suit <- ag[x]
    range <-0
    log_hpd <- mean(model.data$log_hpd)
    log_access <- mean(model.data$log_access)
    log_AREA.PA <- mean(model.data$log_AREA.PA )
    DoP.PA <- mean(model.data$DoP.PA )
    Zone <-"Temperate"
    taxon_of_interest<- "Invertebrates"
    Predominant_habitat <-"Primary Vegetation"
    log_bound_dist_km_PA_neg <- mean(model.data$log_bound_dist_km_PA_neg)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, log_access, log_AREA.PA, DoP.PA, range, log_hpd, ag_suit))
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
  


  plot(ag,z, ylim=c(4.8,6), col = 8,
		bty = "l", #log = "x",
		type = "l",ylab = "range Log10(sq km) ± s.e", xlab="Agricultural suitability (higher = more suitable)")
  points(ag,zu,type="l",lty=2, col = 8)
  points(ag,zl,type="l",lty=2, col = 8)




model.data$Within_PA <- relevel(model.data$Within_PA, "yes")

mam<- lmer(range.model$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))


  y<-vector(mode="numeric",length=length(hpd))
  yu<-vector(mode="numeric",length=length(hpd))
  yl<-vector(mode="numeric",length=length(hpd))


 for (x in 1:L)
  {
    ag_suit <- ag[x]
    range <-0
    log_hpd <- mean(model.data$log_hpd)
    log_access <- mean(model.data$log_access)
    log_AREA.PA <- mean(model.data$log_AREA.PA )
    DoP.PA <- mean(model.data$DoP.PA )
    Zone <-"Temperate"
    taxon_of_interest<- "Invertebrates"
    Predominant_habitat <-"Primary Vegetation"
    log_bound_dist_km_PA_neg <- mean(model.data$log_bound_dist_km_PA_neg)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, log_access, log_AREA.PA, DoP.PA, range, log_hpd, ag_suit))
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



  
  points(ag,y,type="l",lty=1 , col = 1, lwd = 2)
  points(ag,yu,type="l",lty=3, lwd = 2, col = 1)
  points(ag,yl,type="l",lty=3, lwd = 2, col = 1)


legend("topright", c("Protected", "Unprotected") , col = c(1,8), lty = c(1,1), lwd = c(1,2))



dev.off()











#### hpd ##########



tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/range vs hpd.tif",
	width = 12, height = 10, units = "cm", pointsize = 12, res = 300)

par(mfrow = c(1,1))


  par(mar=c(4,4.5,4,1.5))
  par(mgp=c(2.5,1,0))
  
 L = 30

  hpd <-seq(from=min(model.data$log_hpd),to=max(model.data$log_hpd),length=L)
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
    log_hpd <- hpd[x]
    range <-0
    log_access <- mean(model.data$log_access)
    ag_suit <- median(model.data$ag_suit)
    log_AREA.PA <- mean(model.data$log_AREA.PA )
    DoP.PA <- mean(model.data$DoP.PA )
    Zone <-"Temperate"
    taxon_of_interest<- "Invertebrates"
    Predominant_habitat <-"Primary Vegetation"
    log_bound_dist_km_PA_neg <- mean(model.data$log_bound_dist_km_PA_neg)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, log_access, range, log_AREA.PA, DoP.PA, log_hpd, ag_suit))
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
  

  hpd <- exp(hpd) - 1


  plot(hpd,z, ylim=c(4.8,6), col = 8,
		bty = "l", log = "x",
		type = "l",ylab = "range Log10(sq km) ± s.e", xlab="Human population density (per km2)")
   rug(data$hpd, ticksize = 0.03, side = 1,  pos = 4.8, lwd = 0.5, col = 8) #### change to be data used in model, poss slightly fewer points


  points(hpd,zu,type="l",lty=2, col = 8)
  points(hpd,zl,type="l",lty=2, col = 8)







dev.off()




####






####  IUCN CAT and hpd ##########


titles <- c("I and II", "III to VI", "NA or Not Reported", "Unprotected")




tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/range IUCN cat vs hpd.tif",
	width = 20, height = 20, units = "cm", pointsize = 12, res = 300)



par(mfrow = c(2,2))

cats <- levels(matched.landuse$IUCN.PA)
c <- cats[1]

for(c in cats){ 

  par(mar=c(4,4.5,1,1.5))
  par(mgp=c(2.5,1,0))
  
 L = 30

  hpd<-seq(from=min(model.data$log_hpd),to=max(model.data$log_hpd),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  


# it doesnt matter what the values for taxon and zone are, as there is no interaction with these




model.data$IUCN.PA <- relevel(model.data$IUCN.PA, c)
mam <- lmer(range.model$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

 for (x in 1:L)
  {
    log_access <- mean(model.data$log_access)
    range <-0
    log_hpd <- hpd[x]
    ag_suit <- median(model.data$ag_suit)
    log_AREA.PA <- mean(model.data$log_AREA.PA )
    DoP.PA <- mean(model.data$DoP.PA )
    IUCN.PA <- c
    Zone <-"Temperate"
    taxon_of_interest<- "Invertebrates"
    Predominant_habitat <-"Primary Vegetation"
    log_bound_dist_km_PA_neg <- mean(model.data$log_bound_dist_km_PA_neg)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, log_access, log_AREA.PA, DoP.PA, range, log_hpd, ag_suit))
    newdat.f<-data.frame(cbind( Zone,taxon_of_interest, Predominant_habitat, IUCN.PA))
    levels(newdat.f$Zone)<-levels(model.data$Zone)
    levels(newdat.f$taxon_of_interest)<-levels(model.data$taxon_of_interest)
    levels(newdat.f$Predominant_habitat)<-levels(model.data$Predominant_habitat)
    levels(newdat.f$IUCN.PA) <- levels(model.data$IUCN.PA)
    newdat<-data.frame(newdat.f,newdat.n)
    mm<-model.matrix(terms(mam),newdat)
    pvar1 <- diag(mm %*% tcrossprod(vcov(mam),mm))
    z[x]<-mm %*% fixef(mam)
    zu[x]<-z[x]+sqrt(pvar1)
    zl[x]<-z[x]-sqrt(pvar1) 
  }
  
  hpd_plot <- (exp(hpd) -1)

  plot(hpd_plot,z, ylim=c(4.7,7), col = 8, 
		bty = "l", log = "x",
		type = "l",ylab = "range Log10(sq km) ± s.e", xlab="Human population density (per km2)",
		main = titles[which(cats == c)])
   rug(data$hpd[which(data$IUCN.PA == c)]
	, ticksize = 0.03, side = 1,   lwd = 1, col = 8) 

  points(hpd_plot,zu,type="l",lty=2, col = 8)
  points(hpd_plot,zl,type="l",lty=2, col = 8)


}


dev.off()








####  IUCN CAT and ag suit ##########


titles <- c("I and II", "III to VI", "NA or Not Reported", "Unprotected")

tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/range IUCN cat vs ag suit.tif",
	width = 20, height = 20, units = "cm", pointsize = 12, res = 300)



par(mfrow = c(2,2))

cats <- levels(matched.landuse$IUCN.PA)
c <- cats[1]

for(c in cats){ 

  par(mar=c(4,4.5,1,1.5))
  par(mgp=c(2.5,1,0))
  
 L = 30

  ag<-seq(from=min(model.data$ag_suit),to=max(model.data$ag_suit),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  



model.data$IUCN.PA <- relevel(model.data$IUCN.PA, c)
mam <- lmer(range.model$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

 for (x in 1:L)
  {
    log_access <- mean(model.data$log_access)
    range <-0
    ag_suit <- ag[x]
    log_hpd<- mean(model.data$log_hpd)
    log_AREA.PA <- mean(model.data$log_AREA.PA )
    DoP.PA <- mean(model.data$DoP.PA )
    IUCN.PA <- c
    Zone <-"Temperate"
    taxon_of_interest<- "Invertebrates"
    Predominant_habitat <-"Primary Vegetation"
    log_bound_dist_km_PA_neg <- mean(model.data$log_bound_dist_km_PA_neg)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, log_access, log_AREA.PA, DoP.PA, range, log_hpd, ag_suit))
    newdat.f<-data.frame(cbind( Zone,taxon_of_interest, Predominant_habitat, IUCN.PA))
    levels(newdat.f$Zone)<-levels(model.data$Zone)
    levels(newdat.f$taxon_of_interest)<-levels(model.data$taxon_of_interest)
    levels(newdat.f$Predominant_habitat)<-levels(model.data$Predominant_habitat)
    levels(newdat.f$IUCN.PA) <- levels(model.data$IUCN.PA)
    newdat<-data.frame(newdat.f,newdat.n)
    mm<-model.matrix(terms(mam),newdat)
    pvar1 <- diag(mm %*% tcrossprod(vcov(mam),mm))
    z[x]<-mm %*% fixef(mam)
    zu[x]<-z[x]+sqrt(pvar1)
    zl[x]<-z[x]-sqrt(pvar1) 
  }
  


  plot(ag,z, ylim=c(4.7,7), col = 8, 
		bty = "l", 
		type = "l",ylab = "range Log10(sq km) ± s.e", xlab="Agricultural suitability",
		main = titles[which(cats == c)])
   rug(data$ag_suit[which(data$IUCN.PA == c)]
	, ticksize = 0.03, side = 1,   lwd = 1, col = 8) 

  points(ag,zu,type="l",lty=2, col = 8)
  points(ag,zl,type="l",lty=2, col = 8)


}


dev.off()




#### diff IUCN CATS ####




IUCN <- vector()
estimate <- vector()
se <- vector()



model.data$IUCN.PA <- relevel(model.data$IUCN.PA, "1.5")
mam <- lmer (range.model$final.call,  data = model.data)
estimate[1] <- fixef(mam)[1]
se[1] <- se.fixef(mam)[1]
IUCN [1] <- "I & II"



model.data$IUCN.PA  <- relevel(model.data$IUCN.PA , "4.5")
mam <- lmer (range.model$final.call,  data = model.data)
estimate[2] <- fixef(mam)[1]
se[2] <- se.fixef(mam)[1]
IUCN [2] <- "III - VI"


model.data$IUCN.PA  <- relevel(model.data$IUCN.PA , "7")
mam <- lmer (range.model$final.call,  data = model.data)
estimate[3] <- fixef(mam)[1]
se[3] <- se.fixef(mam)[1]
IUCN [3] <- "Unknown"


model.data$IUCN.PA  <- relevel(model.data$IUCN.PA , "0")
mam <- lmer (range.model$final.call,  data = model.data)
estimate[4] <- fixef(mam)[1]
se[4] <- se.fixef(mam)[1]
IUCN [4] <- "Unprotected"


plot.data <- data.frame(IUCN = IUCN, se = se, estimate = estimate)
plot.data$upper <- plot.data$estimate + 1.96*plot.data$se
plot.data$lower <- plot.data$estimate - 1.96*plot.data$se
pos <- c(1:4)

cats <- c("I & II", "III - VI", "Unknown", "Unprotected")


tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/range IUCN.tif",
	width = 15, height = 12, units = "cm", pointsize = 12, res = 300)


plot(plot.data$estimate ~ pos, cex = 1, col = 1,
	ylim = c(5,7), 
	bty = "l", xaxt = "n", xlab = "", ylab = "range Log10(sq km) ± s.e")
arrows(pos, plot.data$lower, pos, plot.data$upper, angle = 90, code = 3, length = 0.05)
axis(1, pos, cats)


dev.off()





####  IUCN CAT and PA size ##########


titles <- c("I and II", "III to VI", "NA or Not Reported", "Unprotected")

tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/range vs PA size.tif",
	width = 23, height = 10, units = "cm", pointsize = 12, res = 300)



  par(mar=c(4,4.5,1,1.5))
  par(mgp=c(2.5,1,0))
  
 L = 30

  gis <-seq(from=min(model.data$log_AREA.PA),to=max(model.data$log_AREA.PA),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  



mam <- lmer(range.model$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

 for (x in 1:L)
  {
    log_access <- mean(model.data$log_access)
    range <-0
    ag_suit <- median(model.data$ag_suit)
    log_hpd<- mean(model.data$log_hpd)
    log_AREA.PA <- gis[x]
    DoP.PA <- mean(model.data$DoP.PA )
    IUCN.PA <- c
    Zone <-"Temperate"
    taxon_of_interest<- "Invertebrates"
    Predominant_habitat <-"Primary Vegetation"
    log_bound_dist_km_PA_neg <- mean(model.data$log_bound_dist_km_PA_neg)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, log_access, log_AREA.PA, DoP.PA, range, log_hpd, ag_suit))
    newdat.f<-data.frame(cbind( Zone,taxon_of_interest, Predominant_habitat, IUCN.PA))
    levels(newdat.f$Zone)<-levels(model.data$Zone)
    levels(newdat.f$taxon_of_interest)<-levels(model.data$taxon_of_interest)
    levels(newdat.f$Predominant_habitat)<-levels(model.data$Predominant_habitat)
    levels(newdat.f$IUCN.PA) <- levels(model.data$IUCN.PA)
    newdat<-data.frame(newdat.f,newdat.n)
    mm<-model.matrix(terms(mam),newdat)
    pvar1 <- diag(mm %*% tcrossprod(vcov(mam),mm))
    z[x]<-mm %*% fixef(mam)
    zu[x]<-z[x]+sqrt(pvar1)
    zl[x]<-z[x]-sqrt(pvar1) 
  }
  
 gis <- exp(gis) -1

  plot(gis,z, ylim=c(4.7,7), col = 8, 
		bty = "l",  log = "x",
		type = "l",ylab = "range Log10(sq km) ± s.e", xlab="Agricultural suitability",
		main = titles[which(cats == c)])
   rug(data$AREA.PA	, ticksize = 0.03, side = 1,   lwd = 1, col = 8) 

  points(gis,zu,type="l",lty=2, col = 8)
  points(gis,zl,type="l",lty=2, col = 8)


}


dev.off()




