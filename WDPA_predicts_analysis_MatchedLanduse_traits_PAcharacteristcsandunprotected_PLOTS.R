# PLOTS



###
# mass
####






model.data <- mass.model$data
nrow(model.data)
length(unique(model.data$SS))

data <- matched.landuse[,c("mass", "Zone", "log_bound_dist_km_PA_neg", "bound_dist_km_PA_neg", "DoP.PA", 
		"log_AREA.PA", "AREA.PA",
		"access", "log_access", "hpd", "log_hpd", "ag_suit")]
data <- na.omit(data)
nrow(data)




### mass vs  dist to boundary ####





tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/mass vs dist to PA boundary.tif",
	width = 15, height = 10, units = "cm", pointsize = 12, res = 300)

par(mfrow = c(1,1))


  par(mar=c(4,4.5,1,1.5))
  par(mgp=c(2.5,1,0))
  
 L = 100

  lbd <-seq(from=min(model.data$log_bound_dist_km_PA_neg),to=max(model.data$log_bound_dist_km_PA_neg),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  


# it doesnt matter what the values for taxon and zone are, as there is no interaction with these



model.data$Within_PA <- relevel(model.data$Within_PA, "no")
model.data$Zone<- relevel(model.data$Zone , "Tropical")
mam <- lmer(mass.model$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

 for (x in 1:L)
  {
    ag_suit <- median(model.data$ag_suit)
    log_hpd <- mean(model.data$log_hpd)
    log_AREA.PA <- mean(model.data$log_AREA.PA)
    mass <-0
    log_bound_dist_km_PA_neg <- lbd[x]
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg,mass, log_hpd, ag_suit, log_AREA.PA))
    newdat.f<-data.frame(cbind(Within_PA))
    levels(newdat.f$Within_PA) <- levels(model.data$Within_PA)
    newdat<-data.frame(newdat.f,newdat.n)
    mm<-model.matrix(terms(mam),newdat)
    pvar1 <- diag(mm %*% tcrossprod(vcov(mam),mm))
    z[x]<-mm %*% fixef(mam)
    zu[x]<-z[x]+sqrt(pvar1)
    zl[x]<-z[x]-sqrt(pvar1) 
  }
  


  z  <- 10^(z) 
  zu <- 10^(zu)
  zl <- 10^(zl)

 lbd1 <- exp(abs(lbd)) -1
 inside <- which(lbd < 0)
 lbd1[inside] <- lbd1[inside]*-1

  plot(lbd1,z, ylim=c(10,80), col = 1,
		bty = "l", 
		type = "l",ylab = "CWM Mass (g) ± s.e", xlab="Distance to protected area boundary (km)")
  points(lbd1,zu,type="l",lty=2, col = 1)
  points(lbd1,zl,type="l",lty=2, col = 1)
  rug(data$bound_dist_km_PA_neg
	, ticksize = 0.03, side = 1,  lwd = 0.5, col = 1) #### change to be data used in model, poss slightly fewer points

abline(v = 0, lty = 2, col = 8)



dev.off()







### mass vs ag suit ####



tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/mass vs ag_suit.tif",
	width = 15, height = 10, units = "cm", pointsize = 12, res = 300)

par(mfrow = c(1,1))


  par(mar=c(4,4.5,1,1.5))
  par(mgp=c(2.5,1,0))
  
 L = 100

  ag <-seq(from=0,to=max(model.data$ag_suit),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  


# it doesnt matter what the values for taxon and zone are, as there is no interaction with these



model.data$Within_PA <- relevel(model.data$Within_PA, "no")
model.data$Zone<- relevel(model.data$Zone , "Tropical")
mam <- lmer(mass.model$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

 for (x in 1:L)
  {
 ag_suit <- ag[x]
    mass <-0
    log_hpd <- mean(model.data$log_hpd)
    log_AREA.PA <- mean(model.data$log_AREA.PA)
    log_bound_dist_km_PA_neg <- mean(model.data$log_bound_dist_km_PA_neg)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg,mass, log_hpd, log_AREA.PA, ag_suit))
    newdat.f<-data.frame(cbind(Within_PA))
    levels(newdat.f$Within_PA) <- levels(model.data$Within_PA)
    newdat<-data.frame(newdat.f,newdat.n)
    mm<-model.matrix(terms(mam),newdat)
    pvar1 <- diag(mm %*% tcrossprod(vcov(mam),mm))
    z[x]<-mm %*% fixef(mam)
    zu[x]<-z[x]+sqrt(pvar1)
    zl[x]<-z[x]-sqrt(pvar1) 
  }
  


  z  <- 10^(z) 
  zu <- 10^(zu)
  zl <- 10^(zl)

  plot(ag,z, ylim=c(10,80), col = 1,
		bty = "l", #log = "x",
		type = "l",ylab = "CWM Mass (g) ± s.e ± s.e", xlab="Agricultural suitability (higher = more suitable)")
  points(ag,zu,type="l",lty=2, col = 1)
  points(ag,zl,type="l",lty=2, col = 1)





dev.off()









#### mass vs hpd #####



tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/mass vs hpd.tif",
	width = 15, height = 10, units = "cm", pointsize = 12, res = 300)

par(mfrow = c(1,1))


  par(mar=c(4,4.5,1,1.5))
  par(mgp=c(2.5,1,0))
  
 L = 30

  hpd <-seq(from=0,to=max(model.data$log_hpd),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  


# it doesnt matter what the values for taxon and zone are, as there is no interaction with these



model.data$Within_PA <- relevel(model.data$Within_PA, "no")
model.data$Zone<- relevel(model.data$Zone , "Tropical")
mam <- lmer(mass.model$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

 for (x in 1:L)
  {
    ag_suit <- median(model.data$ag_suit)
    log_hpd <- hpd[x]
    mass <-0
    log_bound_dist_km_PA_neg <- mean(model.data$log_bound_dist_km_PA_neg)
    log_AREA.PA <- mean(model.data$log_AREA.PA)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg,mass, log_hpd, log_AREA.PA, ag_suit))
    newdat.f<-data.frame(cbind(Within_PA))
    levels(newdat.f$Within_PA) <- levels(model.data$Within_PA)
    newdat<-data.frame(newdat.f,newdat.n)
    mm<-model.matrix(terms(mam),newdat)
    pvar1 <- diag(mm %*% tcrossprod(vcov(mam),mm))
    z[x]<-mm %*% fixef(mam)
    zu[x]<-z[x]+sqrt(pvar1)
    zl[x]<-z[x]-sqrt(pvar1) 
  }
  


  z  <- 10^(z) 
  zu <- 10^(zu)
  zl <- 10^(zl)

  hpd <- exp(hpd)

  plot(hpd,z, ylim=c(10,80), col = 1,
		bty = "l", log = "x",
		type = "l",ylab = "CWM Mass (g) ± s.e", xlab="Human population density (per km2)")
  points(hpd,zu,type="l",lty=2, col = 1)
  points(hpd,zl,type="l",lty=2, col = 1)
  rug(data$hpd
	, ticksize = 0.03, side = 1,  lwd = 0.5, col = 1) #### change to be data used in model, poss slightly fewer points





dev.off()






#### mass vs area #####



tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/mass vs PA size.tif",
	width = 15, height = 10, units = "cm", pointsize = 12, res = 300)

par(mfrow = c(1,1))


  par(mar=c(4,4.5,1,1.5))
  par(mgp=c(2.5,1,0))
  
 L = 30

  gis <-seq(from=min(model.data$log_AREA.PA),to=max(model.data$log_AREA.PA),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  


# it doesnt matter what the values for taxon and zone are, as there is no interaction with these



model.data$Within_PA <- relevel(model.data$Within_PA, "no")
model.data$Zone<- relevel(model.data$Zone , "Tropical")
mam <- lmer(mass.model$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

 for (x in 1:L)
  {
    ag_suit <- median(model.data$ag_suit)
    log_hpd <- mean(model.data$log_hpd)
    mass <-0
    log_bound_dist_km_PA_neg <- mean(model.data$log_bound_dist_km_PA_neg)
    log_AREA.PA <- gis[x]
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg,mass, log_hpd, log_AREA.PA, ag_suit))
    newdat.f<-data.frame(cbind(Within_PA))
    levels(newdat.f$Within_PA) <- levels(model.data$Within_PA)
    newdat<-data.frame(newdat.f,newdat.n)
    mm<-model.matrix(terms(mam),newdat)
    pvar1 <- diag(mm %*% tcrossprod(vcov(mam),mm))
    z[x]<-mm %*% fixef(mam)
    zu[x]<-z[x]+sqrt(pvar1)
    zl[x]<-z[x]-sqrt(pvar1) 
  }
  


  z  <- 10^(z) 
  zu <- 10^(zu)
  zl <- 10^(zl)

  gis <- exp(gis) -1

  plot(gis,z, ylim=c(10,80), col = 1,
		bty = "l", log = "x",
		type = "l",ylab = "CWM Mass (g) ± s.e", xlab="PA size (km2)")
  points(gis,zu,type="l",lty=2, col = 1)
  points(gis,zl,type="l",lty=2, col = 1)
  rug(data$AREA.PA
	, ticksize = 0.03, side = 1,  lwd = 0.5, col = 1) #### change to be data used in model, poss slightly fewer points





dev.off()












###
#veg height
#####




model.data <- veg.model$data

nrow(model.data)
length(unique(model.data$SS))

data <- matched.landuse[,c("veg", "Zone", "log_bound_dist_km_PA_neg", "bound_dist_km_PA_neg", "DoP.PA", 
		"log_AREA.PA", "AREA.PA",
		"access", "log_access", "hpd", "log_hpd", "ag_suit")]
data <- na.omit(data)
nrow(data)




#### ag suit ##########



tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/veg height vs ag_suit.tif",
	width = 15, height = 10, units = "cm", pointsize = 12, res = 300)

par(mfrow = c(1,1))


  par(mar=c(4,4.5,4,1.5))
  par(mgp=c(2.5,1,0))
  
 L = 30

  ag <-seq(from=0,to=max(model.data$ag_suit),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  


# it doesnt matter what the values for taxon and zone are, as there is no interaction with these



model.data$Within_PA <- relevel(model.data$Within_PA, "no")
model.data$Zone<- relevel(model.data$Zone , "Tropical")
mam <- lmer(veg.model$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

 for (x in 1:L)
  {
    ag_suit <- ag[x]
    veg <-0
    log_hpd <- mean(model.data$log_hpd)
    Zone <-"Temperate"
    log_AREA.PA <- mean(model.data$log_AREA.PA)
    Predominant_habitat <-"Primary Vegetation"
    log_bound_dist_km_PA_neg <- mean(model.data$log_bound_dist_km_PA_neg)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, veg, log_hpd, log_AREA.PA, ag_suit))
    newdat.f<-data.frame(cbind( Zone, Predominant_habitat, Within_PA))
    levels(newdat.f$Zone)<-levels(model.data$Zone)
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
  


  z  <- 10^(z) 
  zu <- 10^(zu)
  zl <- 10^(zl)

  plot(ag,z, ylim=c(0.1,20), col = 1,
		bty = "l", #log = "x",
		type = "l",ylab = "CWM Vegetative height (m) ± s.e", xlab="Agricultural suitability (higher = more suitable)")
  points(ag,zu,type="l",lty=2, col = 1)
  points(ag,zl,type="l",lty=2, col = 1)


dev.off()





### PA size ##########



tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/veg height vs ag_suit.tif",
	width = 15, height = 10, units = "cm", pointsize = 12, res = 300)

par(mfrow = c(1,1))


  par(mar=c(4,4.5,4,1.5))
  par(mgp=c(2.5,1,0))
  
 L = 30

  gis <-seq(from=min(model.data$log_AREA.PA),to=max(model.data$log_AREA.PA),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  


# it doesnt matter what the values for taxon and zone are, as there is no interaction with these



model.data$Within_PA <- relevel(model.data$Within_PA, "no")
model.data$Zone<- relevel(model.data$Zone , "Tropical")
mam <- lmer(veg.model$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

 for (x in 1:L)
  {
    ag_suit <- mean(model.data$ag_suit)
    veg <-0
    log_hpd <- mean(model.data$log_hpd)
    Zone <-"Temperate"
    log_AREA.PA <- gis[x]
    Predominant_habitat <-"Primary Vegetation"
    log_bound_dist_km_PA_neg <- mean(model.data$log_bound_dist_km_PA_neg)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, veg, log_hpd, log_AREA.PA, ag_suit))
    newdat.f<-data.frame(cbind( Zone, Predominant_habitat, Within_PA))
    levels(newdat.f$Zone)<-levels(model.data$Zone)
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
  


  z  <- 10^(z) 
  zu <- 10^(zu)
  zl <- 10^(zl)

  gis <- exp(gis) -1

  plot(gis,z, ylim=c(0.1,20), col = 1,
		bty = "l", log = "x",
		type = "l",ylab = "CWM Vegetative height (m) ± s.e", xlab="Agricultural suitability (higher = more suitable)")
  points(gis,zu,type="l",lty=2, col = 1)
  points(gis,zl,type="l",lty=2, col = 1)


dev.off()












####
# inverts vol
#####



invert.vol.data <- inverts.data[,c("vol","access","ag_suit")]

invert.vol.data <- na.omit(invert.vol.data)




#### ag suit ##########


model.data <- vol.model$data

tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/invert vol vs ag suit.tif",
	width = 15, height = 10, units = "cm", pointsize = 12, res = 300)

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
model.data$Zone<- relevel(model.data$Zone , "Tropical")
mam <- lmer(vol.model$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

 for (x in 1:L)
  {
    ag_suit <- ag[x]
    vol <-0
    log_access <- mean(model.data$log_access)
    log_hpd <- mean(model.data$log_hpd)
    log_AREA.PA <- mean(model.data$log_AREA.PA)
    Within_PA <- "no"
    Zone <- "Tropical"
    newdat.n<-data.frame(cbind(vol, log_access, ag_suit, log_AREA.PA, log_hpd))
    newdat.f<-data.frame(cbind(Within_PA, Zone))
    levels(newdat.f$Within_PA) <- levels(model.data$Within_PA)
    levels(newdat.f$Zone) <- levels(model.data$Zone)
    newdat<-data.frame(newdat.f,newdat.n)
    mm<-model.matrix(terms(mam),newdat)
    pvar1 <- diag(mm %*% tcrossprod(vcov(mam),mm))
    z[x]<-mm %*% fixef(mam)
    zu[x]<-z[x]+sqrt(pvar1)
    zl[x]<-z[x]-sqrt(pvar1) 
  }
  


  z  <- 10^(z) 
  zu <- 10^(zu)
  zl <- 10^(zl)


  plot(ag,z, ylim=c(0,7000), col = 1,
		bty = "l", #log = "x",
		type = "l",ylab = "CWM Length derived vol(mm3) ± s.e", xlab="Agricultural suitability")
  points(ag,zu,type="l",lty=2, col = 1)
  points(ag,zl,type="l",lty=2, col = 1)





dev.off()





#### invert vol vs access ##########


model.data <- vol.model$data

tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/invert vol vs access.tif",
	width = 15, height = 10, units = "cm", pointsize = 12, res = 300)

par(mfrow = c(1,1))


  par(mar=c(4,4.5,4,1.5))
  par(mgp=c(2.5,1,0))
  
 L = 30

  acc <-seq(from=min(model.data$log_access),to=max(model.data$log_access),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  


# it doesnt matter what the values for taxon and zone are, as there is no interaction with these



model.data$Within_PA <- relevel(model.data$Within_PA, "no")
model.data$Zone<- relevel(model.data$Zone , "Tropical")
mam <- lmer(vol.model$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

 for (x in 1:L)
  {
    log_access <- acc[x]
    vol <-0
    ag_suit<- median(model.data$ag_suit)
    log_hpd <- mean(model.data$log_hpd)
    log_AREA.PA <- mean(model.data$log_AREA.PA)
    Within_PA <- "no"
    Zone <- "Tropical"
    newdat.n<-data.frame(cbind(vol, log_access, ag_suit, log_AREA.PA, log_hpd))
    newdat.f<-data.frame(cbind(Within_PA, Zone))
    levels(newdat.f$Within_PA) <- levels(model.data$Within_PA)
    levels(newdat.f$Zone) <- levels(model.data$Zone)
    newdat<-data.frame(newdat.f,newdat.n)
    mm<-model.matrix(terms(mam),newdat)
    pvar1 <- diag(mm %*% tcrossprod(vcov(mam),mm))
    z[x]<-mm %*% fixef(mam)
    zu[x]<-z[x]+sqrt(pvar1)
    zl[x]<-z[x]-sqrt(pvar1) 
  }



  acc <- exp(acc)/60

  z  <- 10^(z) 
  zu <- 10^(zu)
  zl <- 10^(zl)


  plot(acc,z, ylim=c(0,7000), col = 1,
		bty = "l", log = "x",
		type = "l",ylab = "CWM Length derived vol(mm3) ± s.e", xlab="Accessibility (hours to city >50 000)")
  points(acc,zu,type="l",lty=2, col = 1)
  points(acc,zl,type="l",lty=2, col = 1)

  rug(invert.vol.data$access/60
	, ticksize = 0.03, side = 1,  lwd = 0.5, col = 1) 


dev.off()








#### invert vol vs hpd ##########


model.data <- vol.model$data

tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/invert vol vs hpd.tif",
	width = 15, height = 10, units = "cm", pointsize = 12, res = 300)

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
mam <- lmer(vol.model$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

 for (x in 1:L)
  {
    log_access <-mean(model.data$log_access)
    vol <-0
    ag_suit<- median(model.data$ag_suit)
    log_hpd <- hpd[x]
    log_AREA.PA <- mean(model.data$log_AREA.PA)
    Within_PA <- "no"
    Zone <- "Tropical"
    newdat.n<-data.frame(cbind(vol, log_access, ag_suit, log_AREA.PA, log_hpd))
    newdat.f<-data.frame(cbind(Within_PA, Zone))
    levels(newdat.f$Within_PA) <- levels(model.data$Within_PA)
    levels(newdat.f$Zone) <- levels(model.data$Zone)
    newdat<-data.frame(newdat.f,newdat.n)
    mm<-model.matrix(terms(mam),newdat)
    pvar1 <- diag(mm %*% tcrossprod(vcov(mam),mm))
    z[x]<-mm %*% fixef(mam)
    zu[x]<-z[x]+sqrt(pvar1)
    zl[x]<-z[x]-sqrt(pvar1) 
  }



  hpd <- exp(hpd)-1

  z  <- 10^(z) 
  zu <- 10^(zu)
  zl <- 10^(zl)


  plot(hpd,z, ylim=c(0,7000), col = 1,
		bty = "l", log = "x",
		type = "l",ylab = "CWM Length derived vol(mm3) ± s.e", xlab="Human population density (per km2)")
  points(hpd,zu,type="l",lty=2, col = 1)
  points(hpd,zl,type="l",lty=2, col = 1)

  rug(invert.vol.data$hpd
	, ticksize = 0.03, side = 1,  lwd = 0.5, col = 1) 


dev.off()


















#### invert vol vs PA size vs zone ##########


model.data <- vol.model$data

tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/invert vol vs p asize vs zone.tif",
	width = 15, height = 10, units = "cm", pointsize = 12, res = 300)

par(mfrow = c(1,2))


  par(mar=c(4,4.5,4,1.5))
  par(mgp=c(2.5,1,0))
  
 L = 30

  gis <-seq(from=min(model.data$log_AREA.PA),to=max(model.data$log_AREA.PA),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  


# it doesnt matter what the values for taxon and zone are, as there is no interaction with these




model.data$Zone<- relevel(model.data$Zone , "Tropical")
mam <- lmer(vol.model$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

 for (x in 1:L)
  {
    log_access <-mean(model.data$log_access)
    vol <-0
    ag_suit<- median(model.data$ag_suit)
    log_hpd <- mean(model.data$log_hpd)
    log_AREA.PA <- gis[x]
    Within_PA <- "no"
    Zone <- "Tropical"
    newdat.n<-data.frame(cbind(vol, log_access, ag_suit, log_AREA.PA, log_hpd))
    newdat.f<-data.frame(cbind(Within_PA, Zone))
    levels(newdat.f$Within_PA) <- levels(model.data$Within_PA)
    levels(newdat.f$Zone) <- levels(model.data$Zone)
    newdat<-data.frame(newdat.f,newdat.n)
    mm<-model.matrix(terms(mam),newdat)
    pvar1 <- diag(mm %*% tcrossprod(vcov(mam),mm))
    z[x]<-mm %*% fixef(mam)
    zu[x]<-z[x]+sqrt(pvar1)
    zl[x]<-z[x]-sqrt(pvar1) 
  }



  gis <- exp(gis)-1

  z  <- 10^(z) 
  zu <- 10^(zu)
  zl <- 10^(zl)


  plot(gis,z, ylim=c(0,7000), col = 1,
		bty = "l", log = "x", main = "Tropical",
		type = "l",ylab = "CWM Length derived vol(mm3) ± s.e", xlab="PA Size (km2)")
  points(gis,zu,type="l",lty=2, col = 1)
  points(gis,zl,type="l",lty=2, col = 1)

  rug(invert.vol.data$AREA.PA
	, ticksize = 0.03, side = 1,  lwd = 0.5, col = 1) 


  gis <-seq(from=min(model.data$log_AREA.PA),to=max(model.data$log_AREA.PA),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  


model.data$Zone<- relevel(model.data$Zone , "Temperate")
mam <- lmer(vol.model$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

 for (x in 1:L)
  {
    log_access <-mean(model.data$log_access)
    vol <-0
    ag_suit<- median(model.data$ag_suit)
    log_hpd <- mean(model.data$log_hpd)
    log_AREA.PA <- gis[x]
    Within_PA <- "no"
    Zone <- "Tropical"
    newdat.n<-data.frame(cbind(vol, log_access, ag_suit, log_AREA.PA, log_hpd))
    newdat.f<-data.frame(cbind(Within_PA, Zone))
    levels(newdat.f$Within_PA) <- levels(model.data$Within_PA)
    levels(newdat.f$Zone) <- levels(model.data$Zone)
    newdat<-data.frame(newdat.f,newdat.n)
    mm<-model.matrix(terms(mam),newdat)
    pvar1 <- diag(mm %*% tcrossprod(vcov(mam),mm))
    z[x]<-mm %*% fixef(mam)
    zu[x]<-z[x]+sqrt(pvar1)
    zl[x]<-z[x]-sqrt(pvar1) 
  }



  gis <- exp(gis)-1


  z  <- 10^(z) 
  zu <- 10^(zu)
  zl <- 10^(zl)



  plot(gis,z, ylim=c(0,7000), col = 1,
		bty = "l", log = "x", main = "Temperate",
		type = "l",ylab = "CWM Length derived vol(mm3) ± s.e", xlab="PA Size (km2)")
  points(gis,zu,type="l",lty=2, col = 1)
  points(gis,zl,type="l",lty=2, col = 1)

  rug(invert.vol.data$AREA.PA
	, ticksize = 0.03, side = 1,  lwd = 0.5, col = 1) 



dev.off()

















####
#verts vol
####





#### ag suit ##########


model.data <- vol.model$data

tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/amphib vol vs ag suit.tif",
	width = 15, height = 10, units = "cm", pointsize = 12, res = 300)

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
model.data$Zone<- relevel(model.data$Zone , "Tropical")
mam <- lmer(vol.model$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

 for (x in 1:L)
  {
    ag_suit <- ag[x]
    vol <-0
    log_bound_dist_km_PA_neg <- mean(model.data$log_bound_dist_km_PA_neg)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(vol, log_bound_dist_km_PA_neg, ag_suit))
    newdat.f<-data.frame(cbind(Within_PA))
    levels(newdat.f$Within_PA) <- levels(model.data$Within_PA)
    newdat<-data.frame(newdat.f,newdat.n)
    mm<-model.matrix(terms(mam),newdat)
    pvar1 <- diag(mm %*% tcrossprod(vcov(mam),mm))
    z[x]<-mm %*% fixef(mam)
    zu[x]<-z[x]+sqrt(pvar1)
    zl[x]<-z[x]-sqrt(pvar1) 
  }
  


  z  <- 10^(z) 
  zu <- 10^(zu)
  zl <- 10^(zl)


  plot(ag,z/1000, ylim=c(0,1700), col = 1,
		bty = "l", #log = "x",
		type = "l",ylab = "CWM Length derived vol(cm3) ± s.e", xlab="Agricultural suitability")
  points(ag,zu/1000,type="l",lty=2, col = 1)
  points(ag,zl/1000,type="l",lty=2, col = 1)





dev.off()







#### lbd ##########


model.data <- vol.model$data

tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/amphib vol vs log bound dist.tif",
	width = 15, height = 10, units = "cm", pointsize = 12, res = 300)

par(mfrow = c(1,1))


  par(mar=c(4,4.5,4,1.5))
  par(mgp=c(2.5,1,0))
  
 L = 30

  lbd <-seq(from=min(model.data$log_bound_dist_km_PA_neg),to=max(model.data$log_bound_dist_km_PA_neg),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  


# it doesnt matter what the values for taxon and zone are, as there is no interaction with these



model.data$Within_PA <- relevel(model.data$Within_PA, "no")
model.data$Zone<- relevel(model.data$Zone , "Tropical")
mam <- lmer(vol.model$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

 for (x in 1:L)
  {
    log_bound_dist_km_PA_neg <- lbd[x]
    vol <-0
    ag_suit<- median(model.data$ag_suit)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(vol, log_bound_dist_km_PA_neg, ag_suit))
    newdat.f<-data.frame(cbind(Within_PA))
    levels(newdat.f$Within_PA) <- levels(model.data$Within_PA)
    newdat<-data.frame(newdat.f,newdat.n)
    mm<-model.matrix(terms(mam),newdat)
    pvar1 <- diag(mm %*% tcrossprod(vcov(mam),mm))
    z[x]<-mm %*% fixef(mam)
    zu[x]<-z[x]+sqrt(pvar1)
    zl[x]<-z[x]-sqrt(pvar1) 
  }
  

  lbd <- exp(lbd)

  z  <- 10^(z) 
  zu <- 10^(zu)
  zl <- 10^(zl)


  plot(lbd,z/1000, ylim=c(0,600), col = 1,
		bty = "l", #log = "x",
		type = "l",ylab = "CWM Length derived vol(cm3) ± s.e", xlab="Distance to protected area boundary (km)")
  points(lbd,zu/1000,type="l",lty=2, col = 1)
  points(lbd,zl/1000,type="l",lty=2, col = 1)





dev.off()





