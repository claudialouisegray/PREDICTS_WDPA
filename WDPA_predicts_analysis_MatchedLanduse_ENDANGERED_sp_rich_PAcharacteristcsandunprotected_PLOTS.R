# PLOTS




model.data <- Species_richness.model$data

data <- multiple.taxa.matched.landuse_VU_EN_CR[,c("Species_richness", "Zone", "ag_suit", "access", "hpd", "SS")]


data <- na.omit(data) #have to get rid of NA to try poly on ag_suit

names(matched.landuse)

nrow(data)
length(unique(data$SS))



nrow(model.data)
length(unique(model.data$SS))




display.brewer.all()
cols <- brewer.pal(8, "Paired")
display.brewer.pal(8, "Paired")
#drop the red to avoid red-green colorblindness issues

zone.cols <- cols[c(6,5)]


inside.col <- rgb(0.1,0.8,0.2)

display.brewer.all()
cols <- brewer.pal(8, "Paired")
display.brewer.pal(8, "Paired")
#drop the red to avoid red-green colorblindness issues

cols <- cols[c(2,4,8)]




#### ag suit ##########

inside.col <- rgb(0.3,0.7,0.3)

tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/11_14/VU_EN_CR Species_richness within pa vs  ag_suit.tif",
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
  


# it doesnt matter what the values for taxon and zone are, as there is no interaction with these



model.data$Within_PA <- relevel(model.data$Within_PA, "no")
model.data$Zone<- relevel(model.data$Zone , "Tropical")
mam <- glmer(Species_richness.model$final.call, model.data, family = "poisson")

 for (x in 1:L)
  {
    ag_suit <- ag[x]
    Species_richness <-0
    log_hpd <- mean(model.data$log_hpd)
    IUCN.PA <- "1.5"
    Zone <-"Tropical"
    Predominant_habitat <-"Primary Vegetation"
    log_access <- mean(model.data$log_access)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(Species_richness, log_hpd, log_access, ag_suit))
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
  
  #backtransform from poisson log-link
  z <- exp(z)
  zu <- exp(zu)
  zl <- exp(zl)


  plot(ag,z, ylim=c(0,1), col = 1,
		bty = "l", #log = "x",
		type = "l",ylab = "Species richness per site ± s.e", 
		xlab="Agricultural suitability (higher = more suitable)")
  points(ag,zu,type="l",lty=2, 1)
  points(ag,zl,type="l",lty=2,1)




dev.off()







#### hpd ##########



tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/Species_richness vs hpd.tif",
	width = 12, height = 10, units = "cm", pointsize = 12, res = 300)



par(mfrow = c(1,1))


  par(mar=c(4,4.5,4,1.5))
  par(mgp=c(2.5,1,0))
  
 L = 30

  hpd <-seq(from=min(model.data$log_hpd),
		to=max(model.data$log_hpd),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  


# it doesnt matter what the values for taxon and zone are, as there is no interaction with these



model.data$Within_PA <- relevel(model.data$Within_PA, "no")
model.data$Zone<- relevel(model.data$Zone , "Tropical")
#mam <- glmer(Species_richness.model$final.call, model.data, family = "poisson")

 for (x in 1:L)
  {
    ag_suit <- median(model.data$ag_suit)
    Species_richness <-0
    log_hpd <- hpd[x]
    IUCN.PA <- "1.5"
    Zone <-"Tropical"
    Predominant_habitat <-"Primary Vegetation"
    log_access <- mean(model.data$log_access)
    newdat.n<-data.frame(cbind(Species_richness, log_hpd, log_access, ag_suit))
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
  

  hpd <- exp(hpd) -1

  #backtransform from poisson log-link
  z <- exp(z)
  zu <- exp(zu)
  zl <- exp(zl)


  plot(hpd,z, ylim=c(0, 1), col = 1,
		bty = "l", log = "x",
		type = "l",ylab = "Species richness per site ± s.e", 
		xlab="Human population density (per km2)")
  points(hpd,zu,type="l",lty=2, col = 1)
  points(hpd,zl,type="l",lty=2, col = 1)
  rug(data$hpd,
		col = 1, lwd = 0.8, pos = 0)


dev.off()




#### acc ##########



tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/11_14/Species_richness vs access.tif",
	width = 12, height = 10, units = "cm", pointsize = 12, res = 300)



par(mfrow = c(1,1))


  par(mar=c(4,4.5,4,1.5))
  par(mgp=c(2.5,1,0))
  
 L = 30

  acc <-seq(from=min(model.data$log_access),
		to=max(model.data$log_access),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  


# it doesnt matter what the values for taxon and zone are, as there is no interaction with these



model.data$Within_PA <- relevel(model.data$Within_PA, "no")
model.data$Zone<- relevel(model.data$Zone , "Tropical")
#mam <- glmer(Species_richness.model$final.call, model.data, family = "poisson")

 for (x in 1:L)
  {
    ag_suit <- median(model.data$ag_suit)
    Species_richness <-0
    log_access <- acc[x]
    IUCN.PA <- "1.5"
    Zone <-"Tropical"
    Predominant_habitat <-"Primary Vegetation"
    log_hpd <- mean(model.data$log_hpd)
    newdat.n<-data.frame(cbind(Species_richness, log_access, log_hpd, ag_suit))
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
  

  acc <- exp(acc) -1

  #backtransform from poisson log-link
  z <- exp(z)
  zu <- exp(zu)
  zl <- exp(zl)


  plot(acc,z, ylim=c(0, 1), col = 1,
		bty = "l", log = "x",
		type = "l",ylab = "Species richness per site ± s.e", 
		xlab="Human population density (per km2)")
  points(acc,zu,type="l",lty=2, col = 1)
  points(acc,zl,type="l",lty=2, col = 1)
  rug(data$acc,
		col = 1, lwd = 0.8, pos = 0)


dev.off()









