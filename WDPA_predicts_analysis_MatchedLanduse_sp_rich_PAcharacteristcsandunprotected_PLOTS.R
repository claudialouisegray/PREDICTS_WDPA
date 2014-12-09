# PLOTS



model.data <- Species_richness.model$data

data <- multiple.taxa.matched.landuse[,c("Species_richness", "Zone", "taxon_of_interest", "ag_suit", "log_access", "access", "DoP.PA", "AREA.PA",
	 "log_hpd", "hpd", "log_bound_dist_km_PA_neg", "bound_dist_km_PA_neg", "Within_PA", "Predominant_habitat", "SS", "SSBS")]

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



display.brewer.all()
cols <- brewer.pal(8, "Paired")
display.brewer.pal(8, "Paired")
#drop the red to avoid red-green colorblindness issues

cols <- cols[c(2,4,8)]

#green one
#inside.col <- rgb(0.1,0.8,0.2)



inside.col <- 1
outside.col <- rgb(0.7,0.7,0.7)


lu <- c("Primary Vegetation", "Secondary Vegetation", "Plantation forest", "Cropland", "Pasture", "Urban")

lu.cols = c("#5B8A3B", "#1B9E77", "#7570B3", "#E6AB02", "#D95F02", "#E7298A")




###### dist to boundary  #############
###### for different taxa ##################



tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/11_14/Species_richness vs taxon vs dist to boundary.tif",
	width = 10, height = 25, units = "cm", pointsize = 12, res = 300)

par(mfrow = c(3,1))
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
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, Species_richness, log_hpd,log_AREA.PA, DoP.PA, log_access, ag_suit))
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

 lbd1 <- exp(abs(lbd)) -1
 inside <- which(lbd < 0)
 lbd1[inside] <- lbd1[inside]*-1

  plot(lbd1,z, ylim=c(0,50), xlim = c(-80, 200), col = cols[i],
		bty = "l", #log = "x", #yaxt = "n", 
		type = "l",ylab = "Species richness per site ± s.e", xlab="Distance to PA boundary (km)", main = t)
  rug(data$bound_dist_km_PA_neg[which(data$taxon_of_interest == t)]
	, ticksize = 0.03, side = 1, lwd = 0.5, pos = 0, col = cols[i])
  points(lbd1,zu,type="l",lty=2, col = cols[i])
  points(lbd1,zl,type="l",lty=2, col = cols[i])

 abline(v = 0, lty = 2, col = 8)

}



dev.off()

















#### within PA and ag suit ##########



tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/11_14/Species_richness within pa vs  ag_suit.tif",
	width = 12, height = 15, units = "cm", pointsize = 12, res = 300)

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
model.data$taxon_of_interest <- relevel(model.data$taxon_of_interest , "Invertebrates")
mam <- glmer(Species_richness.model$final.call, model.data, family = "poisson")

 for (x in 1:L)
  {
    ag_suit <- ag[x]
    log_bound_dist_km_PA_neg <-  mean(model.data$log_bound_dist_km_PA_neg)
    Species_richness <-0
    log_hpd <- mean(model.data$log_hpd)
    log_AREA.PA <-mean(model.data$log_AREA.PA)
    DoP.PA <- mean(model.data$DoP.PA)
    IUCN.PA <- "1.5"
    Zone <-"Tropical"
    taxon_of_interest<- "Vertebrates"
    Predominant_habitat <-"Primary Vegetation"
    log_access <- mean(model.data$log_access)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, Species_richness, log_hpd,log_AREA.PA, DoP.PA, log_access, ag_suit))
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


  plot(ag,z, ylim=c(5, 30), col = outside.col,
		bty = "l", #log = "x",
		type = "l",ylab = "Species richness per site ± s.e", 
		xlab="Agricultural suitability (higher = more suitable)")
  points(ag,zu,type="l",lty=2, col = outside.col)
  points(ag,zl,type="l",lty=2, col = outside.col)






model.data$Within_PA <- relevel(model.data$Within_PA, "yes")

mam<- glmer(Species_richness.model$final.call, model.data,  family = "poisson", control= glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

  ag <-seq(from=min(model.data$ag_suit[which(model.data$Within_PA == "yes")]),
		to=max(model.data$ag_suit[which(model.data$Within_PA == "yes")]),length=L)
  y <-vector(mode="numeric",length=length(hpd))
  yu <-vector(mode="numeric",length=length(hpd))
  yl <-vector(mode="numeric",length=length(hpd))


 for (x in 1:L)
  {
    ag_suit <- ag[x]
    log_bound_dist_km_PA_neg <-  mean(model.data$log_bound_dist_km_PA_neg)
    Species_richness <-0
    log_hpd <- mean(model.data$log_hpd)
    log_AREA.PA <-mean(model.data$log_AREA.PA)
    DoP.PA <- mean(model.data$DoP.PA)
    IUCN.PA <- "1.5"
    Zone <-"Tropical"
    taxon_of_interest<- "Vertebrates"
    Predominant_habitat <-"Primary Vegetation"
    log_access <- mean(model.data$log_access)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, Species_richness, log_hpd,log_AREA.PA, DoP.PA, log_access, ag_suit))
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

  
  points(ag,y,type="l",lty=1 , col = inside.col, lwd = 2)
  points(ag,yu,type="l",lty=3, lwd = 2, col = inside.col)
  points(ag,yl,type="l",lty=3, lwd = 2, col = inside.col)



# add on distribution of ag_suit values for landuses
ph <- coef(Species_richness.model$model)$Predominant_habitat[,1]
exp(ph)

#from data 

lu.means <- aggregate(ag_suit ~ Predominant_habitat, model.data, mean)
lu.sd <- aggregate(ag_suit ~ Predominant_habitat, model.data, sd)
lu.n <- aggregate(ag_suit ~ Predominant_habitat, model.data, length)

lu.se <- data.frame(Predominant_habitat = lu.sd$Predominant_habitat, 
		ag_suit = lu.sd$ag_suit/sqrt(lu.n$ag_suit))

for(y in 1:6){
	ylim.min <- 5
	xpoint <- lu.means$ag_suit[which(lu.means$Predominant_habitat == lu[y])]
	xse <- lu.sd$ag_suit[which(lu.sd$Predominant_habitat == lu[y])]
	points(y = ylim.min + y, x = xpoint, pch = 16, col = lu.cols[y])
	xmax <- xpoint + xse
	xmin <- xpoint - xse
	arrows(xmax,ylim.min+y,xmin,ylim.min+y, 
	cex = 0.7, code = 3, angle = 90, length = 0.04, col = lu.cols[y])
}


legend("topright", c("Protected", "Unprotected") , cex = 0.7,
	col = c(inside.col, outside.col), lty = c(1,1), lwd = c(1,2))
legend("topleft", rev(lu), cex = 0.7,
	col = rev(lu.cols), pch = 16)


dev.off()






#### accessibility ##########



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
model.data$taxon_of_interest <- relevel(model.data$taxon_of_interest , "Invertebrates")
mam <- glmer(Species_richness.model$final.call, model.data, family = "poisson")

 for (x in 1:L)
  {
    ag_suit <- median(model.data$ag_suit)
    log_bound_dist_km_PA_neg <-  mean(model.data$log_bound_dist_km_PA_neg)
    Species_richness <-0
    log_access <- acc[x]
    log_AREA.PA <-mean(model.data$log_AREA.PA)
    DoP.PA <- mean(model.data$DoP.PA)
    IUCN.PA <- "1.5"
    Zone <-"Tropical"
    taxon_of_interest<- t
    Predominant_habitat <-"Primary Vegetation"
    log_hpd <- mean(model.data$log_hpd)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, Species_richness, log_hpd,log_AREA.PA, DoP.PA, log_access, ag_suit))
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
  

  acc_plot <- (exp(acc) -1)/60

  #backtransform from poisson log-link
  z <- exp(z)
  zu <- exp(zu)
  zl <- exp(zl)


  plot(acc_plot,z, ylim=c(0, 30), col = 1,
		bty = "l", log = "x",
		type = "l",ylab = "Species richness per site ± s.e", 
		xlab="Accessibility (hours to city >50000)")
  points(acc_plot,zu,type="l",lty=2, col = 1)
  points(acc_plot,zl,type="l",lty=2, col = 1)
   rug((data$access)/60
	, ticksize = 0.03, side = 1,   lwd = 0.2, col = 1) 


dev.off()







#### within PA and hpd ##########



tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/11_14/Species_richness within pa vs hpd.tif",
	width = 12, height = 10, units = "cm", pointsize = 12, res = 300)

par(mfrow = c(1,1))


  par(mar=c(4,4.5,4,1.5))
  par(mgp=c(2.5,1,0))
  
 L = 30

  hpd <-seq(from=min(model.data$log_hpd[which(model.data$Within_PA == "no")]),
		to=max(model.data$log_hpd[which(model.data$Within_PA == "no")]),length=L)
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
    log_hpd <- hpd[x]
    log_AREA.PA <-mean(model.data$log_AREA.PA)
    DoP.PA <- mean(model.data$DoP.PA)
    IUCN.PA <- "1.5"
    Zone <-"Tropical"
    taxon_of_interest<- "Vertebrates"
    Predominant_habitat <-"Primary Vegetation"
    log_access <- mean(model.data$log_access)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, Species_richness, log_hpd,log_AREA.PA, DoP.PA, log_access, ag_suit))
    newdat.f<-data.frame(cbind( Zone,taxon_of_interest, Predominant_habitat, Within_PA))
    levels(newdat.f$Zone)<-levels(model.data$Zone)
    levels(newdat.f$taxon_of_interest)<-levels(model.data$taxon_of_interest)
    levels(newdat.f$Predominant_habitat)<-levels(model.data$Predominant_habitat)
    levels(newdat.f$hpd) <- levels(model.data$hpd)
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


  plot(hpd,z, ylim=c(5, 30), col = outside.col,
		bty = "l", log = "x",
		type = "l",ylab = "Species richness per site ± s.e", 
		xlab="Human population density (per km2)")
  points(hpd,zu,type="l",lty=2, col = outside.col)
  points(hpd,zl,type="l",lty=2, col = outside.col)
  rug(data$hpd[which(data$Within_PA == "no")],
	ticksize = 0.03, side = 1,   lwd = 1, col = outside.col) 





model.data$Within_PA <- relevel(model.data$Within_PA, "yes")

mam<- glmer(Species_richness.model$final.call, model.data,  family = "poisson", control= glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

  hpd <-seq(from=min(model.data$log_hpd[which(model.data$Within_PA == "yes")]),
		to=max(model.data$log_hpd[which(model.data$Within_PA == "yes")]),length=L)
  y <-vector(mode="numeric",length=length(hpd))
  yu <-vector(mode="numeric",length=length(hpd))
  yl <-vector(mode="numeric",length=length(hpd))


 for (x in 1:L)
  {
    ag_suit <- median(model.data$ag_suit)
    log_bound_dist_km_PA_neg <-  mean(model.data$log_bound_dist_km_PA_neg)
    Species_richness <-0
    log_hpd <- hpd[x]
    log_AREA.PA <-mean(model.data$log_AREA.PA)
    DoP.PA <- mean(model.data$DoP.PA)
    IUCN.PA <- "1.5"
    Zone <-"Tropical"
    taxon_of_interest<- "Vertebrates"
    Predominant_habitat <-"Primary Vegetation"
    log_access <- mean(model.data$log_access)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, Species_richness, log_hpd,log_AREA.PA, DoP.PA, log_access, ag_suit))
    newdat.f<-data.frame(cbind( Zone,taxon_of_interest, Predominant_habitat, Within_PA))
    levels(newdat.f$Zone)<-levels(model.data$Zone)
    levels(newdat.f$taxon_of_interest)<-levels(model.data$taxon_of_interest)
    levels(newdat.f$Predominant_habitat)<-levels(model.data$Predominant_habitat)
    levels(newdat.f$hpd) <- levels(model.data$hpd)
    levels(newdat.f$Within_PA) <- levels(model.data$Within_PA)
    newdat<-data.frame(newdat.f,newdat.n)
    mm<-model.matrix(terms(mam),newdat)
    pvar1 <- diag(mm %*% tcrossprod(vcov(mam),mm))

    y[x]<-mm %*% fixef(mam)
    yu[x]<-y[x]+sqrt(pvar1)
    yl[x]<-y[x]-sqrt(pvar1) 
  }


  hpd <- exp(hpd) -1


  #backtransform from poisson log-link
  y <- exp(y)
  yu <- exp(yu)
  yl <- exp(yl)

  
  points(hpd,y,type="l",lty=1 , col = inside.col, lwd = 2)
  points(hpd,yu,type="l",lty=3, lwd = 2, col = inside.col)
  points(hpd,yl,type="l",lty=3, lwd = 2, col = inside.col)
  rug(data$hpd[which(data$Within_PA == "yes")],
	ticksize = 0.03, side = 1,   lwd = 1, col = inside.col, pos = 5) 

legend("topright", c("Protected", "Unprotected") , col = c(inside.col, outside.col), lty = c(1,1), lwd = c(1,2))


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
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, Species_richness, log_hpd,log_AREA.PA, DoP.PA, log_access, ag_suit))
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
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, Species_richness, log_hpd,log_AREA.PA, DoP.PA, log_access, ag_suit))
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
model.data$taxon_of_interest <- relevel(model.data$taxon_of_interest , "Invertebrates")
mam <- glmer(Species_richness.model$final.call, model.data, family = "poisson")

 for (x in 1:L)
  {
    ag_suit <- median(model.data$ag_suit)
    log_bound_dist_km_PA_neg <-  mean(model.data$log_bound_dist_km_PA_neg)
    Species_richness <-0
    log_hpd <- hpd[x]
    log_AREA.PA <-mean(model.data$log_AREA.PA)
    DoP.PA <- mean(model.data$DoP.PA)
    IUCN.PA <- "1.5"
    Zone <-"Tropical"
    taxon_of_interest<- t
    Predominant_habitat <-"Primary Vegetation"
    log_access <- mean(model.data$log_access)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, Species_richness, log_hpd,log_AREA.PA, DoP.PA, log_access, ag_suit))
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
  

  hpd <- exp(hpd) -1

  #backtransform from poisson log-link
  z <- exp(z)
  zu <- exp(zu)
  zl <- exp(zl)


  plot(hpd,z, ylim=c(0, 30), col = 1,
		bty = "l", log = "x",
		type = "l",ylab = "Species richness per site ± s.e", 
		xlab="Human population density (per km2)")
  points(hpd,zu,type="l",lty=2, col = 1)
  points(hpd,zl,type="l",lty=2, col = 1)
  rug(data$hpd,
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
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, Species_richness, log_hpd,log_AREA.PA, DoP.PA, log_access, ag_suit))
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
    log_hpd <- mean(model.data$log_hpd)
    log_AREA.PA <- mean(model.data$log_AREA.PA )
    DoP.PA <- mean(model.data$DoP.PA )
    IUCN.PA <- "1.5"
    Zone <-"Temperate"
    taxon_of_interest<- "Vertebrates"
    Predominant_habitat <-"Primary Vegetation"
    ag_suit <- median(model.data$ag_suit)
    log_access <- mean(model.data$log_access)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, Species_richness, log_hpd,log_AREA.PA,  DoP.PA, log_access, ag_suit))
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







