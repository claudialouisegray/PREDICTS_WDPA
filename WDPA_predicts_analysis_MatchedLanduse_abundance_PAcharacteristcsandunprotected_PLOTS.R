# PLOTS


model.data <- log_abundance.model$data



data <- matched.landuse[,c("log_abundance", "Zone", "taxon_of_interest", "ag_suit", "log_access", "access", 
	"log_AREA.PA", "AREA.PA", "DoP.PA", "IUCN.PA",
	 "log_hpd", "hpd", "log_bound_dist_km_PA_neg", "bound_dist_km_PA_neg", "Within_PA", "Predominant_habitat", "SS", "SSBS")]

data <- na.omit(data) #have to get rid of NA to try poly on ag_suit

nrow(data)
length(unique(data$SS))

nrow(model.data)
length(unique(model.data$SS))




display.brewer.all()
cols <- brewer.pal(8, "Paired")
display.brewer.pal(8, "Paired")
#drop the red to avoid red-green colorblindness issues

taxa.cols <- cols[c(2,4,8)]



cols <- brewer.pal(8, "Paired")
#drop the red to avoid red-green colorblindness issues

zone.cols <- cols[c(6,2)]





###### dist to boundary vs taxon  #############



tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/11_14/abundance vs dist to boundary vs taxon.tif",
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

  lbd <-seq(from=min(model.data$log_bound_dist_km_PA_neg[which(model.data$taxon_of_interest == t)]),
	to=max(model.data$log_bound_dist_km_PA_neg),length=L)
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
    log_hpd <- mean(model.data$log_hpd)
    log_access<- mean(model.data$log_access)
    log_AREA.PA <- mean(model.data$log_AREA.PA)
    IUCN.PA <- "1.5"
    DoP.PA <- mean(model.data$DoP.PA)
    Zone <-"Tropical"
    taxon_of_interest<- "Invertebrates"
    Predominant_habitat <-"Primary Vegetation"
    ag_suit <- median(model.data$ag_suit)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, log_abundance, log_AREA.PA, DoP.PA, log_access, log_hpd, ag_suit))
    newdat.f<-data.frame(cbind(Zone,taxon_of_interest, Predominant_habitat, Within_PA))
    levels(newdat.f$Zone)<-levels(model.data$Zone)
    levels(newdat.f$taxon_of_interest)<-levels(model.data$taxon_of_interest)
    levels(newdat.f$Predominant_habitat)<-levels(model.data$Predominant_habitat)
    levels(newdat.f$ag_suit) <- levels(model.data$ag_suit)
    levels(newdat.f$Within_PA) <- levels(model.data$Within_PA)
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



  plot(lbd1,z, ylim=c(0, 600), xlim = c(-80,200), col = taxa.cols[i], main = t,
		bty = "l", #yaxt = "n", 
		type = "l",ylab = "Abundance per site ± s.e", xlab="Distance to PA boundary (km)")
  rug(data$bound_dist_km_PA_neg[which(data$taxon_of_interest == t)]
	, ticksize = 0.03, side = 1, lwd = 0.5, col = taxa.cols[i]) 
  axis(2,at = c(100,200,500,1000), c(100,200,500,1000))
  points(lbd1,zu,type="l",lty=2, col = taxa.cols[i])
  points(lbd1,zl,type="l",lty=2, col = taxa.cols[i])




#  plot(lbd,z, ylim=c(0, 10), col = 8,
#		bty = "l", #log = "x",
#		type = "l",ylab = "Abundance per site ± s.e", xlab="(Log) Distance to PA boundary (km)")
 # rug(model.data$log_bound_dist_km_PA_neg_PA_neg, ticksize = 0.03, side = 1, lwd = 0.5) #### change to be data used in model, poss slightly fewer points 
 # points(lbd,zu,type="l",lty=2, col =8)
 # points(lbd,zl,type="l",lty=2, col = 8)

 abline(v = 0, lty = 2, col = 8)

}





dev.off()





#### PA age #####


tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/11_14/abundance vs DoP.tif",
	width = 15, height = 10, units = "cm", pointsize = 12, res = 300)




  par(mar=c(4,4.5,4,1.5))
  par(mgp=c(2.5,1,0))
  
 L = 30

  dop <-seq(from=min(model.data$DoP.PA),
		to=max(model.data$DoP.PA),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  



model.data$Zone<- relevel(model.data$Zone , "Tropical")
mam <- lmer(log_abundance.model$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

 for (x in 1:L)
  {
    log_bound_dist_km_PA_neg <-  mean(model.data$log_bound_dist_km_PA_neg)
    log_abundance <-0
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
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, log_abundance, log_hpd,log_AREA.PA, DoP.PA, log_access, ag_suit))
    newdat.f<-data.frame(cbind( Zone,taxon_of_interest, Predominant_habitat,  Within_PA))
    levels(newdat.f$Zone)<-levels(model.data$Zone)
    levels(newdat.f$taxon_of_interest)<-levels(model.data$taxon_of_interest)
    levels(newdat.f$Predominant_habitat)<-levels(model.data$Predominant_habitat)
    levels(newdat.f$ag_suit) <- levels(model.data$ag_suit)
    levels(newdat.f$Within_PA) <- levels(model.data$Within_PA)
  #  levels(newdat.f$IUCN.PA) <- levels(model.data$IUCN.PA)
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

  plot(dop,z, ylim=c(0,150), col = 1, lwd = 1,
		bty = "l", #log = "x", #yaxt = "n", 
		type = "l",ylab = "Abundance per site ± s.e", xlab="Duration of protection (yr)")
  rug(data$DoP.PA
	, ticksize = 0.03, side = 1, lwd = 1, col = 1)
  points(dop,zu,type="l",lty=2,  lwd = 1, col =  1)
  points(dop,zl,type="l",lty=2,  lwd = 1, col =  1)





dev.off()





#### PA age, for different zones #####






tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/11_14/abundance vs zone vs DoP.tif",
	width = 15, height = 10, units = "cm", pointsize = 12, res = 300)




  par(mar=c(4,4.5,4,1.5))
  par(mgp=c(2.5,1,0))
  
 L = 30

  dop <-seq(from=min(model.data$DoP.PA[which(model.data$Zone == "Tropical")]),
		to=max(model.data$DoP.PA[which(model.data$Zone == "Tropical")]),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  



#model.data$IUCN.PA <- relevel(model.data$IUCN.PA, "1.5")
model.data$Zone<- relevel(model.data$Zone , "Tropical")
mam <- lmer(log_abundance.model$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

 for (x in 1:L)
  {
    log_bound_dist_km_PA_neg <-  mean(model.data$log_bound_dist_km_PA_neg)
    log_abundance <-0
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
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, log_abundance, log_hpd,log_AREA.PA, DoP.PA, log_access, ag_suit))
    newdat.f<-data.frame(cbind( Zone,taxon_of_interest, Predominant_habitat,  Within_PA))
    levels(newdat.f$Zone)<-levels(model.data$Zone)
    levels(newdat.f$taxon_of_interest)<-levels(model.data$taxon_of_interest)
    levels(newdat.f$Predominant_habitat)<-levels(model.data$Predominant_habitat)
    levels(newdat.f$ag_suit) <- levels(model.data$ag_suit)
    levels(newdat.f$Within_PA) <- levels(model.data$Within_PA)
  #  levels(newdat.f$IUCN.PA) <- levels(model.data$IUCN.PA)
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

  plot(dop,z, ylim=c(0,150), col = zone.cols[1], lwd = 2,
		bty = "l", #log = "x", #yaxt = "n", 
		type = "l",ylab = "Abundance per site ± s.e", xlab="Duration of protection (yr)")
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
mam <- lmer(log_abundance.model$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

 for (x in 1:L)
  {
    log_bound_dist_km_PA_neg <-  mean(model.data$log_bound_dist_km_PA_neg)
    log_abundance <-0
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
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, log_abundance, log_hpd,log_AREA.PA, DoP.PA, log_access, ag_suit))
    newdat.f<-data.frame(cbind( Zone,taxon_of_interest, Predominant_habitat, Within_PA))
    levels(newdat.f$Zone)<-levels(model.data$Zone)
    levels(newdat.f$taxon_of_interest)<-levels(model.data$taxon_of_interest)
    levels(newdat.f$Predominant_habitat)<-levels(model.data$Predominant_habitat)
    levels(newdat.f$ag_suit) <- levels(model.data$ag_suit)
    levels(newdat.f$Within_PA) <- levels(model.data$Within_PA)
   # levels(newdat.f$IUCN.PA) <- levels(model.data$IUCN.PA)
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
  points(dop,yu,type="l",lty=2, col =  zone.cols[2])
  points(dop,yl,type="l",lty=2, col =  zone.cols[2])
  rug(data$DoP.PA[which(data$Zone == "Temperate")],
		col = zone.cols[2], lwd = 1, pos = 0)



legend("topleft", c("Tropical", "Temperate") , col = zone.cols, lty = 1, lwd = c(2,1))

dev.off()





#### within PA and hpd ##########




tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/11_14/abundance vs within pa vs hpd.tif",
	width = 15, height = 12, units = "cm", pointsize = 12, res = 300)



  par(mar=c(4,4.5,4,1.5))
  par(mgp=c(2.5,1,0))
  
 L = 100

  hpd <-seq(from=min(model.data$log_hpd[which(data$Within_PA == "yes")])
		,to=max(model.data$log_hpd[which(data$Within_PA == "yes")]),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  



model.data$Within_PA <- relevel(model.data$Within_PA, "yes")
model.data$taxon_of_interest <- relevel(model.data$taxon_of_interest , "Invertebrates")
mam <- lmer(log_abundance.model$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

 for (x in 1:L)
  {
    ag_suit <- median(model.data$ag_suit)
    log_abundance <-0
    log_hpd <- hpd[x]
    log_access <- mean(model.data$log_access)
    log_AREA.PA <- mean(model.data$log_AREA.PA)
    DoP.PA <- mean(model.data$DoP.PA)
    Zone <- "Tropical"
    taxon_of_interest<- "Invertebrates"
    Predominant_habitat <-"Primary Vegetation"
    log_bound_dist_km_PA_neg <- mean(model.data$log_bound_dist_km_PA_neg)
    Within_PA <- "yes"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, log_abundance, log_AREA.PA, DoP.PA, log_access, log_hpd, ag_suit))
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

  hpd_plot <- exp(hpd) -1

  z  <- exp(z)
  zu <- exp(zu)
  zl <- exp(zl)


  plot(hpd_plot,z, ylim=c(0, 300), xlim = c(0.05,10000),
		 col =1, 
		bty = "l", log = "x", lwd =2, xaxt = "n",
		type = "l",ylab = "Abundance per site ± s.e", xlab="Human population density (per km2)")
  points(hpd_plot,zu,type="l",lty=2,lwd = 2, col = 1)
  points(hpd_plot,zl,type="l",lty=2, lwd = 2,col = 1)
  axis(1, c(0.1,1,10,100,1000,10000), c(0.1,1,10,100,1000,10000) )
  rug(data$hpd[which(data$Within_PA == "yes")],
		col = 1, lwd = 1)




  hpd <-seq(from=min(model.data$log_hpd[which(data$Within_PA == "no")])
		,to=max(model.data$log_hpd[which(data$Within_PA == "no")]),length=L)
  y<-vector(mode="numeric",length=L)
  yu<-vector(mode="numeric",length=L)
  yl<-vector(mode="numeric",length=L)
  





model.data$Within_PA <- relevel(model.data$Within_PA, "no")
mam <- lmer(log_abundance.model$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

 for (x in 1:L)
  {
    ag_suit <- median(model.data$ag_suit)
    log_abundance <-0
    log_hpd <- hpd[x]
    log_access <- mean(model.data$log_access)
    log_AREA.PA <- mean(model.data$log_AREA.PA)
    DoP.PA <- mean(model.data$DoP.PA)
    Zone <-"Temperate"
    taxon_of_interest<- "Invertebrates"
    Predominant_habitat <-"Primary Vegetation"
    log_bound_dist_km_PA_neg <- mean(model.data$log_bound_dist_km_PA_neg)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, log_abundance, log_AREA.PA, DoP.PA, log_access, log_hpd, ag_suit))
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

  hpd_plot <- exp(hpd) -1

  points(hpd_plot,y,type="l",lty=1, col =  8)
  points(hpd_plot,yu,type="l",lty=2, col =  8)
  points(hpd_plot,yl,type="l",lty=2, col =  8)
  rug(data$hpd[which(data$Within_PA == "no")],
		col = 8, lwd = 1, pos = 0)


legend("topright", c("Protected", "Unprotected") , col = c(1,8), lty = 1, lwd = c(2))



dev.off()






#### within PA and hpd ##########
### and diff zones ####



tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/abundance vs within pa vs hpd vs zone.tif",
	width = 25, height = 12, units = "cm", pointsize = 12, res = 300)


zones <- c("Tropical", "Temperate")
t <- zones[1]

par(mfrow = c(1,2))

for(t in zones) { 


  par(mar=c(4,4.5,4,1.5))
  par(mgp=c(2.5,1,0))
  
 L = 100

  hpd <-seq(from=min(model.data$log_hpd[which(data$Within_PA == "yes" & data$Zone == t)])
		,to=max(model.data$log_hpd[which(data$Within_PA == "yes" & data$Zone == t)]),length=L)
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
    log_hpd <- hpd[x]
    log_access <- mean(model.data$log_access)
    log_AREA.PA <- mean(model.data$log_AREA.PA)
    DoP.PA <- mean(model.data$DoP.PA)
    Zone <-t
    taxon_of_interest<- "Invertebrates"
    Predominant_habitat <-"Primary Vegetation"
    log_bound_dist_km_PA_neg <- mean(model.data$log_bound_dist_km_PA_neg)
    Within_PA <- "yes"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, log_abundance, log_AREA.PA, DoP.PA, log_access, log_hpd, ag_suit))
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

  hpd_plot <- exp(hpd) -1

  z  <- exp(z)
  zu <- exp(zu)
  zl <- exp(zl)


  plot(hpd_plot,z, ylim=c(0, 500), xlim = c(0.05,10000),
		 col =1, main = zones[which(zones == t)],
		bty = "l", log = "x", lwd =2, xaxt = "n",
		type = "l",ylab = "Abundance per site ± s.e", xlab="Human population density (per km2)")
  points(hpd_plot,zu,type="l",lty=2,lwd = 2, col = 1)
  points(hpd_plot,zl,type="l",lty=2, lwd = 2,col = 1)
  axis(1, c(0.1,1,10,100,1000,10000), c(0.1,1,10,100,1000,10000) )
  rug(data$hpd[which(data$Within_PA == "yes" & data$Zone == t)],
		col = 1, lwd = 1)



  hpd <-seq(from=min(model.data$log_hpd[which(data$Within_PA == "no" & data$Zone == t)])
		,to=max(model.data$log_hpd[which(data$Within_PA == "no" & data$Zone == t)]),length=L)
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
    log_hpd <- hpd[x]
    log_access <- mean(model.data$log_access)
    log_AREA.PA <- mean(model.data$log_AREA.PA)
    DoP.PA <- mean(model.data$DoP.PA)
    Zone <-"Temperate"
    taxon_of_interest<- "Invertebrates"
    Predominant_habitat <-"Primary Vegetation"
    log_bound_dist_km_PA_neg <- mean(model.data$log_bound_dist_km_PA_neg)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, log_abundance, log_AREA.PA, DoP.PA, log_access, log_hpd, ag_suit))
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

  hpd_plot <- exp(hpd) -1

  points(hpd_plot,y,type="l",lty=1, col =  8)
  points(hpd_plot,yu,type="l",lty=2, col =  8)
  points(hpd_plot,yl,type="l",lty=2, col =  8)
  rug(data$hpd[which(data$Within_PA == "no" & data$Zone == t)],
		col = 8, lwd = 1, pos = 0)

}



legend("topright", c("Protected", "Unprotected") , col = c(1,8), lty = 1, lwd = c(2))


dev.off()





### PA size ####


tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/11_14/abundance vs PA size.tif",
	width = 15, height = 12, units = "cm", pointsize = 12, res = 300)



  par(mar=c(4,4.5,4,1.5))
  par(mgp=c(2.5,1,0))
  


 L = 30

  gis <-seq(from=min(model.data$log_AREA.PA),to=max(model.data$log_AREA.PA),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  




model.data$Zone<- relevel(model.data$Zone , "Tropical")
mam <- lmer(log_abundance.model$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

 for (x in 1:L)
  {
    log_bound_dist_km_PA_neg <- mean (model.data$log_bound_dist_km_PA_neg)
    log_abundance <-0
    log_hpd <- mean(model.data$log_hpd)
    log_access<- mean(model.data$log_access)
    log_AREA.PA <- gis[x]
    IUCN.PA <- "1.5"
    DoP.PA <- mean(model.data$DoP.PA)
    Zone <-"Tropical"
    taxon_of_interest<- "Invertebrates"
    Predominant_habitat <-"Primary Vegetation"
    ag_suit <- median(model.data$ag_suit)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, log_abundance, log_AREA.PA, DoP.PA, log_access, log_hpd, ag_suit))
    newdat.f<-data.frame(cbind(Zone,taxon_of_interest, Predominant_habitat, Within_PA))
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
  


 # transforming back - take - off first and put back
 gis <- exp(gis) - 1

 	
  z  <- exp(z) 
  zu <- exp(zu)
  zl <- exp(zl)

  plot(gis,z, ylim=c(0, 300), col = 1,
		bty = "l", log = "x", 
		type = "l",ylab = "Abundance per site ± s.e", xlab="Protected area size (km2)")
  rug(data$AREA.PA
	, ticksize = 0.03, side = 1, lwd = 0.5, col = 1) 
  axis(2,at = c(100,200,500,1000), c(100,200,500,1000))
  points(gis,zu,type="l",lty=2, col = 1)
  points(gis,zl,type="l",lty=2, col = 1)


dev.off()




#### within PA and hpd in diff IUCN cats ##########



tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/abundance vs IUCN cat vs hpd.tif",
	width = 20, height = 20, units = "cm", pointsize = 12, res = 300)

par(mfrow = c(2,2))


titles <- c("I and II", "III to VI", "NA or Not Reported", "Unprotected")


  par(mar=c(4,4.5,4,1.5))
  par(mgp=c(2.5,1,0))
  
 L = 30

  hpd <-seq(from=0,to=max(model.data$log_hpd),length=L)
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
    log_hpd <- hpd[x]
    log_access <- mean(model.data$log_access)
    log_AREA.PA <- mean(model.data$log_AREA.PA)
    IUCN.PA <- c
    DoP.PA <- mean(model.data$DoP.PA)
    Zone <-"Temperate"
    taxon_of_interest<- "Invertebrates"
    Predominant_habitat <-"Primary Vegetation"
    log_bound_dist_km_PA_neg <- mean(model.data$log_bound_dist_km_PA_neg)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, log_abundance, log_AREA.PA, DoP.PA, log_access, log_hpd, ag_suit))
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

  hpd_plot <- exp(hpd) -1

  z  <- exp(z)
  zu <- exp(zu)
  zl <- exp(zl)

  plot(hpd_plot,z, ylim=c(0, 600), col = 8, main = titles[which(cats == c)],
		bty = "l", log = "x", lwd =2,
		type = "l",ylab = "Abundance per site ± s.e", xlab="Human population density (per km2)")
  points(hpd_plot,zu,type="l",lty=2,lwd = 2, col = 8)
  points(hpd_plot,zl,type="l",lty=2, lwd = 2,col = 8)
  rug(data$hpd[which(data$IUCN.PA == c)], pos = 0,
	 ticksize = 0.03, side = 1,  lwd = 0.5, col = 8)


}


dev.off()






#### within PA and hpd in diff IUCN cats ##########
### and diff zones ####



tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/abundance vs IUCN cat vs hpd vs zone.tif",
	width = 20, height = 20, units = "cm", pointsize = 12, res = 300)

par(mfrow = c(2,2))


titles <- c("I and II", "III to VI", "NA or Not Reported", "Unprotected")


  par(mar=c(4,4.5,4,1.5))
  par(mgp=c(2.5,1,0))
  
 L = 30

  hpd <-seq(from=0,to=max(model.data$log_hpd),length=L)
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
    log_hpd <- hpd[x]
    log_access <- mean(model.data$log_access)
    log_AREA.PA <- mean(model.data$log_AREA.PA)
    IUCN.PA <- c
    DoP.PA <- mean(model.data$DoP.PA)
    Zone <-"Temperate"
    taxon_of_interest<- "Invertebrates"
    Predominant_habitat <-"Primary Vegetation"
    log_bound_dist_km_PA_neg <- mean(model.data$log_bound_dist_km_PA_neg)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, log_abundance, log_AREA.PA, DoP.PA, log_access, log_hpd, ag_suit))
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

  hpd_plot <- exp(hpd) -1

  z  <- exp(z)
  zu <- exp(zu)
  zl <- exp(zl)

  plot(hpd_plot,z, ylim=c(-50, 600), col = zone.cols[1], main = titles[which(cats == c)],
		bty = "l", log = "x", lwd =2,
		type = "l",ylab = "Abundance per site ± s.e", xlab="Human population density (per km2)")
  points(hpd_plot,zu,type="l",lty=2,lwd = 2, col = zone.cols[1])
  points(hpd_plot,zl,type="l",lty=2, lwd = 2,col = zone.cols[1])
  rug(data$hpd[which(data$IUCN.PA == c & data$Zone == "Tropical")],
	 ticksize = 0.03, side = 1,  lwd = 0.5, col = zone.cols[1])




model.data$Zone<- relevel(model.data$Zone , "Temperate")
mam <- lmer(log_abundance.model$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

 for (x in 1:L)
  {
    ag_suit <- median(model.data$ag_suit)
    log_abundance <-0
    log_hpd <- hpd[x]
    log_access <- mean(model.data$log_access)
    log_AREA.PA <- mean(model.data$log_AREA.PA)
    IUCN.PA <- c
    DoP.PA <- mean(model.data$DoP.PA)
    Zone <-"Temperate"
    taxon_of_interest<- "Invertebrates"
    Predominant_habitat <-"Primary Vegetation"
    log_bound_dist_km_PA_neg <- mean(model.data$log_bound_dist_km_PA_neg)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, log_abundance, log_AREA.PA, DoP.PA, log_access, log_hpd, ag_suit))
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

  points(hpd_plot,y,type="l",lty=1, col =  zone.cols[2])
  points(hpd_plot,yu,type="l",lty=2, col =  zone.cols[2])
  points(hpd_plot,yl,type="l",lty=2, col =  zone.cols[2])
  rug(data$hpd[which(data$IUCN.PA == c & data$Zone == "Temperate")],
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
    log_hpd <- mean(model.data$log_hpd)
    log_access<- mean(model.data$log_access)
    log_AREA.PA <- gis[x]
    IUCN.PA <- "1.5"
    DoP.PA <- mean(model.data$DoP.PA)
    Zone <-"Tropical"
    taxon_of_interest<- "Invertebrates"
    Predominant_habitat <-"Primary Vegetation"
    ag_suit <- median(model.data$ag_suit)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, log_abundance, log_AREA.PA, DoP.PA, log_access, log_hpd, ag_suit))
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












#### within PA and ag suit ##########



tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/abundance within pa vs ag_suit.tif",
	width = 15, height = 10, units = "cm", pointsize = 12, res = 300)

par(mfrow = c(1,1))


  par(mar=c(4,4.5,4,1.5))
  par(mgp=c(2.5,1,0))
  
 L = 30

  ag <-seq(from=0,to=max(model.data$ag_suit),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  


model.data$Within_PA <- relevel(model.data$Within_PA, "no")
model.data$Zone<- relevel(model.data$Zone , "Tropical")
model.data$taxon_of_interest <- relevel(model.data$taxon_of_interest , "Invertebrates")
mam <- lmer(log_abundance.model$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

 for (x in 1:L)
  {
    ag_suit <- ag[x]
    log_bound_dist_km_PA_neg <- mean (model.data$log_bound_dist_km_PA_neg)
    log_abundance <-0
    log_hpd <- mean(model.data$log_hpd)
    log_access<- mean(model.data$log_access)
    log_AREA.PA <- mean(model.data$log_AREA.PA)
    DoP.PA <- mean(model.data$DoP.PA)
    Zone <-"Tropical"
    taxon_of_interest<- "Invertebrates"
    Predominant_habitat <-"Primary Vegetation"
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, log_abundance, log_AREA.PA, DoP.PA, log_access, log_hpd, ag_suit))
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

  plot(ag,z, ylim=c(0, 300), col = 8,
		bty = "l", #log = "x",
		type = "l",ylab = "Abundance per site ± s.e", xlab="Agricultural suitability (higher = more suitable)")
  points(ag,zu,type="l",lty=2, col = 8)
  points(ag,zl,type="l",lty=2, col = 8)







model.data$Within_PA <- relevel(model.data$Within_PA, "yes")

mam<- lmer(log_abundance.model$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))


  y<-vector(mode="numeric",length=L)
  yu<-vector(mode="numeric",length=L)
  yl<-vector(mode="numeric",length=L)


 for (x in 1:L)
  {
    ag_suit <- ag[x]
    log_bound_dist_km_PA_neg <- mean (model.data$log_bound_dist_km_PA_neg)
    log_abundance <-0
    log_hpd <- mean(model.data$log_hpd)
    log_access<- mean(model.data$log_access)
    log_AREA.PA <- mean(model.data$log_AREA.PA)
    IUCN.PA <- "1.5"
    DoP.PA <- mean(model.data$DoP.PA)
    Zone <-"Tropical"
    taxon_of_interest<- "Invertebrates"
    Predominant_habitat <-"Primary Vegetation"
    Within_PA <- "yes"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, log_abundance, log_AREA.PA, DoP.PA, log_access, log_hpd, ag_suit))
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


  
  points(ag,y,type="l",lty=1 , col = 1, lwd = 2)
  points(ag,yu,type="l",lty=3, lwd = 2, col = 1)
  points(ag,yl,type="l",lty=3, lwd = 2, col = 1)





legend("topright", c("Protected", "Unprotected") , col = c(1,8), lty = c(1,1), lwd = c(1,2))


dev.off()










#### accessibility ##########



tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/11_14/abundance vs accessibility.tif",
	width = 19, height = 10, units = "cm", pointsize = 12, res = 300)

par(mfrow = c(1,1))


  par(mar=c(4,4.5,4,1.5))
  par(mgp=c(2.5,1,0))
  
 L = 30

  acc <-seq(from=0,to=max(model.data$log_acc),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  


# it doesnt matter what the values for taxon and zone are, as there is no interaction with these



model.data$Within_PA <- relevel(model.data$Within_PA, "no")
model.data$Zone<- relevel(model.data$Zone , "Tropical")
model.data$taxon_of_interest <- relevel(model.data$taxon_of_interest , "Invertebrates")
mam <- lmer(log_abundance.model$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

 for (x in 1:L)
  {
    ag_suit <- median(model.data$ag_suit)
    log_bound_dist_km_PA_neg <- mean (model.data$log_bound_dist_km_PA_neg)
    log_abundance <-0
    log_hpd <- mean(model.data$log_hpd)
    log_access<- acc[x]
    log_AREA.PA <- mean(model.data$log_AREA.PA)
    IUCN.PA <- "1.5"
    DoP.PA <- mean(model.data$DoP.PA)
    Zone <-"Tropical"
    taxon_of_interest<- "Invertebrates"
    Predominant_habitat <-"Primary Vegetation"
    Within_PA <- "yes"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, log_abundance, log_AREA.PA, DoP.PA, log_access, log_hpd, ag_suit))
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

  acc <- exp(acc)/60

  z  <- exp(z)
  zu <- exp(zu)
  zl <- exp(zl)

  plot(acc,z, ylim=c(0, 300), col = 1,
		bty = "l", log = "x",
		type = "l",ylab = "Abundance per site ± s.e", xlab="Accessibility (hours to city >50 000)")
  points(acc,zu,type="l",lty=2, col = 1)
  points(acc,zl,type="l",lty=2, col = 1)
  rug(data$access/60, ticksize = 0.03, side = 1,  lwd = 0.5, col = 1) #### change to be data used in model, poss slightly fewer points


dev.off()


