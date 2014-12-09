
# PLOTS


model.data <- RLS.model$data

data <- matched.landuse_amp_mam_bir[,c("RLS.y", "Zone", "taxon_of_interest", "ag_suit", "log_access", "access", "DoP.PA", "AREA.PA", "IUCN.PA",
	 "log_hpd", "hpd", "log_bound_dist_km_PA_neg", "bound_dist_km_PA_neg", "Within_PA", "Predominant_habitat", "SS", "SSBS")]

data <- na.omit(data) #have to get rid of NA to try poly on ag_suit



nrow(data)
length(unique(data$SS))



nrow(model.data)
length(unique(model.data$SS))


L = 100

#green one
#inside.col <- rgb(0.1,0.8,0.2)


inside.col <- 1
outside.col <- rgb(0.7,0.7,0.7)




display.brewer.all()
cols <- brewer.pal(8, "Paired")
display.brewer.pal(8, "Paired")
#drop the red to avoid red-green colorblindness issues

zone.cols <- cols[c(6,5)]








####  Within PA and hpd ##########




tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/11_14/prop threatened within pa vs hpd.tif",
	width = 12, height = 10, units = "cm", pointsize = 12, res = 300)



par(mfrow = c(1,1))


  par(mar=c(4,4.5,1,1.5))
  par(mgp=c(2.5,1,0))
  
 L = 100

  hpd <-seq(from=min(model.data$log_hpd[which(model.data$Within_PA == "no")]),
		to=max(model.data$log_hpd[which(model.data$Within_PA == "no")]),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  


# it doesnt matter what the values for taxon and zone are, as there is no interaction with these



model.data$Within_PA <- relevel(model.data$Within_PA, "no")
model.data$Zone<- relevel(model.data$Zone , "Tropical")
mam <- glmer(RLS.model$final.call, model.data, family = "binomial")

 for (x in 1:L)
  {
    log_hpd <- hpd[x]
    RLS.y <- 0
    log_access <- mean(model.data$log_access)
    ag_suit <- median(model.data$ag_suit)
    log_AREA.PA <- mean(model.data$log_AREA.PA)
    DoP.PA <- mean(model.data$DoP.PA )
    Zone <-"Temperate"
    Predominant_habitat <-"Primary Vegetation"
    log_bound_dist_km_PA_neg <- mean(model.data$log_bound_dist_km_PA_neg)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, log_access, log_AREA.PA, DoP.PA, log_hpd, ag_suit, RLS.y))
    newdat.f<-data.frame(cbind(Zone, Predominant_habitat, Within_PA))
    levels(newdat.f$Zone)<-levels(model.data$Zone)
    levels(newdat.f$Predominant_habitat)<-levels(model.data$Predominant_habitat)
    levels(newdat.f$Within_PA) <- levels(model.data$Within_PA)
    newdat<-data.frame(newdat.f,newdat.n)
    mm<-model.matrix(terms(mam),newdat)
    pvar1 <- diag(mm %*% tcrossprod(vcov(mam),mm))
    z[x]<-mm %*% fixef(mam)
    zu[x]<-z[x]+sqrt(pvar1)
    zl[x]<-z[x]-sqrt(pvar1) 
  }
  
  hpd_plot <- exp(hpd) -1

 # original backtransform, for binomial 
  z<-1/(1+exp(-(z)))
  zu<-1/(1+exp(-(zu)))
  zl<-1/(1+exp(-(zl)))

  plot(hpd_plot,z, ylim=c(-0.002,0.01),  col = outside.col,
		bty = "l", log = "x", xaxt = "n",
		type = "l",ylab = "Proportion threatened species ± s.e", 
		xlab="Human population density (per km2)")
   rug(data$hpd[which(data$Within_PA == "no")]
	, ticksize = 0.03, side = 1,   lwd = 1, col = outside.col, pos = -0.002) 
  axis(1, c(0.05, 0.1, 0.2, 0.5, 1, 2, 5, 20, 50,500), c(0.05, 0.1, 0.2, 0.5, 1, 2, 5, 20, 50,500))
  points(hpd_plot,zu,type="l",lty=2, col = outside.col)
  points(hpd_plot,zl,type="l",lty=2, col = outside.col)





model.data$Within_PA <- relevel(model.data$Within_PA, "yes")

mam <- glmer(RLS.model$final.call, model.data, family = "binomial")

hpd.y <-seq(from=min(model.data$log_hpd[which(model.data$Within_PA == "yes")]),
		to=max(model.data$log_hpd[which(model.data$Within_PA == "yes")]),length=L)
  y<-vector(mode="numeric",length=length(hpd))
  yu<-vector(mode="numeric",length=length(hpd))
  yl<-vector(mode="numeric",length=length(hpd))


 for (x in 1:L)
  {
    log_hpd <- hpd[x]
    RLS.y <-0
    log_access <- mean(model.data$log_access)
    ag_suit <- median(model.data$ag_suit)
    log_AREA.PA <- mean(model.data$log_AREA.PA )
    DoP.PA <- mean(model.data$DoP.PA )
    Zone <-"Temperate"
    Predominant_habitat <-"Primary Vegetation"
    log_bound_dist_km_PA_neg <- mean(model.data$log_bound_dist_km_PA_neg)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, log_access, log_AREA.PA, DoP.PA, RLS.y, log_hpd, ag_suit))
    newdat.f<-data.frame(cbind( Zone, Predominant_habitat, Within_PA))
    levels(newdat.f$Zone)<-levels(model.data$Zone)
    levels(newdat.f$Predominant_habitat)<-levels(model.data$Predominant_habitat)
    levels(newdat.f$Within_PA) <- levels(model.data$Within_PA)
    newdat<-data.frame(newdat.f,newdat.n)
    mm<-model.matrix(terms(mam),newdat)
    pvar1 <- diag(mm %*% tcrossprod(vcov(mam),mm))
    y[x]<-mm %*% fixef(mam)
    yu[x]<-y[x]+sqrt(pvar1)
    yl[x]<-y[x]-sqrt(pvar1) 
  }

  hpd.y_plot <- exp(hpd.y) -1
  y<-1/(1+exp(-(y)))
  yu<-1/(1+exp(-(yu)))
  yl<-1/(1+exp(-(yl)))


  points(hpd.y_plot,y,type="l",lty=1 , col = inside.col, lwd = 1.5)
  points(hpd.y_plot,yu,type="l",lty=3, lwd = 1.5, col = inside.col)
  points(hpd.y_plot,yl,type="l",lty=3, lwd = 1.5, col = inside.col)
   rug(data$hpd[which(data$Within_PA == "yes")]
	, ticksize = 0.03, side = 1,   lwd = 2, col = inside.col) 





legend("topright", c("Protected", "Unprotected") , col = c(inside.col,outside.col), lty = c(1,1), lwd = c(2,1))


dev.off()




### Access 


tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/11_14/prop threatened vs access.tif",
	width = 12, height = 10, units = "cm", pointsize = 12, res = 300)




  par(mar=c(4,4.5,1,1.5))
  par(mgp=c(2.5,1,0))
  
 L = 100

  acc <-seq(from=min(model.data$log_access),
		to=max(model.data$log_access),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  


# it doesnt matter what the values for taxon and zone are, as there is no interaction with these



model.data$Within_PA <- relevel(model.data$Within_PA, "no")
model.data$Zone<- relevel(model.data$Zone , "Tropical")
mam <- glmer(RLS.model$final.call, model.data, family = "binomial")

 for (x in 1:L)
  {
    log_access <- acc[x]
    RLS.y <- 0
    log_hpd <- mean(model.data$log_hpd)
    ag_suit <- median(model.data$ag_suit)
    log_AREA.PA <- mean(model.data$log_AREA.PA)
    DoP.PA <- mean(model.data$DoP.PA )
    Zone <-"Temperate"
    Predominant_habitat <-"Primary Vegetation"
    log_bound_dist_km_PA_neg <- mean(model.data$log_bound_dist_km_PA_neg)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, log_access, log_AREA.PA, DoP.PA, log_hpd, ag_suit, RLS.y))
    newdat.f<-data.frame(cbind(Zone, Predominant_habitat, Within_PA))
    levels(newdat.f$Zone)<-levels(model.data$Zone)
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

 # original backtransform, for binomial 
  z<-1/(1+exp(-(z)))
  zu<-1/(1+exp(-(zu)))
  zl<-1/(1+exp(-(zl)))

  plot(acc_plot,z, ylim=c(-0.0005,0.001),  col = 1,
		bty = "l", log = "x", xaxt = "n", #yaxt = "n",
		type = "l",ylab = "Proportion threatened species ± s.e", 
		xlab="Distance to city >50000 (hours)")
   rug((data$acc)/60
	, ticksize = 0.03, side = 1,   lwd = 1, col = 1) 
 # axis(2, c(0, 0.0005, 0.001),  c(0, 0.0005, 0.001))
  axis(1, c(0.05, 0.1, 0.2, 0.5, 1, 2, 5, 20, 50,500), c(0.05, 0.1, 0.2, 0.5, 1, 2, 5, 20, 50,500))
  points(acc_plot,zu,type="l",lty=2, col = 1)
  points(acc_plot,zl,type="l",lty=2, col = 1)




dev.off()





######################################### below here extra code just in case needed 





###### dist to boundary for different zones ##################



tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/prop threatened vs zone  vs dist to boundary.tif",
	width = 15, height = 12, units = "cm", pointsize = 12, res = 300)




  par(mar=c(4,4.5,4,1.5))
  par(mgp=c(2.5,1,0))
  
 L = 30

 lbd <-seq(from=min(model.data$log_bound_dist_km_PA_neg[which(data$Zone == "Tropical")]),
		to=max(model.data$log_bound_dist_km_PA_neg[which(data$Zone == "Tropical")]),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  





model.data$Zone<- relevel(model.data$Zone , "Tropical")
mam <- glmer(RLS.model$final.call, model.data, family = "binomial")

 for (x in 1:L)
  {
    log_bound_dist_km_PA_neg <- lbd[x]
    RLS.y <-0
    log_hpd <- mean(model.data$log_hpd)
    log_AREA.PA <- mean(model.data$log_AREA.PA )
    DoP.PA <- mean(model.data$DoP.PA )
    IUCN.PA <- "1.5"
    Zone <-"Tropical"
    Predominant_habitat <-"Primary Vegetation"
    ag_suit <- median(model.data$ag_suit)
    log_access <- mean(model.data$log_access)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, RLS.y, log_hpd,log_AREA.PA, DoP.PA, log_access, ag_suit))
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
  

  z<-1/(1+exp(-(z)))
  zu<-1/(1+exp(-(zu)))
  zl<-1/(1+exp(-(zl)))


 lbd1 <- exp(abs(lbd)) -1
 inside <- which(lbd < 0)
 lbd1[inside] <- lbd1[inside]*-1

  plot(lbd1,z, ylim=c(-0.0003,0.003), xlim = c(-50, 200), 
		col = zone.cols[1], lwd = 2,
		bty = "l", #log = "x", #yaxt = "n", 
		type = "l",ylab = "Proportion species threatened ± s.e", xlab="Distance to PA boundary (km)")
 rug(data$bound_dist_km_PA_neg[which(data$Zone == "Tropical")]
	, ticksize = 0.03, side = 1, lwd = 0.5, col = zone.cols[1])
  points(lbd1,zu,type="l",lty=2,  lwd = 2, col =  zone.cols[1])
  points(lbd1,zl,type="l",lty=2,  lwd = 2, col =  zone.cols[1])



#  plot(lbd,z, ylim=c(4.5,8), xlim = c(-1*log(50+1), log(200+1)), 
#		col = zone.cols[1], lwd = 2,
#		bty = "l", #log = "x", 
#		xaxt = "n", 
#		type = "l",ylab = "RLS.y size Log10(sq km) ± s.e", xlab="Distance to PA boundary (km)", main = t)
 # rug(data$log_bound_dist_km_PA_neg[which(data$taxon_of_interest == t & data$Zone == "Tropical")]
#	, ticksize = 0.03, side = 1, lwd = 0.5, col = zone.cols[1])
#  points(lbd,zu,type="l",lty=2,  lwd = 2, col =  zone.cols[1])
#  points(lbd,zl,type="l",lty=2,  lwd = 2, col =  zone.cols[1])
#  axis(1, c(-1*log(10+1), 0, log(10+1),log(50+1), log(100+1), log(200+1)), c(-10, 0, 10, 50, 100, 200))

 

  lbd <-seq(from=min(model.data$log_bound_dist_km_PA_neg[which(data$Zone == "Temperate")]),
		to=max(model.data$log_bound_dist_km_PA_neg[which(data$Zone == "Temperate")]),length=L)
  y<-vector(mode="numeric",length=L)
  yu<-vector(mode="numeric",length=L)
  yl<-vector(mode="numeric",length=L)
  



model.data$Zone<- relevel(model.data$Zone , "Temperate")
mam <- glmer(RLS.model$final.call, model.data, family = "binomial")


 for (x in 1:L)
  {
    log_bound_dist_km_PA_neg <- lbd[x]
    RLS.y <-0
    log_hpd <- mean(model.data$log_hpd)
    log_AREA.PA <- mean(model.data$log_AREA.PA )
    DoP.PA <- mean(model.data$DoP.PA )
    IUCN.PA <- "1.5"
    Zone <-"Tropical"
    Predominant_habitat <-"Primary Vegetation"
    ag_suit <- median(model.data$ag_suit)
    log_access <- mean(model.data$log_access)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, RLS.y, log_hpd,log_AREA.PA, DoP.PA, log_access, ag_suit))
    newdat.f<-data.frame(cbind( Zone, Predominant_habitat, Within_PA))
    levels(newdat.f$Zone)<-levels(model.data$Zone)
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
  

 lbd2 <- exp(abs(lbd)) -1
 inside <- which(lbd < 0)
 lbd2[inside] <- lbd2[inside]*-1

  y<-1/(1+exp(-(y)))
  yu<-1/(1+exp(-(yu)))
  yl<-1/(1+exp(-(yl)))

  points(lbd2,y,type="l",lty=1, col =  zone.cols[2])
  points(lbd2,yu,type="l",lty=2, col =  zone.cols[2])
  points(lbd2,yl,type="l",lty=2, col =  zone.cols[2])
  rug(data$bound_dist_km_PA_neg[which(data$Zone == "Temperate")],
		col = zone.cols[2], lwd = 1, pos = -0.0003)



#  points(lbd,y,type="l",lty=2, col =  zone.cols[2])
#  points(lbd,yu,type="l",lty=2, col =  zone.cols[2])
#  points(lbd,yl,type="l",lty=2, col =  zone.cols[2])
#  rug(data$log_bound_dist_km_PA_neg[which(data$Zone == "Temperate")],
#		col = zone.cols[2], lwd = 1, pos = 0.1)


 abline(v = 0, lty = 2, col = 8)



legend("topright", c("Tropical", "Temperate") , col = zone.cols, lty = c(1,1), lwd = c(2,1.2))





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
mam <- lmer(RLS.model$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

 for (x in 1:L)
  {
    log_bound_dist_km_PA_neg <-  mean(model.data$log_bound_dist_km_PA_neg)
    RLS <-0
    log_hpd <- mean(model.data$log_hpd)
    log_AREA.PA <-gis[x]
    DoP.PA <- mean(model.data$DoP.PA)
    IUCN.PA <- "1.5"
    Zone <-"Tropical"
    Predominant_habitat <-"Primary Vegetation"
    ag_suit <- median(model.data$ag_suit)
    log_access <- mean(model.data$log_access)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, RLS, log_hpd,log_AREA.PA, DoP.PA, log_access, ag_suit))
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
  


 gis_plot <- exp(gis) -1


  plot(gis_plot,z, ylim=c(3,7), col = zone.cols[1], lwd = 2,
		bty = "l", log = "x", #yaxt = "n", 
		type = "l",ylab = "CWM IUCN Red List Score ± s.e", xlab="PA size (km2)")
  rug(data$AREA.PA[which(data$Zone == "Tropical")]
	, ticksize = 0.03, side = 1, lwd = 2, col = zone.cols[1])
  points(gis_plot,zu,type="l",lty=2,  lwd = 2, col =  zone.cols[1])
  points(gis_plot,zl,type="l",lty=2,  lwd = 2.5, col =  zone.cols[1])


 gis <-seq(from=min(model.data$log_AREA.PA[which(model.data$Zone == "Temperate")]),
		to=max(model.data$log_AREA.PA[which(model.data$Zone == "Temperate")]),length=L)

  y<-vector(mode="numeric",length=L)
  yu<-vector(mode="numeric",length=L)
  yl<-vector(mode="numeric",length=L)
  


model.data$Zone<- relevel(model.data$Zone , "Temperate")
mam <- lmer(RLS.model$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

 for (x in 1:L)
  {
    log_bound_dist_km_PA_neg <-  mean(model.data$log_bound_dist_km_PA_neg)
    RLS <-0
    log_hpd <- mean(model.data$log_hpd)
    log_AREA.PA <-gis[x]
    DoP.PA <- mean(model.data$DoP.PA)
    IUCN.PA <- "1.5"
    Zone <-"Tropical"
    Predominant_habitat <-"Primary Vegetation"
    ag_suit <- median(model.data$ag_suit)
    log_access <- mean(model.data$log_access)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, RLS, log_hpd,log_AREA.PA, DoP.PA, log_access, ag_suit))
    newdat.f<-data.frame(cbind( Zone, Predominant_habitat, Within_PA))
    levels(newdat.f$Zone)<-levels(model.data$Zone)
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
  

  gis_plot <- exp(gis) -1

  points(gis_plot,y,type="l",lty=2, col =  zone.cols[2])
  points(gis_plot,yu,type="l",lty=2, col =  zone.cols[2])
  points(gis_plot,yl,type="l",lty=2, col =  zone.cols[2])
  rug(data$AREA.PA[which(data$Zone == "Temperate")],
		col = zone.cols[2], lwd = 0.8, pos = 0.14)



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
mam <- lmer(RLS.model$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

 for (x in 1:L)
  {
    log_bound_dist_km_PA_neg <- mean(model.data$log_bound_dist_km_PA_neg)
    RLS <-0
    log_hpd <- mean(model.data$log_hpd)
    log_AREA.PA <- mean(model.data$log_AREA.PA )
    DoP.PA <- dop[x]
    IUCN.PA <- "1.5"
    Zone <-"Tropical"
    Predominant_habitat <-"Primary Vegetation"
    ag_suit <- median(model.data$ag_suit)
    log_access <- mean(model.data$log_access)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, RLS, log_hpd,log_AREA.PA, DoP.PA, log_access, ag_suit))
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
  



  plot(dop,z, ylim=c(2, 7), col = zone.cols[1], lwd = 2,
		bty = "l",
		type = "l",ylab = "CWM IUCN Red List Score ± s.e", xlab="Duration of protection (yr)")
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
mam <- lmer(RLS.model$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

 for (x in 1:L)
  {
    log_bound_dist_km_PA_neg <- mean(model.data$log_bound_dist_km_PA_neg)
    RLS <-0
    log_hpd <- mean(model.data$log_hpd)
    log_AREA.PA <- mean(model.data$log_AREA.PA )
    DoP.PA <- dop[x]
    IUCN.PA <- "1.5"
    Zone <-"Tropical"
    Predominant_habitat <-"Primary Vegetation"
    ag_suit <- median(model.data$ag_suit)
    log_access <- mean(model.data$log_access)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, RLS, log_hpd,log_AREA.PA, DoP.PA, log_access, ag_suit))
    newdat.f<-data.frame(cbind( Zone, Predominant_habitat, Within_PA))
    levels(newdat.f$Zone)<-levels(model.data$Zone)
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
  
  points(dop,y,type="l",lty=1, col =  zone.cols[2])
  points(dop,yu,type="l",lty=2, col =  zone.cols[2])
  points(dop,yl,type="l",lty=2, col =  zone.cols[2])
  rug(data$DoP.PA[which(data$Zone == "Temperate")],
		col = zone.cols[2], lwd = 1, pos = 0.14)



legend("topright", c("Tropical", "Temperate") , col = zone.cols, lty = 1, lwd = c(2,1))

dev.off()






####  Within PA and access ##########

inside.col <- rgb(0.3,0.7,0.3)


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
mam <- lmer(RLS.model$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

 for (x in 1:L)
  {
    log_access <- acc[x]
    RLS <-0
    log_hpd <- mean(model.data$log_hpd)
    ag_suit <- median(model.data$ag_suit)
    log_AREA.PA <- mean(model.data$log_AREA.PA )
    DoP.PA <- mean(model.data$DoP.PA )
    Zone <-"Temperate"
    taxon_of_interest<- "Invertebrates"
    Predominant_habitat <-"Primary Vegetation"
    log_bound_dist_km_PA_neg <- mean(model.data$log_bound_dist_km_PA_neg)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, log_access, log_AREA.PA, DoP.PA, RLS, log_hpd, ag_suit))
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

  plot(acc_plot,z, #ylim=c(0.16,0.21), xlim = c(0.02,60), col = 1,
		bty = "l", log = "x", xaxt = "n",
		type = "l",ylab = "Endemicity (1/CWM RLS) ± s.e", 
		xlab="Accessibility (hours to city >50000)")
   rug((data$access[which(data$Within_PA == "no")])/60
	, ticksize = 0.03, side = 1,   lwd = 1, col = 1) 
  axis(1, c(0.05, 0.1, 0.2, 0.5, 1, 2, 5, 20, 50), c(0.05, 0.1, 0.2, 0.5, 1, 2, 5, 20, 50))
  points(acc_plot,zu,type="l",lty=2, col = 1)
  points(acc_plot,zl,type="l",lty=2, col = 1)





model.data$Within_PA <- relevel(model.data$Within_PA, "yes")

mam<- lmer(RLS.model$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

acc <-seq(from=min(model.data$log_access[which(model.data$Within_PA == "yes")]),
		to=max(model.data$log_access[which(model.data$Within_PA == "yes")]),length=L)
  y<-vector(mode="numeric",length=length(hpd))
  yu<-vector(mode="numeric",length=length(hpd))
  yl<-vector(mode="numeric",length=length(hpd))


 for (x in 1:L)
  {
    log_access <- acc[x]
    RLS <-0
    log_hpd <- mean(model.data$log_hpd)
    ag_suit <- median(model.data$ag_suit)
    log_AREA.PA <- mean(model.data$log_AREA.PA )
    DoP.PA <- mean(model.data$DoP.PA )
    Zone <-"Temperate"
    taxon_of_interest<- "Invertebrates"
    Predominant_habitat <-"Primary Vegetation"
    log_bound_dist_km_PA_neg <- mean(model.data$log_bound_dist_km_PA_neg)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, log_access, log_AREA.PA, DoP.PA, RLS, log_hpd, ag_suit))
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

  points(acc_plot,y,type="l",lty=1 , col = inside.col, lwd = 3)
  points(acc_plot,yu,type="l",lty=3, lwd = 3, col = inside.col)
  points(acc_plot,yl,type="l",lty=3, lwd = 3, col = inside.col)
   rug((data$access[which(data$Within_PA == "yes")])/60
	, ticksize = 0.03, side = 1,   lwd = 3, col = inside.col, pos = 0.16) 





legend("topright", c("Protected", "Unprotected") , col = c(inside.col,1), lty = c(1,1), lwd = c(1,2))


dev.off()






#### within PA and ag suit ##########



tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/endemicity within pa vs  ag_suit.tif",
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
model.data$Zone<- relevel(model.data$Zone, "Tropical")
model.data$taxon_of_interest <- relevel(model.data$taxon_of_interest , "Invertebrates")
mam <- lmer(RLS.model$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

 for (x in 1:L)
  {
    ag_suit <- ag[x]
    RLS <-0
    log_hpd <- mean(model.data$log_hpd)
    log_access <- mean(model.data$log_access)
    log_AREA.PA <- mean(model.data$log_AREA.PA )
    DoP.PA <- mean(model.data$DoP.PA )
    Zone <-"Temperate"
    taxon_of_interest<- "Invertebrates"
    Predominant_habitat <-"Primary Vegetation"
    log_bound_dist_km_PA_neg <- mean(model.data$log_bound_dist_km_PA_neg)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, log_access, log_AREA.PA, DoP.PA, RLS, log_hpd, ag_suit))
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

  plot(ag,z, ylim=c(0.18,0.22), col = 8,
		bty = "l", #log = "x",
		type = "l",ylab = "Endemicity (1/CWM RLS) ± s.e", xlab="Agricultural suitability (higher = more suitable)")
  points(ag,zu,type="l",lty=2, col = 8)
  points(ag,zl,type="l",lty=2, col = 8)




model.data$Within_PA <- relevel(model.data$Within_PA, "yes")

mam<- lmer(RLS.model$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

  ag <-seq(from=min(model.data$ag_suit[which(model.data$Within_PA == "yes")]),
		to=max(model.data$ag_suit[which(model.data$Within_PA == "yes")]),length=L)
  y<-vector(mode="numeric",length=length(hpd))
  yu<-vector(mode="numeric",length=length(hpd))
  yl<-vector(mode="numeric",length=length(hpd))


 for (x in 1:L)
  {
    ag_suit <- ag[x]
    RLS <-0
    log_hpd <- mean(model.data$log_hpd)
    log_access <- mean(model.data$log_access)
    log_AREA.PA <- mean(model.data$log_AREA.PA )
    DoP.PA <- mean(model.data$DoP.PA )
    Zone <-"Temperate"
    taxon_of_interest<- "Invertebrates"
    Predominant_habitat <-"Primary Vegetation"
    log_bound_dist_km_PA_neg <- mean(model.data$log_bound_dist_km_PA_neg)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, log_access, log_AREA.PA, DoP.PA, RLS, log_hpd, ag_suit))
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
  
  points(ag,y,type="l",lty=1 , col = 1, lwd = 2)
  points(ag,yu,type="l",lty=3, lwd = 2, col = 1)
  points(ag,yl,type="l",lty=3, lwd = 2, col = 1)


legend("topright", c("Protected", "Unprotected") , col = c(1,8), lty = c(1,1), lwd = c(1,2))



dev.off()






#### hpd ##########



tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/endemicity vs hpd.tif",
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
mam <- lmer(RLS.model$final.call, model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

 for (x in 1:L)
  {
    log_hpd <- hpd[x]
    RLS <-0
    log_access <- mean(model.data$log_access)
    ag_suit <- median(model.data$ag_suit)
    log_AREA.PA <- mean(model.data$log_AREA.PA )
    DoP.PA <- mean(model.data$DoP.PA )
    Zone <-"Temperate"
    taxon_of_interest<- "Invertebrates"
    Predominant_habitat <-"Primary Vegetation"
    log_bound_dist_km_PA_neg <- mean(model.data$log_bound_dist_km_PA_neg)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, log_access, RLS, log_AREA.PA, DoP.PA, log_hpd, ag_suit))
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

  hpd <- exp(hpd) - 1


  plot(hpd,z, ylim=c(0.16,0.21), col = 1 ,
		bty = "l", log = "x",
		type = "l",ylab = "Endemicity (1/CWM RLS) ± s.e", xlab="Human population density (per km2)")
   rug(data$hpd, ticksize = 0.03, side = 1,  lwd = 0.5, col = 1 ) #### change to be data used in model, poss slightly fewer points


  points(hpd,zu,type="l",lty=2, col = 1 )
  points(hpd,zl,type="l",lty=2, col = 1 )







dev.off()




####




