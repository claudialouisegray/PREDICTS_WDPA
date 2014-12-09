# PLOTS


model.data <- Richness_rarefied.model.int$data

data <- matched.landuse[,c("Richness_rarefied", "Zone", "taxon_of_interest", "ag_suit", "log_access", "access", "DoP.PA", "AREA.PA", "IUCN.PA",
	 "log_hpd", "hpd", "log_bound_dist_km_PA_neg", "bound_dist_km_PA_neg", "Within_PA", "Predominant_habitat", "SS", "SSBS")]

data <- na.omit(data) #have to get rid of NA to try poly on ag_suit

names(matched.landuse)

nrow(data)
length(unique(data$SS))



nrow(model.data)
length(unique(model.data$SS))


table(data$IUCN.PA)


display.brewer.all()
cols <- brewer.pal(8, "Paired")
display.brewer.pal(8, "Paired")
#drop the red to avoid red-green colorblindness issues

zone.cols <- cols[c(6,5)]
taxa.cols <- cols[c(2,4,8)]

inside.col <- 1
outside.col <- 8


#### within PA  vs ag suit ###########



### ag suit ####






tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/11_14/Richness_rarefied vs ag suit.tif",
	width = 12, height = 10, units = "cm", pointsize = 12, res = 300)



  par(mar=c(4,4.5,4,1.5))
  par(mgp=c(2.5,1,0))
  
 L = 30

  ag <-seq(from=min(model.data$ag_suit[which(model.data$Within_PA == "yes")]),
	to=max(model.data$ag_suit[which(model.data$Within_PA == "yes")]),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  





model.data$Zone<- relevel(model.data$Zone , "Tropical")
model.data$Within_PA<- relevel(model.data$Within_PA , "yes")

mam <- glmer(Richness_rarefied.model.int$final.call, model.data, family = "poisson")

 for (x in 1:L)
  {
    log_bound_dist_km_PA_neg <-  mean(model.data$log_bound_dist_km_PA_neg)
    Richness_rarefied <-0
    log_hpd <- mean(model.data$log_hpd)
    log_AREA.PA <- mean(model.data$log_AREA.PA)
    DoP.PA <- mean(model.data$DoP.PA)
    IUCN.PA <- "1.5"
    Zone <-"Tropical"
    taxon_of_interest<- "Vertebrates"
    Predominant_habitat <-"Primary Vegetation"
    ag_suit <- ag[x]
    log_access <- mean(model.data$log_access)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, Richness_rarefied, log_hpd,log_AREA.PA, DoP.PA, log_access, ag_suit))
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


  plot(ag,z, ylim=c(0,10), col =inside.col, lwd = 1, 
		bty = "l", 
		type = "l",ylab = "Rarefied species richness per site ± s.e", xlab="Agricultural suitability (higher = more suitable)")
#  rug(data$ag_suit[which(model.data$Within_PA == "yes")]
#	, ticksize = 0.03, side = 1, lwd = 1.5, col = inside.col)
  points(ag,zu,type="l",lty=2,  lwd = 1, col =  inside.col)
  points(ag,zl,type="l",lty=2,  lwd = 1, col =  inside.col)



  ag <-seq(from=min(model.data$ag_suit[which(model.data$Within_PA == "no")]),
	to=max(model.data$ag_suit[which(model.data$Within_PA == "no")]),length=L)
  y<-vector(mode="numeric",length=L)
  yu<-vector(mode="numeric",length=L)
  yl<-vector(mode="numeric",length=L)
  





model.data$Zone<- relevel(model.data$Zone , "Tropical")
model.data$Within_PA<- relevel(model.data$Within_PA , "no")

mam <- glmer(Richness_rarefied.model.int$final.call, model.data, family = "poisson")

 for (x in 1:L)
  {
    log_bound_dist_km_PA_neg <-  mean(model.data$log_bound_dist_km_PA_neg)
    Richness_rarefied <-0
    log_hpd <- mean(model.data$log_hpd)
    log_AREA.PA <- mean(model.data$log_AREA.PA)
    DoP.PA <- mean(model.data$DoP.PA)
    IUCN.PA <- "1.5"
    Zone <-"Tropical"
    taxon_of_interest<- "Vertebrates"
    Predominant_habitat <-"Primary Vegetation"
    ag_suit <- ag[x]
    log_access <- mean(model.data$log_access)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, Richness_rarefied, log_hpd,log_AREA.PA, DoP.PA, log_access, ag_suit))
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


#  rug(data$ag_suit[which(model.data$Within_PA == "no")]
#	, ticksize = 0.03, side = 1, lwd = 1.5, col = outside.col)
  points(ag,y,type="l",lty=2,  lwd = 1, col =  outside.col)
  points(ag,yu,type="l",lty=2,  lwd = 1, col =  outside.col)
  points(ag,yl,type="l",lty=2,  lwd = 1, col =  outside.col)

legend("topleft", lty = 1,  c("Protected", "Unprotected"), col = c(inside.col, outside.col))


dev.off()





###### PA area ##########

tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/11_14/Richness_rarefied vs PA size.tif",
	width = 12, height = 10, units = "cm", pointsize = 12, res = 300)



  par(mar=c(4,4.5,4,1.5))
  par(mgp=c(2.5,1,0))
  
 L = 30

  gis <-seq(from=min(model.data$log_AREA.PA),to=max(model.data$log_AREA.PA),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  


model.data$Zone<- relevel(model.data$Zone , "Tropical")
#mam <- glmer(Richness_rarefied.model.int$final.call, model.data, family = "poisson")

 for (x in 1:L)
  {
    log_bound_dist_km_PA_neg <-  mean(model.data$log_bound_dist_km_PA_neg)
    Richness_rarefied <-0
    log_hpd <- mean(model.data$log_hpd)
    log_AREA.PA <- gis[x]
    DoP.PA <- mean(model.data$DoP.PA)
    IUCN.PA <- "1.5"
    Zone <-"Tropical"
    taxon_of_interest<- "Invertebrates"
    Predominant_habitat <-"Primary Vegetation"
    ag_suit <- median(model.data$ag_suit)
    log_access <- mean(model.data$log_access)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, Richness_rarefied, log_hpd,log_AREA.PA, DoP.PA, log_access, ag_suit))
    newdat.f<-data.frame(cbind( Zone,taxon_of_interest, Predominant_habitat, Within_PA))
    levels(newdat.f$Zone)<-levels(model.data$Zone)
    levels(newdat.f$taxon_of_interest)<-levels(model.data$taxon_of_interest)
    levels(newdat.f$Predominant_habitat)<-levels(model.data$Predominant_habitat)
    levels(newdat.f$ag_suit) <- levels(model.data$ag_suit)
   # levels(newdat.f$IUCN.PA) <- levels(model.data$IUCN.PA)
    levels(newdat.f$Within_PA) <- levels(model.data$Within_PA)
    newdat<-data.frame(newdat.f,newdat.n)
    mm<-model.matrix(terms(mam),newdat)
    pvar1 <- diag(mm %*% tcrossprod(vcov(mam),mm))
    z[x]<-mm %*% fixef(mam)
    zu[x]<-z[x]+sqrt(pvar1)
    zl[x]<-z[x]-sqrt(pvar1) 
  }
  
  gis <- exp(gis) - 1

  z <- exp(z)
  zu <- exp(zu)
  zl <- exp(zl)



  plot(gis,z, ylim=c(0,10), col = 1, lwd = 1,
		bty = "l",  log = "x",xaxt = "n",
		type = "l",ylab = "Rarefied species richness per site ± s.e", xlab="PA size (km2)")
  rug(data$AREA.PA, ticksize = 0.03, side = 1, lwd = 1, col = 1)
  axis(1, c(0.1, 1, 10, 100, 1000),  c(0.1, 1, 10, 100, 1000))
  points(gis,zu,type="l",lty=2,  lwd = 1, col =  1)
  points(gis,zl,type="l",lty=2,  lwd = 1, col =  1)



dev.off()





###### DoP vs taxon of interest ##########


tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/11_14/Richness_rarefied vs DoP vs taxon .tif",
	width = 12, height = 25, units = "cm", pointsize = 12, res = 300)


taxa <- levels(model.data$taxon_of_interest)
 t <- taxa[2]


par(mfrow = c(3,1))

for(t in taxa){ 


  par(mar=c(4,4.5,4,1.5))
  par(mgp=c(2.5,1,0))
  
 L = 30

  dop <-seq(from=min(model.data$DoP.PA[which(model.data$taxon_of_interest == t)])
		,to=max(model.data$DoP.PA[which(model.data$taxon_of_interest == t)]),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  


model.data$Zone<- relevel(model.data$Zone , "Tropical")
model.data$taxon_of_interest <- relevel(model.data$taxon_of_interest , t)
mam <- glmer(Richness_rarefied.model.int$final.call, model.data, family = "poisson")

 for (x in 1:L)
  {
    log_bound_dist_km_PA_neg <-  mean(model.data$log_bound_dist_km_PA_neg)
    Richness_rarefied <-0
    log_hpd <- mean(model.data$log_hpd)
    log_AREA.PA <- mean(model.data$log_AREA.PA)
    DoP.PA <- dop[x]
    IUCN.PA <- "1.5"
    Zone <-"Tropical"
    taxon_of_interest<- "Invertebrates"
    Predominant_habitat <-"Primary Vegetation"
    ag_suit <- median(model.data$ag_suit)
    log_access <- mean(model.data$log_access)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, Richness_rarefied, log_hpd,log_AREA.PA, DoP.PA, log_access, ag_suit))
    newdat.f<-data.frame(cbind( Zone,taxon_of_interest, Predominant_habitat, Within_PA))
    levels(newdat.f$Zone)<-levels(model.data$Zone)
    levels(newdat.f$taxon_of_interest)<-levels(model.data$taxon_of_interest)
    levels(newdat.f$Predominant_habitat)<-levels(model.data$Predominant_habitat)
    levels(newdat.f$ag_suit) <- levels(model.data$ag_suit)
   # levels(newdat.f$IUCN.PA) <- levels(model.data$IUCN.PA)
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



  plot(dop,z, ylim=c(0,17), xlim = c(0,81),
		col = taxa.cols[which(taxa == t)], lwd = 2, main =t,
		bty = "l",  
		type = "l",ylab = "Rarefied species richness per site ± s.e", xlab="Duration of protection (yr)")
  rug(model.data$DoP.PA[which(model.data$taxon_of_interest == t)], ticksize = 0.03, side = 1, lwd = 1.5, col = taxa.cols[which(taxa == t)])
  points(dop,zu,type="l",lty=2,  lwd = 2, col =  taxa.cols[which(taxa == t)])
  points(dop,zl,type="l",lty=2,  lwd = 2.5, col =  taxa.cols[which(taxa == t)])


}




dev.off()












####################### OLDER PLOTS ###################

###### DoP ##########


tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/Richness_rarefied vs DoP.tif",
	width = 12, height = 10, units = "cm", pointsize = 12, res = 300)



  par(mar=c(4,4.5,4,1.5))
  par(mgp=c(2.5,1,0))
  
 L = 30

  dop <-seq(from=min(model.data$DoP.PA),to=max(model.data$DoP.PA),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  


model.data$Zone<- relevel(model.data$Zone , "Tropical")
mam <- glmer(Richness_rarefied.model.int$final.call, model.data, family = "poisson")

 for (x in 1:L)
  {
    log_bound_dist_km_PA_neg <-  mean(model.data$log_bound_dist_km_PA_neg)
    Richness_rarefied <-0
    log_hpd <- mean(model.data$log_hpd)
    log_AREA.PA <- mean(model.data$log_AREA.PA)
    DoP.PA <- dop[x]
    IUCN.PA <- "1.5"
    Zone <-"Tropical"
    taxon_of_interest<- "Invertebrates"
    Predominant_habitat <-"Primary Vegetation"
    ag_suit <- median(model.data$ag_suit)
    log_access <- mean(model.data$log_access)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, Richness_rarefied, log_hpd,log_AREA.PA, DoP.PA, log_access, ag_suit))
    newdat.f<-data.frame(cbind( Zone,taxon_of_interest, Predominant_habitat, Within_PA))
    levels(newdat.f$Zone)<-levels(model.data$Zone)
    levels(newdat.f$taxon_of_interest)<-levels(model.data$taxon_of_interest)
    levels(newdat.f$Predominant_habitat)<-levels(model.data$Predominant_habitat)
    levels(newdat.f$ag_suit) <- levels(model.data$ag_suit)
   # levels(newdat.f$IUCN.PA) <- levels(model.data$IUCN.PA)
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



  plot(dop,z, ylim=c(0,17), col = 1, lwd = 2,
		bty = "l",  
		type = "l",ylab = "Rarefied species richness per site ± s.e", xlab="Duration of protection (yr)")
  rug(data$DoP.PA, ticksize = 0.03, side = 1, lwd = 1.5, col = 1)
  points(dop,zu,type="l",lty=2,  lwd = 2, col =  1)
  points(dop,zl,type="l",lty=2,  lwd = 2.5, col =  1)



dev.off()












##### PA and non PA ######
### different taxa #### 



tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/Richness_rarefied vs PA size.tif",
	width = 12, height = 10, units = "cm", pointsize = 12, res = 300)


PA <- vector()
taxon <- vector()
estimate <- vector()
se <- vector()

model.data$taxon_of_interest <- relevel(model.data$taxon_of_interest , "Invertebrates")
model.data$Within_PA <- relevel(model.data$Within_PA, "yes")
mam <- glmer(Richness_rarefied.model.int$final.call, model.data, family = "poisson")
estimate[1] <- fixef(mam)[1]
se[1] <- se.fixef(mam)[1]
PA[1] <- "Protected"
taxon[1] <- "Invertebrates"


model.data$Within_PA <- relevel(model.data$Within_PA, "no")
mam <- glmer(Richness_rarefied.model.int$final.call, model.data, family = "poisson")
estimate[2] <- fixef(mam)[1]
se[2] <- se.fixef(mam)[1]
PA[2] <- "Unprotected"
taxon[2] <- "Invertebrates"

model.data$taxon_of_interest <- relevel(model.data$taxon_of_interest , "Vertebrates")
model.data$Within_PA <- relevel(model.data$Within_PA, "yes")
mam <- glmer(Richness_rarefied.model.int$final.call, model.data, family = "poisson")
estimate[3] <- fixef(mam)[1]
se[3] <- se.fixef(mam)[1]
PA[3] <- "Protected"
taxon[3] <- "Vertebrates"


model.data$Within_PA <- relevel(model.data$Within_PA, "no")
mam <- glmer(Richness_rarefied.model.int$final.call, model.data, family = "poisson")
estimate[4] <- fixef(mam)[1]
se[4] <- se.fixef(mam)[1]
PA[4] <- "Unprotected"
taxon[4] <- "Vertebrates"


model.data$taxon_of_interest <- relevel(model.data$taxon_of_interest , "Plants and Fungi")
model.data$Within_PA <- relevel(model.data$Within_PA, "yes")
mam <- glmer(Richness_rarefied.model.int$final.call, model.data, family = "poisson")
estimate[5] <- fixef(mam)[1]
se[5] <- se.fixef(mam)[1]
PA[5] <- "Protected"
taxon[5] <- "Plants and Fungi"


multiple.taxa.matched.landuse$Within_PA <- relevel(multiple.taxa.matched.landuse$Within_PA, "no")
mam <- glmer(Richness_rarefied.model.int$final.call, model.data, family = "poisson")
estimate[6] <- fixef(mam)[1]
se[6] <- se.fixef(mam)[1]
PA[6] <- "Unprotected"
taxon[6] <- "Plants and Fungi"



plot.data <- data.frame(PA = PA, taxon = taxon, se = se, estimate = estimate)
plot.data$upper <- plot.data$estimate + 1.96*plot.data$se
plot.data$lower <- plot.data$estimate - 1.96*plot.data$se
pos <- c(1:6)

plot(plot.data$estimate ~ pos, pch = c(16,21,16,21), cex = 2,
	ylim = c(0, 6), 
	bty = "l", xaxt = "n", xlab = "", ylab = "Log Species richness per site")
arrows(pos, plot.data$lower, pos, plot.data$upper, angle = 90, code = 3, length = 0.05)
axis(1, c(1.5, 3.5, 5.5), c("Invertebrates", "Vertebrates", "Plants and Fungi"))
legend("topleft", pch = c(16, 21), pt.cex = 2, c("Protected", "Unprotected"), bty = "n")

dev.off()









###### PA size in different zones  #############
###### for different taxa ##################

tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/Richness_rarefied vs zone vs PA size vs taxon.tif",
	width = 12, height = 10, units = "cm", pointsize = 12, res = 300)



  par(mar=c(4,4.5,4,1.5))
  par(mgp=c(2.5,1,0))
  
 L = 30

  gis <-seq(from=min(model.data$log_AREA.PA),to=max(model.data$log_AREA.PA),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  

par(mfrow = c(3,1))
taxa <- levels(matched.landuse$taxon_of_interest)
i <- 0
t <- taxa[1]

for(t in taxa){

model.data$taxon_of_interest <- relevel(model.data$taxon_of_interest, t)
model.data$Zone<- relevel(model.data$Zone , "Tropical")
mam <- glmer(Richness_rarefied.model$final.call, model.data, family = "poisson")

 for (x in 1:L)
  {
    log_bound_dist_km_PA_neg <-  mean(model.data$log_bound_dist_km_PA_neg)
    Richness_rarefied <-0
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
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, Richness_rarefied, log_hpd,log_AREA.PA, DoP.PA, log_access, ag_suit))
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

 gis_plot <- exp(gis) -1


  plot(gis_plot,z, ylim=c(0,25), col = zone.cols[1], lwd = 2, main = t,
		bty = "l", log = "x", #yaxt = "n", 
		type = "l",ylab = "Richness_rarefied size Log10(sq km) ± s.e", xlab="PA size (km2)")
  rug(data$AREA.PA[which(data$Zone == "Tropical" & data$taxon_of_interest == t)]
	, ticksize = 0.03, side = 1, lwd = 1.5, col = zone.cols[1])
  points(gis_plot,zu,type="l",lty=2,  lwd = 2, col =  zone.cols[1])
  points(gis_plot,zl,type="l",lty=2,  lwd = 2.5, col =  zone.cols[1])



  y<-vector(mode="numeric",length=L)
  yu<-vector(mode="numeric",length=L)
  yl<-vector(mode="numeric",length=L)
  

model.data$taxon_of_interest <- relevel(model.data$taxon_of_interest, t)
model.data$Zone<- relevel(model.data$Zone , "Temperate")
mam <- glmer(Richness_rarefied.model$final.call, model.data, family = "poisson")

 for (x in 1:L)
  {
    log_bound_dist_km_PA_neg <-  mean(model.data$log_bound_dist_km_PA_neg)
    Richness_rarefied <-0
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
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, Richness_rarefied, log_hpd,log_AREA.PA, DoP.PA, log_access, ag_suit))
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
    y[x]<-mm %*% fixef(mam)
    yu[x]<-y[x]+sqrt(pvar1)
    yl[x]<-y[x]-sqrt(pvar1) 
  }
  

  y <- exp(y)
  yu <- exp(yu)
  yl <- exp(yl)


  points(gis_plot,y,type="l",lty=1, col =  zone.cols[2])
  points(gis_plot,yu,type="l",lty=2, col =  zone.cols[2])
  points(gis_plot,yl,type="l",lty=2, col =  zone.cols[2])
  rug(data$AREA.PA[which(data$Zone == "Temperate" & data$taxon_of_interest == t)],
		col = zone.cols[2], lwd = 0.8, pos = 0)

}

legend("topright", c("Tropical", "Temperate") , col = zone.cols, lty = 1, lwd = c(2,1))

dev.off()




### ag suit ####






tiff( "N:/Documents/PREDICTS/WDPA analysis/plots/Richness_rarefied vs ag suit.tif",
	width = 12, height = 10, units = "cm", pointsize = 12, res = 300)



  par(mar=c(4,4.5,4,1.5))
  par(mgp=c(2.5,1,0))
  
 L = 30

  ag <-seq(from=min(model.data$ag_suit),to=max(model.data$ag_suit),length=L)
  z<-vector(mode="numeric",length=L)
  zu<-vector(mode="numeric",length=L)
  zl<-vector(mode="numeric",length=L)
  




model.data$Zone<- relevel(model.data$Zone , "Tropical")
mam <- glmer(Richness_rarefied.model.int$final.call, model.data, family = "poisson")

 for (x in 1:L)
  {
    log_bound_dist_km_PA_neg <-  mean(model.data$log_bound_dist_km_PA_neg)
    Richness_rarefied <-0
    log_hpd <- mean(model.data$log_hpd)
    log_AREA.PA <- mean(model.data$log_AREA.PA)
    DoP.PA <- mean(model.data$DoP.PA)
    IUCN.PA <- "1.5"
    Zone <-"Tropical"
    taxon_of_interest<- t
    Predominant_habitat <-"Primary Vegetation"
    ag_suit <- ag[x]
    log_access <- mean(model.data$log_access)
    Within_PA <- "no"
    newdat.n<-data.frame(cbind(log_bound_dist_km_PA_neg, Richness_rarefied, log_hpd,log_AREA.PA, DoP.PA, log_access, ag_suit))
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


  plot(ag,z, ylim=c(0,25), col =1, lwd = 1, 
		bty = "l", 
		type = "l",ylab = "Rarefied species richness per site ± s.e", xlab="Agricultural suitability (higher = more suitable)")
  rug(data$ag_suit
	, ticksize = 0.03, side = 1, lwd = 1.5, col = 1)
  points(ag,zu,type="l",lty=2,  lwd = 1, col =  1)
  points(ag,zl,type="l",lty=2,  lwd = 1, col =  1)


dev.off()

