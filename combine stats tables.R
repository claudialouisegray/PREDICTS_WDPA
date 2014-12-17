
rm(list=ls()) 


setwd("N:/Documents/PREDICTS/WDPA analysis/stats tables")


files <- list.files()

# all the models that have 2 in the title are the within_PA analyses - belong together
dist.files <- files[grep("model.stats",files)]
pa.files <- files[grep("model2",files)]

# read in dist to boundary analyses stats
 dist.stats <- lapply(dist.files, read.csv) 
 names(dist.stats) <- gsub(".csv","",dist.files)

# read in pa analyses stats
 pa.stats <- lapply(pa.files, read.csv) 
 names(pa.stats) <- gsub(".csv","",pa.files)





# get all terms names 
 dist.terms <- lapply(dist.stats, function(dtf)
			unique(dtf$terms)
			)
 dist.terms <- unique(unlist(dist.terms))
 dist.terms <- sort(dist.terms)

 pa.terms <- lapply(pa.stats, function(dtf)
			unique(dtf$terms)
			)
 pa.terms <- unique(unlist(pa.terms))
 pa.terms <- sort(pa.terms)
# zones not together but sort this out later




### dist to boundary results ###

# build dataframe to hold results

dist.res <- data.frame(terms = dist.terms,
				sp.rich.ChiSq = rep(NA, length(dist.terms)),
				sp.rich.Df = rep(NA, length(dist.terms)),
				sp.rich.P = rep(NA, length(dist.terms)),
				rar.ChiSq = rep(NA, length(dist.terms)),
				rar.Df = rep(NA, length(dist.terms)),
				rar.P = rep(NA, length(dist.terms)),
				ab.ChiSq = rep(NA, length(dist.terms)),
				ab.Df = rep(NA, length(dist.terms)),
				ab.P = rep(NA, length(dist.terms)),
				e.ChiSq = rep(NA, length(dist.terms)),
				e.Df = rep(NA, length(dist.terms)),
				e.P = rep(NA, length(dist.terms)),
				pt.ChiSq = rep(NA, length(dist.terms)),
				pt.Df = rep(NA, length(dist.terms)),
				pt.P = rep(NA, length(dist.terms))
	)

 for(x in 1:length(dist.stats)){
	for(i in 1:nrow(dist.stats[[x]])){
		if(length(grep("Species_richness", names(dist.stats[x]))) >0){
			dist.res$sp.rich.ChiSq[which(dist.res$terms == paste(dist.stats[[x]]$terms[i]))]<- dist.stats[[x]]$ChiSq[i]
			dist.res$sp.rich.Df[which(dist.res$terms == paste(dist.stats[[x]]$terms[i]))] <- paste(dist.stats[[x]]$Df[i])
			dist.res$sp.rich.P[which(dist.res$terms == paste(dist.stats[[x]]$terms[i]))] <- dist.stats[[x]]$P[i]
		} else if(length(grep("rarefied",names(dist.stats[x]))) >0){
			dist.res$rar.ChiSq[which(dist.res$terms == paste(dist.stats[[x]]$terms[i]))] <- dist.stats[[x]]$ChiSq[i]
			dist.res$rar.Df[which(dist.res$terms == paste(dist.stats[[x]]$terms[i]))] <- paste(dist.stats[[x]]$Df[i])
			dist.res$rar.P[which(dist.res$terms == paste(dist.stats[[x]]$terms[i]))] <- dist.stats[[x]]$P[i]
		} else if(length(grep("ab",names(dist.stats[x]))) >0){
			dist.res$ab.ChiSq[which(dist.res$terms == paste(dist.stats[[x]]$terms[i]))] <- dist.stats[[x]]$ChiSq[i]
			dist.res$ab.Df[which(dist.res$terms == paste(dist.stats[[x]]$terms[i]))] <- paste(dist.stats[[x]]$Df[i])
			dist.res$ab.P[which(dist.res$terms == paste(dist.stats[[x]]$terms[i]))] <- dist.stats[[x]]$P[i]
		} else if(length(grep("range",names(dist.stats[x]))) >0){
			dist.res$e.ChiSq[which(dist.res$terms == paste(dist.stats[[x]]$terms[i]))] <- dist.stats[[x]]$ChiSq[i]
			dist.res$e.Df[which(dist.res$terms == paste(dist.stats[[x]]$terms[i]))] <- paste(dist.stats[[x]]$Df[i])
			dist.res$e.P[which(dist.res$terms == paste(dist.stats[[x]]$terms[i]))] <- dist.stats[[x]]$P[i]
		} else if(length(grep("RLS",names(dist.stats[x]))) >0){
			dist.res$pt.ChiSq[which(dist.res$terms == paste(dist.stats[[x]]$terms[i]))] <- dist.stats[[x]]$ChiSq[i]
			dist.res$pt.Df[which(dist.res$terms == paste(dist.stats[[x]]$terms[i]))] <- paste(dist.stats[[x]]$Df[i])
			dist.res$pt.P[which(dist.res$terms == paste(dist.stats[[x]]$terms[i]))] <- dist.stats[[x]]$P[i]
		}	
	}
 }




### Within PA results ###

# build dataframe to hold results

pa.res <- data.frame(terms = pa.terms,
				sp.rich.ChiSq = rep(NA, length(pa.terms)),
				sp.rich.Df = rep(NA, length(pa.terms)),
				sp.rich.P = rep(NA, length(pa.terms)),
				rar.ChiSq = rep(NA, length(pa.terms)),
				rar.Df = rep(NA, length(pa.terms)),
				rar.P = rep(NA, length(pa.terms)),
				ab.ChiSq = rep(NA, length(pa.terms)),
				ab.Df = rep(NA, length(pa.terms)),
				ab.P = rep(NA, length(pa.terms)),
				e.ChiSq = rep(NA, length(pa.terms)),
				e.Df = rep(NA, length(pa.terms)),
				e.P = rep(NA, length(pa.terms)),
				pt.ChiSq = rep(NA, length(pa.terms)),
				pt.Df = rep(NA, length(pa.terms)),
				pt.P = rep(NA, length(pa.terms))
	)


x <- 5

 for(x in 1:length(pa.stats)){
	for(i in 1:nrow(pa.stats[[x]])){
		if(length(grep("Species_richness", names(pa.stats[x]))) >0){
			pa.res$sp.rich.ChiSq[which(pa.res$terms == paste(pa.stats[[x]]$terms[i]))]<- pa.stats[[x]]$ChiSq[i]
			pa.res$sp.rich.Df[which(pa.res$terms == paste(pa.stats[[x]]$terms[i]))] <- paste(pa.stats[[x]]$Df[i])
			pa.res$sp.rich.P[which(pa.res$terms == paste(pa.stats[[x]]$terms[i]))] <- pa.stats[[x]]$P[i]
		} else if(length(grep("rarefied",names(pa.stats[x]))) >0){
			pa.res$rar.ChiSq[which(pa.res$terms == paste(pa.stats[[x]]$terms[i]))] <- pa.stats[[x]]$ChiSq[i]
			pa.res$rar.Df[which(pa.res$terms == paste(pa.stats[[x]]$terms[i]))] <- paste(pa.stats[[x]]$Df[i])
			pa.res$rar.P[which(pa.res$terms == paste(pa.stats[[x]]$terms[i]))] <- pa.stats[[x]]$P[i]
		} else if(length(grep("ab",names(pa.stats[x]))) >0){
			pa.res$ab.ChiSq[which(pa.res$terms == paste(pa.stats[[x]]$terms[i]))] <- pa.stats[[x]]$ChiSq[i]
			pa.res$ab.Df[which(pa.res$terms == paste(pa.stats[[x]]$terms[i]))] <- paste(pa.stats[[x]]$Df[i])
			pa.res$ab.P[which(pa.res$terms == paste(pa.stats[[x]]$terms[i]))] <- pa.stats[[x]]$P[i]
		} else if(length(grep("range",names(pa.stats[x]))) >0){
			pa.res$e.ChiSq[which(pa.res$terms == paste(pa.stats[[x]]$terms[i]))] <- pa.stats[[x]]$ChiSq[i]
			pa.res$e.Df[which(pa.res$terms == paste(pa.stats[[x]]$terms[i]))] <- paste(pa.stats[[x]]$Df[i])
			pa.res$e.P[which(pa.res$terms == paste(pa.stats[[x]]$terms[i]))] <- pa.stats[[x]]$P[i]
		} else if(length(grep("RLS",names(pa.stats[x]))) >0){
			pa.res$pt.ChiSq[which(pa.res$terms == paste(pa.stats[[x]]$terms[i]))] <- pa.stats[[x]]$ChiSq[i]
			pa.res$pt.Df[which(pa.res$terms == paste(pa.stats[[x]]$terms[i]))] <- paste(pa.stats[[x]]$Df[i])
			pa.res$pt.P[which(pa.res$terms == paste(pa.stats[[x]]$terms[i]))] <- pa.stats[[x]]$P[i]
		}	
	}
 }


# add in info on which tests could not be done

# all taxon information for prop threatened

#distance to boundary results
dist.res$pt.ChiSq[grep("taxon_of_interest", dist.res$terms)] <- "NOT TESTED"
dist.res$pt.Df[grep("taxon_of_interest", dist.res$terms)] <- "NOT TESTED"
dist.res$pt.P[grep("taxon_of_interest", dist.res$terms)] <- "NOT TESTED"

#Within_PA results
pa.res$pt.ChiSq[grep("taxon_of_interest", pa.res$terms)] <- "NOT TESTED"
pa.res$pt.Df[grep("taxon_of_interest", pa.res$terms)] <- "NOT TESTED"
pa.res$pt.P[grep("taxon_of_interest", pa.res$terms)] <- "NOT TESTED"


# cubic and quadratic polynomials for ag_suit, elevation and slope in proportion threated Within_PA analysis
pa.res$pt.ChiSq[grep("poly(log_elevation,2)", pa.res$terms)] <- "NOT TESTED"
pa.res$pt.Df[grep("poly(log_elevation,2)", pa.res$terms)] <- "NOT TESTED"
pa.res$pt.P[grep("poly(log_elevation,2)", pa.res$terms)] <- "NOT TESTED"
pa.res$pt.ChiSq[grep("poly(log_elevation,3)", pa.res$terms)] <- "NOT TESTED"
pa.res$pt.Df[grep("poly(log_elevation,3)", pa.res$terms)] <- "NOT TESTED"
pa.res$pt.P[grep("poly(log_elevation,3)", pa.res$terms)] <- "NOT TESTED"

pa.res$pt.ChiSq[grep("poly(log_slope,2)", pa.res$terms)] <- "NOT TESTED"
pa.res$pt.Df[grep("poly(log_slope,2)", pa.res$terms)] <- "NOT TESTED"
pa.res$pt.P[grep("poly(log_slope,2)", pa.res$terms)] <- "NOT TESTED"
pa.res$pt.ChiSq[grep("poly(log_slope,3)", pa.res$terms)] <- "NOT TESTED"
pa.res$pt.Df[grep("poly(log_slope,3)", pa.res$terms)] <- "NOT TESTED"
pa.res$pt.P[grep("poly(log_slope,3)", pa.res$terms)] <- "NOT TESTED"

pa.res$pt.ChiSq[grep("poly(ag_suit,2)", pa.res$terms)] <- "NOT TESTED"
pa.res$pt.Df[grep("poly(ag_suit,2)", pa.res$terms)] <- "NOT TESTED"
pa.res$pt.P[grep("poly(ag_suit,2)", pa.res$terms)] <- "NOT TESTED"
pa.res$pt.ChiSq[grep("poly(ag_suit,3)", pa.res$terms)] <- "NOT TESTED"
pa.res$pt.Df[grep("poly(ag_suit,3)", pa.res$terms)] <- "NOT TESTED"
pa.res$pt.P[grep("poly(ag_suit,3)", pa.res$terms)] <- "NOT TESTED"


# interactions that didnt work for prop threatened Within_PA analysis

#"Within_PA:Zone"
pa.res$pt.ChiSq[which(pa.res$terms == "Within_PA:Zone")] <- "NOT TESTED"
pa.res$pt.Df[which(pa.res$terms == "Within_PA:Zone")] <- "NOT TESTED"
pa.res$pt.P[which(pa.res$terms == "Within_PA:Zone")] <- "NOT TESTED"

#"Zone:poly(log_AREA.PA,1)"
pa.res$pt.ChiSq[which(pa.res$terms == "Zone:poly(log_AREA.PA,1)")] <- "NOT TESTED"
pa.res$pt.Df[which(pa.res$terms == "Zone:poly(log_AREA.PA,1)")] <- "NOT TESTED"
pa.res$pt.P[which(pa.res$terms == "Zone:poly(log_AREA.PA,1)")] <- "NOT TESTED"
pa.res$pt.ChiSq[which(pa.res$terms == "Zone:poly(log_AREA.PA,2)")] <- "NOT TESTED"
pa.res$pt.Df[which(pa.res$terms == "Zone:poly(log_AREA.PA,2)")] <- "NOT TESTED"
pa.res$pt.P[which(pa.res$terms == "Zone:poly(log_AREA.PA,2)")] <- "NOT TESTED"
pa.res$pt.ChiSq[which(pa.res$terms == "Zone:poly(log_AREA.PA,3)")] <- "NOT TESTED"
pa.res$pt.Df[which(pa.res$terms == "Zone:poly(log_AREA.PA,3)")] <- "NOT TESTED"
pa.res$pt.P[which(pa.res$terms == "Zone:poly(log_AREA.PA,3)")] <- "NOT TESTED"


# interactions that didnt work for sp rich Within_PA analysis

#"Zone:poly(log_AREA.PA,3)"
pa.res$sp.rich.ChiSq[which(pa.res$terms == "Zone:poly(log_AREA.PA,3)")] <- "NOT TESTED"
pa.res$sp.rich.Df[which(pa.res$terms == "Zone:poly(log_AREA.PA,3)")] <- "NOT TESTED"
pa.res$sp.rich.P[which(pa.res$terms == "Zone:poly(log_AREA.PA,3)")] <- "NOT TESTED"
pa.res$sp.rich.ChiSq[which(pa.res$terms == "Zone:poly(log_AREA.PA,2)")] <- "NOT TESTED"
pa.res$sp.rich.Df[which(pa.res$terms == "Zone:poly(log_AREA.PA,2)")] <- "NOT TESTED"
pa.res$sp.rich.P[which(pa.res$terms == "Zone:poly(log_AREA.PA,2)")] <- "NOT TESTED"

#"taxon_of_interest:poly(log_AREA.PA,3)"
pa.res$sp.rich.ChiSq[which(pa.res$terms == "taxon_of_interest:poly(log_AREA.PA,3)")] <- "NOT TESTED"
pa.res$sp.rich.Df[which(pa.res$terms == "taxon_of_interest:poly(log_AREA.PA,3)")] <- "NOT TESTED"
pa.res$sp.rich.P[which(pa.res$terms == "taxon_of_interest:poly(log_AREA.PA,3)")] <- "NOT TESTED"
pa.res$sp.rich.ChiSq[which(pa.res$terms == "taxon_of_interest:poly(log_AREA.PA,2)")] <- "NOT TESTED"
pa.res$sp.rich.Df[which(pa.res$terms == "taxon_of_interest:poly(log_AREA.PA,2)")] <- "NOT TESTED"
pa.res$sp.rich.P[which(pa.res$terms == "taxon_of_interest:poly(log_AREA.PA,2)")] <- "NOT TESTED"


#check rarefied richness




# sort out the order so Zones are together
 pa.res$terms <- gsub("Z", "z", pa.terms)
 pa.res <- pa.res[order(pa.res$terms),]

 dist.res$terms <- gsub("Z", "z", dist.terms)
 dist.res <- dist.res[order(dist.res$terms),]



# save files

setwd("N:/Documents/PREDICTS/WDPA analysis")

write.csv(pa.res, "pa.results.table.csv")
write.csv(dist.res, "dist.results.table.csv")

