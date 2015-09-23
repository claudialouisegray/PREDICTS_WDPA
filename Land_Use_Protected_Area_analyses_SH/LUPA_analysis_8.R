library(roquefort)
##LUPA analysis 8 level

#abundance
CompareRandoms(dataset=data,responseVar="LogAbund",fitFamily="gaussian",fixedFactors=c("Predominant_habitat","Within_PA", "Zone","taxon_of_interest"),fixedTerms=c(list(ag_suit=3), list(log_slope=3), list(log_elevation=3)),fixedInteractions="Predominant_habitat:Within_PA",siteRandom=F,fitInteractions=F,verbose=T,optimizer="bobyqa",randomSlopes=T)
ModelSelect(all.data=data,responseVar="LogAbund",fitFamily="gaussian", fixedFactors=c("Predominant_habitat","Within_PA", "Zone","taxon_of_interest"), fixedInteractions="Predominant_habitat:Within_PA",randomStruct="(1+Within_PA|SS)+(1|SSB)",siteRandom=FALSE,fitInteractions=F, verbose=T,fixedTerms=c(list(ag_suit=3), list(log_slope=3), list(log_elevation=3)))
FixedModel(all.data=data,responseVar="LogAbund",fitFamily="gaussian", fixedStruct="Predominant_habitat + Zone+taxon_of_interest + Within_PA+Predominant_habitat:Within_PA", randomStruct="(1+Within_PA|SS)+(1|SSB)",REML=T,optimizer="bobyqa")

#range size
CompareRandoms(dataset=data,responseVar="CWM_Geographic_range_log10_square_km",fitFamily="gaussian",fixedFactors=c("Predominant_habitat","Within_PA", "Zone","taxon_of_interest"),fixedTerms=c(list(ag_suit=3), list(log_slope=3), list(log_elevation=3)),fixedInteractions="Predominant_habitat:Within_PA",siteRandom=F,fitInteractions=F,verbose=T,optimizer="bobyqa",randomSlopes=T)
ModelSelect(all.data=data,responseVar="CWM_Geographic_range_log10_square_km",fitFamily="gaussian", fixedFactors=c("Predominant_habitat","Within_PA", "Zone","taxon_of_interest"), fixedInteractions=c("Predominant_habitat:Within_PA"),randomStruct="(1+Within_PA|SS)+(1|SSB)",siteRandom=FALSE,fitInteractions=F, verbose=T,fixedTerms=c(list(ag_suit=3), list(log_slope=3), list(log_elevation=3)))
FixedModel(all.data=data,responseVar="CWM_Geographic_range_log10_square_km",fitFamily="gaussian", fixedStruct="Predominant_habitat + Zone+taxon_of_interest + Within_PA+Predominant_habitat:Within_PA+ poly(log_slope,3) + poly(log_elevation,2) + poly(ag_suit,3)", randomStruct="(1+Within_PA|SS)+(1|SSB)",REML=T,optimizer="bobyqa")

#species richness
CompareRandoms(dataset=data.spr.rr,responseVar="Species_richness",fitFamily="poisson",fixedFactors=c("Predominant_habitat","Within_PA", "Zone","taxon_of_interest"),fixedTerms=c(list(ag_suit=3), list(log_slope=3), list(log_elevation=3)),fixedInteractions="Predominant_habitat:Within_PA",siteRandom=T,fitInteractions=F,verbose=T,optimizer="bobyqa",randomSlopes=T)
ModelSelect(all.data=data.spr.rr,responseVar="Species_richness",fitFamily="poisson", fixedFactors=c("Predominant_habitat","Within_PA", "Zone","taxon_of_interest"), fixedInteractions=c("Predominant_habitat:Within_PA"),randomStruct="(1+Within_PA|SS)+(1|SSB)+(1|SSBS)",siteRandom=T,fitInteractions=F, verbose=T,fixedTerms=c(list(ag_suit=3), list(log_slope=3), list(log_elevation=3)))
FixedModel(all.data=data.spr.rr,responseVar="Species_richness",fitFamily="poisson", fixedStruct="Predominant_habitat+Within_PA+Predominant_habitat:Within_PA+Zone+taxon_of_interest+poly(log_elevation,2) + poly(ag_suit, 1)", randomStruct="(1+Within_PA|SS)+(1|SSB)+(1|SSBS)",REML=TRUE,optimizer="bobyqa")

#richness rarefied
CompareRandoms(dataset=data.spr.rr,responseVar="Richness_rarefied",fitFamily="poisson",fixedFactors=c("Predominant_habitat","Within_PA", "Zone","taxon_of_interest"),fixedTerms=c(list(ag_suit=3), list(log_slope=3), list(log_elevation=3)),fixedInteractions="Predominant_habitat:Within_PA",siteRandom=T,fitInteractions=F,verbose=T,optimizer="bobyqa",randomSlopes=T)
ModelSelect(all.data=data.spr.rr,responseVar="Richness_rarefied",fitFamily="poisson", fixedFactors=c("Predominant_habitat","Within_PA", "Zone","taxon_of_interest"),fixedInteractions=c("Predominant_habitat:Within_PA"),randomStruct="(1+Within_PA|SS)+(1|SSB)+(1|SSBS)",siteRandom=T,fitInteractions=F, verbose=T,fixedTerms=c(list(ag_suit=3), list(log_slope=3), list(log_elevation=3)))
FixedModel(all.data=data.spr.rr,responseVar="Richness_rarefied",fitFamily="poisson", fixedStruct="Predominant_habitat:Within_PA+ Predominant_habitat+Within_PA+Zone+poly(log_elevation,3) + poly(ag_suit, 3)", randomStruct="(1+Within_PA|SS)+(1|SSB)+(1|SSBS)",REML=TRUE,optimizer="bobyqa")

