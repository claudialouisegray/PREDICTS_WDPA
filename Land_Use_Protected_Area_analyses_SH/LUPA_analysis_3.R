library(roquefort)

#LUPA analysis 3 levels

#taxon####

#species richness
CompareRandoms(dataset=data.spr.rr.com,responseVar="Species_richness",fitFamily="poisson",fixedFactors=c("LUPA", "Zone","taxon_of_interest"),fixedTerms=c(list(ag_suit=3), list(log_slope=3), list(log_elevation=3)),fixedInteractions="LUPA:taxon_of_interest",siteRandom=T,fitInteractions=F,verbose=T,optimizer="bobyqa",randomSlopes=T)
ModelSelect_PA(all.data=data.spr.rr.com,responseVar="Species_richness",fitFamily="poisson", fixedFactors=c("Zone","taxon_of_interest", "LUPA"), fixedInteractions=c("LUPA:taxon_of_interest"), randomStruct="(1+Within_PA|SS)+(1|SSB)+(1|SSBS)",siteRandom=TRUE,fitInteractions=F, verbose=T,fixedTerms=c(list(ag_suit=3,log_elevation=3,log_slope=3)))
FixedModel(all.data=data.spr.rr.com,responseVar="Species_richness",fitFamily="poisson", fixedStruct="LUPA+taxon_of_interest+Zone+LUPA:taxon_of_interest+poly(ag_suit,3)+poly(log_elevation,3)", randomStruct="(1+Within_PA|SS)+(1|SSB)+(1|SSBS)",REML=TRUE,optimizer="bobyqa")

#abundance
CompareRandoms(dataset=data.com,responseVar="LogAbund",fitFamily="gaussian",fixedFactors=c("LUPA", "Zone","taxon_of_interest"),fixedTerms=c(list(ag_suit=3), list(log_slope=3), list(log_elevation=3)),fixedInteractions="LUPA:taxon_of_interest",siteRandom=F,fitInteractions=F,verbose=T,optimizer="bobyqa",randomSlopes=T)
ModelSelect_PA(all.data=data.com,responseVar="LogAbund",fitFamily="gaussian", fixedFactors=c("LUPA", "Zone","taxon_of_interest"), fixedInteractions=c("LUPA:taxon_of_interest"),randomStruct="(1+Within_PA|SS)+(1|SSB)",siteRandom=FALSE,fitInteractions=F, verbose=T,fixedTerms=c(list(ag_suit=3), list(log_slope=3), list(log_elevation=3)))
FixedModel(all.data=data.com,responseVar="LogAbund",fitFamily="gaussian", fixedStruct="LUPA + taxon_of_interest + LUPA:taxon_of_interest", randomStruct="(1+Within_PA|SS)+(1|SSB)",REML=T,optimizer="bobyqa")

#rarefied richness
CompareRandoms(dataset=data.spr.rr.com,responseVar="Richness_rarefied",fitFamily="poisson",fixedFactors=c("LUPA", "Zone","taxon_of_interest"),fixedTerms=c(list(ag_suit=3), list(log_slope=3), list(log_elevation=3)),fixedInteractions="LUPA:taxon_of_interest",siteRandom=T,fitInteractions=F,verbose=T,optimizer="bobyqa",randomSlopes=T)
ModelSelect_PA(all.data=data.spr.rr.com,responseVar="Richness_rarefied",fitFamily="poisson", fixedFactors=c("LUPA", "Zone", "taxon_of_interest"), randomStruct="(1+Within_PA|SS)+(1|SSB) +(1|SSBS)",siteRandom=T,fitInteractions=F, fixedInteractions="LUPA:taxon_of_interest", verbose=T,fixedTerms=c(list(ag_suit=3,log_elevation=3,log_slope=3)))
FixedModel(all.data=data.spr.rr.com,responseVar="Richness_rarefied",fitFamily="poisson", fixedStruct="LUPA:taxon_of_interest+LUPA+Zone+taxon_of_interest+poly(log_elevation,3) + poly(ag_suit, 3)", randomStruct="(1+Within_PA|SS)+(1|SSB)+(1|SSBS)",REML=TRUE,optimizer="bobyqa")

#range size
CompareRandoms(dataset=data.com,responseVar="CWM_Geographic_range_log10_square_km",fitFamily="gaussian",fixedFactors=c("LUPA", "Zone","taxon_of_interest"),fixedTerms=c(list(ag_suit=3), list(log_slope=3), list(log_elevation=3)),fixedInteractions="LUPA:taxon_of_interest",siteRandom=F,fitInteractions=F,verbose=T,optimizer="bobyqa",randomSlopes=T)
ModelSelect_PA(all.data=data.com,responseVar="CWM_Geographic_range_log10_square_km",fitFamily="gaussian",fixedFactors=c("LUPA", "Zone","taxon_of_interest"), fixedInteractions="LUPA:taxon_of_interest", randomStruct="(1+Within_PA|SS)+(1|SSB)",siteRandom=FALSE,fitInteractions=F, verbose=T,fixedTerms=c(list(ag_suit=3), list(log_slope=3), list(log_elevation=3)))
FixedModel(all.data=data.com,responseVar="CWM_Geographic_range_log10_square_km",fitFamily="gaussian", fixedStruct="LUPA + taxon_of_interest + Zone + LUPA:taxon_of_interest + poly(log_slope,3) + poly(log_elevation,3) + poly(ag_suit,3)", randomStruct="(1+Within_PA|SS)+(1|SSB)",REML=T,optimizer="bobyqa")

#zone####
#species richness
CompareRandoms(dataset=data.spr.rr.com,responseVar="Species_richness",fitFamily="poisson",fixedFactors=c("LUPA", "Zone","taxon_of_interest"),fixedTerms=c(list(ag_suit=3), list(log_slope=3), list(log_elevation=3)),fixedInteractions="LUPA:Zone",siteRandom=T,fitInteractions=F,verbose=T,optimizer="bobyqa",randomSlopes=T)
ModelSelect_PA(all.data=data.spr.rr.com,responseVar="Species_richness",fitFamily="poisson", fixedFactors=c("Zone","taxon_of_interest", "LUPA"), fixedInteractions=c("LUPA:Zone"), randomStruct="(1+Within_PA|SS)+(1|SSB)+(1|SSBS)",siteRandom=TRUE,fitInteractions=F, verbose=T,fixedTerms=c(list(ag_suit=3,log_elevation=3,log_slope=3)))
FixedModel(all.data=data.spr.rr.com,responseVar="Species_richness",fitFamily="poisson", fixedStruct="LUPA+taxon_of_interest+Zone+LUPA:Zone+poly(ag_suit,3)+poly(log_elevation,3)", randomStruct="(1+Within_PA|SS)+(1|SSB)+(1|SSBS)",REML=TRUE,optimizer="bobyqa")

#abundance
CompareRandoms(dataset=data.com,responseVar="LogAbund",fitFamily="gaussian",fixedFactors=c("LUPA", "Zone","taxon_of_interest"),fixedTerms=c(list(ag_suit=3), list(log_slope=3), list(log_elevation=3)),fixedInteractions="LUPA:Zone",siteRandom=F,fitInteractions=F,verbose=T,optimizer="bobyqa",randomSlopes=T)
ModelSelect_PA(all.data=data.com,responseVar="LogAbund",fitFamily="gaussian", fixedFactors=c("LUPA", "Zone","taxon_of_interest"), fixedInteractions=c("LUPA:Zone"),randomStruct="(1+Within_PA|SS)+(1|SSB)",siteRandom=FALSE,fitInteractions=F, verbose=T,fixedTerms=c(list(ag_suit=3), list(log_slope=3), list(log_elevation=3)))
FixedModel(all.data=data.com,responseVar="LogAbund",fitFamily="gaussian", fixedStruct="LUPA + taxon_of_interest + Zone + LUPA:Zone", randomStruct="(1+Within_PA|SS)+(1|SSB)",REML=T,optimizer="bobyqa")

#rarefied richness
CompareRandoms(dataset=data.spr.rr.com,responseVar="Richness_rarefied",fitFamily="poisson",fixedFactors=c("LUPA", "Zone","taxon_of_interest"),fixedTerms=c(list(ag_suit=3), list(log_slope=3), list(log_elevation=3)),fixedInteractions="LUPA:Zone",siteRandom=T,fitInteractions=F,verbose=T,optimizer="bobyqa",randomSlopes=T)
ModelSelect_PA(all.data=data.spr.rr.com,responseVar="Richness_rarefied",fitFamily="poisson", fixedFactors=c("LUPA", "Zone", "taxon_of_interest"), randomStruct="(1+Within_PA|SS)+(1|SSB) +(1|SSBS)",siteRandom=T,fitInteractions=F, fixedInteractions="LUPA:Zone", verbose=T,fixedTerms=c(list(ag_suit=3,log_elevation=3,log_slope=3)))
FixedModel(all.data=data.spr.rr.com,responseVar="Richness_rarefied",fitFamily="poisson", fixedStruct="LUPA:Zone+LUPA+Zone+poly(log_elevation,3) + poly(ag_suit, 3)", randomStruct="(1+Within_PA|SS)+(1|SSB)+(1|SSBS)",REML=TRUE,optimizer="bobyqa")

#range size
CompareRandoms(dataset=data.com,responseVar="CWM_Geographic_range_log10_square_km",fitFamily="gaussian",fixedFactors=c("LUPA", "Zone","taxon_of_interest"),fixedTerms=c(list(ag_suit=3), list(log_slope=3), list(log_elevation=3)),fixedInteractions="LUPA:Zone",siteRandom=F,fitInteractions=F,verbose=T,optimizer="bobyqa",randomSlopes=T)
ModelSelect_PA(all.data=data.com,responseVar="CWM_Geographic_range_log10_square_km",fitFamily="gaussian",fixedFactors=c("LUPA", "Zone","taxon_of_interest"), fixedInteractions="LUPA:Zone", randomStruct="(1+Within_PA|SS)+(1|SSB)",siteRandom=FALSE,fitInteractions=F, verbose=T,fixedTerms=c(list(ag_suit=3), list(log_slope=3), list(log_elevation=3)))
FixedModel(all.data=data.com,responseVar="CWM_Geographic_range_log10_square_km",fitFamily="gaussian", fixedStruct="LUPA + taxon_of_interest + Zone + LUPA:Zone + poly(log_slope,3) + poly(log_elevation,3) + poly(ag_suit,3)", randomStruct="(1+Within_PA|SS)+(1|SSB)",REML=T,optimizer="bobyqa")

#use intensity####

#species richness
CompareRandoms(dataset=data.spr.rr.com,responseVar="Species_richness",fitFamily="poisson",fixedFactors=c("LUUI", "Zone","taxon_of_interest","Within_PA"),fixedTerms=c(list(ag_suit=3), list(log_slope=3), list(log_elevation=3)),fixedInteractions="LUUI:Within_PA",siteRandom=T,fitInteractions=F,verbose=T,optimizer="bobyqa",randomSlopes=T)
ModelSelect_PA(all.data=data.spr.rr.com,responseVar="Species_richness",fitFamily="poisson", fixedFactors=c("Zone","taxon_of_interest", "LUUI","Within_PA"), fixedInteractions=c("LUUI:Within_PA"), randomStruct="(1+Within_PA|SS)+(1|SSB)+(1|SSBS)",siteRandom=TRUE,fitInteractions=F, verbose=T,fixedTerms=c(list(ag_suit=3,log_elevation=3,log_slope=3)))
FixedModel(all.data=data.spr.rr.com,responseVar="Species_richness",fitFamily="poisson", fixedStruct="LUUI:Within_PA+Zone+taxon_of_interest+LUUI+Within_PA+poly(log_elevation,2) + poly(ag_suit, 1)", randomStruct="(1+Within_PA|SS)+(1|SSB)+(1|SSBS)",REML=TRUE,optimizer="bobyqa")

#abundance
CompareRandoms(dataset=data.com,responseVar="LogAbund",fitFamily="gaussian",fixedFactors=c("LUUI", "Zone","taxon_of_interest","Within_PA"),fixedTerms=c(list(ag_suit=3), list(log_slope=3), list(log_elevation=3)),fixedInteractions="LUUI:Within_PA",siteRandom=F,fitInteractions=F,verbose=T,optimizer="bobyqa",randomSlopes=T)
ModelSelect_PA(all.data=data.com,responseVar="LogAbund",fitFamily="gaussian", fixedFactors=c("LUUI", "Zone","taxon_of_interest","Within_PA"), fixedInteractions=c("LUUI:Within_PA"),randomStruct="(1+Within_PA|SS)+(1|SSB)",siteRandom=FALSE,fitInteractions=F, verbose=T,fixedTerms=c(list(ag_suit=3), list(log_slope=3), list(log_elevation=3)))
FixedModel(all.data=data.com,responseVar="LogAbund",fitFamily="gaussian", fixedStruct="LUUI:PA +LUUI+Within_PA+ taxon_of_interest + Zone", randomStruct="(1+Within_PA|SS)+(1|SSB)",REML=T,optimizer="bobyqa")

#rarefied richness
CompareRandoms(dataset=data.spr.rr.com,responseVar="Richness_rarefied",fitFamily="poisson",fixedFactors=c("LUUI", "Zone","taxon_of_interest","Within_PA"),fixedTerms=c(list(ag_suit=3), list(log_slope=3), list(log_elevation=3)),fixedInteractions="LUUI:Within_PA",siteRandom=T,fitInteractions=F,verbose=T,optimizer="bobyqa",randomSlopes=T)
ModelSelect_PA(all.data=data.spr.rr.com,responseVar="Richness_rarefied",fitFamily="poisson", fixedFactors=c("LUUI", "Zone", "taxon_of_interest","Within_PA"), randomStruct="(1+Within_PA|SS)+(1|SSB) +(1|SSBS)",siteRandom=T,fitInteractions=F, fixedInteractions="LUUI:Within_PA", verbose=T,fixedTerms=c(list(ag_suit=3,log_elevation=3,log_slope=3)))
FixedModel(all.data=data.spr.rr.com,responseVar="Richness_rarefied",fitFamily="poisson", fixedStruct="LUUI:Within_PA+LUUI+Within_PA+Zone+poly(log_elevation,3) + poly(ag_suit, 3)", randomStruct="(1+Within_PA|SS)+(1|SSB)+(1|SSBS)",REML=TRUE,optimizer="bobyqa")

#range size
CompareRandoms(dataset=data.com,responseVar="CWM_Geographic_range_log10_square_km",fitFamily="gaussian",fixedFactors=c("LUUI", "Zone","taxon_of_interest","Within_PA"),fixedTerms=c(list(ag_suit=3), list(log_slope=3), list(log_elevation=3)),fixedInteractions="LUUI:Within_PA",siteRandom=F,fitInteractions=F,verbose=T,optimizer="bobyqa",randomSlopes=T)
ModelSelect_PA(all.data=data.com,responseVar="CWM_Geographic_range_log10_square_km",fitFamily="gaussian",fixedFactors=c("LUUI", "Zone","taxon_of_interest","Within_PA"), fixedInteractions="LUUI:Within_PA", randomStruct="(1+Within_PA|SS)+(1|SSB)",siteRandom=FALSE,fitInteractions=F, verbose=T,fixedTerms=c(list(ag_suit=3), list(log_slope=3), list(log_elevation=3)))
FixedModel(all.data=data.com,responseVar="CWM_Geographic_range_log10_square_km",fitFamily="gaussian", fixedStruct="LUUI:Within_PA +LUUI+Within_PA+taxon_of_interest + Zone + poly(log_slope,3) + poly(log_elevation,3) + poly(ag_suit,2)", randomStruct="(1+Within_PA|SS)+(1|SSB)",REML=T,optimizer="bobyqa")
