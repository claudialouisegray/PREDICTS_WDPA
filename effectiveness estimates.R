


setwd("N:/Documents/PREDICTS/WDPA analysis/effectiveness estimates")

#these text files are exported from the attribute tables of the rasters created in python scripts for effectiveness estimates

sp_in_aes <- read.csv("sp_in_aes.txt", header = T)
sp.mean.val.per.cell.in <-weighted.mean(sp_in_aes$VALUE_,sp_in_aes$COUNT_)
#this is the calculated species richness inside PAs based on PREDICTS data

sp_out_aes <- read.csv("sp_out_aes.txt", header = T)
sp.mean.val.per.cell.out <-weighted.mean(sp_out_aes$VALUE_,sp_out_aes$COUNT_)
#this is the calculated species richness outside PAs based on PREDICTS data

ab_in_aes <- read.csv("ab_in_aes.txt", header = T)
ab.mean.val.per.cell.in <-weighted.mean(ab_in_aes$VALUE_,ab_in_aes$COUNT_)

ab_out_aes <- read.csv("ab_out_aes.txt", header = T)
ab.mean.val.per.cell.out <-weighted.mean(ab_out_aes$VALUE_,ab_out_aes$COUNT_)

plot(exp(c(1,1.1,1.2,1.3,1.4)) ~ c(1,1.1,1.2,1.3,1.4))


#PA effectiveness estimates

# benefit for species richness


b.sp <- sp.mean.val.per.cell.in/sp.mean.val.per.cell.out -1
b.a <-  ab.mean.val.per.cell.in/ab.mean.val.per.cell.out -1



b.vals <- c(b.sp, b.a)



vals <- data.frame(name =  c("b.sp", "b.a"),
			b.vals = b.vals, 
			metric = rep(c("sp.rich", "abundance"), each = 1), 
			NPA.abs = rep(NA,length(b.vals)),
			PA.abs = rep(NA,length(b.vals)),
			est = rep(NA,length(b.vals)))


#By 2005, we estimate that human impacts had reduced local richness by an average of 13.6% (95% CI: 9.1 – 17.8%) and 
#total abundance by 10.7% (95% CI: 3.8% gain – 23.7% reduction) compared with pre-impact times. 
#Approximately 60% of the decline in richness was independent of effects on abundance: 
#average rarefied richness has fallen by 8.1% (95% CI: 3.5 – 12.9%).
			
for(b in vals$b.vals){

benefit <- as.numeric(b)		# percentage increase in metric in PAs
PA.pct <- 13 				# percentage of total land area in PAs
# global loss of biodiversity (from Newbold et al)
if(vals$metric[which(vals$b.vals == b)] == "sp.rich"){
	global.loss <- 0.136			
	}else if (vals$metric[which(vals$b.vals == b)] == "abundance"){
	global.loss <- 0.107
	}else if (vals$metric[which(vals$b.vals == b)] == "rar.rich"){
	global.loss <- 0.081
	}

NPA.rel <- 1-benefit			# relative biodiversity in unprotected sites
PA.rel <- 1					# biodiversity in protected sites
NPA.pct <- 100-PA.pct			# land area unprotected 
global.int <- 1-global.loss		# overall status of biodiversity relative to pristine


# we want NPA.abs and PA.abs - where these are the biodiv metrics in unprotected and protected relative to pristine
# simultaneous equations are

# global.int = NPA.pct/100*NPA.abs + PA.pct/100*PA.abs #overall loss is loss in PAs and nonPAs relative to pristine
# NPA.abs = PA.abs*NPA.rel			     


# so
#global.int = (NPA.pct/100)*PA.abs*NPA.rel + (PA.pct/100)*PA.abs
# which is the same as
#global.int = PA.abs*(NPA.pct/100*(NPA.rel) + PA.pct/100)

# then
PA.abs <- global.int/(PA.pct/100 + (NPA.pct/100 * (NPA.rel)))
NPA.abs <- PA.abs*NPA.rel
vals$NPA.abs[which(vals$b.vals == b)] <- NPA.abs
vals$PA.abs[which(vals$b.vals == b)] <- PA.abs

# if pristine is 1, then where between 1 and NPA.abs does PA.abs fall?
# get difference between PA.abs and NPA.abs as a percentage of NPA.abs
#((1-NPA.abs)-(1-PA.abs))/(1-NPA.abs)

est <- 1-(1-PA.abs)/(1-NPA.abs)
vals$est[which(vals$b.vals == b)] <- est

# ie
#est <- ((1-NPA.abs)-(1-PA.abs))/(1-NPA.abs)
}

vals

write.csv(vals, "simple.effectiveness.estimates.aes.csv")

