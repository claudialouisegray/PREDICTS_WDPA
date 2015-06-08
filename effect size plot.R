
setwd("R:\\ecocon_d\\shared\\PREDICTS PA paper")
x <- read.csv("effect sizes 20150507.csv")


col.key <- data.frame(Predominant_habitat = c("Primary Vegetation",
                                               "Mature secondary vegetation",
                                                "Intermediate secondary vegetation",
                                                "Young secondary vegetation",
								"Plantation forest",
								"Cropland",
								"Pasture",
								"Urban"),
                                           col=c("#66A61E","#147659","#1B9E77","#8ecfbc","#7570B3",
                                                 "#E6AB02","#D95F02","#E7298A"),
                                           col.ci=c("#66A61E90","#14765990","#1B9E7790","#8ecfbc90","#7570B390",
                                                 "#E6AB0290","#D95F0290","#E7298A90"))
col.key[,2] <- as.character(col.key[,2])
col.key[,3] <- as.character(col.key[,3])


### PLOT


tiff("R:/ecocon_d/clg32/PREDICTS/WDPA analysis/plots/02_15/effect size plot 20150527.tif",
	width = 10, height = 20, units = "cm", pointsize = 12, res = 300)


par(mfrow = c(3,1), mar = c(2,6,2,2), oma = c(12,1,1,1))
sz <- 1.5
pos <- c(1:4)
pos2 <- c(5:9)

#species richness
plot(-1,-1, xlim = c(0, nrow(x) + 1), ylim = c(-45,65),
		pch = 0, bty = "l",
		xaxt = "n",
		ylab =  "% difference between sites outside and \ninside a protected area",
		xlab = "",
		main = "Species Richness")

points( x$Species.Richness[pos] ~ pos, pch = 16, cex = sz,
	col = col.key[c(5,6,7,8), 2])
arrows( pos, x$sp_max[pos], pos, x$sp_min[pos], 
	code = 3, length = 0, lwd = 5, 
	col = col.key[c(5,6,7,8), "col.ci"])

points( x$Species.Richness[pos2] ~ pos2, pch = 16, cex = sz,
	col = col.key[c(5,6,7,8,1), 2])
arrows(pos2, x$sp_max[pos2], pos2, x$sp_min[pos2 ],
	code = 3, length = 0, lwd = 5, 
	col = col.key[c(5,6,7,8,1), "col.ci"])

abline( h = 0, lty = 2, col = 8)



#abundance
plot(-1,-1, xlim = c(0, nrow(x) + 1), ylim = c(-45,65),
		pch = 0, bty = "l",
		xaxt = "n",
		ylab =  "% difference between sites outside and \ninside a protected area",
		xlab = "",
		main = "Abundance")

points( x$Abundance[pos] ~ pos, pch = 16, cex = sz,
	col = col.key[c(5,6,7,8), 2])
arrows( pos, x$ab_max[pos], pos, x$ab_min[pos], 
	code = 3, length = 0, lwd = 5, 
	col = col.key[c(5,6,7,8), "col.ci"])

points( x$Abundance[pos2] ~ pos2, pch = 16, cex = sz,
	col = col.key[c(5,6,7,8,1), 2])
arrows(pos2, x$ab_max[pos2], pos2, x$ab_min[pos2 ],
	code = 3, length = 0, lwd = 5, 
	col = col.key[c(5,6,7,8,1), "col.ci"])

abline( h = 0, lty = 2, col = 8)



# endemicity
plot(-1,-1, xlim = c(0, nrow(x) + 1), ylim = c(-45,65),
		pch = 0, bty = "l",
		xaxt = "n",
		ylab =  "% difference between sites outside and \ninside a protected area",
		xlab = "",
		main = "Endemicity")


points( x$Endemicity[pos] ~ pos, pch = 16, cex = sz,
	col = col.key[c(5,6,7,8), 2])
arrows( pos, x$en_max[pos], pos, x$en_min[pos], 
	code = 3, length = 0, lwd = 5, 
	col = col.key[c(5,6,7,8), "col.ci"])

points( x$Endemicity[pos2] ~ pos2, pch = 16, cex = sz,
	col = col.key[c(5,6,7,8,1), 2])
arrows(pos2, x$en_max[pos2], pos2, x$en_min[pos2 ],
	code = 3, length = 0, lwd = 5, 
	col = col.key[c(5,6,7,8,1), "col.ci"])

abline( h = 0, lty = 2, col = 8)


axis(1,1:9, x$description[c(1:9)], las = 2, cex.axis = 1)

dev.off()





