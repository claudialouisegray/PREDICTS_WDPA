addDataDistribution <- function(b = 50, x, var, data, levels, 
	axis.text.pos, axis.text, spacing = 30, legend.spacing = 0, include.lowest = TRUE,
	levels.col, xlim = c(0,1000)){

# b is the number of bins to cut the continuous data into
# var is the variable for which the number of sites to be shown at the base of the plot, given between ""
# x is the variable that is plotted along the x axis, given between ""
# levels must be a list of the unique values of data$var, in the order to be seen from bottom to top
# levels.col is the colors of the categories, in the order corresponding to categories
# axis.text.pos is the positions at which to put x axis labels
# axis.text is the text for the x axis labels
# legend spacing gives the gap between the end of the x axis and the legend, use units of x. 

y <- 0

break.points <- seq(min(data[,paste(x)]), max(data[,paste(x)]), length = b)

# get maximum y value

y.lim.upper <- 10

for(l in levels){
	data.subset <- data[which(data[,paste(var)] == l),]
	breaks <- cut(data.subset[,paste(x)], breaks = break.points, include.lowest = include.lowest)
	d <- data.frame(table(breaks))
	y.lim.upper <- y.lim.upper + max(d$Freq)+ spacing
	}

plot(rep(0, length(data[,paste(x)])) ~ data[,paste(x)], ylim = c(0,y.lim.upper), xlim = xlim, 
	bty = "l", type = "n",
	ylab = "Number of sites", xlab = "" , axes = F)
axis(1, axis.text.pos,axis.text)

for (i in 1:length(levels)){
	data.subset <- data[which(data[,paste(var)] == levels[i]),]
	breaks <- cut(data.subset[,paste(x)], breaks = break.points, include.lowest = include.lowest)
	d <- data.frame(table(breaks))
	d$breaks <- gsub("\\(", "", d$breaks)
	d$breaks <- gsub("\\]", "", d$breaks)
	d$breaks <- gsub("\\[", "", d$breaks)
	d$lower <- lapply(strsplit(as.character(d$breaks), ","), "[",1)
	d$upper <- lapply(strsplit(as.character(d$breaks), ","), "[",2)
	d$lower <- as.numeric(d$lower)
	d$upper <- as.numeric(d$upper)
	d$mean <- rowMeans(d[,c("lower", "upper")])

	polygon(c(d$mean,rev(d$mean)), c(d$Freq+y,rep(y,length(d$mean))), lty = 1, border = levels.col[i], col = levels.col[i])
	y <- y + max(d$Freq) + spacing 
	}

arrows(0,0,0,100, length = 0.015, code = 3, angle = 90)
text(x = 0, y = 50, pos = 2, cex = 0.7, offset = 0.6,
	"100 \n sites")
legend(x = max(axis.text.pos)+ legend.spacing, y = y.lim.upper/2, rev(levels), xjust = 0, yjust = 0.5,
	cex = 1, col = rev(levels.col), pch = 16, bty = "n", ncol = 1, text.width = 0.7)

}

