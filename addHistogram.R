addHistogram  <- function(data, var, x, levels,levels.col,bar.breaks, xlim = c(0,8), axis.labels = seq(1,8,1),
		spacing = 30, legend.spacing = 0.3){

#data = model.data
#var = "Predominant_habitat"
#x =   "ag_suit"
#levels = lu
#levels.col = lu.cols2.ci
# bar.breaks = ag2 

# var is the variable for which the number of sites to be shown at the base of the plot, given between ""
# x is the variable that is plotted along the x axis, given between ""
# levels must be a list of the unique values of data$var, in the order to be seen from bottom to top
# levels.col is the colors of the categories, in the order corresponding to categories
# bar.breaks gives the values of the beginning and end of each bar

struct <- paste(var, x, sep = "~")
for (i in 1:length(levels)){
	t2 <- aggregate(formula(struct), data = model.data[which(model.data[,paste(var)] == levels[i]),], FUN= length)
	names(t2) <- c(paste(x), levels[i])
	#combine into table so missing levels can be replaced with 0 where appropriate
		if(i == 1){
		t1 <- t2
		}else{
		t1 <- merge(t1,t2, by = paste(x), all.x = T, all.y = T)
	}
}

t1[is.na(t1)] <- 0
# get maximum y value
y.lim.upper <- 50
for(i in 2:ncol(t1)){
	y.lim.upper <- y.lim.upper + max(t1[,i])+ spacing
	}

plot(rep(0, length(bar.breaks)) ~ bar.breaks, ylim = c(0,y.lim.upper), xlim = xlim, 
	bty = "l", type = "n",
	ylab = "Number of sites", xlab = "", axes = F)
axis(1, axis.labels,axis.labels)
y <- 0
for(i in 2:ncol(t1)){
	yvals <- t1[,i] + y
	polygon(c(bar.breaks,rev(bar.breaks)), c(rep(yvals, each = 2),rep(min(ylims)+y,length(bar.breaks))), lty = 0, col = levels.col[i-1],
		border = levels.col[i-1])
	y.pos <- y + max(t1[,i])/2
#	text(x = max(bar.breaks), y = y.pos, pos = 4, cex = 0.7, paste(levels[i - 1]))
	y <- y + min(ylims) + max(t1[,i]) + spacing 
	}

arrows(0.45,0,0.45,100, length = 0.015, code = 3, angle = 90)
text(x = 0.5, y = 50, pos = 2, cex = 0.7, offset = 0.6,
	"100 \n sites")
legend(x = max(bar.breaks)+ legend.spacing, y = max(yvals)/2, rev(levels), xjust = 0, yjust = 0.5,
	cex = 1, col = rev(levels.col), pch = 16, bty = "n", ncol = 1, text.width = 0.7)

}


