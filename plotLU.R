

plotLU <- function (model,responseVar, 
				xvar = "Predominant_habitat",
				intVar = "Within_PA",
				level2 = "yes",
				col.key = NULL,
				logLink = "n",
				seMultiplier=1.96,
				forPaper = FALSE,
				cex.txt = 0.5){

# reference level must be set to Primary Vegetation in data, and as model intercept too
# col.key is the dataframe of colours to use, the first 


if(is.null(col.key)){
col.key <- data.frame(Predominant_habitat = c("Primary Vegetation",
                                               "Mature secondary vegetation",
                                                "Intermediate secondary vegetation",
                                                "Young secondary vegetation",
								"Plantation forest",
								"Cropland",
								"Pasture",
								"Urban"),
                                           col=c("#66A61E","#147659","#1B9E77","#8ecfbc","#7570B3",
                                                 "#E6AB02","#D95F02","#E7298A"))
}

col.key[,2] <- as.character(col.key[,2])

coef.labels <- paste(xvar, col.key[,1], sep = "")

#get estimated for interaction with reference level
	ref.in <- fixef(model)[which(names(fixef(model)) == paste(intVar, level2, sep = ""))]
	ref.in.se <- se.fixef(model)[which(names(fixef(model)) == paste(intVar, level2, sep = ""))]

# get estimates and standard error for other landuses given interaction
	var.cov<-vcov(model)
	est <- c(0,ref.in)
	se <- c(0,ref.in.se)
	for(p in 2:length(coef.labels)){
		lu <- coef.labels[p]
		est1 <- fixef(model)[which(names(fixef(model))== lu)]
		se1 <- se.fixef(model)[which(names(fixef(model))== lu)]
		est <- c(est, est1)
		se <- c(se, se1)

		if(length(grep(":", names(fixef(model)))) > 0){
			#get positions of interaction terms
			pos <- grep(paste(intVar, level2, sep = ""), names(fixef(model)))
			#get position with the specified land use in 
			posInt <- grep(lu, names(fixef(model))[pos])
			est2 <- fixef(model)[pos[posInt]]
			est.p <- est1 + est2 + ref.in
			se2 <- se.fixef(model)[pos[posInt]]
    			cov <-var.cov[which(row.names(var.cov)==coef.labels[p]), which(names(var.cov[1,])==names(fixef(model))[pos[posInt]])]
			se.p <- sqrt(se1^2 + se2^2 + 2*cov)
			est <- c(est, est.p)
			se <- c(se, se.p)
			}
		}

  y <- est
  yplus<-y+se*seMultiplier
  yminus<-y-se*seMultiplier

  if(responseVar == "Endemicity"){
   # convert to real values and get reciprocal
  intercept <- fixef(model)[1]
  y1 <- 1/(y+intercept)
  yplus1<-1/(yplus+intercept)
  yminus1<-1/(yminus+intercept)

  # then get difference between this and 1/intercept
  #y <- y - 1/intercept
  #yplus <- yplus - 1/intercept
  #yminus <- yminus - 1/intercept

  # express this as a percentange of intercept, relative to 100
  y <- 100*(y1/(1/intercept)) -100
  yplus <- 100*yplus1/(1/intercept) -100
  yminus <- 100*yminus1/(1/intercept) -100
  }


# backtransform if needed
  
  if (logLink=="e"){
    y<-(exp(y)*100)-100
    yplus<-(exp(yplus)*100)-100
    yminus<-(exp(yminus)*100)-100
  } else if (logLink=="10") {
    y<-(10^(y)*100)-100
    yplus<-(10^(yplus)*100)-100
    yminus<-(10^(yminus)*100)-100
  } else if (logLink=="b"){
    intercept<-fixef(model)['(Intercept)']
    y<-(((1/(1+exp(-(intercept+y))))/(1/(1+exp(-(intercept)))))*100)-100
    yplus<-(((1/(1+exp(-(intercept+yplus))))/(1/(1+exp(-(intercept)))))*100)-100
    yminus<-(((1/(1+exp(-(intercept+yminus))))/(1/(1+exp(-(intercept)))))*100)-100
  } else if (logLink=="n"){
    
  } else {
    stop("Error: the specified log link is not supported")
  }



  if(forPaper){
    par(mar=c(0.2,3.5,0.2,0.2))
    par(cex.lab=1)
    par(cex.axis=0.7)
    cex.pt<-0.5
    par(mgp=c(2,1,0))
    cex.txt<-ifelse(is.null(cex.txt),0.75,cex.txt)
  } else {
    par(mar=c(1,5,1,1))
    par(cex.lab=1.8)
    par(cex.axis=1.5)
    cex.pt<-1
    par(mgp=c(3.5,1,0))
    cex.txt<-ifelse(is.null(cex.txt),0.8,cex.txt)
  }    
  par(las=1)


#set the max and min as the model CI range, or at least +100 or -100
	predRange <- max(yplus) - min(yminus) 
      plotLims <- c(min(c(yplus,yminus)- 0.4 * predRange) , 
			max(c(yminus,yplus) + 0.4* predRange))

#make y axis label dependent on whether the response was log transformed or not
        if (logLink == "n" & responseVar != "Endemicity") {
          ylabel = paste("Relative", responseVar)
        }else {
	    ylabel = paste(responseVar, "(%)")
	 }


if(length(grep(":", names(fixef(model)))) > 0){
 	xlims <- c(0.5, 2*length(coef.labels))
	labels <- rep(coef.labels, each = 2)
	}else{
	xlims = c(0.5, length(coef.labels))
	labels <- coef.labels
}

plot(-1, -1, xlim = xlims, 
		ylim = plotLims, xlab = NA, ylab = ylabel, xaxt = "n")


# set backdrop

	if (!forPaper) {
                rect(min(grep("Primary", labels)) - 0.5, plotLims[1],
                  max(grep("Primary", labels)) + 0.5, plotLims[2], 
                  col = "#66A61E33", border = NA)
            }
            text(mean(grep("Primary", labels)), plotLims[1] +  # get the position as the mean between all points where "Primary" is found in labels
                0.05 * predRange, "Primary", col = "black", cex = cex.txt)

	if (!forPaper) {
                  rect(min(grep("Mature secondary", labels)) - 0.5, plotLims[1], 
				max(grep("Mature secondary", labels)) + 0.5, plotLims[2], 	
				col = "#14765933", border = NA)
                }
                text(mean(grep("Mature secondary", labels)), 
                  plotLims[1] + 0.05 * predRange, "MSV", col = "black", 
                  cex = cex.txt)

	if (!forPaper) {
                  rect(min(grep("Intermediate secondary", labels)) - 0.5, plotLims[1], 
				max(grep("Intermediate secondary", labels)) + 0.5, plotLims[2], 
				col = "#1B9E7733", border = NA)
                }
                text(mean(grep("Intermediate secondary", labels)), 
                  plotLims[1] + 0.05 * predRange, "ISV", col = "black", 
                  cex = cex.txt)

	if (!forPaper) {
                  rect(min(grep("Young secondary", labels)) - 0.5, plotLims[1], 
				max(grep("Young secondary", labels)) + 0.5, plotLims[2], 
				col = "#8ecfbc33", border = NA)
                }
                text(mean(grep("Young secondary", labels)), plotLims[1] + 
                  0.05 * predRange, "YSV", col = "black", cex = cex.txt)


	if (!forPaper) {
                rect(min(grep("Plantation", labels)) - 0.5, plotLims[1], 
                  max(grep("Plantation", labels)) + 0.5, plotLims[2], 
                  col = "#7570B333", border = NA)
            }
            text(mean(grep("Plantation", labels)), plotLims[1] + 
                0.05 * predRange, "Plantation", col = "black", cex = cex.txt)

	if (!forPaper) {
                rect(min(grep("Cropland", labels)) - 0.5, plotLims[1], 
                  max(grep("Cropland", labels)) + 0.5, plotLims[2], 
                  col = "#F7EBD0", border = NA)
            }
            text(mean(grep("Cropland", labels)), plotLims[1] + 
                0.05 * predRange, "Cropland", col = "black", cex = cex.txt)

	if (!forPaper) {
                rect(min(grep("Pasture", labels)) - 0.5, plotLims[1], 
                  max(grep("Pasture", labels)) + 0.5, plotLims[2], 
                  col = "#D95F0233", border = NA)
            }
            text(mean(grep("Pasture", labels)), plotLims[1] + 
                0.05 * predRange, "Pasture", col = "black", cex = cex.txt)

	if (!forPaper) {
                rect(min(grep("Urban", labels)) - 0.5, plotLims[1], 
                  max(grep("Urban", labels)) + 0.5, plotLims[2], 
                  col = "#E7298A33", border = NA)
            }
            text(mean(grep("Urban", labels)), plotLims[1] + 0.05 * 
                predRange, "Urban", col = "black", cex = cex.txt)

shapes <- rep(c(16,1), length(coef.labels))


if (!forPaper) {
	points(1, 0, pch = 20, cex = 1.5 * cex.pt, col = 1)
	errbar(2:length(labels), y[2:length(y)], 
		yplus[2:length(yplus)], yminus[2:length(yminus)], 
		col = 1,  pch = shapes[2:length(shapes)],
            errbar.col = 1, add = T, cex = cex.pt)
	}else{
	cols <- rep(col.key[,2], each = 2)
	points(1, 0, pch = 20, cex = 1.5 * cex.pt, col = cols[1])
	errbar(2:length(labels), y[2:length(y)], 
		yplus[2:length(yplus)], yminus[2:length(yminus)], 
		col = cols[2:length(cols)],  pch = shapes[2:length(shapes)],
            errbar.col = cols[2:length(cols)], add = T, cex = cex.pt)
	}

abline(h = 0, lty = 2, col = 1)


legend(x = max(xlims), y = max(plotLims), c("Protected", "Not protected"), 
	xjust = 1,pch = c(1,16))

}




