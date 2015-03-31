

plotLUPA <- function (model,responseVar, data,
				xvar = "Predominant_habitat",
				xvar.order = levels(data[,xvar]),
				col.key = NULL,
				labels = NULL,
				logLink = "n",
				seMultiplier=1.96,
				forPaper = FALSE,
				cex.txt = 0.5,
				offSet = 0,
				ylims = c(0,200)){

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


if(is.null(labels)){
	labels <- rep(col.key$Predominant_habitat, each = 2)}

ests <- NULL
se <- NULL
#get estimates for all levels
for(i in 1:length(xvar.order)){
	ests[i] <- fixef(model)[which(names(fixef(model)) == paste("LUPA", xvar.order[i], sep = ""))]
	se[i] <- se.fixef(model)[which(names(fixef(model)) == paste("LUPA", xvar.order[i], sep = ""))]
	}

  y <- ests
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
    cex.pt<- 0.8
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

#make y axis label dependent on whether the response was log transformed or not
        if (logLink == "n" & responseVar != "Endemicity") {
          ylabel = paste("Relative", responseVar)
        }else {
	    ylabel = paste(responseVar, "(%)")
	 }


plot(-1, -1, xlim = c(0,length(xvar.order)+1), bty = "n",
		ylim = ylims, xlab = NA, ylab = ylabel, xaxt = "n")



# set backdrop

	if (!forPaper) {
                rect(min(grep("Primary", labels)) - 0.5, ylims[1],
                  max(grep("Primary", labels)) + 0.5, ylims[2], 
                  col = "#66A61E33", border = NA)
            }
            text(mean(grep("Primary", labels)), ylims[1] +  # get the position as the mean between all points where "Primary" is found in labels
                0.05 * predRange, "Primary", col = "black", cex = cex.txt)
            text(min(grep("Primary", labels)), ylims[1] + 
                0.1 * predRange, paste(nrow(data[which(data$Predominant_habitat == "Primary Vegetation" & data$Within_PA == "no"),])), 
			col = "black", cex = cex.txt)		
            text(max(grep("Primary", labels)), ylims[1] + 
                0.1 * predRange, paste(nrow(data[which(data$Predominant_habitat == "Primary Vegetation" & data$Within_PA == "yes"),])), 
			col = "black", cex = cex.txt)	


	if (!forPaper) {
                  rect(min(grep("Mature secondary", labels)) - 0.5, ylims[1], 
				max(grep("Mature secondary", labels)) + 0.5, ylims[2], 	
				col = "#14765933", border = NA)
                }
                text(mean(grep("Mature secondary", labels)), 
                  ylims[1] + 0.05 * predRange, "MSV", col = "black", 
                  cex = cex.txt)
            text(min(grep("Mature secondary", labels)), ylims[1] + 
                0.1 * predRange, paste(nrow(data[which(data$Predominant_habitat == "Mature secondary vegetation" & data$Within_PA == "no"),])), 
			col = "black", cex = cex.txt)		
            text(max(grep("Mature secondary", labels)), ylims[1] + 
                0.1 * predRange, paste(nrow(data[which(data$Predominant_habitat == "Mature secondary vegetation" & data$Within_PA == "yes"),])), 
			col = "black", cex = cex.txt)	

	if (!forPaper) {
                  rect(min(grep("Intermediate secondary", labels)) - 0.5, ylims[1], 
				max(grep("Intermediate secondary", labels)) + 0.5, ylims[2], 
				col = "#1B9E7733", border = NA)
                }
                text(mean(grep("Intermediate secondary", labels)), 
                  ylims[1] + 0.05 * predRange, "ISV", col = "black", 
                  cex = cex.txt)
            text(min(grep("Intermediate secondary", labels)), ylims[1] + 
                0.1 * predRange, paste(nrow(data[which(data$Predominant_habitat == "Intermediate secondary vegetation" & data$Within_PA == "no"),])), 
			col = "black", cex = cex.txt)		
            text(max(grep("Intermediate secondary", labels)), ylims[1] + 
                0.1 * predRange, paste(nrow(data[which(data$Predominant_habitat == "Intermediate secondary vegetation" & data$Within_PA == "yes"),])), 
			col = "black", cex = cex.txt)	

	if (!forPaper) {
                  rect(min(grep("Young secondary", labels)) - 0.5, ylims[1], 
				max(grep("Young secondary", labels)) + 0.5, ylims[2], 
				col = "#8ecfbc33", border = NA)
                }
                text(mean(grep("Young secondary", labels)), ylims[1] + 
                  0.05 * predRange, "YSV", col = "black", cex = cex.txt)
            text(min(grep("Young secondary", labels)), ylims[1] + 
                0.1 * predRange, paste(nrow(data[which(data$Predominant_habitat == "Young secondary vegetation" & data$Within_PA == "no"),])), 
			col = "black", cex = cex.txt)		
            text(max(grep("Young secondary", labels)), ylims[1] + 
                0.1 * predRange, paste(nrow(data[which(data$Predominant_habitat == "Young secondary vegetation" & data$Within_PA == "yes"),])), 
			col = "black", cex = cex.txt)	

	if (!forPaper) {
                rect(min(grep("Plantation", labels)) - 0.5, ylims[1], 
                  max(grep("Plantation", labels)) + 0.5, ylims[2], 
                  col = "#7570B333", border = NA)
            }
            text(mean(grep("Plantation", labels)), ylims[1] + 
                0.05 * predRange, "Plantation", col = "black", cex = cex.txt)
            text(min(grep("Plantation", labels)), ylims[1] + 
                0.1 * predRange, paste(nrow(data[which(data$Predominant_habitat == "Plantation forest" & data$Within_PA == "no"),])), 
			col = "black", cex = cex.txt)		
            text(max(grep("Plantation", labels)), ylims[1] + 
                0.1 * predRange, paste(nrow(data[which(data$Predominant_habitat == "Plantation forest" & data$Within_PA == "yes"),])), 
			col = "black", cex = cex.txt)	

	if (!forPaper) {
                rect(min(grep("Cropland", labels)) - 0.5, ylims[1], 
                  max(grep("Cropland", labels)) + 0.5, ylims[2], 
                  col = "#F7EBD0", border = NA)
            }
            text(mean(grep("Cropland", labels)), ylims[1] + 
                0.05 * predRange, "Cropland", col = "black", cex = cex.txt)
            text(min(grep("Cropland", labels)), ylims[1] + 
                0.1 * predRange, paste(nrow(data[which(data$Predominant_habitat == "Cropland" & data$Within_PA == "no"),])), 
			col = "black", cex = cex.txt)		
            text(max(grep("Cropland", labels)), ylims[1] + 
                0.1 * predRange, paste(nrow(data[which(data$Predominant_habitat == "Cropland" & data$Within_PA == "yes"),])), 
			col = "black", cex = cex.txt)	

	if (!forPaper) {
                rect(min(grep("Pasture", labels)) - 0.5, ylims[1], 
                  max(grep("Pasture", labels)) + 0.5, ylims[2], 
                  col = "#D95F0233", border = NA)
            }
            text(mean(grep("Pasture", labels)), ylims[1] + 
                0.05 * predRange, "Pasture", col = "black", cex = cex.txt)
            text(min(grep("Pasture", labels)), ylims[1] + 
                0.1 * predRange, paste(nrow(data[which(data$Predominant_habitat == "Pasture" & data$Within_PA == "no"),])), 
			col = "black", cex = cex.txt)		
            text(max(grep("Pasture", labels)), ylims[1] + 
                0.1 * predRange, paste(nrow(data[which(data$Predominant_habitat == "Pasture" & data$Within_PA == "yes"),])), 
			col = "black", cex = cex.txt)	

	if (!forPaper) {
                rect(min(grep("Urban", labels)) - 0.5, ylims[1], 
                  max(grep("Urban", labels)) + 0.5, ylims[2], 
                  col = "#E7298A33", border = NA)
            }
            text(mean(grep("Urban", labels)), ylims[1] + 0.05 * 
                predRange, "Urban", col = "black", cex = cex.txt)
            text(min(grep("Urban", labels)), ylims[1] + 
                0.1 * predRange, paste(nrow(data[which(data$Predominant_habitat == "Urban" & data$Within_PA == "no"),])), 
			col = "black", cex = cex.txt)		
            text(max(grep("Urban", labels)), ylims[1] + 
                0.1 * predRange, paste(nrow(data[which(data$Predominant_habitat == "Urban" & data$Within_PA == "yes"),])), 
			col = "black", cex = cex.txt)	


#set filled and empty points if interaction, else all points the same

#make y axis label dependent on whether the response was log transformed or not

if (logLink == "n" & responseVar != "Endemicity") {
yplot <- c(0,y)
yminus.plot <- c(0, yminus)
yplus.plot <- c(0, yplus)
 
 }else {
yplot <- c(100,y)
yminus.plot <- c(100, yminus)
yplus.plot <- c(100, yplus)
}

shapes <- rep(c(21,16), length(yplot/2))	 

abline(h = 0, lty = 2, col = 1)

if (!forPaper) {
	points(0, 0, pch = 22, cex = 1.5 * cex.pt, col = 1)
	errbar(2:length(yplot)+ offSet, yplot[2:length(yplot)], 
		yplus.plot[2:length(yplus.plot)], yminus.plot[2:length(yminus.plot)], 
		col = 1,  pch = shapes[2:length(shapes)],
            errbar.col = 1, add = T, cex = cex.pt)
	}else{
	cols <- col.key[which(col.key$Predominant_habitat %in% labels),2]
	cols <- rep(cols, each = 2)
	points(1, 0, pch = 21, cex = 1.5 * cex.pt, col = cols[1], bg = "white")
	errbar(2:length(yplot)+ offSet, yplot[2:length(yplot)], 
		yplus.plot[2:length(yplus.plot)], yminus.plot[2:length(yminus.plot)], 
		col = cols[2:length(cols)],  pch = shapes[2:length(shapes)],
            errbar.col = cols[2:length(cols)], add = T, cex = cex.pt, cap = 0)
	#put points back on top
	points(1:length(yplot), yplot, pch = shapes, cex = 1.5 * cex.pt, col = cols, bg = "white")
	}




#legend(x = max(xlims), y = max(ylims), c("Protected", "Not protected"), 
#	xjust = 1,pch = c(21,16))

}




