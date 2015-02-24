

# only where the error for the second fixed effect is not of interest
# ie. only where for taxon or zone is the other factor in the interaction 



plotFactorInteraction <- function (model,responseVar, data,
				xvar = NULL,
				intvar = NULL,
				xvar.order = NULL,
				intvar.order = NULL, 
				logLink = "n",
				seMultiplier=1.96,
				forPaper = FALSE,
				cex.txt = 0.5,
				ref.text = "unprotected",
				xlab = "",
				ylab = ""){

# reference level must be set to age and size of 0
trop.col <- rgb(0.9,0,0)
temp.col <- rgb(0,0.1,0.7)
p.col <- rgb(0.2,0.7,0.2)
i.col <- rgb(0,0.7,0.9)
v.col <- rgb(0.9,0.5,0)


#get estimates for level(s) of xvar, in order required for plotting
pos.x <- NULL 
xests <- expand.grid(x = xvar.order, int = "reference", est = NA, se = NA)
for(i in 1:length(xvar.order)){
	pos.x <- which(names(fixef(model)) == paste(xvar, xvar.order[i], sep = ""))
	xests$est[i] <- fixef(model)[pos.x]
	xests$se[i] <- se.fixef(model)[pos.x]
	}

#get estimates for reference at level(s) of interaction variable, in order required for plotting
if(!is.null(intvar)){
	pos.int.ref <- NULL 
	int.ests <- expand.grid(x = "reference", int = intvar.order, est= 0, se = 0)
# not actually going to use code below as ref for each level of intvar should also be on the 100 line for PREDICTS WDPA plots 
# for(i in 1:length(intvar.order)){
#	pos.int.ref <- which(names(fixef(model)) == paste(intvar, intvar.order[i], sep = ""))
#	int.ests$est[i] <- fixef(model)[pos.int.ref]
#	int.ests$se[i] <- se.fixef(model)[pos.int.ref]
#	}
}

# get estimates and standard error for other levels of x given interaction with intvar
	var.cov<-vcov(model)
	if(length(grep(":", names(fixef(model)))) > 0){
		ests <- expand.grid(x = xvar.order,int = intvar.order,est = NA, se = NA)
		for(i in 1:nrow(ests)){
			#get positions of terms for x and int 
			pos.x <- grep(paste(xvar, ests$x[i], sep = ""), names(fixef(model)))
			pos.int <- grep(paste(intvar, ests$int[i], sep = ""), names(fixef(model)))
			# combine to get positions of interaction terms
			pos <- intersect(pos.x, pos.int)
			#get estimate
			ests$est[i] <- xests$est[which(xests$x == ests$x[i])] +
				# int.ests$est[which(int.ests$int == ests$int[i])]+ # dont need this as reference for all levels of intvar is set to 0
				 as.numeric(fixef(model)[pos])
    			cov <-var.cov[which(row.names(var.cov)== paste(xvar,ests$x[i], sep = "")), which(names(var.cov[1,])==names(fixef(model))[pos])]
			se1 <- xests$se[which(xests$x == ests$x[i])]  #standard error of the xvar term
			se2 <- se.fixef(model)[pos] #standard error of the interaction 
			ests$se[i] <- sqrt(se1^2 + se2^2 + 2*cov) #i.e. not concerned with error associated with estimate of level(s) of intvar
			}
		}

  reference <- data.frame(x = "reference", int = "reference", est = 0, se = 0)

if(length(grep(":", names(fixef(model)))) > 0){
	yData <- rbind(reference, int.ests, xests,  ests)
	}else{
	yData <- rbind(reference, xests)
	}

  y <- yData$est
  yplus<-yData$est + yData$se*seMultiplier
  yminus<-yData$est - yData$se*seMultiplier

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


xlims <- c(1, (length(xvar.order)+1)*(length(intvar.order)+1))
ylims <- c(min(yminus)-10, max(yplus) + 10)
par(mar = c(8,4,1,2), mgp = c(3,0.6,0))
plot(-1, -1, xlim = xlims, bty = "l",
		ylim = ylims, xlab = "", ylab = ylab, cex.lab = 1.2, axes = F)
axis(2,c(-50, 0, 50, 100), c(-50, 0, 50, 100))
mtext(side = 1, text = xlab, line =6, cex = 1.2)

shapes <- NULL
#make shapes for plot
if(intvar == "taxon_of_interest"){
	shapes <- rep(17, length(yData$int))
	}else if (intvar == "Zone"){
	shapes <- rep(15, length(yData$int))
	}else{
	shapes <- as.numeric(yData$int) + 14
	}

colours <- NULL
#make colours for plot
if(intvar == "taxon_of_interest"){
	colours[which(yData$int == "reference")]  <- p.col
	colours[which(yData$int == "Invertebrates")]  <- i.col
	colours[which(yData$int == "Vertebrates")]  <- v.col
	}else if (intvar == "Zone"){
	colours[which(yData$int == "reference")]  <- temp.col
	colours[which(yData$int == "reference")]  <- trop.col
	}else{
	colours <- rep(1,length(yData$int))
	}


p <- 1
for(L in c("reference",intvar.order)){
	xPositions <- p:(p -1 +length(y[which(yData$int == L )]))
	errbar(xPositions, y[which(yData$int == L )], 
		yplus[which(yData$int == L )], 
		yminus[which(yData$int == L )],
		col = colours[which(yData$int == L )],  
		pch = shapes[which(yData$int == L )],
		errbar.col = colours[which(yData$int == L )], 
		add = T, 
		cap = 0)
	if(L == "reference"){
	text(mean(xPositions), max(ylims), levels(data[,c(intvar)])[1])
	}else{
	text(mean(xPositions), max(ylims), L)
	}
	axis(1, p:(p -1 +length(y[which(yData$int == L )])), c(ref.text, xvar.order), 
		las = 2, cex.axis = 1, tick = 0)
	p <- p + length(y[which(yData$int == L )])
	abline(v = p - 0.5, col = 8)
}
if(!is.null(intvar)){
	nlabels <- aggregate(data[,c("SS")], by = list(xvar = data[,xvar], intvar = data[,intvar]), length)
	}else{
	nlabels <- aggregate(data[,c("SS")], by = list(xvar = data[,xvar]), length)
	}

text(1:max(xlims),min(ylims), nlabels$x)
abline(h = 0, lty = 2, col = 1)

}




