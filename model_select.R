model_select<-function(all.data,responseVar,
				alpha = 0.05,
                        fitFamily = "gaussian",
                        fixedFactors= character(0),
                        fixedTerms=character(0),
                       fixedInteractions=character(0),
                       randomStruct,
                       siteRandom=FALSE,
                       fitInteractions=FALSE,verbose=FALSE,
                       otherRandoms=character(0)){
  
  contEffectNames<-names(fixedTerms)
  
  model.data<-subset(all.data,select=c("SS", "SS_PH","SSB","SSBS",fixedFactors,
                                    names(fixedTerms),responseVar,otherRandoms))
    
  model.data<-na.omit(model.data)
  cat<-sapply(model.data,is.factor)
  model.data[cat]<-lapply(model.data[cat],factor)
  
  
  for (fe in fixedFactors){
    eval(substitute(model.data$x<-factor(model.data$x),list(x=fe)))
  }
  
  results<-list(fixef=character(),AIC=numeric())
  all.warnings <- list()
  
  allTerms<-character(0)
  fixedStruct<-""
  for (i in 1:length(fixedFactors)){
    fixedStruct<-paste(fixedStruct,fixedFactors[i],sep="")
    allTerms<-c(allTerms,fixedFactors[i])
    if ((i != length(fixedFactors)) | (length(fixedTerms)>0) | 
          ((length(fixedTerms)==0) & (
            length(fixedInteractions)>0))){
      fixedStruct<-paste(fixedStruct,"+",sep="")
    }
  }
  if (length(fixedTerms)>0){
    for (i in 1:length(fixedTerms)){
      term<-paste("poly(",names(fixedTerms)[i],
                  ",",fixedTerms[i],")",sep="")
      fixedStruct<-paste(fixedStruct,term,sep="")
      allTerms<-c(allTerms,term)
      if ((i != length(fixedTerms)) | (length(fixedInteractions)>0)){
        fixedStruct<-paste(fixedStruct,"+",sep="")
      }
    }
  }
  
  if (fitInteractions){
    fixedStruct<-paste(fixedStruct,"+")
    
    mainTerms<-allTerms
    
    for (i in 1:(length(mainTerms)-1)){
      for (j in (i+1):length(mainTerms)){
        term<-paste(mainTerms[i],mainTerms[j],sep=":")
        fixedStruct<-paste(fixedStruct,term)
        allTerms<-c(allTerms,term)
      }
    }
    
  }
  
  if (length(fixedInteractions)>0){
    for (i in 1:length(fixedInteractions)){
      fixedStruct<-paste(fixedStruct,fixedInteractions[i],sep="")
      allTerms<-c(allTerms,fixedInteractions[i])
      if (i != length(fixedInteractions)){
        fixedStruct<-paste(fixedStruct,"+",sep="")
      }
    }
  }
  
  randomStruct<-gsub(" ","",randomStruct)
  
  call.old<-paste(responseVar,"~",allTerms[1],sep="")
  if (length(allTerms)>1){
    for (t in 2:length(allTerms)){
      call.old<-paste(call.old,"+",allTerms[t],sep="")
    }
  }
  call.old<-paste(call.old,"+ ",randomStruct,sep="")
  
 
  #make stats table
  inters <- allTerms[grep(":", allTerms)]
  if(length(inters) > 0){
  	mainTerms <- allTerms[-grep(":", allTerms)]
	}else{
	mainTerms <- allTerms}
  polyTerms<-mainTerms[grep("poly",mainTerms)]
  c.polyTerms.to.simplify <- polyTerms[grep("3", polyTerms)]
  q.polyTerms.to.simplify <- polyTerms[grep("2", polyTerms)]
  c.polyTerms.to.add1 <-gsub(",3",",2",c.polyTerms.to.simplify)
  c.polyTerms.to.add2 <-gsub(",3",",1",c.polyTerms.to.simplify)
  q.polyTerms.to.add <-gsub(",2",",1",q.polyTerms.to.simplify)
  mainTerms <- c(mainTerms,c.polyTerms.to.add1,c.polyTerms.to.add2, q.polyTerms.to.add)
  mainTerms <- mainTerms[order(mainTerms)]	
  terms <- c(mainTerms, inters)
  stats<-data.frame(terms=terms)
  stats$terms<-paste(stats$terms)
  stats$ChiSq<-NA
  stats$Df<-NA
  stats$P<-NA
  stats$dBIC<-NA
  stats$dAIC<-NA
  
  iter<-1
  
  
  repeat {
    
    print(paste("Performing round ",iter," of interaction-term removal",sep=""))
    
    if (verbose) print(call.old)
      
    if(fitFamily == "gaussian"){
      mOld<-lmer(call.old,data=model.data, REML = F, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
    }else{
      mOld<-glmer(call.old,data=model.data, family = fitFamily, control= glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
    }
    
    # add warnings to list 
    if(is.null(mOld@optinfo$conv$lme4$messages) ==F){
    all.warnings[[paste(iter, "old.int")]] <- list(call = paste(call.old), warnings = mOld@optinfo$conv$lme4$messages)
    print(mOld@optinfo$conv$lme4$messages)}
    
    
    if (iter==1) iTerms <- allTerms[grep(":",allTerms)]
    
    iTerms<-gsub(" ","",iTerms)
    
    pVals<-numeric()
    Chis<-numeric()
    Dfs<-character()
    dBICs<-numeric()
    dAICs<-numeric()
    
    for (t in iTerms){
      
        t1<-gsub("[(]","[(]",t)
        t2<-gsub("[)]","[)]",t1)
        t3<-paste(t2,"[+]",sep="")
        
        call.new<-gsub(t3,"",call.old)
      
        if(fitFamily == "gaussian"){
          mNew <-lmer(call.new,data=model.data, REML = F, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
        }else{
          mNew<-glmer(call.new,data=model.data, family = fitFamily, control= glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
        }
        
        # add warnings to list 
        if(is.null(mNew@optinfo$conv$lme4$messages)==F){
        all.warnings[[paste(iter, "new.int", t)]] <- list(call = paste(call.new), warnings = mNew@optinfo$conv$lme4$messages)
	  print(mNew@optinfo$conv$lme4$messages)}
        
      pVals<-c(pVals,anova(mOld,mNew)$Pr[2])
      Chis<-c(Chis,anova(mOld,mNew)$Chisq[2])
      Dfs<-c(Dfs,paste(anova(mOld,mNew)$'Chi Df'[2],",",anova(mOld,mNew)$Df[2]))
      dBICs<-c(dBICs,BIC(mNew)-BIC(mOld))
      dAICs<-c(dAICs,AIC(mNew)-AIC(mOld))
    }
    
    if (verbose){
      print(iTerms)
      print(pVals)
    }
    
    print(paste(length(which(pVals>alpha))," interaction terms have P-values > ", alpha ,sep=""))
    
    if (length(which(pVals>alpha))==0) break
    
    dropI<-iTerms[order(pVals)[length(order(pVals))]]
    
    stats$ChiSq[which(stats$terms==dropI)]<-Chis[order(pVals)[length(order(pVals))]]
    stats$Df[which(stats$terms==dropI)]<-Dfs[order(pVals)[length(order(pVals))]]
    stats$P[which(stats$terms==dropI)]<-pVals[order(pVals)[length(order(pVals))]]
    stats$dBIC[which(stats$terms==dropI)]<-dBICs[order(pVals)[length(order(pVals))]]
    stats$dAIC[which(stats$terms==dropI)]<-dAICs[order(pVals)[length(order(pVals))]]
    
    print(paste("Dropping ",dropI,sep=""))
    
      t1<-gsub("[(]","[(]",dropI)
      t2<-gsub("[)]","[)]",t1)
      t3<-paste(t2,"[+]",sep="")
      
      call.old<-gsub(t3,"",call.old)
    
    
    iTerms<-iTerms[-order(pVals)[length(order(pVals))]]
    allTerms<-allTerms[-which(allTerms==dropI)]
    
    iter<-iter+1 
    
  }
  
  stats$ChiSq[na.omit(match(iTerms,stats$terms))]<-Chis
  stats$Df[na.omit(match(iTerms,stats$terms))]<-Dfs
  stats$P[na.omit(match(iTerms,stats$terms))]<-pVals
  stats$dBIC[na.omit(match(iTerms,stats$terms))]<-dBICs
  stats$dAIC[na.omit(match(iTerms,stats$terms))]<-dAICs 
  
  #remove all the interaction terms from the call
    for (t in iTerms){
      t1<-gsub("[(]","[(]",t)
      t2<-gsub("[)]","[)]",t1)
      t3<-paste(t2,"[+]",sep="")
      
      call.old<-gsub(t3,"",call.old)
      
    }

  
  # remove remaining interaction terms from allTerms
  itersRemaining <-which(allTerms %in% iTerms)
  if (length(itersRemaining)>0) allTerms<-allTerms[-itersRemaining]
  
  iter<-1
  
  repeat {
    
    print(paste("Performing round ",iter," of main-effect removal",sep=""))
    
    if (verbose) print(call.old)
    
    
    if(fitFamily == "gaussian"){
      mOld<-lmer(call.old,data=model.data, REML = F, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
    }else{
      mOld<-glmer(call.old,data=model.data, family = fitFamily, control= glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
    }
    
    # add warnings to list 
    if(is.null(mOld@optinfo$conv$lme4$messages) == F){
    all.warnings[[paste(iter, "old.main")]] <- list(call = paste(call.old), warnings = mOld@optinfo$conv$lme4$messages)
    print(mOld@optinfo$conv$lme4$messages)}
    
    mTerms<-allTerms
    
    mTerms<-gsub(" ","",mTerms)
    
    pVals<-numeric()
    Chis<-numeric()
    Dfs<-character()
    dBICs<-numeric()
    dAICs <- numeric()
    
    for (t in mTerms){

    #remove the term from the call
   
	if ((grepl("poly",t)) & grepl(",3",t)){
        t1<-gsub("[(]","[(]",t)
        t2<-gsub("[)]","[)]",t1)
        t3<-gsub(",3",",2",t)
        
        call.new<-gsub(t2,t3,call.old)
        
      } else if ((grepl("poly",t)) & grepl(",2",t)){
        t1<-gsub("[(]","[(]",t)
        t2<-gsub("[)]","[)]",t1)
        t3<-gsub(",2",",1",t)
        
        call.new<-gsub(t2,t3,call.old)
        
      } else {

      t1<-gsub("[(]","[(]",t)
      t2<-gsub("[)]","[)]",t1)
      t3<-paste(t2,"[+]",sep="")

      #drop the term, but avoid losing it from random structure
      # this only works because there is a space after the "+" so it splits the call in two
      # the "+ " was added specifically for this in the model call building stage
      struct.split <- strsplit( call.old, "+ ")
      struct.new <-gsub(t3,"",struct.split[[1]][1])
      call.new <- paste(struct.new," ",struct.split[[1]][2], sep = "")
    
      }

	if(verbose){print(call.new)}


    # fit the new model
      if(fitFamily == "gaussian"){
        mNew <-lmer(call.new,data=model.data, REML = F, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
      }else{
        mNew<-glmer(call.new,data=model.data, family = fitFamily, control= glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
      }
      
      # add warnings to list 
    if(is.null(mNew@optinfo$conv$lme4$messages)==F){
      all.warnings[[paste(iter, "new.main",t)]] <- list(call = paste(call.new), warnings = mNew@optinfo$conv$lme4$messages)
	print(mNew@optinfo$conv$lme4$messages)}
     
      #add the stats from dropping each term to a vector (each iteration adds a new value) 
      pVals<-c(pVals,anova(mOld,mNew)$Pr[2])
      Chis<-c(Chis,anova(mOld,mNew)$Chisq[2])
      Dfs<-c(Dfs,paste(anova(mOld,mNew)$'Chi Df'[2],",",anova(mOld,mNew)$Df[2]))
      dBICs<-c(dBICs,BIC(mNew)-BIC(mOld))
      dAICs<-c(dAICs,AIC(mNew)-AIC(mOld))
    }
    
    print(paste(length(which(pVals>alpha))," candidate main effects have P-values > ", alpha,sep=""))
    
    if (verbose){
      print(mTerms)
      print(pVals)
    }
    
    if (length(which(pVals>alpha))==0) break
    
    dropM<-mTerms[order(pVals)[length(order(pVals))]]
    
    stats$ChiSq[which(stats$terms==dropM)]<-Chis[order(pVals)[length(order(pVals))]]
    stats$Df[which(stats$terms==dropM)]<-Dfs[order(pVals)[length(order(pVals))]]
    stats$P[which(stats$terms==dropM)]<-pVals[order(pVals)[length(order(pVals))]]
    stats$dBIC[which(stats$terms==dropM)]<-dBICs[order(pVals)[length(order(pVals))]]
    stats$dAIC[which(stats$terms==dropM)]<-dAICs[order(pVals)[length(order(pVals))]]
    

  if ((grepl("poly",dropM)) & grepl(",3",dropM)){
      print(paste("Simplifying ",dropM,sep=""))
      
      d1<-gsub("[(]","[(]",dropM)
      d2<-gsub("[)]","[)]",d1)
      d3<-gsub(",3",",2",dropM)
      
      call.old<-gsub(d2,d3,call.old)
      
      mTerms<-gsub(d2,d3,mTerms)
      allTerms<-gsub(d2,d3,allTerms)
      
    } else if ((grepl("poly",dropM)) & grepl(",2",dropM)){
      print(paste("Simplifying ",dropM,sep=""))
      
      d1<-gsub("[(]","[(]",dropM)
      d2<-gsub("[)]","[)]",d1)
      d3<-gsub(",2",",1",dropM)
      
      call.old<-gsub(d2,d3,call.old)
      
      mTerms<-gsub(d2,d3,mTerms)
      allTerms<-gsub(d2,d3,allTerms)
      
    } else {

      print(paste("Dropping ",dropM,sep=""))
      t1<-gsub("[(]","[(]",dropM)
      t2<-gsub("[)]","[)]",t1)
      t3<-paste(t2,"[+]",sep="")

    #drop the term, but avoid losing it from random structure
    struct.split <- strsplit( call.old, "+ ")
    struct.new <-gsub(t3,"",struct.split[[1]][1])
    call.old <- paste(struct.new," ",struct.split[[1]][2], sep = "")
    

    # take the dropped term out of mTerms and allTerms
    mTerms<-mTerms[-order(pVals)[length(order(pVals))]]
    allTerms<-allTerms[-which(allTerms==dropM)]
      
    }

      
    iter<-iter+1 
    
    
  }
  
  stats$ChiSq[na.omit(match(mTerms,stats$terms))]<-Chis
  stats$Df[na.omit(match(mTerms,stats$terms))]<-Dfs
  stats$P[na.omit(match(mTerms,stats$terms))]<-pVals
  stats$dBIC[na.omit(match(mTerms,stats$terms))]<-dBICs
  stats$dAIC[na.omit(match(mTerms,stats$terms))]<-dAICs
  
  
  fixedStruct<-""
  
  sig.terms<-stats[stats$P<alpha,]
  
  
  sig.terms<-na.omit(sig.terms)
  
  if (dim(sig.terms)[1]>0){
    
    	sig.terms<-paste(sig.terms$terms)
    	sig.inter<-sig.terms[grepl(":",sig.terms)]

	# add main terms that are in an interaction but not otherwise present
    	inter.mains<-unique(unlist(strsplit(sig.inter,":")))
    	sig.terms<-c(sig.terms,inter.mains[!(inter.mains %in% sig.terms)]) 

	#take out any lower order polynomials that are in the main terms...
    for (t in names(fixedTerms)){
      mainMatches<-sig.terms[which((grepl(paste("poly[(]",t,",[0-9]{1}[)]",sep=""),sig.terms)) & 
                                     !(grepl(":",sig.terms)))]
      mainMatchPosits<-which((grepl(paste("poly[(]",t,",[0-9]{1}[)]",sep=""),sig.terms)) & 
                               !(grepl(":",sig.terms)))
      if (length(mainMatches)>1){
        sig.terms<-sig.terms[-mainMatchPosits[order(mainMatches,decreasing=TRUE)][-1]]
      }
    }

    for (i in 1:length(sig.terms)){
      fixedStruct<-paste(fixedStruct,sig.terms[i],sep="")
      if (i != length(sig.terms)) fixedStruct<-paste(fixedStruct,"+",sep="")
    }
    call.best<-construct_call(responseVar,fixedStruct,randomStruct)
    if (verbose) print(call.best)
    
    if(fitFamily == "gaussian"){
      mBest<-lmer(call.best,data=model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
    }else{
      mBest<-glmer(call.best,data=model.data, family = fitFamily, control= glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
    }
    
    # cat("Estimating the influence of different studies in the model\n")
    # infl<-influence(model=mBest,group="SS")
    # cook<-cooks.distance.estex(infl,sort=TRUE)

    return(list(model=mBest,data=model.data,stats=stats,final.call=call.best, fixedStruct = fixedStruct, randomStruct = randomStruct, warnings = all.warnings))

  } else {

    print("Warning: all terms were dropped from the model")
    call.best<-construct_call(responseVar,"1",randomStruct)
    if (verbose) print(call.best)
    
    if(fitFamily == "gaussian"){
      mBest<-lmer(call.best,data=model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
    }else{
      mBest<-glmer(call.best,data=model.data, family = fitFamily, control= glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
    }
    
    return(list(model=mBest,data=model.data,stats=stats,final.call=call.best, warnings = all.warnings))  
  }
  
}