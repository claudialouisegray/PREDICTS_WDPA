#names(pairwise_log_abundance)
#pairwise_log_abundance$IUCN_CAT_number <- factor(pairwise_log_abundance$IUCN_CAT_number)

#dataset <- matched.landuse.s
#responseVar <- "log_abundance"
#fitFamily <- "gaussian"
#fixedFactors= fF
#fixedTerms = fT
#fixedInteractions= character(0) # c("Predominant_habitat:Realm")
#otherRandoms= "SS_PH"
#fitInteractions=FALSE
#fixed_RandomSlopes = c("log_GIS_AREA", "DoP_yr") # give the terms that should be fitted as a random slope, i.e. the ones that a study level gradient could be calculated for
#verbose=TRUE



# NB should really make all either be glmer or lmer, but this runs fine for now


compare_randoms <-function(dataset,responseVar,fixedFactors=character(0),
                          fixedTerms=list(),fixedInteractions=character(0),
				  fitFamily = "gaussian",
                          otherRandoms=character(0),
				  siteRandom = FALSE,
				  fixed_RandomSlopes = character(0),
                          fitInteractions=FALSE,verbose=FALSE,
				  keepVars = list()){
  
  if ((length(fixedInteractions)>0) & (fitInteractions)){
    stop("Error: specifying particular interactions and all two-way interactions will not work!")
  }

  
  
  inter.split<-unlist(strsplit(fixedInteractions,":"))
  inter.split<-gsub("poly[(]","",inter.split)
  inter.mains<-gsub("[,][:0-9:][)]","",inter.split)
  
  # Subset the data frame, retaining just the relevant columns
  # include Source_ID and remove blocks as not relevant
  model.data<-subset(dataset,select=c("Source_ID","SS","SS_PH","SSB","SSBS",fixedFactors,
							names(fixedTerms),responseVar,otherRandoms, names(keepVars)))

  model.data<-na.omit(model.data)
  cat<-sapply(model.data,is.factor)
  model.data[cat]<-lapply(model.data[cat],factor)
  
  # Create a list to store the results
  results<-list(ranef=character(),AIC=numeric())

  # Create a list to store warnings
  all.warnings <- list()
  
  # Construct the fixed-effects part of the original model call
  fixedStruct<-""
  for (i in 1:length(fixedFactors)){
    fixedStruct<-paste(fixedStruct,fixedFactors[i],sep="")
    if ((i != length(fixedFactors)) | (length(fixedTerms)>0) | 
          ((length(fixedTerms)==0) & (length(fixedInteractions)>0))){
      fixedStruct<-paste(fixedStruct,"+",sep="")
    }
  }
  if (length(fixedTerms)>0){
    for (i in 1:length(fixedTerms)){
      fixedStruct<-paste(fixedStruct,"poly(",names(fixedTerms)[i],
                         ",",fixedTerms[i],")",sep="")
      if ((i != length(fixedTerms)) | (length(fixedInteractions)>0)){
        fixedStruct<-paste(fixedStruct,"+",sep="")
      }
    }
  }
 
 if (length(keepVars)>0){
    for (i in 1:length(keepVars)){
      term<-paste("poly(",names(keepVars)[i],
                  ",",keepVars[i],")",sep="")
      fixedStruct<-paste(fixedStruct,term,sep="")
      allTerms<-c(allTerms,term)
	keepTerms <- c(keepTerms, term)
      if ((i != length(keepVars)) | (length(fixedInteractions)>0)){
        fixedStruct<-paste(fixedStruct,"+",sep="")
      }
    }
  }

  if (fitInteractions) fixedStruct<-paste("(",fixedStruct,")^2",sep="")
  
  if (length(fixedInteractions)>0){
    for (i in 1:length(fixedInteractions)){
      fixedStruct<-paste(fixedStruct,fixedInteractions[i],sep="")
      if (i != length(fixedInteractions)){
        fixedStruct<-paste(fixedStruct,"+",sep="")
      }
    }
  }
  
  
 
  print(paste("Fixed structure:",fixedStruct))
  
  # Set the original value for the best model AIC
  best.aic<-Inf
  
  # First, try fitting just landuse group within study as a random effect
  print(paste("Comparing random structure 1 of ",length(fixed_RandomSlopes)+
                length(otherRandoms)+
                3,sep=""))
  
  new.random<-"(1|SS)" #switched back from SS_PH
  best.random<-new.random
  
  # Construct the complete model call
  new.call<-construct_call(responseVar,fixedStruct,new.random)
  if (verbose) print(new.call)
  
  # Run the model
  # using try() allows for the model to fail and gives back any error messages
  if(fitFamily == "gaussian"){
    mod<-try(lmer(new.call,data=model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=1000000))))
	}else{
    mod<-try(glmer(new.call,data=model.data, family = fitFamily, control= glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=1000000))))
	}  

  
  # add warnings to the list
  all.warnings[[1]] <- list(call = paste(new.call), warnings = mod@optinfo$conv$lme4$messages)
  

  # Check whether this model has the best AIC value and update the results
  if ((class(mod)!="try-error")){
    if (!is.na(AIC(mod)) & !is.null(AIC(mod))){
      results$ranef<-c(results$ranef,new.random)
      results$AIC<-c(results$AIC,AIC(mod))
      
      if(AIC(mod)<best.aic){
        best.aic<-AIC(mod)
        best.random<-new.random
      }
    }
  }
  
  
  
  
#    # add SS to random structure
#  
#    print(paste("Comparing random structure 2 of ",length(fixed_RandomSlopes)+
#                  length(otherRandoms)+
#                  4,sep=""))
#  
#    new.random <- gsub("1[|]","1|SS/",best.random)
#  
#    # Construct the complete model call
#    new.call<-construct_call(responseVar,fixedStruct,new.random)
#    if (verbose) print(new.call)
#    
#    # Run the model
#    # using try() allows for the model to fail and gives back any error messages
#    mod<-try(lmer(new.call,data=model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=1000000))))
#  
#  # add warnings to the list
#  all.warnings[[2]] <- list(call = paste(new.call), warnings = mod@optinfo$conv$lme4$messages)
#  
#  
#    
#    # Check whether this model has the best AIC value and update the results
#    if ((class(mod)!="try-error")){
#      if (!is.na(AIC(mod)) & !is.null(AIC(mod))){
#        results$ranef<-c(results$ranef,new.random)
#        results$AIC<-c(results$AIC,AIC(mod))
#        
#        if(AIC(mod)<best.aic){
#          best.aic<-AIC(mod)
#          best.random<-new.random
#        }
#      }
#    }
#  
  


#  # Try adding Source_ID to the random structure
#  print(paste("Comparing random structure 3 of ",length(fixed_RandomSlopes)+
#                length(otherRandoms)+
#                4,sep=""))
#  
#  new.random <- gsub("1[|]","1|Source_ID/",best.random)
#  
#  # Construct the complete model call
#  new.call<-construct_call(responseVar,fixedStruct,new.random)
#  if (verbose) print(new.call)
#  
#  # Run the model
#  mod<-try(lmer(new.call,data=model.data, control= lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000))))
#
#  # add warnings to the list
#  all.warnings[[3]] <- list(call = paste(new.call), warnings = mod@optinfo$conv$lme4$messages)
#  
#
#  # Check whether this model has the best AIC value and update the results
#
#  if ((class(mod)!="try-error")){
#    if (!is.na(AIC(mod)) & !is.null(AIC(mod))){
#      results$ranef<-c(results$ranef,new.random)
#      results$AIC<-c(results$AIC,AIC(mod))
#      
#      if(AIC(mod)<best.aic){
#        best.aic<-AIC(mod)
#        best.random<-new.random
#      }
#    }
#  }
#  

    
  if(siteRandom){
    print(paste("Comparing random structure 2 of " ,length(fixed_RandomSlopes)+
                length(otherRandoms) + 3,sep=""))
    
    new.random <- paste(best.random, "+ (1|SSBS)",sep="")
   
    
    # Construct the complete model call
    new.call<-construct_call(responseVar,fixedStruct,new.random)
    if (verbose) print(new.call)
    
  # Run the model
   if(fitFamily == "gaussian"){
    mod<-try(lmer(new.call,data=model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=1000000))))
	}else{
    mod<-try(glmer(new.call,data=model.data, family = fitFamily, control= glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=1000000))))
  }

  # add warnings to the list
  all.warnings[[2]] <- list(call = paste(new.call), warnings = mod@optinfo$conv$lme4$messages)
  
    
    # Check whether this model has the best AIC value and update the results
    if ((class(mod)!="try-error")){
      if (!is.na(AIC(mod)) & !is.null(AIC(mod))){
        results$ranef<-c(results$ranef,new.random)
        results$AIC<-c(results$AIC,AIC(mod))
        
        if(AIC(mod)<best.aic){
          best.aic<-AIC(mod)
          best.random<-new.random
        }
      }
    }
    
    
  }





# For each specified term, try adding the factor in a random slopes-intercept model

  i<-1
  for (f in fixed_RandomSlopes){
    print(paste("Comparing random structure ",2+i," of ",
                length(fixed_RandomSlopes)+length(otherRandoms)+3,sep=""))

    
   
if(siteRandom){
	randomStruct <- strsplit(best.random, "[+] ")
	new.random1 <- gsub("[|]SS",paste("+",f,"|SS",sep=""),randomStruct[[1]][1])
	new.random <- paste(new.random1, "+ (1|SSBS)", sep="")
	}else{
    new.random <- gsub("[|]",paste("+",f,"|",sep=""),best.random)
	}
     if (verbose) print(new.random)  

   
   
    # Construct the complete model call
    new.call<-construct_call(responseVar,fixedStruct,new.random)
    if (verbose) print(new.call)


    # Run the model
     if(fitFamily == "gaussian"){
    mod<-try(lmer(new.call,data=model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=1000000))))
	}else{
    mod<-try(glmer(new.call,data=model.data, family = fitFamily, control= glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=1000000))))
  }

  # add warnings to the list
  all.warnings[[2+i]] <- list(call = paste(new.call), warnings = mod@optinfo$conv$lme4$messages)
  
    
    # Check whether this model has the best AIC value and update the results
    if ((class(mod)!="try-error")){
      if (!is.na(AIC(mod)) & !is.null(AIC(mod))){
        results$ranef<-c(results$ranef,new.random)
        results$AIC<-c(results$AIC,AIC(mod))
        
        if(AIC(mod)<best.aic){
          best.aic<-AIC(mod)
          best.random<-new.random
        }
      }
    }
    
    i <- i+1
  }
 



   # Try adding block to the random structure
  print(paste("Comparing random structure ", 2 + length(fixed_RandomSlopes) + 1, " of ",length(fixed_RandomSlopes)+
                length(otherRandoms)+
                3,sep=""))

    new.random <- paste(best.random, "+ (1|SSB)",sep="")
    
  # Construct the complete model call
  new.call<-construct_call(responseVar,fixedStruct,new.random)
  if (verbose) print(new.call)
  
  # Run the model
  mod<-try(lmer(new.call,data=model.data, control= lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000))))

  # add warnings to the list
  all.warnings[[4]] <- list(call = paste(new.call), warnings = mod@optinfo$conv$lme4$messages)
 

  # Check whether this model has the best AIC value and update the results
  if ((class(mod)!="try-error")){
    if (!is.na(AIC(mod)) & !is.null(AIC(mod))){
      results$ranef<-c(results$ranef,new.random)
      results$AIC<-c(results$AIC,AIC(mod))
      
      if(AIC(mod)<best.aic){
        best.aic<-AIC(mod)
        best.random<-new.random
      }
    }
  }





  # Try adding any other specified randoms

  i <- 1

  for (re in otherRandoms){

    print(paste("Comparing random structure ",3+ length(fixed_RandomSlopes)+i," of ",
                length(fixed_RandomSlopes)+
                  length(otherRandoms)+3,sep=""))


    new.random<-paste(best.random,"+(1|",re,")",sep="")
        
    # Construct the complete model call
    new.call<-construct_call(responseVar,fixedStruct,new.random)
    if (verbose) print(new.call)

  # Run the model
   if(fitFamily == "gaussian"){
    mod<-try(lmer(new.call,data=model.data, control= lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=1000000))))
	}else{
    mod<-try(glmer(new.call,data=model.data, family = fitFamily, control= glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=1000000))))
  }

  # add warnings to the list
  all.warnings[[2 + length(otherRandoms)+ i]] <- list(call = paste(new.call), warnings = mod@optinfo$conv$lme4$messages)
  

    # Check whether this model has the best AIC value and update the results
    if ((class(mod)!="try-error")){
      if (!is.na(AIC(mod)) & !is.null(AIC(mod))){
        results$ranef<-c(results$ranef,new.random)
        results$AIC<-c(results$AIC,AIC(mod))
        
        if(AIC(mod)<best.aic){
          best.aic<-AIC(mod)
          best.random<-new.random
        }
      }
    }
    
    i<-i+1
  }
  
 
  
  return(list(full.results=results,best.random=best.random, warnings = all.warnings))
  
}