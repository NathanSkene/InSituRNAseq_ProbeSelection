# Applies constraints to proportions and scores to obtain list of candidate genes for a given cell class

# INPUTS:
# target: name of cell group, i.e. "Neurons"
# groupMeans:  means of cell groups for the appropriate level of annotation
# proportions: proportion of expression within the cell groups of that level
# annot:       annotations for level of target
# annotBelow:  annotations for the level below target

# EXAMPLE:
# SelectCandidates(target="Neurons",groupMeans=levels[[1]]$stats$groupMedians,proportions=levels[[1]]$stats$proportions)
SelectCandidatesForPairs <- function(target,groupMeans,proportions,annot,annotBelow,numCandidates=20){
    if(!is.na(annotBelow)){
        targetBelow=unique(annotBelow[annot==target])
        cellsFromBelow = annotBelow==targetBelow
    }else{
        cellsFromBelow = rep(TRUE,length(annot))
    }
    useGroups      = unique(annot[cellsFromBelow])  
    offTargGroups  = setdiff(useGroups,target)
    
    candidates=data.frame()
    for(i in 1:length(offTargGroups)){
        pairGroups = c(target,offTargGroups[i])
        
        # If this level of definition does not provide any further information... then ignore it
        if(length(useGroups)==1){return(NA)}
        propAdj        = proportions[,pairGroups]/(apply(proportions[,pairGroups],1,sum)+0.000000001)
        
        # Check for errors
        if(!is.na(annotBelow)){
            if(length(targetBelow)>1){stop("targetBelow should have length of 1")}
        }
        
        # Get score: (Multiply proportion by group mean)
        # score = propAdj * groupMeans[,useGroups]    
        
        # Use instead the scoring metric that Kenneth suggested
        #g="Tspan12"
        #offTargV = groupMeans[g,setdiff(useGroups,target)]
        #targV = groupMeans[g,target]    
        #score = (targV-sum(offTargV))/(sum(groupMeans[g,])+10)    
        
        score=rep(0,dim(expDat$readCounts)[1])
        vAdj=0
        nC = numCandidates
        while(nC>0){
            vAdj = vAdj + 0.1
            
            # Get 90th percentile highest expression in off-target cell group
            offtarg_exp = expDat$readCounts[,annot==offTargGroups[i]]
            offtarg_q = apply(offtarg_exp,1,quantile,probs=0.90)
            # Get 10th percentile lowest expression in target group
            targ_exp = expDat$readCounts[,annot==target]
            targ_q = apply(targ_exp,1,quantile,probs=vAdj)     
    
            #offTargV = groupMeans[,offTargGroups[i]]
            #targV = groupMeans[,target]    
            numerator = (targ_q-offtarg_q) 
            #numerator[numerator<0]=0
            denom = (apply(cbind(targ_q,offtarg_q),1,sum)+10)
            score = numerator / denom 
            #sort(score,decreasing=TRUE)[1:10]
            
            scoreOut = score[!(names(score) %in% as.character(candidates$genes))]
            
            candsToUse = nC
            if(sum(scoreOut>0)<nC){candsToUse=sum(scoreOut>0)}
            
            if(sum(scoreOut>0)>0){
                candidates = rbind(candidates,data.frame(target=target,offTarg=offTargGroups[i],vAdj=vAdj,genes=names(sort(scoreOut,decreasing=TRUE)[1:candsToUse])))
            }
            
            nC = nC - candsToUse
        }
        
        #nC = numCandidates
        #if(sum(score>0)<nC){
        #    nC = sum(score>0)
        #}
        
        # Return the results

    }
    return(candidates)
}