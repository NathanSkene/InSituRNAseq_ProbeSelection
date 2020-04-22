# Applies constraints to proportions and scores to obtain list of candidate genes for a given cell class

# INPUTS:
# target: name of cell group, i.e. "Neurons"
# groupMeans:  means of cell groups for the appropriate level of annotation
# proportions: proportion of expression within the cell groups of that level
# annot:       annotations for level of target
# annotBelow:  annotations for the level below target

# EXAMPLE:
# SelectCandidates(target="Neurons",groupMeans=levels[[1]]$stats$groupMedians,proportions=levels[[1]]$stats$proportions)
SelectCandidates <- function(target,groupMeans,proportions,annot,annotBelow,numCandidates=20){
    if(!is.na(annotBelow)){
        targetBelow=unique(annotBelow[annot==target])
        cellsFromBelow = annotBelow==targetBelow
    }else{
        cellsFromBelow = rep(TRUE,length(annot))
    }
    useGroups      = unique(annot[cellsFromBelow])  
    # If this level of definition does not provide any further information... then ignore it
    if(length(useGroups)==1){return(NA)}
    propAdj        = proportions[,useGroups]/apply(proportions[,useGroups],1,sum)
    
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
    
    offTargV = groupMeans[,setdiff(useGroups,target)]
    targV = groupMeans[,target]    
    if(length(setdiff(useGroups,target))>1){
        numerator = (targV-apply(offTargV,1,sum))
    }else{
        numerator = targV-offTargV
    }
    denom = (apply(groupMeans,1,sum)+10)
    score =  numerator / denom 
    #sort(score,decreasing=TRUE)[1:10]
    
    # Extract key statistics
    scoreOut = score#[,target]
    targ_Proportions = propAdj[,target]
    numOffTarg = sum(colnames(score) != target)
    if(sum(colnames(score) != target)>1){
        min_OffTarg_Proportions = apply(propAdj[,colnames(score) != target],1,min)
    }else{
        min_OffTarg_Proportions = propAdj[,colnames(score) != target]
    }
    
    # Apply selection critera
    #scoreOut[targ_Proportions<(0.50/((numOffTarg+1)/2))]=0
    #scoreOut[min_OffTarg_Proportions>(0.2/numOffTarg)]=0
    #scoreOut[groupMeans[,target]<=2]=0
    
    # Check if numCandidates is higher than possible candidates
    if(numCandidates>sum(scoreOut>0)){numCandidates=sum(scoreOut>0)}
    
    # Return the results
    candidates = data.frame(target=target,genes=names(sort(scoreOut,decreasing=TRUE)[1:numCandidates]))
    return(candidates)
}