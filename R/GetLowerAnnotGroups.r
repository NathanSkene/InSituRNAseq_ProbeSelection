GetLowerAnnotGroups <- function(target,targLevel,levels){
    # Get annotations
    annot=levels[[targLevel]]$annot
    if(targLevel==1){a0 = NA}else{a0 = levels[[targLevel-1]]$annot}
    
    if(!is.na(a0)){
        targetBelow=unique(a0[annot==target])
        cellsFromBelow = a0==targetBelow
    }else{
        cellsFromBelow = rep(TRUE,length(annot))
    }    
    useGroups      = setdiff(unique(annot[cellsFromBelow]),target)
    return(list(otherGroups=useGroups,annotBelow=a0))
}