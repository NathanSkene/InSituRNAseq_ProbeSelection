# DESCRIPTION:
# - This function checks whether the proportions (of target and off-target cell types) obay certain constraints
# - It is intended as a method for limiting the number of gene set that need testing with LM
# - Basic idea is that average proportion should be high in target and low elsewhere

# EXAMPLE:
# - CheckGeneSetProportions(target=g,GeneSetExp=propAdj[curGenes,])

CheckGeneSetProportions <- function(target,GeneSetExp){
    # Multiply the proportions together...
    # - You want a probe to be approximately zero in every off target cell type
    iProps = apply(GeneSetExp*100,2,prod)#apply(GeneSetExp-apply(GeneSetExp,1,mean)/2,2,sum)
    validI = TRUE
    
    # Get the proportion data
    howMuchBigger = iProps[target] / max(iProps[setdiff(names(iProps),target)])
    mean_iProps = mean(iProps[setdiff(names(iProps),target)])
    targ_iProps = iProps[target]
    
    # Test the constraints
    if(iProps[target]!=max(iProps)){validI=FALSE}
    if(howMuchBigger<1){validI=FALSE} # 2 would be reasonable
    #if(mean_iProps < iProps[target]){validI=FALSE} 
    if(targ_iProps < 5*mean_iProps){validI=FALSE} # 5 would be reasonable
    
    # Return findings
    return(data.frame(validI=validI,howMuchBigger=howMuchBigger,mean_iProps=mean_iProps,targ_iProps=targ_iProps))
}