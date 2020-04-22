GetGeneSetOptions <- function(candidates,target){
    cands   = candidates
    candsG  = as.character(cands[cands$target==target,2])
    options = combn(candsG,3)
    return(options)
}