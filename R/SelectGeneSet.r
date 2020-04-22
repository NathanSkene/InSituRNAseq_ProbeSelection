

SelectGeneSet <- function(target,targLevel,levels,readCounts){
    # Find candidate gene triplets for this celltype
    options = GetGeneSetOptions(target=target,candidates=levels[[targLevel]]$candidates)
    
    # Initial variables
    GeneSetSummary=data.frame()
    
    # Get group proportions
    proportions = levels[[targLevel]]$stats$proportions
    annotGroups=GetLowerAnnotGroups(target,targLevel,levels)
    allGroups = c(target,annotGroups$otherGroups)
    propAdj = proportions[,allGroups]/apply(proportions[,allGroups],1,sum)        
    
    for(i in 1:dim(options)[2]){
    #for(i in 1:5){
        if(i%%50==0){print(i/dim(options)[2])}
        curGenes=options[,i]
        
        # Adjust read counts for lower efficiency of In Situ Sequencing
        if(targLevel<=2){curExp  = round(readCounts[curGenes,]/100)}
        if(targLevel==3){curExp  = round(readCounts[curGenes,]/20)}
        if(targLevel>=4){curExp  = round(readCounts[curGenes,]/7)}
        
        # Initial check on proportions
        propChecks = CheckGeneSetProportions(target=target,GeneSetExp=propAdj[curGenes,])
        #propChecks$validI=TRUE
        
        if(propChecks$validI==TRUE){
        #if(TRUE==TRUE){
            # Run linear model and calculate prediction error
            if(targLevel==1){a0 = NA}else{a0 = levels[[targLevel-1]]$annot}
            CellSets = GetTrainingAndTestSets(target=target,annot=levels[[targLevel]]$annot,annotBelow=a0,expI=curExp)
            mod        = lm(annot~g1+g2+g3+g1*g2*g3,data=CellSets$trainDat)
            PredError = GetPredictionError(model=mod,data=CellSets,offTargetCells=unique(annotGroups$otherGroups))
        }else{ PredError = data.frame(FPos=1,FNeg=1,max_OffTarget_FPos=1) }
        
        GeneSetSummary = rbind(GeneSetSummary,data.frame(Group=target,t(curGenes),PredError,propChecks,i))
    }
    colnames(GeneSetSummary)[2:4] = c("g1","g2","g3")
    GeneSetSummary = GeneSetSummary[order(as.numeric(GeneSetSummary[,5])+as.numeric(GeneSetSummary[,6])),]
    
    selectedMarkers = data.frame(Group=target,Markers=as.vector(GeneSetSummary[1,2:4]))
    print(selectedMarkers)
    
    pdf(file=sprintf("MarkerPlots/MarkerPlot-Level%s-%s.pdf",targLevel,target),width=8,height=4)
    print(plotOptI(i=as.numeric(GeneSetSummary[1,"i"]),target=target,readCounts=expDat$readCounts,options=options,annot=levels[[targLevel]]$annot,annotBelow=annotGroups$annotBelow))
    dev.off()
    
    return(list(selectedMarkers=selectedMarkers,GeneSetSummary=GeneSetSummary))
}