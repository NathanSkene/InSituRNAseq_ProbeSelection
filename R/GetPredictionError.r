# EXAMPLES:
# - PredError = GetPredictionError(model=mod,data=CellSets,offTargetCells=unique(annotGroups$otherGroups))

GetPredictionError <- function(model,data,offTargetCells){
    # Get false positives+negatives
    #output_training  = predict(data$trainDat,object=model)
    output_test      = predict(data$testDat,object=model)
    FPos = sum(output_test[data$TestAnnot==FALSE]>0.5)/sum(data$TestAnnot==FALSE)
    FNeg = sum(output_test[data$TestAnnot==TRUE]<0.5)/sum(data$TestAnnot==TRUE)
    
    # Find the FPos for each cell group
    G_FPos = rep(1,length(offTargetCells))
    for(jjj in 1:length(offTargetCells)){
        G_FPos[jjj] = sum(output_test[data$TestGroups==offTargetCells[jjj]]>0.5)/sum(data$TestGroups==offTargetCells[jjj])
    }
    max_Group_FPos = max(G_FPos)
    
    return(data.frame(FPos=FPos,FNeg=FNeg,max_OffTarget_FPos=max_Group_FPos))
}