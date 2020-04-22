# INPUTS:
# target: name of the target group, i.e. "NEURONS"
# annot: annotation at the level of group
# annotBelow:  annotations for the level below target

GetTrainingAndTestSets <- function(target,annot,annotBelow,expI){
    if(!is.na(annotBelow)){
        targetBelow=unique(annotBelow[annot==target])
        cellsFromBelow = annotBelow==targetBelow
    }else{
        cellsFromBelow = rep(TRUE,length(annot))
    }
    useGroups      = unique(annot[cellsFromBelow])  
    # If this level of definition does not provide any further information... then ignore it
    if(length(useGroups)==1){return(NA)}    
    
    posIdx = which(annot==target & cellsFromBelow)
    negIdx = which(annot!=target & cellsFromBelow)
    posTrainIdx = sample(posIdx,floor(length(posIdx)/2))
    posTestIdx  = setdiff(posIdx,posTrainIdx)
    negTrainIdx = sample(negIdx,floor(length(negIdx)/2))
    negTestIdx  = setdiff(negIdx,negTrainIdx)      
    TrainIdx    = c(posTrainIdx,negTrainIdx)
    TrainAnnot  = levels[[targLevel]]$annot[TrainIdx]==target
    TrainGroups = levels[[targLevel]]$annot[TrainIdx]
    TestIdx    = c(posTestIdx,negTestIdx)
    TestAnnot  = levels[[targLevel]]$annot[TestIdx]==target
    TestGroups = levels[[targLevel]]$annot[TestIdx]
    
    TrainExp = t(expI[,TrainIdx])
    trainDat = data.frame(g1=TrainExp[,1],g2=TrainExp[,2],g3=TrainExp[,3],annot=TrainAnnot)
    colnames(trainDat) = c("g1","g2","g3","annot")
    
    TestExp = t(expI[,TestIdx])
    testDat = data.frame(g1=TestExp[,1],g2=TestExp[,2],g3=TestExp[,3],annot=TestAnnot)
    colnames(testDat) = c("g1","g2","g3","annot")
    
    return(list(TrainIdx=TrainIdx,TestIdx=TestIdx,TrainAnnot=TrainAnnot,TestAnnot=TestAnnot,trainDat=trainDat,testDat=testDat,TestGroups=TestGroups,TrainGroups=TrainGroups))
}