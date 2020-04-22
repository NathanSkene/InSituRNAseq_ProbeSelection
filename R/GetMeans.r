GetMeans <- function(expDat,annot){
    count=0
    for(cT in unique(annot)){
        count=count+1
        
        reads= expDat$readCounts[,annot==cT]
        
        medianCT=apply(expDat$readCounts[,annot==cT],1,median)
        meanCT=apply(expDat$readCounts[,annot==cT],1,mean,trim=0.33)
        #meanCT=apply(expDat$readCounts[,annot==cT],1,mean,trim=0)
        sdCT=apply(expDat$readCounts[,annot==cT],1,sd)
        medianDat = matrix(medianCT)
        meanDat = matrix(meanCT)
        sdDat = matrix(sdCT)
        colnames(meanDat)=colnames(sdDat)=colnames(medianDat)=cT
        rownames(meanDat)=rownames(sdDat)=rownames(medianDat)=names(medianCT)
        
        if(count==1){
            allMedianDat=medianDat
            allMeanDat=meanDat
            allSdDat=sdDat
        }else{
            allMedianDat=cbind(allMedianDat,medianDat)
            allMeanDat=cbind(allMeanDat,meanDat)
            allSdDat=cbind(allSdDat,sdDat)
        }
    }
    mean_cell_molecules=allMeanDat
    median_cell_molecules=allMedianDat
    sd_cell_molecules=allSdDat
    effectSize=mean_cell_molecules/sd_cell_molecules
    #proportions = mean_cell_molecules/apply(mean_cell_molecules,1,sum)
    proportions = median_cell_molecules/(apply(median_cell_molecules,1,sum)+0.000001)
    
    return(list(groupMeans=mean_cell_molecules,groupMedians=median_cell_molecules,groupSd=sd_cell_molecules,groupEffectSize=effectSize,proportions=proportions))
}