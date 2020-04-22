# EXAMPLES:
# - plotOptI(i=1,target=target,readCounts=expDat$readCounts,options=options,annot=levels[[targLevel]]$annot,annotBelow=levels[[targLevel-1]]$annot)

plotOptI <- function(i,target,readCounts,options,annot,annotBelow){
    if(!is.na(annotBelow)){
        targetBelow=unique(annotBelow[annot==target])
        cellsFromBelow = annotBelow==targetBelow
    }else{
        cellsFromBelow = rep(TRUE,length(annot))
    }    
    useGroups      = unique(annot[cellsFromBelow])  
    
    curGenes=options[,i]
    curExp  = round(readCounts[curGenes,])
    plotDat = melt(data.frame(g1=t(curExp[1,]),g2=t(curExp[2,]),g3=t(curExp[3,]),annot=annot))
    plotDat  = plotDat[plotDat$annot %in% useGroups,]
    upperLim=quantile(plotDat$value,0.95)
    ggplot(plotDat)+geom_boxplot(aes(x=annot,y=value,fill=variable))+ylab("Expression Level")+xlab("Cell Type")+ylim(0,upperLim)
}