#plotGeneSet(geneSet=c("Sst","Pvalb","Gad1"),target=target,readCounts=expDat$readCounts,annot=levels[[targLevel]]$annot,annotBelow=annotGroups$annotBelow)

plotGeneSet <- function(geneSet,target,readCounts,annot,annotBelow){
    library(reshape2)
    if(!is.na(annotBelow)){
        targetBelow=unique(annotBelow[annot==target])
        cellsFromBelow = annotBelow==targetBelow
    }else{
        cellsFromBelow = rep(TRUE,length(annot))
    }    
    useGroups      = unique(annot[cellsFromBelow])  
    
    curGenes=geneSet
    curExp  = round(readCounts[curGenes,])
    plotDat = melt(data.frame(g1=t(curExp[1,]),g2=t(curExp[2,]),g3=t(curExp[3,]),annot=annot))
    plotDat  = plotDat[as.character(plotDat$annot) %in% useGroups,]

    upperLim=quantile(plotDat$value,0.999999)
    if(upperLim>200){
        print(ggplot(plotDat)+geom_boxplot(aes(x=annot,y=value,fill=variable))+ylab("Expression Level")+xlab("Cell Type")+scale_y_log10(limits = c(1,upperLim)))
    }else{
        print(ggplot(plotDat)+geom_boxplot(aes(x=annot,y=value,fill=variable))+ylab("Expression Level")+xlab("Cell Type")+ylim(0,upperLim))
    }
}