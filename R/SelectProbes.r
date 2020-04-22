
## GET MEAN NUMBER OF READS IN A CELL TYPE
load("/Users/ns9/Google Drive/G2C.Transcriptome.Analysis/Data/Merged celltype_data (pyramidal as one group).rda")
exp_dat = read.csv("/Users/ns9/Desktop/Chat with Jens about In Situ Seq/expression_mRNA_17-Aug-2014_Reformed.csv")
exp_dat2 = exp_dat[!duplicated(exp_dat[,1]),]
rownames(exp_dat2)=exp_dat2[,1]
exp_dat2 = exp_dat2[,-1]
colN = gsub("\\..*","",colnames(exp_dat2))
count=0
for(cT in unique(colN)){
    count=count+1
    medianCT=apply(exp_dat2[,colN==cT],1,median)
    meanCT=apply(exp_dat2[,colN==cT],1,mean)
    sdCT=apply(exp_dat2[,colN==cT],1,sd)
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
allMedianDat=allMedianDat[,colnames(allMedianDat)!="X"]
allMeanDat=allMeanDat[,colnames(allMeanDat)!="X"]
allSdDat=allSdDat[,colnames(allSdDat)!="X"]
mean_cell_molecules=allMeanDat
median_cell_molecules=allMedianDat
sd_cell_molecules=allSdDat
save(mean_cell_molecules,file="mean_cell_molecules.rda")
save(median_cell_molecules,file="median_cell_molecules.rda")
save(sd_cell_molecules,file="sd_cell_molecules.rda")

## GET INTERNEURON SPECIFIC CELLTYPE PROPORTIONS
#  - Using Mean
mean_int_molecules = mean_cell_molecules[,grep("^Int",colnames(mean_cell_molecules))]
sum_int_molecules = apply(mean_int_molecules,1,sum)
proportion_int_molecules_MEAN = mean_int_molecules/sum_int_molecules
#  - Using Median
median_int_molecules = median_cell_molecules[,grep("^Int",colnames(median_cell_molecules))]
sum_int_molecules_Median = apply(median_int_molecules,1,sum)
proportion_int_molecules_MEDIAN = median_int_molecules/sum_int_molecules_Median

# CONSIDER AS MARKERS ALL GENES WITH AT LEAST 50% OF MEAN EXPRESSION IN ONE CELL TYPE
# - Using Mean
potential_markers_MEAN = proportion_int_molecules_MEAN[apply(proportion_int_molecules_MEAN,1,max)>0.5,]
potential_markers_MEAN = potential_markers_MEAN[!is.na(rownames(potential_markers_MEAN)),]
# - Using median
potential_markers_MEDIAN = proportion_int_molecules_MEDIAN[apply(proportion_int_molecules_MEDIAN,1,max)>0.5,]
potential_markers_MEDIAN = potential_markers_MEDIAN[!is.na(rownames(potential_markers_MEDIAN)),]

proportion_int_molecules_MEDIAN[is.nan(proportion_int_molecules_MEDIAN)]=0

count=0
for(gg in unique(rownames(potential_markers_MEAN),rownames(potential_markers_MEDIAN))){
    count=count+1
    maxProp_mean=max(proportion_int_molecules_MEAN[gg,])
    maxProp_median=max(proportion_int_molecules_MEDIAN[gg,])
    if(is.nan(maxProp_median)){maxProp_median=0}
    
    # Find the cell which most expresses that gene
    targetCell_mean=colnames(proportion_int_molecules_MEAN)[proportion_int_molecules_MEAN[gg,]==maxProp_mean]
    targetCell_median=colnames(proportion_int_molecules_MEDIAN)[proportion_int_molecules_MEDIAN[gg,]==maxProp_median][1]
    
    # Get mean/median expression in that cell type
    target_exp_mean = mean_cell_molecules[gg,targetCell_mean]
    target_exp_median = median_cell_molecules[gg,targetCell_median]
    
    off_target_interneurons = setdiff(colnames(mean_int_molecules)[mean_int_molecules[gg,]>1],targetCell_mean)
    off_target_interneurons_txt = paste0(off_target_interneurons,collapse=" & ")
    off_target_cells = setdiff(colnames(mean_cell_molecules)[mean_cell_molecules[gg,]>1],c(targetCell_mean,off_target_interneurons))
    off_target_cells_txt = paste0(off_target_cells,collapse=" & ")
    
    gene_dat = data.frame(gene=gg,maxProportion_mean=maxProp_mean,maxProportion_median=maxProp_median,targetCell=targetCell,expLevel_target_mean=target_exp_mean,expLevel_target_median=target_exp_median,offTargetInts=off_target_interneurons_txt,offTargetCells=off_target_cells_txt)
    if(count==1){
        allGenes=gene_dat
    }else{
        allGenes=rbind(allGenes,gene_dat)
    }
}
write.csv(allGenes,file="allGenes.csv")

# Where are Kenneth's targets amongst these?
kenneths = read.csv("kenneths_targets.csv")
allGenes$isKenneths = allGenes$gene %in% kenneths$GENE
#ggplot(allGenes)+geom_jitter(aes(x=maxProportion,y=expLevel_target,colour=isKenneths))+scale_y_log10()
pdf(file="MeanExp vs InterneuronSpecificity.pdf",width=16,height=14)
ggplot(allGenes)+geom_text(aes(x=maxProportion_mean,y=expLevel_target_mean,colour=isKenneths,label=gene))+scale_y_log10()
dev.off()

## Try finding the cells with highest proportion of expression in Int16
int16=data.frame(gene=rownames(mean_cell_molecules),mean_exp=mean_cell_molecules[,"Int16"],prop=proportion_int_molecules_MEAN[,"Int16"])
int16_sub = int16[int16$mean_exp>1 & int16$prop>0.1,]
ggplot(int16_sub)+geom_text(aes(x=prop,y=mean_exp,label=gene))+scale_y_log10()

## Find at most 10 genes, for each int type, with highest proportion (over 0.2) and mean_exp>10 in at least one interneuron cell type
geneList=c()
for(intT in 1:16){
    intName = sprintf("Int%s",intT)
    int16=data.frame(gene=rownames(mean_cell_molecules),mean_exp=mean_cell_molecules[,intName],prop=proportion_int_molecules_MEAN[,intName])
    int16_sub = int16[int16$mean_exp>10 & int16$prop>0.2,]
    int16_sub = int16_sub[order(int16_sub$prop*int16_sub$mean_exp,decreasing=TRUE),]
    geneList=c(geneList,as.character(int16_sub$gene[1:15]))
    #ggplot(int16_sub)+geom_text(aes(x=prop,y=mean_exp,label=gene))+scale_y_log10()
}

unique(geneList)[unique(geneList) %in% kenneths$GENE]
unique(geneList)[!unique(geneList) %in% kenneths$GENE]

kenneths$GENE[!kenneths$GENE %in% unique(geneList)]


ggplot(allGenes)+geom_text(aes(x=maxProportion_median,y=expLevel_target_median,colour=isKenneths,label=gene))+scale_y_log10()

ggplot(allGenes)+geom_text(aes(x=maxProportion_median,y=expLevel_target_median,colour=isKenneths,label=gene))+scale_y_log10()

Look at "Agtr2" to understand what is wrong

exp_dat3 = exp_dat2
colnames(exp_dat3)=
exp_dat_medians = 