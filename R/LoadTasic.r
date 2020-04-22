LoadTasic <- function(){
    path = "/Users/ns9/Datasets that are too large to store elsewhere/Tasic et al Cell Type dataset/"
    gene_rpkm = read.csv(sprintf("%sgenes_rpkm.csv",path))
    rownames(gene_rpkm) = gene_rpkm$X
    gene_rpkm = gene_rpkm[,-1]
    cell_class  = read.csv(sprintf("%scell_classification.csv",path))
    rownames(cell_class) = cell_class$X
    #cell_class = cell_class[,-1]    
    clusters    = read.csv(sprintf("%scluster_metadata.csv",path))
    rownames(clusters) = clusters$cluster_id
    cell_meta   = read.csv(sprintf("%scell_metadata.csv",path))
    rownames(cell_meta) = cell_meta$short_name
    cell_class$long_name = cell_meta$long_name[cell_class$X]
    rownames(cell_class) = cell_class$long_name
    return(list(rpkm=gene_rpkm,cells=cell_class,clusters=clusters,cell_meta=cell_meta))
}

PlotTasicGene <- function(data,gSym){
    library(ggplot2)
    exp = tasic$rpkm[gSym,]
    cell_cluster = as.character(data$cells[names(exp),]$primary)
    cell_type    = as.character(data$clusters[cell_cluster,]$vignette_label)
    tmp = data.frame(exp=as.numeric(exp),cell_type=cell_type)
    ggplot(tmp) + geom_boxplot(aes(x=cell_type,y=exp))+scale_y_log10(limits=c(50,max(exp)))+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
}

tasic = LoadTasic()