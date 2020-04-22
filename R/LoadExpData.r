LoadExpData <- function(){
    #load("/Users/ns9/Google Drive/G2C.Transcriptome.Analysis/Data/Merged celltype_data (pyramidal as one group).rda")
    exp_dat = read.csv("Data/expression_mRNA_17-Aug-2014_Reformed.csv")
    exp_dat2 = exp_dat[!duplicated(exp_dat[,1]),]
    rownames(exp_dat2)=exp_dat2[,1]
    exp_dat2 = exp_dat2[,-1]
    colN = gsub("\\..*","",colnames(exp_dat2))
    
    # Drop cells annotated as "X"
    exp_dat3 = exp_dat2[,colN!="X"]
    colN = colN[colN!="X"]
    
    allDat = list(cellNames=colN,readCounts=exp_dat3)
    return(allDat)
}