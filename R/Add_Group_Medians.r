dropThese = c("Hpcal1","Ptprm","Chat","Tacr1")

pdlist = read.csv("/Users/ns9/Google Drive/In Situ Probe Selection/Post_Discussion_Probelist_v2.csv",stringsAsFactors=FALSE)
pdlist = pdlist[as.character(pdlist$symbols)!="",]
pdlist = pdlist[!as.character(pdlist$symbols) %in% dropThese,]

# First, add median expression
needsMedian = as.character(pdlist$symbols[is.na(pdlist$MedianExpression_inTarget)])
for(sym in needsMedian){
    targ = as.character(pdlist[pdlist$symbols==sym,"target"])
    isFound = FALSE
    for(lev in 1:length(levels)){
        if(targ %in% colnames(levels[[lev]]$stats$groupMedians)){
            isFound = TRUE
            break()
        }
    }
    if(isFound){
        pdlist[pdlist$symbol==sym,"MedianExpression_inTarget"] = levels[[lev]]$stats$groupMedians[sym,targ]
    }else{
        targN = sprintf("%s\\.",targ)
        med = median(as.numeric(expDat$readCounts[sym,grep(targN,colnames(expDat$readCounts))]))
        pdlist[pdlist$symbol==sym,"MedianExpression_inTarget"] = med
    }
    if(pdlist[pdlist$symbol==sym,"MedianExpression_inTarget"]==0){
        pdlist[pdlist$symbol==sym,"MedianExpression_inTarget"]=0.1
    }
}
# Force Cort to have higher median expression
pdlist[pdlist$symbols=="Cort","MedianExpression_inTarget"]=30

# Set number of curProbes
load(file="/Users/ns9/Google Drive/In Situ Analysis/res2.rda")
colnames(res2)[1] = "symbols"
pdlist = pdlist[,c("symbols","MedianExpression_inTarget","Origin","target")]
pdlist = merge(pdlist,res2,by="symbols",all.x=TRUE)
pdlist$nProbes[is.na(pdlist$nProbes)]=0
mSets = pdlist[,c("symbols","nProbes","MedianExpression_inTarget","Origin","target","reads","totalReads_inCells","maxCellularReads","numCells_withAtLeast2reads")]
mSets$curProbes=mSets$nProbes

# Drop duplicates, keeping the one with lowest median expression
mSets = mSets[order(mSets$MedianExpression_inTarget),]
mSets = mSets[!duplicated(mSets$symbols),]

## Set all other genes to get two additional probes by default
#mSets[is.na(mSets[,"nProbes"]),"nProbes"] = 0

## Selective add probes where I want them
mSets$num_probes_to_order = 2
mSets[mSets$symbols=="Gad1","num_probes_to_order"]=5
mSets[mSets$symbols=="Slc6a1","num_probes_to_order"]=5
mSets[mSets$symbols=="Pvalb","num_probes_to_order"]=4
mSets[mSets$symbols=="Sst","num_probes_to_order"]=3
mSets[mSets$symbols=="Vip","num_probes_to_order"]=6
mSets[mSets$symbols=="Htr3a","num_probes_to_order"]=0
mSets[mSets$symbols=="Reln","num_probes_to_order"]=0
mSets$nProbes = mSets$curProbes + mSets$num_probes_to_order

## Add new probes in a greedy manner
numPlates = 2
mSets$num_probes_to_order = mSets$nProbes - mSets$curProbes
while(sum(mSets$num_probes_to_order)<(96*numPlates)){
    mSets$multiple = mSets$MedianExpression_inTarget*mSets$nProbes
    mSets=mSets[order(mSets$multiple),]
    # Restict so only five genes are ordered at most for any one gene
    addToGene=1
    for(i in 1:dim(combined)[1]){
        banned_genes = c("Nrn1","3110035E14Rik","Atp1b1","Npy")
        is_banned_gene = mSets$symbols[i] %in% banned_genes
        if(mSets$num_probes_to_order[i]<5 & !is_banned_gene){
            addToGene=i
            break()
        }
    }
    mSets$nProbes[addToGene] = mSets$nProbes[addToGene]+1
    mSets$num_probes_to_order = mSets$nProbes - mSets$curProbes
}

## Force Nrn1, 3110035E14Rik and Atp1b1 to have at most 2 probes

write.csv(mSets,file="UpdatedListOfProbes_forXiaoyan.csv")