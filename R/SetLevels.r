SetLevels <- function(cellNames){
    # Level0---Neurons vs Oligodendrocytes vs Astrocytes vs Other
    level0 = rep("Other",length(cellNames))
    level0[grep("(^Int)|(^CA1Pyr)|(^S1Pyr)|(^CA2Pyr)|(^ClauPyr)|(^SubPyr)",cellNames)]="Neurons"
    #level0[grep("^Oligo",cellNames)]="Oligodendrocytes"
    #level0[grep("^Astro",cellNames)]="Astrocytes"
    #level0[grep("(^Mgl1$)|(^Mgl2$)|(^Pvm1$)|(^Pvm2$)",cellNames)]="Microglia"
    
    # Level1---Pyramidal vs Interneurons
    level1 = rep("-",length(cellNames))
    level1[grep("(^Int)",cellNames)] = "Interneurons"
    level1[grep("(^CA1Pyr)|(^S1Pyr)|(^CA2Pyr)|(^ClauPyr)|(^SubPyr)",cellNames)] = "Pyramidal"
    
    # Level2--- 
    # -- Interneurons (Pvalb, Sst, 5HT3a (VIP+), 5HT3a (VIP-))
    # -- Pyramdal (S1PyrL6, S1PyrL6b, S1PyrL5, S1PyrL5a)
    level2 = rep("-",length(cellNames))
    level2[grep("(^Int1$)|(^Int2$)|(^Int4$)|(^Int3$)",cellNames)] = "Sst Pvalb"
    level2[intersect(grep("(^Sst Pvalb$)",level2,invert=TRUE),grep("Interneurons",level1))] = "5HT3a"  
    level2[level1=="Pyramidal"]="-"
    level2[grep("(^CA1Pyr)|(^CA2Pyr)",cellNames)]="Hippocampal Pyramidal"
    level2[grep("(^S1Pyr)|(^ClauPyr)|(^SubPyr)",cellNames)]="Cortical Pyramidal"
    
    # Level 3: specific cell types
    level3 = rep("-",length(cellNames))
    level3[grep("(^Int3$)",cellNames)]="Pvalb - Int3"
    level3[grep("(^Int1$)|(^Int2$)|(^Int4$)",cellNames)]="Sst - Int1 + Int2 + Int4"
    level3[grep("(^Int6$)|(^Int9$)|(^Int10$)",cellNames)]="Htr3a Vip+ - Int6 + Int9 + Int10"
    vipNeg1 = "(^Int5$)|(^Int7$)|(^Int8$)|(^Int11$)|(^Int12$)"
    vipNeg2 = "(^Int13$)|(^Int14$)|(^Int15$)|(^Int16$)"
    vipNegBoth = sprintf("%s|%s",vipNeg1,vipNeg2)
    level3[grep(vipNegBoth,cellNames)]="Htr3a Vip-"
    
    # Level 4: specific cell types
    level4 = rep("-",length(cellNames))
    level4[grep("(^Int1$)",cellNames)]="Sst - Int1"
    level4[grep("(^Int2$)|(^Int4$)",cellNames)]="Sst - Int2 + Int4"
    level4[grep("(^Int6$)",cellNames)]="Htr3a Vip+ - Int6"
    level4[grep("(^Int9$)",cellNames)]="Htr3a Vip+ - Int9"
    level4[grep("(^Int10$)",cellNames)]="Htr3a Vip+ - Int10"
    level4[grep(vipNeg1,cellNames)]="Htr3a Vip- Int5 + 7 + 8 + 11 + 12" 
    level4[grep(vipNeg2,cellNames)]="Htr3a Vip- Ivy+NGF Int13 + 14 + 15 + 16" 
    
    
    # Level 5: specific cell types
    level5 = rep("-",length(cellNames))
    level5[grep("(^Int13$)",cellNames)]="Htr3a Vip- Ivy NGF Int13"
    level5[grep("(^Int14$)",cellNames)]="Htr3a Vip- Ivy NGF Int14"
    level5[grep("(^Int15$)|(^Int16$)",cellNames)]="Htr3a Vip- Ivy NGF Int15 + 16"
    
    return(list(list(number=0,annot=level0),list(number=1,annot=level1),list(number=2,annot=level2),list(number=3,annot=level3),list(number=4,annot=level4),list(number=5,annot=level5)))
}

# 
# 
# SetLevels <- function(cellNames){
#     # Level0---Neurons vs Oligodendrocytes vs Astrocytes vs Other
#     level0 = rep("Other",length(cellNames))
#     level0[grep("(^Int)|(^CA1Pyr)|(^S1Pyr)|(^CA2Pyr)|(^ClauPyr)|(^SubPyr)",cellNames)]="Neurons"
#     level0[grep("^Oligo",cellNames)]="Oligodendrocytes"
#     level0[grep("^Astro",cellNames)]="Astrocytes"
#     
#     # Level1---Pyramidal vs Interneurons
#     level1 = rep("-",length(cellNames))
#     level1[grep("(^Int)",cellNames)] = "Interneurons"
#     level1[grep("(^CA1Pyr)|(^S1Pyr)|(^CA2Pyr)|(^ClauPyr)|(^SubPyr)",cellNames)] = "Pyramidal"
#     
#     # Level2--- 
#     # -- Interneurons (Pvalb, Sst, 5HT3a (VIP+), 5HT3a (VIP-))
#     # -- Pyramdal (S1PyrL6, S1PyrL6b, S1PyrL5, S1PyrL5a)
#     level2 = rep("-",length(cellNames))
#     level2[grep("(^Int1$)|(^Int2$)|(^Int4$)",cellNames)] = "Sst"
#     level2[grep("(^Int3$)",cellNames)] = "Pvalb"    
#     level2[grep("(^Int6$)|(^Int10$)|(^Int9$)",cellNames)] = "5HT3a (VIP+)"    
#     level2[intersect(grep("(^Sst$)|(^Pvalb$)|(^5HT3a)",level2,invert=TRUE),grep("Interneurons",level1))] = "5HT3a (VIP-)"  
#     # NOTE: This applies a conservative definition of VIP+... Int7 and Int8 are also kinda high
#     level2[level1=="Pyramidal"]="Pyramidal (other)"
#     level2[grep("(^CA1Pyr)",cellNames)]="CA1 Pyramidal"
#     level2[grep("(^S1PyrL6)",cellNames)]="S1 Pyr Layer 6"
#     level2[grep("(^S1PyrL5)",cellNames)]="S1 Pyr Layer 5"
#     level2[grep("(^S1PyrDL$)",cellNames)]="S1 Pyr Layer 6"
#     #level2[grep("(^ClauPyr$)",cellNames)]="ClauPyr"
#     #level2[grep("(^SubPyr$)",cellNames)]="SubPyr"
#     
#     # Level 3: specific cell types
#     level3 = rep("-",length(cellNames))
#     level3[level1=="Interneurons"] = cellNames[level1=="Interneurons"]
#     level3[grep("(^S1PyrL6b$)",cellNames)]="S1PyrL6b"
#     level3[grep("(^S1PyrL6$)",cellNames)]="S1PyrL6"
#     level3[grep("(^S1PyrL5$)",cellNames)]="S1PyrL5"
#     level3[grep("(^S1PyrL5a$)",cellNames)]="S1PyrL5a"
#     level3[grep("(^S1PyrL23$)",cellNames)]="S1PyrL23"
#     level3[grep("(^S1PyrL4$)",cellNames)]="S1PyrL4"
#     level3[grep("(^S1PyrDL$)",cellNames)]="S1PyrDL"
#     #level2[grep("(^ClauPyr$)",cellNames)]="ClauPyr"
#     #level2[grep("(^SubPyr$)",cellNames)]="SubPyr"    
#     
#     return(list(list(number=0,annot=level0),list(number=1,annot=level1),list(number=2,annot=level2),list(number=3,annot=level3)))
# }
# 
# 
# ##########################################################################################
# ##########################################################################################
# ##########################################################################################
# 
# SetLevels <- function(cellNames){
#     # Level0---Neurons vs Oligodendrocytes vs Astrocytes vs Other
#     level0 = rep("Other",length(cellNames))
#     level0[grep("(^Int)|(^CA1Pyr)|(^S1Pyr)|(^CA2Pyr)|(^ClauPyr)|(^SubPyr)",cellNames)]="Neurons"
#     level0[grep("^Oligo",cellNames)]="Oligodendrocytes"
#     level0[grep("^Astro",cellNames)]="Astrocytes"
#     
#     # Level1---Pyramidal vs Interneurons
#     level1 = rep("-",length(cellNames))
#     level1[grep("(^Int)",cellNames)] = "Interneurons"
#     level1[grep("(^CA1Pyr)|(^S1Pyr)|(^CA2Pyr)|(^ClauPyr)|(^SubPyr)",cellNames)] = "Pyramidal"
#     
#     # Level2--- 
#     # -- Interneurons (Pvalb, Sst, 5HT3a (VIP+), 5HT3a (VIP-))
#     # -- Pyramdal (S1PyrL6, S1PyrL6b, S1PyrL5, S1PyrL5a)
#     level2 = rep("-",length(cellNames))
#     level2[grep("(^Int1$)|(^Int2$)|(^Int4$)",cellNames)] = "Sst"
#     level2[grep("(^Int3$)",cellNames)] = "Pvalb"    
#     level2[grep("(^Int6$)|(^Int10$)|(^Int9$)",cellNames)] = "5HT3a (VIP+)"    
#     level2[intersect(grep("(^Sst$)|(^Pvalb$)|(^5HT3a)",level2,invert=TRUE),grep("Interneurons",level1))] = "5HT3a (VIP-)"  
#     # NOTE: This applies a conservative definition of VIP+... Int7 and Int8 are also kinda high
#     level2[level1=="Pyramidal"]="Pyramidal (other)"
#     level2[grep("(^CA1Pyr)",cellNames)]="CA1 Pyramidal"
#     level2[grep("(^S1PyrL6)",cellNames)]="S1 Pyr Layer 6"
#     level2[grep("(^S1PyrL5)",cellNames)]="S1 Pyr Layer 5"
#     level2[grep("(^S1PyrDL$)",cellNames)]="S1 Pyr Layer 6"
#     #level2[grep("(^ClauPyr$)",cellNames)]="ClauPyr"
#     #level2[grep("(^SubPyr$)",cellNames)]="SubPyr"
#     
#     # Level 3: specific cell types
#     level3 = rep("-",length(cellNames))
#     level3[level1=="Interneurons"] = cellNames[level1=="Interneurons"]
#     level3[grep("(^S1PyrL6b$)",cellNames)]="S1PyrL6b"
#     level3[grep("(^S1PyrL6$)",cellNames)]="S1PyrL6"
#     level3[grep("(^S1PyrL5$)",cellNames)]="S1PyrL5"
#     level3[grep("(^S1PyrL5a$)",cellNames)]="S1PyrL5a"
#     level3[grep("(^S1PyrL23$)",cellNames)]="S1PyrL23"
#     level3[grep("(^S1PyrL4$)",cellNames)]="S1PyrL4"
#     level2[grep("(^S1PyrDL$)",cellNames)]="S1PyrDL"
#     #level2[grep("(^ClauPyr$)",cellNames)]="ClauPyr"
#     #level2[grep("(^SubPyr$)",cellNames)]="SubPyr"    
#     
#     return(list(list(number=0,annot=level0),list(number=1,annot=level1),list(number=2,annot=level2),list(number=3,annot=level3)))
# }
# 
# 

SetLevelsOld <- function(cellNames){
    # Level0---Neurons vs Oligodendrocytes vs Astrocytes vs Other
    level0 = rep("Other",length(cellNames))
    level0[grep("(^Int)|(^CA1Pyr)|(^S1Pyr)|(^CA2Pyr)|(^ClauPyr)|(^SubPyr)",cellNames)]="Neurons"
    level0[grep("^Oligo",cellNames)]="Oligodendrocytes"
    level0[grep("^Astro",cellNames)]="Astrocytes"
    
    # Level1---Pyramidal vs Interneurons
    level1 = rep("-",length(cellNames))
    level1[grep("(^Int)",cellNames)] = "Interneurons"
    level1[grep("(^CA1Pyr)|(^S1Pyr)|(^CA2Pyr)|(^ClauPyr)|(^SubPyr)",cellNames)] = "Pyramidal"
    
    # Level2--- 
    # -- Interneurons (Pvalb, Sst, 5HT3a (VIP+), 5HT3a (VIP-))
    # -- Pyramdal (S1PyrL6, S1PyrL6b, S1PyrL5, S1PyrL5a)
    level2 = rep("-",length(cellNames))
    level2[grep("(^Int1$)|(^Int2$)|(^Int4$)",cellNames)] = "Sst"
    level2[grep("(^Int3$)",cellNames)] = "Pvalb"    
    level2[grep("(^Int6$)|(^Int10$)|(^Int9$)",cellNames)] = "5HT3a (VIP+)"    
    level2[intersect(grep("(^Sst$)|(^Pvalb$)|(^5HT3a)",level2,invert=TRUE),grep("Interneurons",level1))] = "5HT3a (VIP-)"  
    # NOTE: This applies a conservative definition of VIP+... Int7 and Int8 are also kinda high
    level2[level1=="Pyramidal"]="Pyramidal (other)"
    level2[grep("(^CA1Pyr)",cellNames)]="CA1 Pyramidal"
    level2[grep("(^S1PyrL6)",cellNames)]="S1 Pyr Layer 6"
    level2[grep("(^S1PyrL5)",cellNames)]="S1 Pyr Layer 5"
    level2[grep("(^S1PyrDL$)",cellNames)]="S1 Pyr Layer 6"
    #level2[grep("(^ClauPyr$)",cellNames)]="ClauPyr"
    #level2[grep("(^SubPyr$)",cellNames)]="SubPyr"
    
    # Level 3: specific cell types
    level3 = rep("-",length(cellNames))
    level3[level1=="Interneurons"] = cellNames[level1=="Interneurons"]
    level3[grep("(^S1PyrL6b$)",cellNames)]="S1PyrL6b"
    level3[grep("(^S1PyrL6$)",cellNames)]="S1PyrL6"
    level3[grep("(^S1PyrL5$)",cellNames)]="S1PyrL5"
    level3[grep("(^S1PyrL5a$)",cellNames)]="S1PyrL5a"
    level3[grep("(^S1PyrL23$)",cellNames)]="S1PyrL23"
    level3[grep("(^S1PyrL4$)",cellNames)]="S1PyrL4"
    level2[grep("(^S1PyrDL$)",cellNames)]="S1PyrDL"
    #level2[grep("(^ClauPyr$)",cellNames)]="ClauPyr"
    #level2[grep("(^SubPyr$)",cellNames)]="SubPyr"    
    
    return(list(list(number=0,annot=level0),list(number=1,annot=level1),list(number=2,annot=level2),list(number=3,annot=level3)))
}