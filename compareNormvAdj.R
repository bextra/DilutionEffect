# compareNormvsAdj.R
# K. Beck
# 5-31-2013



# Workflow for above functions
# BeckGeneU = read.table("~/Work/1_Milk/ProteomicsWork/Human-BeckvGao/FullDatabase/BeckAllGeneNames3.txt", 
#                        header=TRUE, stringsAsFactors=FALSE)
# GaoGeneU = read.table("~/Work/1_Milk/ProteomicsWork/ProteinLists-Other's/GaoAllGeneList.txt", 
#                       header=TRUE, stringsAsFactors=FALSE)
# byUniProtAccession = compare(set1= BeckGeneU, set2= GaoGeneU, compareColumn=1)
# writeParts(compareList = byUniProtAccession, set1.desc = "Beck", set2.desc= "Gao")


## BOVINE ##
normUp = read.table("~/Work/1_Milk/DilutionEffect/DE-Genes/Bovine-FINAL/Bt_norm_upall.txt", 
                        header=TRUE, stringsAsFactors=FALSE)
adjUp = read.table("~/Work/1_Milk/DilutionEffect/DE-Genes/Bovine-FINAL/Bt_adj_upall.txt", 
                       header=TRUE, stringsAsFactors=FALSE)

byRefSeq = compare(set1= normUp, set2= adjUp, compareColumn=1)

writeParts(compareList = byRefSeq, set1.desc = "upNorm", set2.desc= "upAdj")