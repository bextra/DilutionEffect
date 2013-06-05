# compareNormvsAdj.R
# K. Beck
# 5-31-2013


## BOVINE ##
bt.normUp = read.table("~/Work/1_Milk/DilutionEffect/DE-Genes/Bovine-FINAL/Bt_norm_upall.txt", 
                        header=TRUE, stringsAsFactors=FALSE)
bt.adjUp = read.table("~/Work/1_Milk/DilutionEffect/DE-Genes/Bovine-FINAL/Bt_adj_upall.txt", 
                       header=TRUE, stringsAsFactors=FALSE)

bt.byRefSeq = compare(set1= bt.normUp, set2= bt.adjUp, compareColumn=1)

writeParts(compareList = bt.byRefSeq, set1.desc = "upNorm", set2.desc= "upAdj")



## HUMAN ##
hs.normUp = read.table("~/Work/1_Milk/DilutionEffect/DE-Genes/Human-FINAL/Hs_norm_upall.txt", 
                    header=TRUE, stringsAsFactors=FALSE)
hs.adjUp = read.table("~/Work/1_Milk/DilutionEffect/DE-Genes/Human-FINAL/Hs_adj_upall.txt", 
                   header=TRUE, stringsAsFactors=FALSE)

hs.byRefSeq = compare(set1= hs.normUp, set2= hs.adjUp, compareColumn=1)

writeParts(compareList = hs.byRefSeq, set1.desc = "upNorm", set2.desc= "upAdj")