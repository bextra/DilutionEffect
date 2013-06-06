# compareNormvsAdj.R
# K. Beck
# 5-31-2013

## This script requires functions in fileDifference.R
## This script requires functions in fileDifference.R

### Determine genes that are uniquely upregulated after the dilution adjustment ###
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


### Determine the regulation status of genes prior to the dilution adjustment ###
# Load table with all genes and their differential expression status in the unadjusted list
allBovineDE = read.table("~/Work/1_Milk/DilutionEffect/DE-Genes/Bovine-FINAL/Bt_norm_DEgenes16K.txt", header=TRUE, stringsAsFactors=FALSE)
allHumanDE = read.table("~/Work/1_Milk/DilutionEffect/DE-Genes/Human-FINAL/Hs_norm_DEgenes60K.txt", header = TRUE, stringsAsFactors=FALSE)

# Load gene lists that are the switchers
bovineSwitchers = readLines("~/Work/1_Milk/DilutionEffect/SwitcherAnalysis/Bovine/REFSEQ-upAdj_unique.txt")
humanSwitchers = readLines("~/Work/1_Milk/DilutionEffect/SwitcherAnalysis/Human/IDs-upAdj_unique.txt")

# Get the amount of differential expression from the unadjusted set for the genes that were switchers
BtSwitcherTable = allBovineDE[which(allBovineDE$id %in% bovineSwitchers),]
HsSwitcherTable = allHumanDE[which(allHumanDE$id %in% humanSwitchers),]

# Count how many of those have and adj log2 < 1 or not changed










## Confirm switchers are present in the adjusted set, but not in the unadjusted set ##
# Load the adjusted DE genes that have adj log2 > |1| and p < 0.05
# adjBovineDE = read.table("~/Work/1_Milk/DilutionEffect/DE-Genes/Bovine-FINAL/Bt_adj_DEgenes.txt", header=TRUE, stringsAsFactors=FALSE)
# adjHumanDE = read.table("~/Work/1_Milk/DilutionEffect/DE-Genes/Human-FINAL/Hs_adj_DEgenes.txt", header = TRUE, stringsAsFactors=FALSE)
# Make a table of the switcher genes in their original DE table
# BtSwitcherCheck = adjBovineDE[which(adjBovineDE$id %in% bovineSwitchers),]
# Observations: these genes have adj log2 value near 1 or p-values close to 0.05. Also, max is expreesion for Base Mean B is ~10K
# HsSwitcherCheck = adjHumanDE[which(adjHumanDE$id %in% humanSwitchers),]
# Observations: p-value close to 0.05, fold changes typically higher than in bovine though, but human overall has higher maxima