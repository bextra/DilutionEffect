# compareNormvsAdj.R
# K. Beck
# 5-31-2013

## This script requires functions (writeCompareOutput) from fileDifference.R
## This script requires functions (writeCompareOutput) from fileDifference.R

# # # # # # # # # # # # # # # # # 
#
# Not DE to Up Reg Switchers 
#
# # # # # # # # # # # # # # # # #

### Determine genes that are uniquely upregulated after the dilution adjustment ###
## BOVINE ##
bt.normUp = read.table("~/Work/1_Milk/DilutionEffect/DE-Genes/Bovine-FINAL/Bt_norm_upall.txt", 
                        header=TRUE, stringsAsFactors=FALSE)
bt.adjUp  = read.table("~/Work/1_Milk/DilutionEffect/DE-Genes/Bovine-FINAL/Bt_adj_upall.txt", 
                       header=TRUE, stringsAsFactors=FALSE)
bt.byRefSeq = compare(set1= bt.normUp, set2= bt.adjUp, compareColumn=1)
writeCompareOutput(compareList = bt.byRefSeq, set1.desc = "upNorm", set2.desc= "upAdj")

## HUMAN ##
hs.normUp = read.table("~/Work/1_Milk/DilutionEffect/DE-Genes/Human-FINAL/Hs_norm_upall.txt", 
                    header=TRUE, stringsAsFactors=FALSE)
hs.adjUp  = read.table("~/Work/1_Milk/DilutionEffect/DE-Genes/Human-FINAL/Hs_adj_upall.txt", 
                   header=TRUE, stringsAsFactors=FALSE)
hs.byRefSeq = compare(set1= hs.normUp, set2= hs.adjUp, compareColumn=1)
writeCompareOutput(compareList = hs.byRefSeq, set1.desc = "upNorm", set2.desc= "upAdj")


### Determine the regulation status of genes prior to the dilution adjustment ###
# Load table with all genes and their differential expression status in the unadjusted list
allBovineDE = read.table("~/Work/1_Milk/DilutionEffect/DE-Genes/Bovine-FINAL/Bt_norm_DEgenes16K.txt", header=TRUE, stringsAsFactors=FALSE)
allHumanDE = read.table("~/Work/1_Milk/DilutionEffect/DE-Genes/Human-FINAL/Hs_norm_DEgenes60K.txt", header = TRUE, stringsAsFactors=FALSE)

## Load gene lists that are the switchers
bovineSwitchers = readLines("~/Work/1_Milk/DilutionEffect/SwitcherAnalysis/Bovine/REFSEQ-upAdj_unique.txt")
humanSwitchers = readLines("~/Work/1_Milk/DilutionEffect/SwitcherAnalysis/Human/IDs-upAdj_unique.txt")

## Get the amount of differential expression from the unadjusted set for the genes that were switchers
BtSwitcherTable = allBovineDE[which(allBovineDE$id %in% bovineSwitchers),]
HsSwitcherTable = allHumanDE[which(allHumanDE$id %in% humanSwitchers),]


## Determine what the switchers were in the unadjusted data set
# Switchers are defined as not present in unadjusted up, but they *are* present in adjusted up
# Were they all not previously no change? or were some down regulated? or were some upregulated but excluded for p-val?
## BOVINE ##
length(which((BtSwitcherTable$padj < 0.05) & (BtSwitcherTable$adjLog2FoldChange > 0) & (BtSwitcherTable$adjLog2FoldChange < 1)))
    # 666 genes were statistically significant prior to dilution adjustment
    # all of those were also not differentially expressed (adj log2 < |1|) in the unadjusted set
length(which((BtSwitcherTable$padj > 0.05) & (BtSwitcherTable$adjLog2FoldChange > -1) & (BtSwitcherTable$adjLog2FoldChange < 1))) 
length(which((BtSwitcherTable$padj > 0.05) & (BtSwitcherTable$adjLog2FoldChange > 1)))
    # 305 were not statistically significant
        # broken down into 238 that were not stat sig and were also not differentially expressed
        # and 67 that were not statistically significant but were upregulated

## HUMAN ##
length(which((HsSwitcherTable$padj < 0.05) & (HsSwitcherTable$adjLog2FoldChange > 0) & (HsSwitcherTable$adjLog2FoldChange < 1)))
    # 45 genes were statistically significant prior to dilution adjustment
    # all of those were also not differentially expressed (adj log2 < |1|) in the unadjusted set
length(which((HsSwitcherTable$padj > 0.05) & (HsSwitcherTable$adjLog2FoldChange > -1) & (HsSwitcherTable$adjLog2FoldChange < 1))) 
length(which((HsSwitcherTable$padj > 0.05) & (HsSwitcherTable$adjLog2FoldChange > 1)))
    # 2110 were not statistically significant
        # broken down into 1716 that were not stat sig and were also not differentially expressed
        # and 394 that were not statistically significant but were upregulated


## Confirm switchers are present in the adjusted set, but not in the unadjusted set ##
# Load the adjusted DE genes that have adj log2 > |1| and p < 0.05
# adjBovineDE = read.table("~/Work/1_Milk/DilutionEffect/DE-Genes/Bovine-FINAL/Bt_adj_DEgenes.txt", header=TRUE, stringsAsFactors=FALSE)
# adjHumanDE = read.table("~/Work/1_Milk/DilutionEffect/DE-Genes/Human-FINAL/Hs_adj_DEgenes.txt", header = TRUE, stringsAsFactors=FALSE)
# Make a table of the switcher genes in their original DE table
# BtSwitcherCheck = adjBovineDE[which(adjBovineDE$id %in% bovineSwitchers),]
# Observations: these genes have adj log2 value near 1 or p-values close to 0.05. Also, max is expreesion for Base Mean B is ~10K
# HsSwitcherCheck = adjHumanDE[which(adjHumanDE$id %in% humanSwitchers),]
# Observations: p-value close to 0.05, fold changes typically higher than in bovine though, but human overall has higher maxima

# unadj not DE 1788
# adj not DE 3334
#outputFile="~/Work/1_Milk/DilutionEffect/DE-Genes/Bt_adj_notDEgenes.txt"


# # # # # # # # # # # # # # # # # 
#
# Down Reg to Not DE Switchers 
#
# # # # # # # # # # # # # # # # #

### Determine genes that go from downregulated to not-regulated after the dilution adjustment ###
## BOVINE ##
bt.normDown = read.table("~/Work/1_Milk/DilutionEffect/DE-Genes/Bovine-FINAL/Bt_norm_downall.txt", 
                       header=TRUE, stringsAsFactors=FALSE)
bt.adjDown  = read.table("~/Work/1_Milk/DilutionEffect/DE-Genes/Bovine-FINAL/Bt_adj_downall.txt", 
                       header=TRUE, stringsAsFactors=FALSE)
bt.byRefSeqDwn = compare(set1= bt.normDown, set2= bt.adjDown, compareColumn=1)
writeCompareOutput(compareList = bt.byRefSeqDwn, set1.desc = "downNorm", set2.desc= "downAdj")

## HUMAN ##
hs.normDown = read.table("~/Work/1_Milk/DilutionEffect/DE-Genes/Human-FINAL/Hs_norm_downall.txt", 
                       header=TRUE, stringsAsFactors=FALSE)
hs.adjDown  = read.table("~/Work/1_Milk/DilutionEffect/DE-Genes/Human-FINAL/Hs_adj_downall.txt", 
                       header=TRUE, stringsAsFactors=FALSE)
hs.byRefSeqDwn = compare(set1= hs.normDown, set2= hs.adjDown, compareColumn=1)
writeCompareOutput(compareList = hs.byRefSeqDwn, set1.desc = "downNorm", set2.desc= "downAdj")

### Determine the regulation status of genes prior to the dilution adjustment ###
# Load table with all genes and their differential expression status in the unadjusted list
allBovineDE = read.table("~/Work/1_Milk/DilutionEffect/DE-Genes/Bovine-FINAL/Bt_adj_DEgenes16K.txt", header=TRUE, stringsAsFactors=FALSE)
allHumanDE = read.table("~/Work/1_Milk/DilutionEffect/DE-Genes/Human-FINAL/Hs_adj_DEgenes60K.txt", header = TRUE, stringsAsFactors=FALSE)

## Load gene lists that are the switchers
## Prior to this line take just the IDs in a separate file
bovineSwitchers = readLines("~/Work/1_Milk/DilutionEffect/SwitcherAnalysis/Bovine/DowntoNoReg/REFSEQ-downNorm_unique.txt")
humanSwitchers = readLines("~/Work/1_Milk/DilutionEffect/SwitcherAnalysis/Human/DowntoNoReg/IDs-downNorm_unique.txt")

## Get the amount of differential expression from the adjusted set for the genes that were switchers
BtSwitcherTable = allBovineDE[which(allBovineDE$id %in% bovineSwitchers),]
HsSwitcherTable = allHumanDE[which(allHumanDE$id %in% humanSwitchers),]


## Determine what the down to no reg switchers were in the adjusted data set
# Switchers are defined as *present* in unadjusted down, but they are *not present* in adjusted down
# i.e. dilution adjusting them changed their regulation status
# Are they newly no change? or were some down regulated? but not stat-sig?
## BOVINE ##
length(which((BtSwitcherTable$padj < 0.05) & (BtSwitcherTable$adjLog2FoldChange > -1) & (BtSwitcherTable$adjLog2FoldChange < 1)))
    # 1729 genes were statistically significant after dilution adjustment
    # all of those were also not differentially expressed (adj log2 < |1|) in the adjusted set
length(which((BtSwitcherTable$padj > 0.05) & (BtSwitcherTable$adjLog2FoldChange < 0) & (BtSwitcherTable$adjLog2FoldChange > -1)))
    # 3 were not statistically significant but were still down regulated 
length(which(BtSwitcherTable$adjLog2FoldChange > 1))
    # none were upregulated after adjustment (regardless of statistical significance)


## HUMAN ##
length(which((HsSwitcherTable$padj < 0.05)))
    # 0 genes were statistically significant after dilution adjustment
length(which((HsSwitcherTable$padj > 0.05) & (HsSwitcherTable$adjLog2FoldChange >= -1) & (HsSwitcherTable$adjLog2FoldChange <= 1))) 
    # 3111 were not statistically significant
        # 2526 were not statistically significant and were not differentially expressed
length(which((HsSwitcherTable$padj > 0.05) & (HsSwitcherTable$adjLog2FoldChange < -1)))
        # 585 were not statistically significant and were downregulated
length(which((HsSwitcherTable$padj > 0.05) & (HsSwitcherTable$adjLog2FoldChange > 1)))
     # 0 were not stat sig and up regulated










