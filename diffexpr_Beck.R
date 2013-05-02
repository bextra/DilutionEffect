#########################################################################
# diffexpr_Beck.R
#
# Adapted from diffexpr. R by D. Lemay 
# 5/1/2013
# R script to generate determine differentially expressed genes given raw
# count or FPKM data
##########################################################################

# Source loadCounts function in exprCountsLoad.R
library(DESeq) # Load required package



# Create a table of counts
# NOTE: this assumes all input files have been SORTED by geneID !!!!!!
# NOTE: this assumes all input files have been SORTED by geneID !!!!!!
pflist = list.files("~/Work/1_Milk/DilutionEffect/DE-Genes/", pattern = "^p.+[1-6].txt") # Get a list of the data files required
lflist = list.files("~/Work/1_Milk/DilutionEffect/DE-Genes/", pattern = "^l.+A[1-6].txt") # adjusted
setwd("~/Work/1_Milk/DilutionEffect/DE-Genes/")
nonlac = loadCounts(pflist, n = 6, computeMean=FALSE) # n = number of replicates
lac    = loadCounts(lflist, n = 6, computeMean=FALSE)
countTable = cbind(nonlac[,-1], lac[,-1])
countTable = sapply(countTable[, c(1:ncol(countTable))], as.numeric)
row.names(countTable) = lac[,1]

# Describe the conditions (i.e. treatment vs non-treatment) of each sample
cond = factor(c(rep("prepuberty", 6), rep("lactation", 6))) # number indicated is the number of replicates in each treatment group

## DESeq ##

# Create a CountDataSet which is the central data structure in DESeq
cds = newCountDataSet(countTable, cond)

# Estimate effective library size
cds = estimateSizeFactors(cds)

# Estimate dispersions (should be zero since we generated the data using a Poisson distribution
# cds = estimateVarianceFunctions(cds)

## if using a newer version of DESeq, comment out previous line, uncomment this:
cds = estimateDispersions(cds)

# test for differential expression
res = nbinomTest(cds, "prepuberty", "lactation")

# if the Mean of one of the two comparisons is zero, the log2FoldChange will be INF/-INF so calculating an adjusted log2 fold change
res$adjLog2FoldChange = log2 ((res$baseMeanB + 1)/(res$baseMeanA + 1))

# get transcripts with an adjusted log2FoldChange of greater than 1 or 
# less than -1 and with an adjusted p-value < 0.05
res.filtered = res[which((abs(res$adjLog2FoldChange) > 1) & (res$padj < 0.05 | is.na(res$padj))),]

# order transcripts by highest fold change
res.ordered <- res.filtered[order(res.filtered$adjLog2FoldChange, decreasing=TRUE),]

# save results to a file
write.table(res.ordered, file="./DEgenesAllTesting.txt", quote=FALSE, row.names = FALSE, sep="\t")



#########################################################################

# # Get the raw count data
# p_rep1 <- read.table(file="pseuct_pcountsRep1.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
# p_rep2 <- read.table(file="pseuct_pcountsRep2.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
# p_rep3 <- read.table(file="pseuct_pcountsRep3.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
# p_rep4 <- read.table(file="pseuct_pcountsRep4.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
# p_rep5 <- read.table(file="pseuct_pcountsRep5.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
# p_rep6 <- read.table(file="pseuct_pcountsRep6.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
# #p_rep7 <- read.table(file="pseuct_pcountsRep7.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
# #p_rep8 <- read.table(file="pseuct_pcountsRep8.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
# #p_rep9 <- read.table(file="pseuct_pcountsRep9.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
# #p_rep10 <- read.table(file="pseuct_pcountsRep10.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
# l_rep1 <- read.table(file="lcountsA1.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
# l_rep2 <- read.table(file="lcountsA2.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
# l_rep3 <- read.table(file="lcountsA3.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
# l_rep4 <- read.table(file="lcountsA4.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
# l_rep5 <- read.table(file="lcountsA5.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
# l_rep6 <- read.table(file="lcountsA6.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
# #l_rep7 <- read.table(file="lcountsN7.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
# #l_rep8 <- read.table(file="lcountsN8.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
# #l_rep9 <- read.table(file="lcountsN9.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
# #l_rep10 <- read.table(file="lcountsN10.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
# 
# #Use these lines for 10 reps (below)
# #countsTable <- as.data.frame(cbind(p_rep1$Rep1, p_rep2$Rep2, p_rep3$Rep3, p_rep4$Rep4, p_rep5$Rep5, p_rep6$Rep6, p_rep7$Rep7, p_rep8$Rep8, p_rep9$Rep9, p_rep10$Rep10, l_rep1$Rep1, l_rep2$Rep2, l_rep3$Rep3, l_rep4$Rep4, l_rep5$Rep5, l_rep6$Rep6, l_rep7$Rep7, l_rep8$Rep8, l_rep9$Rep9, l_rep10$Rep10))
# #colnames(countsTable) <- c("p_rep1", "p_rep2", "p_rep3","p_rep4", "p_rep5", "p_rep6", "p_rep7", "p_rep8", "p_rep9", "p_rep10", "l_rep1", "l_rep2", "l_rep3","l_rep4", "l_rep5", "l_rep6", "l_rep7", "l_rep8", "l_rep9", "l_rep10")
# 
# #Use these lines for 6 reps (below)
# countsTable <- as.data.frame(cbind(p_rep1$Rep1, p_rep2$Rep2, p_rep3$Rep3, p_rep4$Rep4, p_rep5$Rep5, p_rep6$Rep6, l_rep1$Rep1, l_rep2$Rep2, l_rep3$Rep3, l_rep4$Rep4, l_rep5$Rep5, l_rep6$Rep6))
# colnames(countsTable) <- c("p_rep1", "p_rep2", "p_rep3","p_rep4", "p_rep5", "p_rep6", "l_rep1", "l_rep2", "l_rep3","l_rep4", "l_rep5", "l_rep6")
# 
# # Use this line for 10 reps
# #cond <- factor(c("prepuberty", "prepuberty", "prepuberty", "prepuberty", "prepuberty", "prepuberty", "prepuberty", "prepuberty", "prepuberty", "prepuberty", "lactation", "lactation", "lactation", "lactation", "lactation", "lactation", "lactation", "lactation", "lactation", "lactation"))
# # Use this line for 6 reps
# cond <- factor(c("prepuberty", "prepuberty", "prepuberty", "prepuberty", "prepuberty", "prepuberty", "lactation", "lactation", "lactation", "lactation", "lactation", "lactation"))
# 
# # modified to test for NA from p-values
# #res.filtered <- res[which(abs(res$adjLog2FoldChange) > 1),]
# 
# # modified to capture not statistically significant or no fold change
# #res.nostatsig <- res[which(abs(res$adjLog2FoldChange) > 1 & res$padj == 5e-02),]
# 
# # save results to a file that are not DE or statistically significant
# #write.table(res.nostatsig, file="./DE_genes_Reps6_EqStatSig.txt", quote=FALSE, row.names = FALSE, sep="\t")
