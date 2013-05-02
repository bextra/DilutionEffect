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
milkTable = cbind(nonlac[,-1], lac[,-1])
milkTable = sapply(milkTable[, c(1:ncol(milkTable))], as.numeric)
row.names(milkTable) = lac[,1]

computeDEgenes(milkTable, outputFile="./atyourmamasfile.txt")

# Data structure should be columns of one sample condition bound to columns of another sample condition

## DESeq ##
computeDEgenes =
    function (countTable, n = 6, outputFile = "./DEgenes.txt") {
        
        # Describe the conditions (i.e. treatment vs non-treatment) of each sample
        cond = factor(c(rep("baseline", n), rep("lactation", n))) # number indicated is the number of replicates in each treatment group
        
        # Create a CountDataSet which is the central data structure in DESeq
        cds = newCountDataSet(countTable, cond)
        cat("Count data set made.\n")
        
        # Estimate effective library size
        cds = estimateSizeFactors(cds)
        
        # Estimate dispersions (should be zero since we generated the data using a Poisson distribution
        # cds = estimateVarianceFunctions(cds)
        
        ## if using a newer version of DESeq, comment out previous line, uncomment this:
        cds = estimateDispersions(cds)
        
        # test for differential expression
        cat("Testing for differential expression.\n")
        res = nbinomTest(cds, "baseline", "lactation")
        
        # if the Mean of one of the two comparisons is zero, the log2FoldChange will be INF/-INF so calculating an adjusted log2 fold change
        res$adjLog2FoldChange = log2 ((res$baseMeanB + 1)/(res$baseMeanA + 1))
        
        # get transcripts with an adjusted log2FoldChange of greater than 1 or 
        # less than -1 and with an adjusted p-value < 0.05
        res.filtered = res[which((abs(res$adjLog2FoldChange) > 1) & (res$padj < 0.05 | is.na(res$padj))),]
        
        # order transcripts by highest fold change
        res.ordered <- res.filtered[order(res.filtered$adjLog2FoldChange, decreasing=TRUE),]
        
        # save results to a file
        cat("Making output file.\n")
        write.table(res.ordered, file=outputFile, quote=FALSE, row.names = FALSE, sep="\t")
        
    }

