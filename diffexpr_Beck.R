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


# # # # # # # # # # # # # # # # # 
#
# Functions 
#
# # # # # # # # # # # # # # # # #
combineSamples = 
    function (data1, data2) {
        # assumes gene IDs listed in first column followed by corresponding gene expression data, one replicate per column  
        # assumes both dataframes are sorted by gene ID
        myTable = cbind(data1[,-1], data2[,-1])
        myTable = sapply(myTable[, c(1:ncol(myTable))], as.numeric)
        row.names(myTable) = data1[,1]        
        return(myTable)
    }


## DESeq ##
computeDEgenes =
    function (countTable, repsCondition1 = 6, repsCondition2 = 6, outputFile = "./DEgenes.txt") {
        
        # Describe the conditions (i.e. treatment vs non-treatment) of each sample
        cond = factor(c(rep("condition1", repsCondition1), rep("condition2", repsCondition2))) # number indicated is the number of replicates in each treatment group
        
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
        res = nbinomTest(cds, "condition1", "condition2")
        
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

# # # # # # # # # # # # # # # # # # # # # # # # 
#
# Compute differentially expressed genes 
#
# # # # # # # # # # # # # # # # # # # # # # # # 

### Source loadCounts function in exprCountsLoad.R
# NOTE: this assumes all input files have been SORTED by geneID !!!!!!
# List baseline/control condition *first* in the counts table and replicate number


# # # BOVINE # # #
# Pre-puberty to lactation comparison of count data from Harhay et al.
pflist = list.files("~/Work/1_Milk/DilutionEffect/DE-Genes/", pattern = "^p.+[1-6].txt") # list data files required - prepuberty
lflist = list.files("~/Work/1_Milk/DilutionEffect/DE-Genes/", pattern = "^l.+A[1-6].txt") # ditto for adjusted set  - lactation
setwd("~/Work/1_Milk/DilutionEffect/DE-Genes/") # change working directy to desired location for output
nonlac = loadCounts(pflist, n = 6, computeMean=FALSE) # n = number of replicates
lac    = loadCounts(lflist, n = 6, computeMean=FALSE)
bovineAdjusted = combineSamples(nonlac, lac)

computeDEgenes(bovineAdjusted, outputFile="./atyourmamasfile2.txt")


