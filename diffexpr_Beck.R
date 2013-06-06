#########################################################################
# diffexpr_Beck.R
#
# Adapted from diffexpr.R by D. Lemay 
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
        res.filtered = res[which((abs(res$adjLog2FoldChange) > 1) & (res$padj < 0.05 | is.na(res$padj))),] # default settings
        # res.filtered = res # use this line to get all genes and their DE values
        
        # order transcripts by highest fold change
        res.ordered = res.filtered[order(res.filtered$adjLog2FoldChange, decreasing=TRUE),]
        
        cat(nrow(res.ordered), "differentially expressed genes were found.\n")
        
        # save results to a file
        write.table(res.ordered, file=outputFile, quote=FALSE, row.names = FALSE, sep="\t")
        
    }

# # # # # # # # # # # # # # # # # # # # # # # # 
#
# Compute differentially expressed genes 
#
# # # # # # # # # # # # # # # # # # # # # # # # 

### Source loadCounts function in exprCountsLoad.R
# NOTE: this assumes all input files have been SORTED by geneID !!!!!!
# List baseline/control condition *first* in the counts table
# Therfore, adj log2 = 1.5 indicates an increase of 1.5 from condition1 to condition2

# # # BOVINE PRE-PUBERTY TO UN-ADJUSTED LACTATION # # #
# Pre-puberty to lactation comparison of count data from Harhay et al.
setwd("~/Work/1_Milk/DilutionEffect/Dilution_Outputs/DilutionAdj-BovineFINAL/") # change working directory to input file location
norm_Bt_pflist = list.files("~/Work/1_Milk/DilutionEffect/Dilution_Outputs/DilutionAdj-BovineFINAL/", pattern = "^psct_p.+[1-6].txt") # list data files required - prepuberty
norm_Bt_lflist = list.files("~/Work/1_Milk/DilutionEffect/Dilution_Outputs/DilutionAdj-BovineFINAL/", pattern = "^norm.+[1-6].txt") # unadjusted lactation
nonlac = loadCounts(norm_Bt_pflist, nreps= 6, computeMean=FALSE) # function in exprCountsLoad.R
lac    = loadCounts(norm_Bt_lflist, nreps= 6, computeMean=FALSE) # function in exprCountsLoad.R
bovineUnadjusted = combineSamples(nonlac, lac)
computeDEgenes(bovineUnadjusted, outputFile="~/Work/1_Milk/DilutionEffect/DE-Genes/Bt_norm_DEgenes.txt")



# # # BOVINE PRE-PUBERTY TO ADJUSTED LACTATION # # #
# Pre-puberty to lactation comparison of count data from Harhay et al.
norm_Bt_pflist = list.files("~/Work/1_Milk/DilutionEffect/Dilution_Outputs/DilutionAdj-BovineFINAL/", pattern = "^psct_p.+[1-6].txt") # list data files required - prepuberty
adj_Bt_lflist  = list.files("~/Work/1_Milk/DilutionEffect/Dilution_Outputs/DilutionAdj-BovineFINAL/", pattern = "^adj.+[1-6].txt") # dilution adjusted lactation
nonlac = loadCounts(norm_Bt_pflist, nreps = 6, computeMean=FALSE)
lac    = loadCounts(adj_Bt_lflist,  nreps = 6, computeMean=FALSE)
bovineAdjusted = combineSamples(nonlac, lac)
computeDEgenes(bovineAdjusted, outputFile="~/Work/1_Milk/DilutionEffect/DE-Genes/Bt_adj_DEgenes.txt")



# # # HUMAN COLOSTRUM TO UN-ADJUSTED MATURE LACTATION # # #
# Colostrum to mature lactation comparison of FPKM data from Nommsen-Rivers transcriptome manuscript
## NOTE: Prior to loading files replace header line from bottom of sort back to top of file where header should be
## NOTE: Prior to loading files replace header line from bottom of sort back to top of file where header should be
setwd("~/Work/1_Milk/DilutionEffect/Dilution_Outputs/DilutionAdj-HumanFINAL/") # change working directory to input file location
norm_Hs_cflist = list.files("~/Work/1_Milk/DilutionEffect/Dilution_Outputs/DilutionAdj-HumanFINAL/", pattern = "^psct_Col.+[1-2].txt")   # list data files required - colostrum
norm_Hs_mflist = list.files("~/Work/1_Milk/DilutionEffect/Dilution_Outputs/DilutionAdj-HumanFINAL/", pattern = "^norm.+[1-6].txt")  # ditto for unadjusted set - mature
colostrum = loadCounts(norm_Hs_cflist, nreps = 2, computeMean=FALSE)
mature    = loadCounts(norm_Hs_mflist, nreps = 6, computeMean=FALSE)
humanUnadjusted = combineSamples(colostrum, mature)
computeDEgenes(humanUnadjusted, repsCondition1= 2, repsCondition2= 6, outputFile="~/Work/1_Milk/DilutionEffect/DE-Genes/Hs_norm_DEgenes.txt")



# # # HUMAN COLOSTRUM TO ADJUSTED MATURE LACTATION # # #
# Colostrum to mature lactation comparison of FPKM data from Nommsen-Rivers transcriptome manuscript
## NOTE: Prior to loading files replace header line from bottom of sort back to top of file where header should be
## NOTE: Prior to loading files replace header line from bottom of sort back to top of file where header should be
norm_Hs_cflist = list.files("~/Work/1_Milk/DilutionEffect/Dilution_Outputs/DilutionAdj-HumanFINAL/", pattern = "^psct_Col.+[1-2].txt") # list data files required - colostrum
adj_Hs_mflist  = list.files("~/Work/1_Milk/DilutionEffect/Dilution_Outputs/DilutionAdj-HumanFINAL/", pattern = "^adj.+[1-6].txt") # ditto for adjusted set   - mature
colostrum = loadCounts(norm_Hs_cflist, nreps = 2, computeMean=FALSE)
mature    = loadCounts(adj_Hs_mflist,  nreps = 6, computeMean=FALSE)
humanAdjusted = combineSamples(colostrum, mature)
computeDEgenes(humanAdjusted, repsCondition1= 2, repsCondition2= 6, outputFile="~/Work/1_Milk/DilutionEffect/DE-Genes/Hs_adj_DEgenes.txt")


# # # GET ALL DIFFERENTIAL EXPR VALUES FOR ALL GENES # # #
# # # BOVINE PRE-PUBERTY TO UN-ADJUSTED LACTATION # # #
# Pre-puberty to lactation comparison of count data from Harhay et al.
setwd("~/Work/1_Milk/DilutionEffect/Dilution_Outputs/DilutionAdj-BovineFINAL/") # change working directory to input file location
norm_Bt_pflist = list.files("~/Work/1_Milk/DilutionEffect/Dilution_Outputs/DilutionAdj-BovineFINAL/", pattern = "^psct_p.+[1-6].txt") # list data files required - prepuberty
norm_Bt_lflist = list.files("~/Work/1_Milk/DilutionEffect/Dilution_Outputs/DilutionAdj-BovineFINAL/", pattern = "^norm.+[1-6].txt") # unadjusted lactation
nonlac = loadCounts(norm_Bt_pflist, nreps= 6, computeMean=FALSE) # function in exprCountsLoad.R
lac    = loadCounts(norm_Bt_lflist, nreps= 6, computeMean=FALSE) # function in exprCountsLoad.R
bovineUnadjusted = combineSamples(nonlac, lac)
# prior to running next line uncomment the res.filtered line to remove filtering for only genes differentially expr > |1|
computeDEgenes(bovineUnadjusted, outputFile="~/Work/1_Milk/DilutionEffect/DE-Genes/Bovine-FINAL/Bt_norm_DEgenes16K.txt")

# # # HUMAN COLOSTRUM TO UN-ADJUSTED MATURE LACTATION # # #
# Colostrum to mature lactation comparison of FPKM data from Nommsen-Rivers transcriptome manuscript
## NOTE: Prior to loading files replace header line from bottom of sort back to top of file where header should be
## NOTE: Prior to loading files replace header line from bottom of sort back to top of file where header should be
setwd("~/Work/1_Milk/DilutionEffect/Dilution_Outputs/DilutionAdj-HumanFINAL/") # change working directory to input file location
norm_Hs_cflist = list.files("~/Work/1_Milk/DilutionEffect/Dilution_Outputs/DilutionAdj-HumanFINAL/", pattern = "^psct_Col.+[1-2].txt")   # list data files required - colostrum
norm_Hs_mflist = list.files("~/Work/1_Milk/DilutionEffect/Dilution_Outputs/DilutionAdj-HumanFINAL/", pattern = "^norm.+[1-6].txt")  # ditto for unadjusted set - mature
colostrum = loadCounts(norm_Hs_cflist, nreps = 2, computeMean=FALSE)
mature    = loadCounts(norm_Hs_mflist, nreps = 6, computeMean=FALSE)
humanUnadjusted = combineSamples(colostrum, mature)
# prior to running next line uncomment the res.filtered line to remove filtering for only genes differentially expr > |1|
computeDEgenes(humanUnadjusted, repsCondition1= 2, repsCondition2= 6, outputFile="~/Work/1_Milk/DilutionEffect/DE-Genes/Human-FINAL/Hs_norm_DEgenes60K.txt")


