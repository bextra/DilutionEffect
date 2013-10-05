# compareDElists.R
# K. Beck
# September 25, 2013


# # # # # # # # # # # # # # #
#
# FUNCTIONS
#
# # # # # # # # # # # # # # #

# Number of genes in each regulation status when dealing with unadjusted data
nDEunadjusted = 
    function (df) {
        # count the number of genes in each regulation status
        up = length(df$id[which(df$adjLog2FoldChange > 1 & (is.na(df$padj) == TRUE | df$padj < 0.05 ))])
        noChange = length(df$id[which(abs(df$adjLog2FoldChange) < 1)])
        down = length(df$id[which((df$adjLog2FoldChange < -1) & (is.na(df$padj) == TRUE | df$padj < 0.05 ))])
        
        nDE = list(up, noChange, down) # store as a list for easy export
        
        # print the values
        cat ("Up, No Change, Down\n")
        cat(unlist(nDE))
        
        return(unlist(nDE)) # return the list of genes in each set
        
    }


# Number of genes in each regulation status when dealing with adjusted data
nDEadjusted = 
    function (unadjDF, adjDF) {
        
        # Create a data frame with adjusted log2 fold change values from unadjusted and adjusted data sets with adj pval
        tt = merge(unadjDF, adjDF, by='id') # combine two data frames
        tt = data.frame(tt$id, tt$adjLog2FoldChange.x, tt$padj.x, tt$adjLog2FoldChange.y, tt$padj.y) # extract desired columns
        names(tt) = c("id", "unadj_FoldChange", "unadj_P", "adj_FoldChange", "adj_P") # rename columns for easy handling
        
        ## Determine the number of genes in each regulation status
        # DOWN
        down = length(which((tt$adj_FoldChange   < -1) & (is.na(tt$adj_P) == TRUE | tt$adj_P < 0.05)))
        
        # DOWN TO NO CHANGE
        down2NoChangeIDs = droplevels(tt$id[(which(
            ((tt$unadj_FoldChange < -1) & (is.na(tt$unadj_P) == TRUE | tt$unadj_P < 0.05)) & # down and significant
                (abs(tt$adj_FoldChange)  <  1)))]) # no change
        down2NoChange = length(down2NoChangeIDs)
        # Note: utilized this calculation when determining number of genes from down regulated to no change.
        # It is important to filter downreg unadjusted genes by p-val because you can state that you are confident that this
        # is their accurate regualation status and that it truly changed after adjustment.
        # This was used opposed to accepting all downregulated genes regardless of statistical significance that changed to no change
        # 2524 in humans # matches previous calculation - safe, # 1732 in bovine
        
        # NO CHANGE
        allNoChangeIDs = droplevels(tt$id[which(abs(tt$adj_FoldChange) < 1)])
        noChange = length(which((allNoChangeIDs %in% down2NoChangeIDs == FALSE))) # how many no change genes are there that weren't switched from down? 
        
        # NO CHANGE TO UP
        noChange2upIDs = droplevels(tt$id[which(
            (abs(tt$unadj_FoldChange) < 1) &
                (tt$adj_FoldChange   > 1 & (is.na(tt$adj_P) == TRUE | tt$adj_P < 0.05)))])
        noChange2up = length(noChange2upIDs)
        # Note: There are genes that were upregulated in the unadjusted set that were not statistically significant
        # but ended up becoming statistically significant and still up in the adjusted set.
        # These were not included here because they are not truly a gene whose regulation status changed
        # These are however captured in the general upregulated list
        # 1761 in human, considerable less than before (formerly 2155)
        # 904 in bovine (formerly 971 which was actually just the uniquely up genes and not accurate)
        
        # UP
        allUpIDs = droplevels(tt$id[which(tt$adj_FoldChange > 1 & (is.na(tt$adj_P) == TRUE | tt$adj_P < 0.05))]) # total up genes matches
        up = length(which((allUpIDs %in% noChange2upIDs == FALSE))) # how many no change genes are there that weren't switched from down? 
        
        # Print out number of genes in each regulation status
        nDE = list(up, noChange2up, noChange, down2NoChange, down)
        cat ("Up, No2Up, No Change, Down2No, Down\n", unlist(nDE))
        
        return(unlist(nDE)) # return the list of genes in each set
        
    }




# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Get number of genes per regulation status
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


# TODO put those numbers into figure function

### --- BOVINE --- ###
## Load data for all 16K genes
unadjAllBt = read.table("~/Work/1_Milk/DilutionEffect/DE-Genes/Bovine-FINAL/Bt_norm_allDEgenes16K.txt", header=TRUE)
adjAllBt   = read.table("~/Work/1_Milk/DilutionEffect/DE-Genes/Bovine-FINAL/Bt_adj_allDEgenes16K.txt" , header=TRUE)

# Get the number of genes for each regulation status
BtUnadjStats = nDEunadjusted(unadjAllBt) # matches
BtAdjStats = nDEadjusted(unadjAllBt, adjAllBt)

# In the unadjusted Bovine Data Set, there were 3,754 downregulated and 2,767 upregulated genes between the puberty and lactation states. 
# In the adjusted Bovine Data Set, there were 2,475 downregulated and 3,701 upregulated genes. 



### --- HUMAN --- ##
unadjAllHs = read.table("~/Work/1_Milk/DilutionEffect/DE-Genes/Human-FINAL/Hs_norm_allDEgenes60K.txt", header=TRUE)
adjAllHs = read.table("~/Work/1_Milk/DilutionEffect/DE-Genes/Human-FINAL/Hs_adj_allDEgenes60K.txt", header = TRUE)

# Get the number of genes for each regulation status
HsUnadjStats = nDEunadjusted(unadjAllHs) # matches
HsAdjStats = nDEadjusted(unadjAllHs, adjAllHs) 

# In the Human Data Set, there were 6,323 downregulated and 1,654 upregulated genes between the colostrum and mature milk states. 
# In the adjusted Human Data Set, there were 3,212 downregulated and 3,808 upregulated genes. 



