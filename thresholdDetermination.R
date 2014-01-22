# thresholdDetermination.R
# K. Beck
# PhD Candidate, UC Davis Genome Center
# www.korflab.ucdavis.edu

# # # # # # 
# 
# OBJECTIVES
#
# # # # # # 
# 1. Determine if dilution effect adjustment is necessary based on Kolmogorov-Smirnov and Wilcox test
# 2. If so, calculate an appropriate threshold for high abundance genes
# Note: The threshold determined for each replicate and subsequently used as input for dilutioneffect.pl


# # # # # # 
# 
# FUNCTIONS
#
# # # # # # 

options(scipen= 10) # Remove scientific notation from output

determineThreshold = 
    function (data, column = 2, q = 0.9995) {
        # Store the column of interest in a temporary variable
        tmp = as.numeric(data[,column])
        
        # Determine the value along the cumulative distribution function at the quantile specified
        qq = quantile(tmp, probs = q, names = FALSE)
        totalExpr = sum(tmp) # get total expression
        threshold = qq/totalExpr # calculate the threshold relating to the quantile
        
        nHighAbundance = length(tmp[tmp > qq])
        highGenes = tmp[tmp > qq]
        lowestValHighAbundance = min(highGenes)
        
        
        # Estimate the dilution adjustment
        aboveQuantile = sum(tmp[tmp > qq]) 
        adjFactor = 1/(1-(aboveQuantile/totalExpr))
        # cat("The estimated dilution adjustment is", adjFactor, "\n")
        
        # cat("The quantile determined threshold is", threshold, "\n")
        
        myVec = c(threshold, lowestValHighAbundance, nHighAbundance, adjFactor) ## TODO get this to export as a matrix or write to a file
        return(myVec)
    }


generateThresholdObject =
    function (data) {
        threshBatch = data.frame(t(sapply(2:ncol(data), function(x) determineThreshold(data, column = x, q= q_defined))))
        colnames(threshBatch) = c("Threshold", "ThresholdIndex", "nHighAbundanceGenes", "estAdjFactor")
        rownames(threshBatch) = colnames(data)[2:ncol(data)]
        
        return(threshBatch)
    }


## 0. Load expression data ---------------------
# Check to see if files are declared on command line or console
if (exists("file1") & exists("file2")) {
    cat("File names specified in console.\n")
} else {
    args = commandArgs(trailingOnly=TRUE)
    if(length(args) < 2) stop("Usage: Files undefined, please specify expression data. 
    If working from command line use
          $ Rscript thresholdDetermination.R <file1.txt> <file2.txt>
    If working in Rstudio console use
          > file1 = \"<file1.txt>\"
          > file2 = \"<file2.txt>\"" )
    file1 = args[1]
    file2 = args[2]
}

cat(  "Expression data 1 specified: ", file1,
    "\nExpression data 2 specified: ", file2, "\n")
expr1 = read.table(file1, stringsAsFactors=FALSE, header=TRUE)
expr2 = read.table(file2, stringsAsFactors=FALSE, header=TRUE)


## 1. Determine if distributions of two transcriptomes are similar ---------------------
# Determine the difference between distributions of both samples using statistical tests
# If p < 0.05 in both Wilcoxon and KS test then you can confidently reject the null hypothesis and 
# conclude that the distributions are different enough to warrant a dilution adjustment

# Mann-Whitney-Wilcoxon Test
w = wilcox.test(x=expr1$mean, y=expr2$mean)
# Null hypothesis: expression data from stage 1 and stage 2 are identical populations
# if p < 0.05 we can reject the null hypothesis, i.e. the two samples are from nonidentical distributions
# ref: http://www.r-tutor.com/elementary-statistics/non-parametric-methods/mann-whitney-wilcoxon-test

# Kolmogorov-Smirnov Test
k = ks.test(x=expr1$mean,     y=expr2$mean)
# interpret p-value as above i.e. p-value < 0.05 means the two samples are from nonidentical distributions
# D = K-S test statistic aka the maximum difference between the x & y cumulative distribution function, higher value means larger difference

# Check if distributions are different and exit if they are not statistically different
if (k$p.value < 0.05 & w$p.value < 0.05) {
    cat("Distributions are statistically different.\nProceeding to calculate threshold.\n")
} else {
    warning("Distributions are not statistically different and dilution adjustment is not required.\nWilcoxon Test p-value = ", w$p.value,
        "\nKS Test p-value = ", k$p.value)
    stop("Exiting script")
}


## 1. Determine thresholds using quantile for top expression ---------------------
# Check for user defined quantile
if (exists("q_defined")) {
    cat("Quantile specified in console.\n")
    if(is.numeric(q_defined) == FALSE) {
        warning("q_defined must be numeric\nPlease specify floating point number i.e. q_defined = 0.9995\n")
    }
} else if (length(args) > 2) {
    q_defined = as.numeric(args[3])
    cat("User defined quantile: ", q_defined, "\n")
} else {
    q_defined = 0.9995
    cat("Using default quantile:", q_defined, "\n")
}



thresholdsFile1 = generateThresholdObject(expr1)
thresholdsFile2 = generateThresholdObject(expr2)

write.table(thresholdsFile1, file="Results/thresholdsFile1.txt", sep = "\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
write.table(thresholdsFile2, file="Results/thresholdsFile2.txt", sep = "\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

