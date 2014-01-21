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
        
        cat("There are", length(tmp[tmp > qq]), "high abundance genes.\n")
        highGenes = tmp[tmp > qq]
        cat("The smallest of which is", min(highGenes), "\n")
        
        # Estimate the dilution adjustment
        aboveQuantile = sum(tmp[tmp > qq]) 
        adjFactor = 1/(1-(aboveQuantile/totalExpr))
        cat("The estimated dilution adjustment is", adjFactor, "\n")
        
        cat("The quantile determined threshold is", threshold, "\n")
        
        return(threshold)
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
threshExpr1 = sapply(2:ncol(expr1), function(x) determineThreshold(expr1, column = x, q= 0.9995))
threshExpr2 = sapply(2:ncol(expr2), function(x) determineThreshold(expr2, column = x, q= 0.9995))
# These thresholds are input for dilution_effect.pl and correspond to each replicate


# # # # # # 
# 
# TESTING
#
# # # # # # 

# TODO maybe report the high abundance genes for comparison between human lactation
intersect(highGenesHs$GeneID, highGenesHsC$GeneID) # 11 shared high abundance genes between colostrum and mature


########## Quantile and other statistics research ##########

# Explore quantile values
# returns the x value at this portion of the cumultive distribution curve
# quantile(lactation$mean, probs = c(0.05, 0.95, 0.999, .999908, 1), names=TRUE) # cow
# quantile(NRmature$mean, probs = c(0.05, 0.95, 0.999, .999908, 1), names=TRUE)  # human


# Get just the p-value like this (or store the wilcoxon results in a variable and use the brackets after that)
# wilcox.test(x=NRcolostrum$mean,    y=NRmature$mean)["p.value"] 
# This needs m and n values, but how do you determine those?
# plot(dwilcox(lactation$mean))


# Description of statistical modeling
# d = probability density values (pdf = probability distribution function)
# returns height of the distribution curve at specified (input) x value(s)
# in some functions you can specify the mean and sd to tailor it to your curve
# otherwise it treats x as the z score
# p = cumulative probabilities (cdf = cumulative distribution function)
# returns the area below the curve up until the value specified
# i.e. pnorm(1) returns the cumulative probabilities up until 1 in the curve specified
# defaults to mean = 0 and sd = 1
# to find the area above that point pnorm(1, lower.tail=FALSE)
# q = quantile values
# gives you quantiles or critical values
# qnorm(0.95) # p = 0.05, one-tailed (upper) will return 1.644854
# qnorm(0.95, mean = 5, sd =1) # can shift the curve to fit your data like so
# r = random numbers from the distribution
# ref: http://ww2.coastal.edu/kingw/statistics/R-tutorials/prob.html
