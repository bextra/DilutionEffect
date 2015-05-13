# thresholdDetermination.R
# K. Beck
# PhD Candidate, UC Davis Genome Center
# www.korflab.ucdavis.edu

# # # # # # 
# 
# OBJECTIVES
#
# # # # # # 

# Determine if dilution effect adjustment is necessary based on Kolmogorov-Smirnov and Wilcox test


options(scipen = 10) # Remove scientific notation from output



## 0. Load expression data ---------------------
args = commandArgs(trailingOnly=TRUE)
file1 = args[1]
file2 = args[2]
expr1 = read.table(file1, stringsAsFactors=FALSE, header=TRUE)
expr2 = read.table(file2, stringsAsFactors=FALSE, header=TRUE)


# Determine if average level of expression is different
# Determine the difference between distributions of both samples using statistical tests
# If p < 0.05 in both Wilcoxon and KS test then you can confidently reject the null hypothesis and 
# conclude that the distributions are different enough to warrant a dilution adjustment

# Mann-Whitney-Wilcoxon Test, comparing means
w = suppressWarnings(wilcox.test(x=expr1$mean, y=expr2$mean))

# Null hypothesis: expression data from stage 1 and stage 2 are identical populations
# if p < 0.05 we can reject the null hypothesis, i.e. the two samples are from nonidentical distributions
# ref: http://www.r-tutor.com/elementary-statistics/non-parametric-methods/mann-whitney-wilcoxon-test

# Kolmogorov-Smirnov Test — goodness of fit test, comparing cumulative distributions
k = suppressWarnings(ks.test(x=expr1$mean,     y=expr2$mean))
# interpret p-value as above i.e. p-value < 0.05 means the two samples are from nonidentical distributions
# D = K-S test statistic aka the maximum difference between the x & y cumulative distribution function, higher value means larger difference

cat("\tRunning Mann-Whitney-Wilcoxon test to compare means: P = ",  w$p.value, "\n")
cat("\tRunning Kolmogorov-Smirnov test to compare distributions: P = ",     k$p.value, "\n")



if (w$p.value < 0.05 & k$p.value < 0.05){
    cat("\tAverage expression and underlying distributions are statistically different.\n")
    cat("\tThe dilution adjustment should be applied.\n")
}else if (w$p.value < 0.05){
    cat("\tAverage expression is statistically different but distributions are not.\n")
    cat("\tThere may be some value in applying dilution adjustment.\n")
} else if (k$p.value < 0.05){
    cat("\tUnderlying distributions are significantly different but mean expression is not.\n")
    cat("\tThere may be some value in applying dilution adjustment.\n")
} else{
    cat("\tAverage expression is not statistically different (P > 0.05) and distributions are not signififcantly different.\n")
    cat ("The dilution adjustment is probably not required.\n")
    stop()
}
