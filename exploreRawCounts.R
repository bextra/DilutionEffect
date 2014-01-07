# exploreRawCounts.R
# K. Beck
# Jan. 15, 2013

# Load data with exprCountsLoad.R
options(scipen= 10)

## Plot expression data for various milk stages
# Scatter plots
plot(sort(NRcolostrum$mean, decreasing=FALSE))
plot(sort(NRtransitional$mean, decreasing=FALSE))
plot(sort(NRmature$mean, decreasing=FALSE))
plot(NRmature$mean)

# Histograms
hist(NRtransitional$mean)
hist(NRmature$mean)

# Density plots
plot(density(NRcolostrum$mean))
plot(density(log10(NRtransitional$mean)))
plot(density(log10(NRmature$mean)))
plot(density(NRmature$mean))

# Box plot
boxplot(log10(NRcolostrum$mean), log10(NRtransitional$mean), log10(NRmature$mean),
        pars = list(boxwex = 1.5, staplewex = 0.5, outwex = 0.5)
)

# Multi-dimensional Scaling Plot
source("http://www.bioconductor.org/biocLite.R")
biocLite("edgeR")
# run this once you get the CDS
# see if there is a similar step in the DESeq package



## Get basic descriptive statistics
# Mean
mean(NRcolostrum$mean)
mean(NRtransitional$mean)
mean(NRmature$mean)

# Median
median(NRcolostrum$mean)
median(NRtransitional$mean)
median(NRmature$mean)

## Subset data to remove lowly abundant genes
shortC = subset(NRcolostrum$mean, NRcolostrum$mean > 3)
shortT = subset(NRtransitional$mean, NRtransitional$mean > 3)
shortM = subset(NRmature$mean, NRmature$mean > 10)

## Explore subsetted data
hist(shortC)
plot(density(shortM))
plot()

## Other descriptive plots to help determine distribution
hist(lactation$mean, main = "Histogram  of observed data")
plot(density(lactation$mean), main = "Density estimate of data")
plot(ecdf(lactation$mean), main = "Empirical cumulative distribution function")

z.norm = (lactation$mean-mean(lactation$mean))/sd(lactation$mean)
qqnorm(z.norm)
abline(0,1)

## Model data to determine the distribution model analytically using:
# - mean
# - variability
# - skewness (type of curve)
# - kurtosis (type of curve)
library("fBasics")

skewness(lactation$mean)
# [1] 90.54249
# attr(,"method")
# [1] "moment"
kurtosis(lactation$mean)
# [1] 9226.221
# attr(,"method")
# [1] "excess"
skewness(NRmature$mean)
# [1] 158.325
# attr(,"method")
# [1] "moment"
kurtosis(NRmature$mean)
# [1] 27416.5
# attr(,"method")
# [1] "excess"

# Skewness: indicator used in distribution analysis as a sign of asymmetry and deviation from a normal distribution. 
# Interpretation: 
# Skewness > 0 - Right skewed distribution - most values are concentrated on left of the mean, with extreme values to the right.
# Skewness < 0 - Left skewed distribution - most values are concentrated on the right of the mean, with extreme values to the left.
# Skewness = 0 - mean = median, the distribution is symmetrical around the mean.
# 
# Kurtosis - indicator used in distribution analysis as a sign of flattening or "peakedness" of a distribution. 
# Interpretation: 
# Kurtosis > 3 - Leptokurtic distribution, sharper than a normal distribution, with values concentrated around the mean and thicker tails. This means high probability for extreme values.
# Kurtosis < 3 - Platykurtic distribution, flatter than a normal distribution with a wider peak. The probability for extreme values is less than for a normal distribution, and the values are wider spread around the mean.
# Kurtosis = 3 - Mesokurtic distribution - normal distribution for example.
# *High values indicate a high and sharp peak indicating that a few extreme differences from the mean drive an increase in variability
# Ref: http://www.tc3.edu/instruct/sbrown/stat/shape.htm

