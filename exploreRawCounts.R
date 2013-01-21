# exploreRawCounts.R
# K. Beck
# Jan. 15, 2012

# Load data with exprCountsLoad.R
options(scipen= 10)
plot(density(NRcolostrum$mean))

plot(sort(NRcolostrum$mean, decreasing=FALSE))
plot(sort(NRtransitional$mean, decreasing=FALSE))
plot(sort(NRmature$mean, decreasing=FALSE))

hist(shortC)
hist(NRtransitional$mean)
hist(NRmature$mean)

mean(NRcolostrum$mean)
mean(NRtransitional$mean)
mean(NRmature$mean)
median(NRcolostrum$mean)
median(NRtransitional$mean)
median(NRmature$mean)

shortC = subset(NRcolostrum$mean, NRcolostrum$mean > 3)
shortT = subset(NRtransitional$mean, NRtransitional$mean > 3)
shortM = subset(NRmature$mean, NRmature$mean > 3)

plot(density(shortC))
plot(density(log10(NRtransitional$mean)))
plot(density(log10(NRmature$mean)))
boxplot(log10(NRcolostrum$mean), log10(NRtransitional$mean), log10(NRmature$mean),
        pars = list(boxwex = 1.5, staplewex = 0.5, outwex = 0.5)
        )
?boxplot
plot(NRmature$mean)
mean(NRcolostrum$mean, NRtransitional$mean, NRmature$mean)

###########
# Mann-Whitney-Wilcoxon Test
# Null hypothesis: non-lactating and lactating are identical populations
# if p < 0.05 we can reject the null hypothesis

wilcox.test(x=NRcolostrum$mean, y=NRmature$mean)
# returns p < 2.2e-16
# therefore non-lac and lac are nonidentical populations
wilcox.test(x=prepuberty$mean, y=lactation$mean)
# returns p < 2.2e-16
# therefore colostrum and mature are nonidentical populations
# ref: http://www.r-tutor.com/elementary-statistics/non-parametric-methods/kruskal-wallis-test


ks.test(x=prepuberty$mean, y=lactation$mean)
ks.test(x=NRcolostrum$mean, y=NRmature$mean)
?ks.test
plot(dwilcox())


hist(lactation$mean, main = "Histogram  of observed data")
plot(density(lactation$mean), main = "Density estimate of data")
plot(ecdf(lactation$mean), main = "Empirical cumulative distribution function")

z.norm = (lactation$mean-mean(lactation$mean))/sd(lactation$mean)
qqnorm(z.norm)
abline(0,1)

### Modeling data
# determine the distribution model analytically using
# - mean
# - variability
# - skewness (type of curve)
# - kurtosis (type of curve)

install.packages("fBasics")
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
# 
# 
# 

