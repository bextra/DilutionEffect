# quantileDetermination.R
# K. Beck
# Jan. 15, 2013

# Load data with exprCountsLoad.R
options(scipen= 10)

## How similar are the distributions of non-lactating and lactating?
## Mann-Whitney-Wilcoxon Test
# Null hypothesis: non-lactating and lactating are identical populations
# if p < 0.05 we can reject the null hypothesis
wilcox.test(x=NRcolostrum$mean,    y=NRmature$mean) # returns p < 2.2e-16
wilcox.test(x=NRtransitional$mean, y=NRmature$mean) # returns p < 2.2e-16
wilcox.test(x=NRcolostrum$mean,    y=NRtransitional$mean) # p-value = 0.01624....ask Danielle about this case???
# therefore non-lac and lac are nonidentical populations

wilcox.test(x=prepuberty$mean,     y=lactation$mean) # returns p < 2.2e-16
# therefore colostrum and mature are nonidentical populations
# ref: http://www.r-tutor.com/elementary-statistics/non-parametric-methods/mann-whitney-wilcoxon-test

# This needs m and n values, but how do you determine those?
plot(dwilcox(lactation$mean))

## Kolmogorov-Smirnov Tests
ks.test(x=prepuberty$mean, y=lactation$mean) # p < 2.2e-16
ks.test(x=NRcolostrum$mean, y=NRmature$mean) # p < 2.2e-16


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

# Explore quantile values
# returns the x value at this portion of the curve

quantile(lactation$mean, probs = c(0.05, 0.95, 0.999, .999908, 1), names=TRUE) # cow
quantile(NRmature$mean, probs = c(0.05, 0.95, 0.999, .999908, 1), names=TRUE)  # human

# ref: http://ww2.coastal.edu/kingw/statistics/R-tutorials/prob.html
?quantile

## Determine what the dilution adjustment factor would be at a specified quantile
## COW ##
qq = quantile(lactation$mean, probs = 0.999, names = FALSE)
totalExpr = sum(lactation$mean)
aboveQQ = sum(subset(lactation$mean, lactation$mean > qq)) # 66 high abundance genes
adjFactorBt = 1/(1-(aboveQQ/totalExpr))


## HUMAN ##
qq = quantile(NRmature$mean, probs = 0.999, names = FALSE)
totalExpr = sum(NRmature$mean)
aboveQQ = sum(subset(NRmature$mean, NRmature$mean > qq)) # 57 high abundance genes
adjFactorHs = 1/(1-(aboveQQ/totalExpr))

cat("Cow dilution adj. factor:", adjFactorBt, "\nHuman dilution adj. factor:", adjFactorHs, "\n")